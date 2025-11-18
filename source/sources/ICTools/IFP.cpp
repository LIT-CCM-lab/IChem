#include "headers/ICTools/switch.h"
#include "headers/ICCalcs/interaction.h"
#include "headers/ICMole/similarity.h"

using namespace std;
using namespace ICMole;

void IChemSwitch::helpIFP() const {
    cout
        << "IFP - Interaction FingerPrint\n"
        << "Usage: IChem [options] IFP protein ligand\n"
        << "       IChem [options] IFP protein ligand ligand_ref\n"
        << "\n"
        << "Description:\n"
        << "  Compute interaction fingerprints between a protein and a ligand\n"
        << "  Compare a reference ligand and detect various interaction types\n"
        << "\n"
        << "Options:\n"
        << "  -name N (LIG)  Name of the fingerprint (default: ligand name)\n"
        << "  --polar        Detect and output only polar interactions\n"
        << "  --metal        Detect and output only metal interactions\n"
        << "  --extended     Include in the fingerprint:\n"
        << "                    * Metal/Acceptor interactions\n"
        << "                    * Weak Hydrogen bonds\n"
        << "                    * Pi-Cation interactions\n"
        << "  --ligD         Print all possible ligand interactions\n"
        << "  --bitstringOFF Do not output the bitstring representation\n"
        << "\n"
        << "Thresholds and angles:\n"
        << "  -D_Hb    N (3.5)  Hbond length (Å)\n"
        << "  -D_Hyd   N (4.5)  Hydrophobic length (Å)\n"
        << "  -D_Io    N (4.0)  Ionic length (Å)\n"
        << "  -D_Me    N (2.8)  Metal/Acceptor length (Å)\n"
        << "  -D_Ar    N (5.0)  Aromatic interaction length (Å)\n"
        << "  -D_Pic   N (5.0)  Pi-Cation interaction length (Å)\n"
        << "  -a_H     N (Pi)   Hbond angle (rad)\n"
        << "  -at_H    N (Pi/3) Hbond tolerance (rad)\n"
        << "  -a_ArFF  N (Pi)   Aromatic Face-to-Face angle (rad)\n"
        << "  -at_ArFF N (Pi/6) Aromatic FF tolerance (rad)\n"
        << "  -a_ArEF  N (Pi/2) Aromatic Edge-to-Face angle (rad)\n"
        << "  -at_ArEF N (Pi/3) Aromatic EF tolerance (rad)\n"
        << "  -a_Pic   N (Pi)   Pi-Cation angle (rad)\n"
        << "  -at_Pic  N (Pi/6) Pi-Cation tolerance (rad)\n"
        << "\n"
        << "Example:\n"
        << "  IChem IFP protein.mol2 ligand.mol2\n"
        << "  IChem IFP protein.mol2 ligand.mol2 ligand_ref.mol2\n"
        << "\n"
        << "---------------------------------------------------------------------------\n\n";
}


namespace {

inline unsigned chooseIFPType(bool polarOnly, bool extended, bool metalOnly) {
    // 0 - standard
    // 1 - polar only
    // 2 - extended
    // 3 - extended polar
    // 4 - metal only
    if (!polarOnly && !extended && !metalOnly) return 0;
    if ( polarOnly && !extended)               return 1;
    if (!polarOnly &&  extended)               return 2;
    if ( metalOnly)                            return 4;
    return 3;
}

struct InteractionOverrides {
    bool has_D_Hb   = false; double D_Hb = 0.0;
    bool has_D_Hyd  = false; double D_Hyd = 0.0;
    bool has_D_Io   = false; double D_Io = 0.0;
    bool has_D_Me   = false; double D_Me = 0.0;
    bool has_D_Ar   = false; double D_Ar = 0.0;
    bool has_D_Pic  = false; double D_Pic = 0.0;

    bool has_a_H    = false; double a_H = 0.0;
    bool has_at_H   = false; double at_H = 0.0;
    bool has_a_ArFF = false; double a_ArFF = 0.0;
    bool has_at_ArFF= false; double at_ArFF = 0.0;
    bool has_a_ArEF = false; double a_ArEF = 0.0;
    bool has_at_ArEF= false; double at_ArEF = 0.0;
    bool has_a_Pic  = false; double a_Pic = 0.0;
    bool has_at_Pic = false; double at_Pic = 0.0;
};

inline void applyOverrides(Interactions& interactions, const InteractionOverrides& o) {
    // Distances
    if (o.has_D_Hb)   interactions.setDist_H(o.D_Hb);
    if (o.has_D_Hyd)  interactions.setDist_Hyd(o.D_Hyd);
    if (o.has_D_Io)   interactions.setDist_Ionic(o.D_Io);
    if (o.has_D_Me)   interactions.setDist_Metal(o.D_Me);
    if (o.has_D_Ar)   interactions.setDist_Arom(o.D_Ar);
    if (o.has_D_Pic)  interactions.setDist_PICation(o.D_Pic);

    // Angles
    if (o.has_a_H)      interactions.setAngl_H(o.a_H);
    if (o.has_at_H)     interactions.setAngl_Tol_H(o.at_H);
    if (o.has_a_ArFF)   interactions.setAngl_AromFF(o.a_ArFF);
    if (o.has_at_ArFF)  interactions.setAngl_Tol_AromFF(o.at_ArFF);
    if (o.has_a_ArEF)   interactions.setAngl_AromEF(o.a_ArEF);
    if (o.has_at_ArEF)  interactions.setAngl_Tol_AromEF(o.at_ArEF);
    if (o.has_a_Pic)    interactions.setAngl_PICation(o.a_Pic);
    if (o.has_at_Pic)   interactions.setAngl_Tol_PICation(o.at_Pic);
}

struct IFPOptions {
    std::string fingerprintName;
    bool polarOnly = false;
    bool metalOnly = false;
    bool extended = false;
    bool includeSolvent = true;
    bool includeCofactor = true;
    bool oldHydrophobic  = true;
    bool ligandDebug = false; // I should remove
    bool outputBitstring = true;
    unsigned layoutIndex = 0;
    InteractionOverrides overrides;
};

struct IFPEntry {
    std::string name;
    Fingerprint fp;
    std::string fpString;
};

using OptionMap = std::map<std::string, std::vector<std::string>>;

// Parse Opt_Values into an IFPOptions struct
IFPOptions parseIFPOptions(const OptionMap& optionsValues) {
    IFPOptions options;

    for (const auto& [key, vals] : optionsValues) {
        const std::string value = vals.empty() ? std::string{} : vals.front();

        if      (key == "-name")          options.fingerprintName = value;
        else if (key == "--polar")        options.polarOnly       = true;
        else if (key == "--metal")        options.metalOnly       = true;
        else if (key == "--extended")     options.extended        = true;
        else if (key == "--solvent")      options.includeSolvent  = false;
        else if (key == "--cofactor")     options.includeCofactor = false;
        else if (key == "--newH")         options.oldHydrophobic  = false;
        else if (key == "--ligD")         options.ligandDebug     = true;   // Currently unused
        else if (key == "--bitstringOFF") options.outputBitstring = false;

        // distances
        else if (key == "-D_Hb")   { options.overrides.has_D_Hb   = true; options.overrides.D_Hb   = std::stod(value); }
        else if (key == "-D_Hyd")  { options.overrides.has_D_Hyd  = true; options.overrides.D_Hyd  = std::stod(value); }
        else if (key == "-D_Io")   { options.overrides.has_D_Io   = true; options.overrides.D_Io   = std::stod(value); }
        else if (key == "-D_Me")   { options.overrides.has_D_Me   = true; options.overrides.D_Me   = std::stod(value); }
        else if (key == "-D_Ar")   { options.overrides.has_D_Ar   = true; options.overrides.D_Ar   = std::stod(value); }
        else if (key == "-D_Pic")  { options.overrides.has_D_Pic  = true; options.overrides.D_Pic  = std::stod(value); }

        // angles
        else if (key == "-a_H")      { options.overrides.has_a_H    = true; options.overrides.a_H    = std::stod(value); }
        else if (key == "-at_H")     { options.overrides.has_at_H   = true; options.overrides.at_H   = std::stod(value); }
        else if (key == "-a_ArFF")   { options.overrides.has_a_ArFF = true; options.overrides.a_ArFF = std::stod(value); }
        else if (key == "-at_ArFF")  { options.overrides.has_at_ArFF= true; options.overrides.at_ArFF= std::stod(value); }
        else if (key == "-a_ArEF")   { options.overrides.has_a_ArEF = true; options.overrides.a_ArEF = std::stod(value); }
        else if (key == "-at_ArEF")  { options.overrides.has_at_ArEF= true; options.overrides.at_ArEF= std::stod(value); }
        else if (key == "-a_Pic")    { options.overrides.has_a_Pic  = true; options.overrides.a_Pic  = std::stod(value); }
        else if (key == "-at_Pic")   { options.overrides.has_at_Pic = true; options.overrides.at_Pic = std::stod(value); }
    }

    options.layoutIndex = chooseIFPType(options.polarOnly, options.extended, options.metalOnly);
    return options;
}

// Apply residue rules once
inline void configureResidueRules(const IFPOptions& opt) {
    if (opt.includeSolvent) {
        Residu::Rules[MoleType::PROTEIN][ResType::WATER] = MoleType::PROTEIN;
    }
    if (opt.includeCofactor) {
        Residu::Rules[MoleType::PROTEIN][ResType::COFACTOR] = MoleType::PROTEIN;
    }
}

void computeIFPForLigand(Interactions& interactions, Molecule& ligand, const IFPOptions& options, InterResults& output) {
   
    // RingDetection not run yet
    // if (ligand.firstCycle() == ligand.lastCycle()) {
    //     ligand.ringPerception();
    // }

    interactions.detectInteractions(ligand, output, true, options.oldHydrophobic);
    interactions.genIFP(output, options.layoutIndex);
}

// 3 arguments mode (protein, ligand, ligand_ref)
// We compute IFPs for all reference ligands first, then all query ligands
void processReferenceMode(const std::string& ligandFile, const std::string& refLigandFile, Interactions& interactions, const IFPOptions& options) {
    std::vector<IFPEntry> entries;
    entries.reserve(128);

    MoleReader reader;
    std::size_t numReferences = 0;

    // Reference ligands
    reader.loadNewFile(refLigandFile);
    while (!reader.isEOF()) {
        Molecule ligand;
        reader.loadNextMolecule(ligand, MoleType::LIGAND);
        //delete
        ligand.ringPerception();

        InterResults result;
        computeIFPForLigand(interactions, ligand, options, result);

        entries.push_back(IFPEntry{
            ligand.getName(),
            result.IFP,
            result.IFPString
        });

        ++numReferences;
    }

    // Query ligands
    reader.loadNewFile(ligandFile);
    while (!reader.isEOF()) {
        Molecule ligand;
        reader.loadNextMolecule(ligand, MoleType::LIGAND);
        //delete
        ligand.ringPerception();

        InterResults result;
        computeIFPForLigand(interactions, ligand, options, result);

        entries.push_back(IFPEntry{
            ligand.getName(),
            result.IFP,
            result.IFPString
        });
    }

    const std::size_t total = entries.size();

    // print per-ligand IFP bitstrings if requested
    if (options.outputBitstring) {
        for (const auto& entry : entries) {
            cout << entry.name << '\t' << entry.fpString << '\n'
                 << entry.name << '\t' << entry.fp.toString() << '\n';
        }
    }

    // similarity matrix:
    // refs = [0..numReferences-1], compare each ref i with j from i..total-1 using Tanimoto
    Similarity sims(false);
    for (std::size_t i = 0; i < numReferences; ++i) {
        sims.setRef(entries[i].fp);
        for (std::size_t j = i; j < total; ++j) {
            sims.setComp(entries[j].fp);
            cout << entries[i].name << '\t'
                 << entries[j].name << '\t'
                 << sims.Tanimoto() << '\n';
        }
    }
}

// 2 arguments mode: protein + ligand file (single or multi ligand)
void processLigandMode(const std::string& ligandFile, Complex& complex, Interactions& interactions, const IFPOptions& options) {
    
    MoleReader reader;
    reader.loadNewFile(ligandFile);
    const std::size_t numLigands = reader.getNumMolecules();

    if (numLigands == 1) {
        // Single ligand, loaded in the complex
        reader.loadInComplex(complex, MoleType::LIGAND);

        if (complex.getMole(MoleType::LIGAND) == nullptr)
            throw MoleExcept(9020102, "IChem::IFP", "No ligand found in " + ligandFile);

        Molecule& ligand = *complex.getMole(MoleType::LIGAND);

        InterResults interRes;
        computeIFPForLigand(interactions, ligand, options, interRes);

        cout << interactions.toString(interRes) << '\n';
        if (options.outputBitstring) {
            cout << interRes.IFPString << '\n'
                 << interRes.IFP.toString() << '\n';
        }
    }
    else {
        // Multi-ligand file, not loaded in the complex
        reader.loadNewFile(ligandFile);

        std::vector<IFPEntry> entries;
        entries.reserve(numLigands);

        while (!reader.isEOF()) {
            Molecule ligand;
            reader.loadNextMolecule(ligand, MoleType::LIGAND);

            InterResults res;
            computeIFPForLigand(interactions, ligand, options, res);

            entries.push_back(IFPEntry{
                ligand.getName(),
                res.IFP,
                res.IFPString
            });
        }

        if (options.outputBitstring) {
            for (const auto& entry : entries) {
                cout << entry.name << '\t' << entry.fpString << '\n'
                     << entry.name << '\t' << entry.fp.toString() << '\n';
            }
        }
    }
}

}

void IChemSwitch::IFP() const {

    const std::size_t inputSize = Input_Values.size();
    const bool withReference = (inputSize == 3);

    if (inputSize != 2 && !withReference)
        throw MoleExcept(9020101, "IChem::IFP", "Number of parameters must be 2 or 3");

    const std::string& proteinFile = Input_Values.at(0);
    const std::string& ligandFile  = Input_Values.at(1);

    // Parse CLI options once into a config struct
    const IFPOptions options = parseIFPOptions(Opt_Values);

    try {
        // set residue rules
        configureResidueRules(options);

        // load protein into complex
        Complex complex;
        MoleReader proteinReader;

        proteinReader.loadNewFile(proteinFile);
        proteinReader.get_format_file();
        proteinReader.loadInComplex(complex, MoleType::PROTEIN);

        if (complex.getMole(MoleType::PROTEIN) == nullptr)
            throw MoleExcept(9020102, "IChem::IFP", "No protein found in " + proteinFile);

        // One interactions engine per protein
        Interactions interactions(complex);
        applyOverrides(interactions, options.overrides);

        if (withReference) {
            const std::string& refLigandFile = Input_Values.at(2);
            processReferenceMode(ligandFile, refLigandFile, interactions, options);
        } else {
            processLigandMode(ligandFile, complex, interactions, options);
        }
    }
    catch (MoleExcept& exception) {
        exception.addTrace("IChem::IFP");
        throw;
    }
}
