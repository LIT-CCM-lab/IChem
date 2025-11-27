#include "headers/ICTools/ifp_module.hpp"
#include "headers/ICCalcs/interaction.h"
#include "headers/ICMole/similarity.h"
#include "ichemo.h"

using namespace std;
using namespace ICMole;

namespace {

using OptionMap = IFPModule::OptionMap;

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

inline void applyOverrides(Interactions& interactions, const InteractionOverrides& overrideInteractionParams) {
    // Distances
    if (overrideInteractionParams.has_D_Hb)   interactions.setDist_H(overrideInteractionParams.D_Hb);
    if (overrideInteractionParams.has_D_Hyd)  interactions.setDist_Hyd(overrideInteractionParams.D_Hyd);
    if (overrideInteractionParams.has_D_Io)   interactions.setDist_Ionic(overrideInteractionParams.D_Io);
    if (overrideInteractionParams.has_D_Me)   interactions.setDist_Metal(overrideInteractionParams.D_Me);
    if (overrideInteractionParams.has_D_Ar)   interactions.setDist_Arom(overrideInteractionParams.D_Ar);
    if (overrideInteractionParams.has_D_Pic)  interactions.setDist_PICation(overrideInteractionParams.D_Pic);

    // Angles
    if (overrideInteractionParams.has_a_H)      interactions.setAngl_H(overrideInteractionParams.a_H);
    if (overrideInteractionParams.has_at_H)     interactions.setAngl_Tol_H(overrideInteractionParams.at_H);
    if (overrideInteractionParams.has_a_ArFF)   interactions.setAngl_AromFF(overrideInteractionParams.a_ArFF);
    if (overrideInteractionParams.has_at_ArFF)  interactions.setAngl_Tol_AromFF(overrideInteractionParams.at_ArFF);
    if (overrideInteractionParams.has_a_ArEF)   interactions.setAngl_AromEF(overrideInteractionParams.a_ArEF);
    if (overrideInteractionParams.has_at_ArEF)  interactions.setAngl_Tol_AromEF(overrideInteractionParams.at_ArEF);
    if (overrideInteractionParams.has_a_Pic)    interactions.setAngl_PICation(overrideInteractionParams.a_Pic);
    if (overrideInteractionParams.has_at_Pic)   interactions.setAngl_Tol_PICation(overrideInteractionParams.at_Pic);
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
inline void configureResidueRules(const IFPOptions& options) {
    if (options.includeSolvent) {
        Residu::Rules[MoleType::PROTEIN][ResType::WATER] = MoleType::PROTEIN;
    }
    if (options.includeCofactor) {
        Residu::Rules[MoleType::PROTEIN][ResType::COFACTOR] = MoleType::PROTEIN;
    }
}

void computeIFPForLigand(Interactions& interactions, Molecule& ligand, const IFPOptions& options, InterResults& output) {
   
    // RingDetection not run yet
    if (ligand.firstCycle() == ligand.lastCycle()) {
        ligand.ringPerception();
    }

    interactions.detectInteractions(ligand, output, true, options.oldHydrophobic);
    interactions.genIFP(output, options.layoutIndex);
}

// 3 arguments mode (protein, ligand, ligand_ref)
// We compute IFPs for all reference ligands first, then all query ligands
void processReferenceMode(const std::string& ligandFile, const std::string& refLigandFile, Interactions& interactions, const IFPOptions& options)
{
    std::vector<IFPEntry> docked; // ligandFile
    std::vector<IFPEntry> refs;   // refLigandFile

    MoleReader reader;

    // docked ligands
    reader.loadNewFile(ligandFile);
    while (!reader.isEOF()) {
        Molecule ligand;
        reader.loadNextMolecule(ligand, MoleType::LIGAND);
        ligand.checkMOL2();

        InterResults result;
        computeIFPForLigand(interactions, ligand, options, result);

        docked.push_back(IFPEntry{
            ligand.getName(),   // e.g. "REF" in your current example
            result.IFP,
            result.IFPString
        });
    }

    // reference ligands
    reader.loadNewFile(refLigandFile);
    while (!reader.isEOF()) {
        Molecule ligand;
        reader.loadNextMolecule(ligand, MoleType::LIGAND);
        ligand.checkMOL2();

        InterResults result;
        computeIFPForLigand(interactions, ligand, options, result);

        refs.push_back(IFPEntry{
            ligand.getName(),
            result.IFP,
            result.IFPString
        });
    }

    // 3) Print all IFPs if requested
    if (options.outputBitstring) {
        // docked ligands
        for (const auto& d : docked) {
            cout << d.name << '\t' << d.fpString     << '\n'
                 << d.name << '\t' << d.fp.toString() << '\n';
        }
        // all references
        for (const auto& r : refs) {
            cout << r.name << '\t' << r.fpString     << '\n'
                 << r.name << '\t' << r.fp.toString() << '\n';
        }
    }

    // 4) Similarities: each docked vs ref
    Similarity sims(false);
    for (const auto& d : docked) {
        sims.setRef(d.fp);
        for (const auto& r : refs) {
            sims.setComp(r.fp);
            cout << d.name << '\t'
                 << r.name << '\t'
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
        ligand.checkMOL2();

        InterResults interRes;
        computeIFPForLigand(interactions, ligand, options, interRes);

        cout << interactions.toString(interRes) << '\n';
        if (options.outputBitstring) {
            cout << interRes.IFPString << '\n' // Residue list
                 << interRes.IFP.toString() << '\n'; // IFP
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
            ligand.checkMOL2();

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

std::vector<IFPEntry> computeIFPsFromFiles(const std::string& proteinFile, const std::string& ligandFile, const IFPOptions&  options) {
    
    std::vector<IFPEntry> entries;

    // Protein
    Complex complex;
    MoleReader proteinReader;

    proteinReader.loadNewFile(proteinFile);
    proteinReader.get_format_file();
    proteinReader.loadInComplex(complex, MoleType::PROTEIN);

    if (complex.getMole(MoleType::PROTEIN) == nullptr)
        throw MoleExcept(9020102, "IChem::IFP", "No protein found in " + proteinFile);
    
    // Interactions object for this protein
    Interactions interactions(complex);
    applyOverrides(interactions, options.overrides);

    // 1 ou many ligands
    MoleReader ligandReader;
    ligandReader.loadNewFile(ligandFile);

    while (!ligandReader.isEOF()) {
        
        Molecule ligand;
        ligandReader.loadNextMolecule(ligand, MoleType::LIGAND);
        ligand.checkMOL2();

        InterResults result;
        computeIFPForLigand(interactions, ligand, options, result);

        IFPEntry entry;
        entry.name = ligand.getName();
        entry.fp = result.IFP;
        entry.fpString = result.IFPString;

        entries.push_back(std::move(entry));
    }

    if (entries.empty()) {
        throw MoleExcept(9020102, "IChem::IFP", "No ligand found in " + ligandFile);
    }

    return entries;
}


void processFourArgMode(const std::string& protein1File, const std::string& ligand1File, const std::string& protein2File, const std::string& ligand2File, const IFPOptions&  options) {
    
    // Convention:
    // (protein1, ligand1) = DOCKED (query)
    // (protein2, ligand2) = REF (reference)

    std::vector<IFPEntry> dockedEntries = computeIFPsFromFiles(protein1File, ligand1File, options);
    std::vector<IFPEntry> refEntries = computeIFPsFromFiles(protein2File, ligand2File, options);

    const std::size_t numDocked = dockedEntries.size();
    const std::size_t numRef = refEntries.size();

    if (options.outputBitstring) {
        // DOCKED block
        for (const auto& entry : dockedEntries) {
            std::cout << entry.name << '\t' << entry.fpString      << '\n'
                      << entry.name << '\t' << entry.fp.toString() << '\n';
        }

        // REF block
        for (const auto& entry : refEntries) {
            std::cout << entry.name << '\t' << entry.fpString      << '\n'
                      << entry.name << '\t' << entry.fp.toString() << '\n';
        }
    }

    // Similarity
    Similarity sims(false);

    for (std::size_t i = 0; i < numDocked; ++i) {
        const auto& docked = dockedEntries[i];

        for (std::size_t j = 0; j < numRef; ++j) {
            const auto& ref = refEntries[j];

            // Fingerprint lengths must match
            if (docked.fpString.size() != ref.fpString.size())
                throw MoleExcept( 9020105, "IChem::IFP", "Cannot compare fingerprints with different sizes in 4-argument mode");

            sims.setRef(ref.fp);
            sims.setComp(docked.fp);

            std::cout << docked.name << '\t' << ref.name    << '\t' << sims.Tanimoto() << '\n';
        }
    }
}



} // anonymous namespace

unsigned chooseIFPType_forTests(bool polarOnly, bool extended, bool metalOnly) {
    return chooseIFPType(polarOnly, extended, metalOnly);
}

namespace IFPModule {

void runIFP(const std::vector<std::string> &Input_Values, const OptionMap &Opt_Values) {
    
    const std::size_t inputSize = Input_Values.size();

    if (inputSize < 2 || inputSize > 4)
        throw MoleExcept(9020101, "IChem::IFP", "Number of parameters must be 2, 3 or 4");

    const IFPOptions options = parseIFPOptions(Opt_Values);

    try {
        // Global rules, once per call
        Residu::loadRules();
        Molecule::loadRules();
        configureResidueRules(options);

        // 4 argument mode
        if (inputSize == 4) {
            const std::string& protein1File = Input_Values[0];
            const std::string& ligand1File  = Input_Values[1];
            const std::string& protein2File = Input_Values[2];
            const std::string& ligand2File  = Input_Values[3];

            processFourArgMode(protein1File, ligand1File, protein2File, ligand2File, options);
            return;
        }

        // 2 and 3 argument modes
        const std::string& proteinFile = Input_Values.at(0);
        const std::string& ligandFile = Input_Values.at(1);
        const bool withReference = (inputSize == 3);

        Complex complex;
        MoleReader proteinReader;
        proteinReader.loadNewFile(proteinFile);
        proteinReader.get_format_file();
        proteinReader.loadInComplex(complex, MoleType::PROTEIN);

        if (complex.getMole(MoleType::PROTEIN) == nullptr)
            throw MoleExcept(9020102, "IChem::IFP", "No protein found in " + proteinFile);

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

} // namespace IFPModule