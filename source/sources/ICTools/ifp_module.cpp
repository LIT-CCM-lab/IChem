#include "headers/ICTools/ifp_module.hpp"
#include "headers/ICCalcs/interaction.h"
#include "headers/ICMole/similarity.h"
#include "headers/ICTools/ifp_engine.hpp"
#include "ichemo.h"

using namespace std;
using namespace ICMole;
using namespace IFPInternal;

namespace {

// 3 argument mode: protein, ligand, ligand_ref
void processReferenceMode(const std::string& ligandFile, const std::string& refLigandFile, Interactions& interactions, const IFPOptions& options) {
    
    std::vector<IFPEntry> docked;
    std::vector<IFPEntry> refs;

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
            ligand.getName(),
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

    // Print all IFPs
    for (const auto& d : docked) {
        cout << d.name << '\t' << d.fpString      << '\n'
                << d.name << '\t' << d.fp.toString() << '\n';
    }
    for (const auto& r : refs) {
        cout << r.name << '\t' << r.fpString      << '\n'
                << r.name << '\t' << r.fp.toString() << '\n';
    }

    // Similarities docked vs ref
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

// 2 argument mode: protein + ligand file (single or multi-ligand)
void processLigandMode(const std::string& ligandFile, Complex& complex, Interactions& interactions, const IFPOptions& options) {
    
    MoleReader reader;
    reader.loadNewFile(ligandFile);
    const std::size_t numLigands = reader.getNumMolecules();

    if (numLigands == 1) {
        reader.loadInComplex(complex, MoleType::LIGAND);

        if (complex.getMole(MoleType::LIGAND) == nullptr) {
            throw MoleExcept(9020102, "IChem::IFP", "No ligand found in " + ligandFile);
        }

        Molecule& ligand = *complex.getMole(MoleType::LIGAND);
        ligand.checkMOL2();

        InterResults interRes;
        computeIFPForLigand(interactions, ligand, options, interRes);

        cout << interactions.toString(interRes) << '\n';
        cout << interRes.IFPString << '\n'
             << interRes.IFP.toString() << '\n';
    }
    else {
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

        for (const auto& entry : entries) {
            cout << entry.name << '\t' << entry.fpString      << '\n'
                    << entry.name << '\t' << entry.fp.toString() << '\n';
        }
    }
}

// 4 argument mode: (protein1, ligand1) vs (protein2, ligand2)
void processFourArgMode(const std::string& protein1File, const std::string& ligand1File, const std::string& protein2File, const std::string& ligand2File, const IFPOptions& options) {
    
    std::vector<IFPEntry> dockedEntries = computeIFPsFromFiles(protein1File, ligand1File, options);
    std::vector<IFPEntry> refEntries = computeIFPsFromFiles(protein2File, ligand2File, options);

    const std::size_t numDocked = dockedEntries.size();
    const std::size_t numRef = refEntries.size();

    for (const auto& entry : dockedEntries) {
        cout << entry.name << '\t' << entry.fpString << '\n'
             << entry.name << '\t' << entry.fp.toString() << '\n';
    }
    for (const auto& entry : refEntries) {
        cout << entry.name << '\t' << entry.fpString << '\n'
             << entry.name << '\t' << entry.fp.toString() << '\n';
    }

    Similarity sims(false);

    for (std::size_t i = 0; i < numDocked; ++i) {
        const auto& docked = dockedEntries[i];

        for (std::size_t j = 0; j < numRef; ++j) {
            const auto& ref = refEntries[j];

            if (docked.fpString.size() != ref.fpString.size()) {
                throw MoleExcept(9020105, "IChem::IFP", "Cannot compare fingerprints with different sizes in 4-argument mode");
            }

            sims.setRef(ref.fp);
            sims.setComp(docked.fp);

            cout << docked.name << '\t'
                 << ref.name    << '\t'
                 << sims.Tanimoto() << '\n';
        }
    }
}

} // anonymous namespace

namespace IFPModule {

void runIFP(const std::vector<std::string> &Input_Values, const OptionMap &Opt_Values) {
    const std::size_t inputSize = Input_Values.size();

    if (inputSize < 2 || inputSize > 4) {
        throw MoleExcept(9020101, "IChem::IFP", "Number of parameters must be 2, 3 or 4");
    }

    IFPOptions options = parseIFPOptions(Opt_Values);

    try {
        Residu::loadRules();
        Molecule::loadRules();
        configureResidueRules(options);

        // 4 argument mode: protein1 ligand1 protein2 ligand2
        if (inputSize == 4) {
            const std::string& protein1File = Input_Values[0];
            const std::string& ligand1File  = Input_Values[1];
            const std::string& protein2File = Input_Values[2];
            const std::string& ligand2File  = Input_Values[3];

            processFourArgMode(protein1File, ligand1File, protein2File, ligand2File, options);
            return;
        }

        const std::string& proteinFile = Input_Values[0];
        const std::string& ligandFile = Input_Values[1];
        const bool withReference = (inputSize == 3);

        Complex complex;
        MoleReader proteinReader;

        proteinReader.loadNewFile(proteinFile);
        proteinReader.get_format_file();
        proteinReader.loadInComplex(complex, MoleType::PROTEIN);

        if (complex.getMole(MoleType::PROTEIN) == nullptr) {
            throw MoleExcept(9020102, "IChem::IFP", "No protein found in " + proteinFile);
        }

        Interactions interactions(complex);
        applyOverrides(interactions, options.overrides);

        if (withReference) {
            const std::string& refLigandFile = Input_Values[2];
            processReferenceMode(ligandFile, refLigandFile, interactions, options);
        }
        else {
            processLigandMode(ligandFile, complex, interactions, options);
        }
    }
    catch (MoleExcept& e) {
        e.addTrace("IChem::IFP");
        throw;
    }
}

} // namespace IFPModule
