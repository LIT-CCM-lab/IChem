#include "headers/ICTools/ifp_api.hpp"
#include "headers/ICTools/ifp_engine.hpp"
#include "headers/ICMole/moledata.h"

#include "headers/ICCalcs/interaction.h"
#include "headers/ICMole/similarity.h"
#include "ichemo.h"

#include <stdexcept> // std::runtime_error
#include <sstream>

using namespace std;
using namespace ICMole;

namespace {

// Convert IFPConfig -> IFPOptions (via the same parser used by the CLI)
IFPInternal::IFPOptions makeIFPOptionsFromConfig(const IFPAPI::IFPConfig& cfg) {
    using namespace IFPInternal;

    OptionMap optMap;

    // One at the time
    if (cfg.all) {
        optMap["--all"] = {};
    }
    if (cfg.basic) {
        optMap["--basic"] = {};
    }
    if (cfg.weakh) {
        optMap["--weakh"] = {};
    }
    if (cfg.picat) {
        optMap["--picat"] = {};
    }
    if (cfg.metal) {
        optMap["--metal"] = {};
    }
    if (cfg.old_layout) {
        optMap["--old"] = {};
    }
   

    // Rules
    if (!cfg.includeSolvent) {
        optMap["--solvent"] = {};
    }
    if (!cfg.includeCofactor) {
        optMap["--cofactor"] = {};
    }
    if (!cfg.oldHydrophobic) {
        optMap["--newH"] = {};
    }

    // Helper to add numeric options if set
    auto setDoubleOpt = [&](const std::string& key, const std::optional<double>& val) {
        if (val) {
            std::ostringstream oss;
            oss.setf(std::ios::fixed, std::ios::floatfield);
            oss.precision(3);
            oss << *val;
            optMap[key] = {oss.str()};
        }
    };

    // Distances max
    setDoubleOpt("-D_Hb",  cfg.D_Hb);
    setDoubleOpt("-D_Hyd", cfg.D_Hyd);
    setDoubleOpt("-D_Io",  cfg.D_Io);
    setDoubleOpt("-D_Me",  cfg.D_Me);
    setDoubleOpt("-D_Ar",  cfg.D_Ar);
    setDoubleOpt("-D_Pic", cfg.D_Pic);
    setDoubleOpt("-D_WHb", cfg.D_WHb);

    // Distances min
    setDoubleOpt("-d_Hb",  cfg.d_Hb);
    setDoubleOpt("-d_Hyd", cfg.d_Hyd);
    setDoubleOpt("-d_Io",  cfg.d_Io);
    setDoubleOpt("-d_Me",  cfg.d_Me);
    setDoubleOpt("-d_Ar",  cfg.d_Ar);
    setDoubleOpt("-d_Pic", cfg.d_Pic);
    setDoubleOpt("-d_WHb", cfg.d_WHb);

    // Angles
    setDoubleOpt("-a_H",     cfg.a_H);
    setDoubleOpt("-at_H",    cfg.at_H);
    setDoubleOpt("-a_ArFF",  cfg.a_ArFF);
    setDoubleOpt("-at_ArFF", cfg.at_ArFF);
    setDoubleOpt("-a_ArEF",  cfg.a_ArEF);
    setDoubleOpt("-at_ArEF", cfg.at_ArEF);
    setDoubleOpt("-a_Pic",   cfg.a_Pic);
    setDoubleOpt("-at_Pic",  cfg.at_Pic);

    // We let the parser build IFPOptions and do all checks
    IFPOptions options = parseIFPOptions(optMap);
    return options;
}


static std::string interactionTypeToString(unsigned int interactionType)
{
    switch (interactionType) {
        case InterType::UNDEFINED:
            return "Undefined";

        case InterType::HBOND_PROT:
            return "HBond_protein";
        case InterType::HBOND_LIG:
            return "HBond_ligand";

        case InterType::IONIC_PROT:
            return "Ionic_protein";
        case InterType::IONIC_LIG:
            return "Ionic_ligand";

        case InterType::HYDROPHOBIC:
            return "Hydrophobic";

        case InterType::METAL:
            return "Metal";

        case InterType::ARFACEFACE:
            return "Aromatic_FaceToFace";
        case InterType::AREDGEFACE:
            return "Aromatic_EdgeToFace";

        case InterType::PICATION:
            return "PiCation";

        case InterType::WHBOND_LIG:
            return "WeakHBond_ligand";
        case InterType::WHBOND_PROT:
            return "WeakHBond_protein";

        default: {
            std::ostringstream oss;
            oss << "Unknown(" << interactionType << ")";
            return oss.str();
        }
    }
}

// Extract all interactions for one ligand into API records
static IFPAPI::LigandInteractions buildLigandInteractions(const ICMole::Molecule& ligand, const ICMole::InterResults& interRes) {
    
    IFPAPI::LigandInteractions out;
    out.ligand_name = ligand.getName();

    const std::size_t n = interRes.listInters.size();
    out.interactions.reserve(n);

    for (std::size_t i = 0; i < n; ++i) {
        const ICMole::InterPoint& ip = interRes.listInters[i];
        IFPAPI::InteractionRecord rec;

        // Type of interaction
        rec.type_interaction = interactionTypeToString(ip.interaction);

        // Protein side
        if (ip.Prot_Ref) {
            rec.atom_prot = ip.Prot_Ref->getName();
            rec.id_atom_prot = ip.Prot_Ref->getNum();
            rec.residue = ip.Prot_Ref->getResiduName();
            rec.chain = ip.Prot_Ref->getChainName();

        } else {
            rec.atom_prot.clear();
            rec.id_atom_prot = -1;
            rec.residue.clear();
            rec.chain.clear();
        }

        // Ligand side
        if (ip.Lig_Ref) {
            rec.atom_lig = ip.Lig_Ref->getName();
            rec.id_atom_lig = ip.Lig_Ref->getNum();
        } else {
            rec.atom_lig.clear();
            rec.id_atom_lig = -1;
        }

        // distance
        rec.distance = ip.dist;

        out.interactions.push_back(std::move(rec));
    }

    return out;
}


// For 2 argument mode, with possible multi-ligand file
static std::vector<IFPAPI::LigandInteractions> computeInteractionsForLigandFile(const std::string& proteinFile, const std::string& ligandFile, const IFPInternal::IFPOptions& options) {
    
    using namespace IFPInternal;

    std::vector<IFPAPI::LigandInteractions> results;

    // Rules
    Residu::loadRules();
    Molecule::loadRules();
    configureResidueRules(options);

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

    MoleReader ligandReader;
    ligandReader.loadNewFile(ligandFile);

    while (!ligandReader.isEOF()) {
        Molecule ligand;
        ligandReader.loadNextMolecule(ligand, MoleType::LIGAND);
        ligand.checkMOL2();

        InterResults interRes;
        computeIFPForLigand(interactions, ligand, options, interRes);

        results.push_back(buildLigandInteractions(ligand, interRes));
    }

    if (results.empty()) {
        throw MoleExcept(9020102, "IChem::IFP", "No ligand found in " + ligandFile);
    }

    return results;
}

// For fingerprints: reuse computeIFPsFromFiles
static std::vector<IFPAPI::LigandFingerprint> buildLigandFingerprints(const std::vector<IFPInternal::IFPEntry>& entries) {
    
    std::vector<IFPAPI::LigandFingerprint> out;
    out.reserve(entries.size());

    for (const auto& entry : entries) {
        IFPAPI::LigandFingerprint ligandFingerprint;
        ligandFingerprint.ligand_name = entry.name;
        ligandFingerprint.residues = entry.fpString;
        ligandFingerprint.bitstring = entry.fp.toString();
        out.push_back(std::move(ligandFingerprint));
    }

    return out;
}

} // anonymous namespace

namespace IFPAPI {

std::vector<LigandInteractions> compute_ifp_interactions(const std::string& protein_file, const std::string& ligand_file, const IFPConfig& cfg) {
    
    try {
        auto options = makeIFPOptionsFromConfig(cfg);
        return computeInteractionsForLigandFile(protein_file, ligand_file, options);
    }
    catch (MoleExcept& e) {
        e.addTrace("IFPAPI::compute_ifp_interactions");
        std::ostringstream oss;
        oss << "IChem::IFP error " << e.getCode() << " - " << e.getData();
        throw std::runtime_error(oss.str());
    }
}

std::vector<LigandFingerprint> compute_ifp_fingerprints(const std::string& protein_file, const std::string& ligand_file, const IFPConfig& cfg) {
    
    using namespace IFPInternal;

    try {
        auto options = makeIFPOptionsFromConfig(cfg);

        Residu::loadRules();
        Molecule::loadRules();
        configureResidueRules(options);

        auto entries = computeIFPsFromFiles(protein_file, ligand_file, options);
        return buildLigandFingerprints(entries);
    }
    catch (MoleExcept& e) {
        e.addTrace("IFPAPI::compute_ifp_fingerprints");
        std::ostringstream oss;
        oss << "IChem::IFP error " << e.getCode() << " - " << e.getData();
        throw std::runtime_error(oss.str());
    }
}

std::vector<TanimotoScore> compute_ifp_tanimoto(const std::string& protein_file, const std::string& ligand_file, const std::string& reference_file, const IFPConfig& cfg) {
    
    using namespace IFPInternal;

    try {
        auto options = makeIFPOptionsFromConfig(cfg);

        Residu::loadRules();
        Molecule::loadRules();
        configureResidueRules(options);

        // One protein, two ligand files
        std::vector<IFPEntry> dockedEntries = computeIFPsFromFiles(protein_file, ligand_file, options);
        std::vector<IFPEntry> refEntries = computeIFPsFromFiles(protein_file, reference_file, options);

        std::vector<TanimotoScore> scores;
        scores.reserve(dockedEntries.size() * refEntries.size());

        Similarity sims(false);

        for (const auto& d : dockedEntries) {
            for (const auto& r : refEntries) {

                // Safety: ensure fingerprints are compatible
                if (d.fpString.size() != r.fpString.size()) {
                    throw MoleExcept(9020105,"IChem::IFP","Cannot compare fingerprints with different sizes");
                }

                sims.setRef(r.fp);   // REF
                sims.setComp(d.fp);  // DOCKED

                TanimotoScore ts;
                // Docked side
                ts.ligand = d.name;
                ts.ligand_residues = d.fp.toString();
                ts.ligand_bitstring = d.fpString;

                // Reference side
                ts.reference = r.name;
                ts.reference_residues = r.fp.toString();
                ts.reference_bitstring = r.fpString;

                ts.tanimoto = sims.Tanimoto();

                scores.push_back(std::move(ts));
            }
        }

        return scores;
    }
    catch (MoleExcept& e) {
        e.addTrace("IFPAPI::compute_ifp_tanimoto");
        std::ostringstream oss;
        oss << "IChem::IFP error " << e.getCode() << " - " << e.getData();
        throw std::runtime_error(oss.str());
    }
}


std::vector<TanimotoScore> compute_ifp_tanimoto_ensembles(const std::string& protein1_file, const std::string& ligand1_file, const std::string& protein2_file, const std::string& ligand2_file, const IFPConfig& cfg) {
    
    using namespace IFPInternal;

    try {
        auto options = makeIFPOptionsFromConfig(cfg);

        Residu::loadRules();
        Molecule::loadRules();
        configureResidueRules(options);

        std::vector<IFPEntry> dockedEntries =
            computeIFPsFromFiles(protein1_file, ligand1_file, options);
        std::vector<IFPEntry> refEntries =
            computeIFPsFromFiles(protein2_file, ligand2_file, options);

        if (!dockedEntries.empty() && !refEntries.empty()) {
            if (dockedEntries.front().fpString.size() !=
                refEntries.front().fpString.size()) {
                throw MoleExcept(9020105,"IChem::IFP","Cannot compare fingerprints with different sizes in ensemble mode"
                );
            }
        }

        std::vector<TanimotoScore> scores;
        scores.reserve(dockedEntries.size() * refEntries.size());

        Similarity sims(false);

        for (const auto& d : dockedEntries) {
            for (const auto& r : refEntries) {

                if (d.fpString.size() != r.fpString.size()) {
                    throw MoleExcept(9020105, "IChem::IFP", "Cannot compare fingerprints with different sizes in ensemble mode");
                }

                sims.setRef(r.fp);   // ensemble2 as "reference"
                sims.setComp(d.fp);  // ensemble1 as "ligand"

                TanimotoScore ts;
                ts.ligand           = d.name;
                ts.ligand_residues  = d.fp.toString();
                ts.ligand_bitstring = d.fpString;

                ts.reference = r.name;
                ts.reference_residues = r.fp.toString();
                ts.reference_bitstring = r.fpString;

                ts.tanimoto = sims.Tanimoto();

                scores.push_back(std::move(ts));
            }
        }

        return scores;
    }
    catch (MoleExcept& e) {
        e.addTrace("IFPAPI::compute_ifp_tanimoto_ensembles");
        std::ostringstream oss;
        oss << "IChem::IFP error " << e.getCode() << " - " << e.getData();
        throw std::runtime_error(oss.str());
    }
}


} // namespace IFPAPI
