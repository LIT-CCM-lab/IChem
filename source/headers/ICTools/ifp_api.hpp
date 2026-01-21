#ifndef IFP_API_HPP
#define IFP_API_HPP

#include <string>
#include <vector>
#include <optional>

namespace IFPAPI {

// One interaction between protein and ligand
struct InteractionRecord {
    std::string type_interaction;  // Hydrophobic, Hbond
    std::string atom_prot;
    int id_atom_prot;
    std::string residue_identifier;
    std::string chain;
    std::string atom_lig;
    int id_atom_lig;
    double distance;
};

// One ligand's interactions (for multi-ligand files)
struct LigandInteractions {
    std::string ligand_name;
    std::vector<InteractionRecord> interactions;
};

// One ligand's fingerprint
struct LigandFingerprint {
    std::string ligand_name;
    std::string residues;
    std::string bitstring;
};

// Result of comparing two ligand sets
struct TanimotoScore {
    // Docked ligand
    std::string ligand;
    std::string ligand_residues;
    std::string ligand_bitstring;

    // Reference ligand
    std::string reference;
    std::string reference_residues;
    std::string reference_bitstring;

    // Similarity
    double tanimoto;
};

// Configuration => mirrors CLI options in a structured way
struct IFPConfig {
    // Which bits are active:
    bool all = false;
    bool basic = false;
    bool weakh = false;
    bool picat = false;
    bool metal = false;
    bool old_layout = false; // 7 bit layout

    // Rules
    bool includeSolvent  = true;
    bool includeCofactor = true;
    bool oldHydrophobic  = true; // --newH for false

    // Distances (optional)
    std::optional<double> D_Hb, D_Hyd, D_Io, D_Me, D_Ar, D_Pic, D_WHb;
    std::optional<double> d_Hb, d_Hyd, d_Io, d_Me, d_Ar, d_Pic, d_WHb;

    // Angles (optional)
    std::optional<double> a_H, at_H;
    std::optional<double> a_ArFF, at_ArFF;
    std::optional<double> a_ArEF, at_ArEF;
    std::optional<double> a_Pic, at_Pic;
};

// 1) IChem IFP protein.mol2 ligand.mol2 -> list of interactions
std::vector<LigandInteractions> compute_ifp_interactions(const std::string& protein_file, const std::string& ligand_file, const IFPConfig& cfg = IFPConfig());

// 2) IChem IFP protein.mol2 ligand.mol2 -> bitstrings
std::vector<LigandFingerprint> compute_ifp_fingerprints(const std::string& protein_file, const std::string& ligand_file, const IFPConfig& cfg = IFPConfig());

// 3) IChem IFP protein.mol2 ligand.mol2 reference.mol2 -> Tanimoto
std::vector<TanimotoScore> compute_ifp_tanimoto(const std::string& protein_file, const std::string& ligand_file, const std::string& reference_file, const IFPConfig& cfg = IFPConfig());

// 4) IChem IFP protein1 ligand1 protein2 ligand2 -> Tanimoto for two ensembles
std::vector<TanimotoScore> compute_ifp_tanimoto_ensembles(const std::string& protein1_file, const std::string& ligand1_file, const std::string& protein2_file, const std::string& ligand2_file, const IFPConfig& cfg = IFPConfig());

} // namespace IFPAPI

#endif // IFP_API_HPP
