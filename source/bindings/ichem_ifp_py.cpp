#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "headers/ICTools/ifp_api.hpp"

namespace py=pybind11;

using IFPAPI::IFPConfig;
using IFPAPI::InteractionRecord;
using IFPAPI::LigandInteractions;
using IFPAPI::LigandFingerprint;
using IFPAPI::TanimotoScore;

PYBIND11_MODULE(ichem_ifp,m)
{
    m.doc()=R"pbdoc(
###################################################################################################
#                                                                                                 #
#                    ~~~===>>>----====::::://///\\\\\\::::====----<<<===~~~                       #
#           F R A N Ç O I S   -   T H E   L A S T   T O R S I O N   B E N D E R   A P I           #
#                    ~~~===>>>----====::::://///\\\\\\::::====----<<<===~~~                       #
#                                                                                                 #
###################################################################################################


Python bindings for IChem IFP module:
- configure IFP profiles (--basic, --weakh etc) and thresholds (minimum maximum distances)
- list protein ligand interactions
- compute fingerprints and Tanimoto scores
)pbdoc";

    // IFPConfig
    py::class_<IFPConfig>(m,"IFPConfig","IFP configuration: profile and thresholds")
        .def(py::init<>())

        // profiles (pick exactly one)
        .def_readwrite("basic",&IFPConfig::basic, "Basic profile: hydrophobic + H-bond + ionic + aromatic")
        .def_readwrite("weakh",&IFPConfig::weakh, "Weak H-bonds only")
        .def_readwrite("picat",&IFPConfig::picat, "Pi-cation only")
        .def_readwrite("metal",&IFPConfig::metal, "Metal only")
        .def_readwrite("old",&IFPConfig::old_layout, "Basic profile for old version, coded on 7 bits")

        // rules
        .def_readwrite("includeSolvent",&IFPConfig::includeSolvent, "Remove water residues")
        .def_readwrite("includeCofactor",&IFPConfig::includeCofactor, "Remove cofactors")
        .def_readwrite("oldHydrophobic",&IFPConfig::oldHydrophobic, "Use legacy hydrophobic definition")

        // max distances (Å)
        .def_readwrite("D_Hb",&IFPConfig::D_Hb, "Max H-bond distance (Å)")
        .def_readwrite("D_Hyd",&IFPConfig::D_Hyd, "Max hydrophobic distance (Å)")
        .def_readwrite("D_Io",&IFPConfig::D_Io, "Max ionic distance (Å)")
        .def_readwrite("D_Me",&IFPConfig::D_Me, "Max metal distance (Å)")
        .def_readwrite("D_Ar",&IFPConfig::D_Ar, "Max aromatic distance (Å)")
        .def_readwrite("D_Pic",&IFPConfig::D_Pic, "Max Pi-cation distance (Å)")
        .def_readwrite("D_WHb",&IFPConfig::D_WHb, "Max weak H-bond distance (Å)")

        // min distances (Å)
        .def_readwrite("d_Hb",&IFPConfig::d_Hb,"Min H-bond distance (Å)")
        .def_readwrite("d_Hyd",&IFPConfig::d_Hyd,"Min hydrophobic distance (Å)")
        .def_readwrite("d_Io",&IFPConfig::d_Io,"Min ionic distance (Å)")
        .def_readwrite("d_Me",&IFPConfig::d_Me,"Min metal distance (Å)")
        .def_readwrite("d_Ar",&IFPConfig::d_Ar,"Min aromatic distance (Å)")
        .def_readwrite("d_Pic",&IFPConfig::d_Pic,"Min Pi-cation distance (Å)")
        .def_readwrite("d_WHb",&IFPConfig::d_WHb,"Min weak H-bond distance (Å)")

        // angles (degrees)
        .def_readwrite("a_H",&IFPConfig::a_H,"H-bond angle cutoff (°)")
        .def_readwrite("at_H",&IFPConfig::at_H,"H-bond angle tolerance (°)")
        .def_readwrite("a_ArFF",&IFPConfig::a_ArFF,"Aromatic face-face angle cutoff (°)")
        .def_readwrite("at_ArFF",&IFPConfig::at_ArFF,"Aromatic face-face angle tolerance (°)")
        .def_readwrite("a_ArEF",&IFPConfig::a_ArEF,"Aromatic edge-face angle cutoff (°)")
        .def_readwrite("at_ArEF",&IFPConfig::at_ArEF,"Aromatic edge-face angle tolerance (°)")
        .def_readwrite("a_Pic",&IFPConfig::a_Pic,"Pi-cation angle cutoff (°)")
        .def_readwrite("at_Pic",&IFPConfig::at_Pic,"Pi-cation angle tolerance (°)")
        ;

    // InteractionRecord
    py::class_<InteractionRecord>(m,"InteractionRecord","Protein-ligand interaction")
        .def_readonly("type_interaction",&InteractionRecord::type_interaction, "Interaction type")
        .def_readonly("atom_prot",&InteractionRecord::atom_prot, "Protein atom name")
        .def_readonly("id_atom_prot",&InteractionRecord::id_atom_prot, "Protein atom index")
        .def_readonly("residue",&InteractionRecord::residue, "Residue identifier")
        .def_readonly("chain",&InteractionRecord::chain, "Chain id")
        .def_readonly("atom_lig",&InteractionRecord::atom_lig, "Ligand atom name")
        .def_readonly("id_atom_lig",&InteractionRecord::id_atom_lig, "Ligand atom index")
        .def_readonly("distance",&InteractionRecord::distance,"Atom-atom distance (Å)");

    // LigandInteractions
    py::class_<LigandInteractions>(m,"LigandInteractions","All interactions for one ligand")
        .def_readonly("ligand_name",&LigandInteractions::ligand_name, "Ligand name")
        .def_readonly("interactions",&LigandInteractions::interactions, "List of InteractionRecord");

    // LigandFingerprint
    py::class_<LigandFingerprint>(m,"LigandFingerprint","Interaction fingerprint for one ligand")
        .def_readonly("ligand_name",&LigandFingerprint::ligand_name, "Ligand name")
        .def_readonly("residues",&LigandFingerprint::residues, "Residues used in the fingerprint")
        .def_readonly("bitstring",&LigandFingerprint::bitstring, "Binary fingerprint string");

    // TanimotoScore
    py::class_<TanimotoScore>(m,"TanimotoScore","Tanimoto similarity between two ligands")
        .def_readonly("ligand",&TanimotoScore::ligand, "Docked ligand name")
        .def_readonly("ligand_residues",&TanimotoScore::ligand_residues, "Residues for the docked fingerprint")
        .def_readonly("ligand_bitstring",&TanimotoScore::ligand_bitstring, "Fingerprint for the docked ligand")
        .def_readonly("reference",&TanimotoScore::reference, "Reference ligand name")
        .def_readonly("reference_residues",&TanimotoScore::reference_residues, "Residues for the reference fingerprint")
        .def_readonly("reference_bitstring",&TanimotoScore::reference_bitstring, "Fingerprint for the reference ligand")
        .def_readonly("tanimoto",&TanimotoScore::tanimoto, "Tanimoto coefficient");

    // Functions
    m.def("compute_ifp_interactions",
          &IFPAPI::compute_ifp_interactions,
          py::arg("protein_file"),
          py::arg("ligand_file"),
          py::arg("cfg")=IFPConfig{},
          "Compute all protein-ligand interactions for a protein and one/many ligands");

    m.def("compute_ifp_fingerprints",
          &IFPAPI::compute_ifp_fingerprints,
          py::arg("protein_file"),
          py::arg("ligand_file"),
          py::arg("cfg")=IFPConfig{},
          "Compute interaction fingerprints for a protein and one/many ligands");

    m.def("compute_ifp_tanimoto",
          &IFPAPI::compute_ifp_tanimoto,
          py::arg("protein_file"),
          py::arg("ligand_file"),
          py::arg("reference_file"),
          py::arg("cfg")=IFPConfig{},
          "Compute Tanimoto similarities: ligands vs reference ligands (same protein)");

    m.def("compute_ifp_tanimoto_ensembles",
          &IFPAPI::compute_ifp_tanimoto_ensembles,
          py::arg("protein1_file"),
          py::arg("ligand1_file"),
          py::arg("protein2_file"),
          py::arg("ligand2_file"),
          py::arg("cfg")=IFPConfig{},
          "Compute Tanimoto similarities between two protein/ligand ensembles");
}