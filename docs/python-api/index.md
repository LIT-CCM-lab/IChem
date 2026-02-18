# Python API - Overview

The `ichem_ifp` module is the Python interface to the IChem Interaction Fingerprint engine

It allows users to compute protein-ligand interactions, generate interaction fingerprints, and evaluate similarity scores directly from Python, while relying on the high-performance native IChem core


## Architecture Overview
### IChem IFP Engine

The Interaction Fingerprint functionality is part of the IChem software suite and is implemented in native C++.
The core engine performs:

- Detection of protein-ligand interactions
- Construction of interaction fingerprints
- Computation of Tanimoto similarity scores

This layer contains the optimized computational logic

### Python Bindings Layer

The ichem_ifp module is built using pybind11, which exposes the C++ IChem IFP engine to Python.
This design provides:

- Native performance
- Direct access to core IFP functionality
- A clean and Pythonic API interface
- All computationally intensive operations are executed within the IChem engine, while Python acts as the orchestration and user interface layer

#### Exposed Functionality

The Python API exposes the following high-level capabilities:

- configure IFP profiles (**basic / weakh / metal / picat / old**) and **geometric thresholds**
- list **protein-ligand interactions**
- compute **interaction fingerprints** (bitstrings)
- compute **Tanimoto similarities** (ligands vs references, or ensemble (protein and ligand) vs ensemble)

> **Input formats**
> Functions expect **MOL2** files (site/protein and ligand(s)).
> Ligand input may consist of:
    - A single ligand MOL2 file
    - A multi-MOL2 file containing multiple ligands

---

## Quick install

See: **Python API section Install**

---

## Typical workflow

1. Create an `IFPConfig()` and enable a profile
2. Call a function (interactions / fingerprints / tanimoto / tanimoto_ensembles)
3. Consume returned objects: `LigandInteractions`, `LigandFingerprint`, `TanimotoScore`

---

## Main functions

| Goal | Function | Returns |
|---|---|---|
| List all interactions | `compute_ifp_interactions(protein_file, ligand_file, cfg)` | `list[LigandInteractions]` |
| Compute fingerprints | `compute_ifp_fingerprints(protein_file, ligand_file, cfg)` | `list[LigandFingerprint]` |
| Tanimoto (ligands vs references, same protein) | `compute_ifp_tanimoto(protein_file, ligand_file, reference_file, cfg)` | `list[TanimotoScore]` |
| Tanimoto between two ensembles | `compute_ifp_tanimoto_ensembles(protein1, ligand1, protein2, ligand2, cfg)` | `list[TanimotoScore]` |

---

## Data model - returned objects

- **`LigandInteractions`**
    - `ligand_name`: ligand name
    - `interactions`: list of `InteractionRecord`

- **`InteractionRecord`**
    - `type_interaction`: interaction type
    - `atom_prot`: protein atom name
    - `id_atom_prot`: protein atom index
    - `residue_identifier`: protein residue identifier
    - `chain`: chain id
    - `atom_lig`: ligand atom name
    - `id_atom_lig`: ligand atom index
    - `distance`: distance (Å)
    - `angle` : angle (°) *(mainly meaningful for HBond / Aromatic / Pi-cation. It returns `None` for other interaction types)*

- **`LigandFingerprint`**
    - `ligand_name` : ligand name
    - `residues`: residue list used to build the fingerprint layout
    - `bitstring`: binary fingerprint string

- **`TanimotoScore`**
    - `ligand` : ligand
    - `ligand_residues` : ligand residues
    - `ligand_bitstring` : ligand binary fingerprint string
    - `reference` : reference
    - `reference_residues` : reference residues
    - `reference_bitstring` : reference binary fingerprint string
    - `tanimoto` : similarity

---

## Important notes

- "Fingerprint compatibility": Fingerprints are **site-dependent**, the residue order (`residues`) defines the fingerprint layout
    Compare fingerprints only when generated from compatible sites


---

## Next
- **Installation**: How to install the Python API
- **Quickstart**: minimal scripts (interactions / fingerprints / tanimoto)
- **Configuration**: `IFPConfig` explained
