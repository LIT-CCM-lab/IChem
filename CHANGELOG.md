# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Refactoring] 


### Added

- Interactions::detectInteractions() replaces the old monolithic function calcInteractions() 
- Interactions::addInteractions(): helper function to centralize creation of InterPoint entries and to update InterResults counters (N, O, Zn, CA, CZ, NZ, OD1)
- Interactions::processAromaticInteractions(): refactored aromatic face/edge detection logic into its own function
- Interactions::processHydrophobicInteraction(): refactored hydrophobic contact handling into its own function
- Interactions::checkAromaticHydrophobicInteractions(): separate function that handles hydrophobic interactions between aromatic rings
- A set of specialized helpers functions called from detectInteractions: 
    checkMetalNitrogenSulfonamideCase(): detects metal interactions with a sulfonamide nitrogen bound to sulfur
    checkMetalInteractions(): identifies metal ligand interactions based on the distance
    checkHydrogenBondLigandAcceptor(): finds hydrogen bonds where the ligand acts as acceptor and protein as donor
    checkWeakHydrogenBondLigandAcceptor(): detects weak hydrogen bonds with ligand as acceptor and protein as weak donor
    checkWeakHydrogenBondLigandWeakAcceptor(): identifies weak hydrogen bonds when the ligand is a weak acceptor and protein is weak donor or just donor
    checkHydrogenBondLigandDonor(): finds hydrogen bonds where the ligand acts as donor and protein as acceptor
    checkWeakHydrogenBondLigandDonorProteinWeakAcceptor(): detects weak hydrogen bonds where the ligand is donor and protein is weak acceptor
    checkWeakHydrogenBondLigandWeakDonorProteinAcceptor(): identifies weak hydrogen bonds with ligand as weak donor and protein as weak acceptor or just acceptor
    checkIonicProteinInteractions(): detects ionic interactions between anionic ligand atoms and cationic protein atoms
    checkIonicLigandInteractions(): detects ionic interactions between cationic ligand atoms and anionic protein atoms
    checkPiCationInteraction(): identifies pi-cation interactions between ligand atoms and aromatic protein cycles
- These helper functions break the previous long function (calcInteractions()) into focused and testable units

- Spatial search / algorithmic data structures:
    - NeighborSearch KD-tree style index used to build a point cloud of protein atoms and perform fast neighborhood queries
    - The KD-tree can handle multiple ligand atoms in a single call returning a list of nearby protein atoms for each ligand atom efficiently and returning contact structs
    - Use of NeighborSearch::Contact (with ligand_idx, protein_idx, distance_squared) to efficiently iterate contacts
- Minor data-structure additions:
    - resbest usage extended now used with both ligand and protein atom pointers in hydlist entries to store best hydrophobic candidate per residue
        - In the original struct with only the protein atom, we only track which protein atom was closest and its distance. But when we need to add the interaction we also need to know which ligand atom was involved. Without atmL we lose the information. it’s essential for being able to reconstruct the ligand-protein pair when flushing hydrophobic interactions
    - addInteraction increments NInter and appends to interResult.listInters using emplace_back for efficiency

### Changed

- Major refactor: monolithic calcInteractions() reworked into a modular pipeline:
    - heavy nested loops over boxes/adjacent boxes replaced by building a KD-tree of protein points and performing nearest-neighbour queries for all ligand atoms (big performance and clarity improvement)
    - grid adjacency based scanning (grid.getAdjacency() per ligand atom) removed in favor of NeighborSearch queries
    - Many responsibilities split into focused functions (see Added), this improves readability testability and future maintenance

- Performance and  complexity:

    - Spatial complexity reduced by using KD-tree queries and performing a single query for all ligand atoms at once instead of iterating many boxes/atoms per ligand atom
    - Sorting by ligand atom groups contacts, allowing hydrophobic contacts to be finalized per atom before processing the next one
    - addInteraction centralizes InterPoint creation and result counters reducing duplication and chances for inconsistent updates

- Behavioural / algorithmic changes:

    - Hydrophobic detection:
        - Old: used grid-adjacency and grid.getAdjacentAtoms then computed local neighbor density by iterating atoms in boxes
        - New: for a candidate protein atom atomP the hydrophobic density check uses neighborSearch.query around atomP and computes nbhyd/nbatm directly from the KD results
        - Hydrophobic candidates are accumulated per protein residue in hydlist and flushed per ligand atom, keeping the "best" protein atom per residue (smallest distance)
    - Aromatic interactions:
        - Old: iterated ligand cycles, got adjacency using grid.getAdjacency of the ligand cycle center; scanned protein atoms to find DuAr centers and then scanned full cycles comparing per atom
        - New: for each ligand aromatic cycle the KD-tree is queried around the ring center (neighborSearch.query({q}, searchRadius)), candidate protein ring centers are collected and made unique then perform per ring pair checks
        - Uses Cycle::calcVector() and getNormVector similarly but the candidate selection is KD-tree based
    - Pi-cation detection:
        -Old: pi-cation checks were done inside the big box loop and inside the aromatic pair loops. the new code isolates these checks into checkPiCationInteraction() and also uses KD-tree hits to find protein cation centers near ligand cycle centers
    - Metal / special cases:
        - Metal interactions are moved to checkMetalInteractions()
        - The special N-linked to sulfonamide is moved to checkMetalNitrogenSulfonamideCase()
    - Partial charges / PLP:
        - All the PLP logic is completely removed from detectInteractions()
    - Iteration and ordering:
        - Old: used grid boxes and nested iteration over boxlist and their firstAtom() loops, rotating molecules into the grid as necessary. 
        - New: builds proteinPoints and ligandPoints vectors and processes pairs returned by KD-tree. This changes processing order, hydrophobic flush per ligand atom uses sorted contact list keyed by ligand index

- API changes:
    - Signature changed: calcInteractions(Molecule& ligand, InterResults& interResult, bool wMerge, bool oldh, bool mono_prop, bool out_lig, bool stdoutt) to detectInteractions(Molecule& ligand, InterResults& interResult, bool wMerge, bool oldh)
    - Removed parameters: mono_prop, out_lig, stdoutt
    - Internal global/static thresholds (Interactions::Dist_H, Dist_Hyd, dist_H etc) have been migrated to a parameter object. All these thresholds are now isolated and organised inside a struct.
    - Interactions::setStdRules() is removed as it is not necessary anymore

- Code quality modernizations:
    - Use of C++-style containers (std::vector, std::array) and modern APIs:
        - emplace_back used to construct InterPoint
        - std::sort and std::unique used for deduplication of cycle lists (replacing SortUniqueVector)
        - std::max with initializer list used to compute max_allowed_dist
        - Batch processing (vectorized queries) and lambda helpers (eg push or flushHydList)

- Robustness and correctness:

    - Multiple small defensive checks added:
        - skipping same residue hydrophobic entries
            - Same residue hydrophobic skip is applied earlier and centrally: instead of checking only at the final add the new code puts the guard into the flushHydList() helper that’s called for every ligand atoms batch. That guarantees that hydrophobic interactions within the same residue should not be added to the output list
        - Deduplication of protein cycles before evaluation to avoid redundant checks
            - Now the code removes duplicate cycles after collecting them from the KD-tree search by sorting and removing duplicates, we make sure each ring is only counted once

### Removed

- Old monolithic implementation calcInteractions() with nested grid.getAdjacency() loops and box scanning removed and replaced by detectInteractions() with neighborsearch algorithm and helpers
- The out_lig/stdoutt logging logic has been removed
- The old setStdRules() initialization routine and the static variables such as Interactions::Dist_H, Interactions::Dist_Hyd, Interactions::Angl_H has been replaced by a params style configuration.
- All PLP logic in the ligand loop are removed
- grid.getAdjacentAtoms() for hydrophobic density calculations are removed in favour of NeighborSearch

## [IChem_5.3.5] - 2025-09-05

### Added

- Feature IFP calculating the similarity of two IFPs (4 arguments)
- We display all the apolar interactions of type C.ar atoms between cycles
- Query optimized in neighborsearch, replace std::pair with a struct
- Nanoflann well integrated with all the interactions  

### Changed

- Remove the best hydrophobic interaction in aromatic interaction: we display them all now.
- Function name calcInteractions() into detectInteractions()
- nanoflann well integrated to processAromaticInteractions()
- Refactoring processAromaticInteractions()

### Removed

- the if block flags for pication/aromatic interactions in calcInteractionsppi() is unecessary
     + (We need to add the change for everything we added as functions and how they are named)
- CalcInteractionPpi() removed from interaction.cpp
- The function mergeSpeInts() removed as it is used within calcInteractionsppi() function
 
## [IChem_5.3.3 - IChem_5.3.4] - 2025-06-02

### Added

- Repository ci-containers for different OS (ubuntu 22.04 / ubi8.8 / debian 12)
- gitlab-ci.yml, pipeline added 
- add nanoflann header into the project
- refactoring calcInteractions: Isolating ionic interactions into its own function 
- Display ligand aromatic cycle atoms after detecting an aromatic interaction
- weakAcceptor and weakDonor properties added in physprop.cpp
- functional tests for grim module 
- functional tests for ints module
- functional tests for IFP
- AtomProps functional tests
- Add data for AtomProps functional tests
- AtomProps feature (atom_properties.cpp) added, included in switch.h and switch.cpp
- Functional tests for sims module
- Adding checking block for the files (PROD/TEST) 
- gitignore added

### Fixed

-  Repetitive aromatic interaction removed 
-  Use IChem indexation to solve the homodimer bug, map data structure do not allow repetitive index 
-  in mergeInteractions() bug invalid memory access pointers
- uninitialized value in physprop.cpp  

### Changed

- Separate metal nitrogen sulfonamide case into its own function
- IChem menu: AtomProps menu description changed
- Refactoring processHydrophobicInteraction() function
- Dist_WHBond value (max dist for weak hydrogen) to 3.5 from 2.8
- Refactoring: remove static variables from interaction.h
- Refactoring: encapsulate the distances into a struct
- Refactoring: move resbest struct into interaction.h
- Change parameters for calcInteractions() - update in all files using it (IFP-Grim-Ints)
- minimum distance of ionic and hydrogen bond is set to 2.3 from 2.5
- props.setWeakAcceptor(true) changed to props.setWeakDonor(true) in atom.cpp for aromatic
- Sending the warnings to stderr: to prevent the pollution of warning messages while parsing mol2 files. 
- Change menu ints for aromatic and Pi-cation interaction length set to 5.0
- Dist_Arom and Dist_PiCation set to 5.0 from 4.0


### Removed

- ppi related files deleted (detectPPI.cpp, ppi.h), removed from CMakeLists.txt and from modifications in switch.h and switch.cpp 