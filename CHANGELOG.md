# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Not released yet
  * Keep the name of the molecule extended (limited to 40 chacraters before)
  * Adding helper function to remove whitespace in molecule's name

## [IChem_5.3.8] - 2025-12-08

### Added

- High-level C++ API for IFP (`ICTools/ifp_api.cpp`, `headers/ICTools/ifp_api.hpp`):
  - New namespace `IFPAPI` providing a structured, programmatic API to the IFP engine, decoupled from the CLI
  - New configuration and data structures:
    - `IFPConfig`  
      - Mirrors the IFP CLI options in a single struct:
        - profile flags: `basic`, `weakh`, `picat`, `metal`, `old_layout`  
        - rules: `includeSolvent`, `includeCofactor`, `oldHydrophobic`  
        - optional distance overrides: `D_Hb`, `D_Hyd`, `D_Io`, `D_Me`, `D_Ar`, `D_Pic`, `D_WHb`  
        - optional minimum distances: `d_Hb`, `d_Hyd`, `d_Io`, `d_Me`, `d_Ar`, `d_Pic`, `d_WHb`  
        - optional angles and tolerances: `a_H`, `at_H`, `a_ArFF`, `at_ArFF`, `a_ArEF`, `at_ArEF`, `a_Pic`, `at_Pic`
    - `InteractionRecord`  
      - One interaction between protein and ligand:
        - `type_interaction`  
        - `atom_prot`, `id_atom_prot`  
        - `residue`, `chain`  
        - `atom_lig`, `id_atom_lig`  
        - `distance`
    - `LigandInteractions`  
      - Per-ligand container:
        - `ligand_name`  
        - `interactions`: `std::vector<InteractionRecord>`
    - `LigandFingerprint`  
      - Per-ligand fingerprint representation:
        - `ligand_name`  
        - `residues`: string with the list of interacting residues  
        - `bitstring`: IFP bitstring corresponding to the active profile
    - `TanimotoScore`  
      - Pairwise similarity container:
        - `ligand`, `ligand_residues`, `ligand_bitstring`  
        - `reference`, `reference_residues`, `reference_bitstring`  
        - `tanimoto` (double)

- New C++ entry points built on top of the IFP engine:
  - `std::vector<LigandInteractions> compute_ifp_interactions(const std::string& protein_file, const std::string& ligand_file, const IFPConfig& cfg = IFPConfig());`  
    - Compute all protein–ligand interactions for a protein and a (possibly multi-ligand) file  
    - Returns one `LigandInteractions` for each ligand in `ligand_file`.
  - `std::vector<LigandFingerprint> compute_ifp_fingerprints(const std::string& protein_file, const std::string& ligand_file, const IFPConfig& cfg = IFPConfig());`  
    - Compute IFP bitstrings for each ligand in `ligand_file` against the same protein
    - Each entry exposes the residues string and the bitstring, mirroring the IChem CLI output
  - `std::vector<TanimotoScore> compute_ifp_tanimoto(const std::string& protein_file, const std::string& ligand_file, const std::string& reference_file, const IFPConfig& cfg = IFPConfig());`  
    - One-protein, two-ligand-files mode:
      - `ligand_file`   : docked ligands  
      - `reference_file`: reference ligands  
    - Returns all pairwise similarities with full residue/bitstring context
  - `std::vector<TanimotoScore> compute_ifp_tanimoto_ensembles(const std::string& protein1_file, const std::string& ligand1_file, const std::string& protein2_file, const std::string& ligand2_file, const IFPConfig& cfg = IFPConfig());`  
    - Two-protein, two-ligand-files “ensemble” mode, matching the 4 argument CLI  
    - Returns all docked vs reference scores with associated residues and bitstrings

- Python bindings for the IFP API (`source/bindings/ichem_ifp_py.cpp`):
  - New Python extension module `ichem_ifp` built with pybind11, wrapping the C++ `IFPAPI` layer
  - Exposed Python classes:
    - `IFPConfig`  
      - Same fields as the C++ struct, accessible as attributes:
        - `basic`, `weakh`, `picat`, `metal`, `old`  
        - `includeSolvent`, `includeCofactor`, `oldHydrophobic`  
        - `D_Hb`, `D_Hyd`, `D_Io`, `D_Me`, `D_Ar`, `D_Pic`, `D_WHb`  
        - `d_Hb`, `d_Hyd`, `d_Io`, `d_Me`, `d_Ar`, `d_Pic`, `d_WHb`  
        - `a_H`, `at_H`, `a_ArFF`, `at_ArFF`, `a_ArEF`, `at_ArEF`, `a_Pic`, `at_Pic`
    - `InteractionRecord`  
      - Read-only view of a single interaction:
        - `type_interaction`, `atom_prot`, `id_atom_prot`,  
          `residue`, `chain`,  
          `atom_lig`, `id_atom_lig`,  
          `distance`
    - `LigandInteractions`  
      - Read-only container:
        - `ligand_name`  
        - `interactions` (list of `InteractionRecord`)
    - `LigandFingerprint`  
      - Read-only container:
        - `ligand_name`  
        - `residues`  
        - `bitstring`
    - `TanimotoScore`  
      - Read-only container:
        - `ligand`, `ligand_residues`, `ligand_bitstring`  
        - `reference`, `reference_residues`, `reference_bitstring`  
        - `tanimoto`
  - Exposed Python functions:
    - `ichem_ifp.compute_ifp_interactions(protein_file, ligand_file, cfg=IFPConfig())`  
      - Returns a list of `LigandInteractions`
    - `ichem_ifp.compute_ifp_fingerprints(protein_file, ligand_file, cfg=IFPConfig())`  
      - Returns a list of `LigandFingerprint`
    - `ichem_ifp.compute_ifp_tanimoto(protein_file, ligand_file, reference_file, cfg=IFPConfig())`  
      - Returns a list of `TanimotoScore` for a single protein and two ligand sets
    - `ichem_ifp.compute_ifp_tanimoto_ensembles(protein1_file, ligand1_file, protein2_file, ligand2_file, cfg=IFPConfig())`  
      - Returns a list of `TanimotoScore` for two protein/ligand ensembles

### Changed

- IFP engine and module layering (`ICTools/ifp_module.cpp`, `ICTools/ifp_engine.cpp`, `headers/ICTools/ifp_engine.hpp`, `headers/ICTools/ifp_module.hpp`):
  - The IFP code is now split into three logical layers:
    - IFP front-end (`IFPModule::runIFP` in `ifp_module.cpp`):
      - Responsible only for:
        - CLI argument dispatch (2/3/4 argument modes)
        - printing interaction/fingerprint/similarity results to stdout
      - Reuses the shared engine functions instead of embedding the full logic
    - IFP engine (`IFPInternal` namespace in `ifp_engine.cpp` / `ifp_engine.hpp`):
      - Hosts the core reusable building blocks:
        - `IFPOptions` and `InteractionOverrides`  
        - `parseIFPOptions()`  
        - `parseIFPNumeric()`  
        - `applyOverrides()`  
        - `configureResidueRules()`  
        - `computeIFPForLigand()`  
        - `computeIFPsFromFiles()`
      - This layer is used by both the CLI (`IFPModule::runIFP`) and the new C++/Python APIs
    - IFP high-level API (`IFPAPI` namespace in `ifp_api.cpp` / `ifp_api.hpp`):
      - Bridges the engine with external code (C++ or Python):
        - converts `IFPConfig` into `IFPOptions` via `makeIFPOptionsFromConfig()`  
        - builds `InteractionRecord`, `LigandInteractions`, `LigandFingerprint` and `TanimotoScore` structures  
        - keeps the CLI behaviour (profiles, thresholds, residue rules) consistent in all entry points

- Error reporting for API consumers:
  - All `IFPAPI` entry points (`compute_ifp_interactions`, `compute_ifp_fingerprints`, `compute_ifp_tanimoto`, `compute_ifp_tanimoto_ensembles`) now:
    - catch `MoleExcept`,  
    - append an API-specific trace (e.g: `"IFPAPI::compute_ifp_interactions"`),  
    - and rethrow as `std::runtime_error` with a message of the form: `"IChem::IFP error <code> - <details>"`
  - This preserves the IChem error codes and messages while making the failures visible as standard exceptions to C++ callers and Python users


## [IChem_5.3.7] - 2025-12-03

### Added

- IFP numeric threshold overrides (in `ICTools/ifp_module.cpp`):
  - Distances (maximum):
    - `-D_WHb`  : max weak H-bond distance
  - Distances (minimum):
    - `-d_Hb`, `-d_Hyd`, `-d_Io`, `-d_Me`, `-d_Ar`, `-d_Pic`, `-d_WHb`

- IFP “profiles” and 11 bit fingerprint layout:
  - New mandatory, mutually exclusive profile options:
    - `--all`   : activate all 11 bits (basic + pi-cation + metal + weak H-bonds)
    - `--basic` : only the basic 7 bits (hydrophobic, aromatic FF/EF, H-bond P/L, ionic P/L)
    - `--picat` : only the pi-cation bit
    - `--metal` : only the metal/acceptor bit
    - `--weakh` : only the weak H-bond bits
    - `--old`   : basic interactions, but using the legacy 7 bit layout
  - Internally, `parseIFPOptions()` computes an 11-bit mask:
    - bit 0 : hydrophobic  
    - bit 1 : aromatic face-to-face  
    - bit 2 : aromatic edge-to-face  
    - bit 3 : H-bond protein  
    - bit 4 : H-bond ligand  
    - bit 5 : ionic protein  
    - bit 6 : ionic ligand  
    - bit 7 : Pi-cation  
    - bit 8 : metal  
    - bit 9 : weak H-bond protein  
    - bit 10: weak H-bond ligand  
  - A dedicated flag bit (`FLAG_OLD_LAYOUT = 1u << 31`) is used to request the legacy 7 bit (basic interactions)

- Strict numeric parsing helper for IFP thresholds:
  - New function `parseIFPNumeric(const std::string& value, const std::string& optName)`:
    - Enforces at most **3 digits after the decimal point**
    - Requires at least one digit before the decimal point
    - Rejects inputs such as `"3."`, `".5"`, `"3.12345"`, `"abc"`, `"+"`, `"-"`, or garbage values
  - All distance and angle options in `parseIFPOptions()` now go through `parseIFPNumeric()`, and invalid values throw an exception

### Changed

- IFP help and documentation (`ICTools/IFP.cpp`):
  - `IChemSwitch::helpIFP()` has been fully rewritten to match the new IFP interface:
    - Documents all three modes:
      - `IChem [options] IFP protein ligand`
      - `IChem [options] IFP protein ligand ligand_ref`
      - `IChem [options] IFP protein1 ligand1 protein2 ligand2`
    - States that **all options must appear before the keyword `IFP`**
    - Documents that **distance thresholds accept at most 3 decimals**, matching `parseIFPNumeric()`
    - Describes the 11 bit vs legacy 7 bit formats and their relation to `--old`
    - Lists the new profile options (`--all`, `--basic`, `--picat`, `--metal`, `--weakh`, `--old`) and their effects
    - Lists default values for all distance and angle overrides, grouped by:
      - maximum distances (`-D_*`),
      - minimum distances (`-d_*`),
      - angle / tolerance options (`-a_*`, `-at_*`).
    - Adds concrete usage examples for 2-, 3- and 4argument modes with thresholds and profiles combined

- IFP options parsing (`ICTools/ifp_module.cpp`):
  - `parseIFPOptions()` now:
    - Requires **exactly one** of `--all`, `--basic`, `--weakh`, `--picat`, `--metal`, or `--old`
      - If none is given, it throws `MoleExcept(9020101, "IChem::IFP", "You must specify one of ...")` 
      - If more than one profile is specified, it throws a “mutually exclusive” error.
    - Populates `IFPOptions::bitMask` according to the chosen profile (basic/all/weakH/Pi-cation/metal/old) using the 11 bit layout and `FLAG_OLD_LAYOUT` for legacy output
    - Validates every min/max distance pair via a `checkInterval()` helper:
      - If both `-d_*` and `-D_*` are provided and `min > max`, an error is raised:
        - `"For H-bond (-d_Hb / -D_Hb) minimal distance (...) is greater than maximal distance (...)"`
      - This is applied for H-bond, hydrophobic, ionic, metal, aromatic, Pi-cation, and weak H-bond distances
  - Unknown IFP options are no longer silently ignored:
    - Any option key not matched by the explicit `if/else if` chain now triggers:
      - `MoleExcept(9020101, "IChem::IFP", "Unknown IFP option: ... Allowed options include: --all, --basic, --weakh, --picat, --metal, --old, --solvent, --cofactor, --newH, -name, -D_*, -d_*, -a_*, -at_*")`

- IChem front-end CLI parsing (`ICTools/switch.cpp`):
  - Both `IChemSwitch` constructors (`(int argc, char** argv)` and `(int argc, const std::vector<std::string>& argv)`) now treat **unknown tokens before the tool name** as hard errors:
    - After checking:
      - long options (`--...`),
      - short options (`-...` with value),
      - and known tool names (`IFP`, `grim`, `volsite`, etc.),
    - any remaining token in the “options zone” (before the tool name is seen) triggers:
      ```cpp
      std::ostringstream oss;
      oss << "Unknown option or misplaced token before tool name: '" << tmpStr << "'";
      throw MoleExcept(9010106, "IChem::CONSTRUCTOR", oss.str());
      ```
  - Once a valid tool name is encountered, `opts` is set to `false` and all following tokens are treated strictly as input parameters (`Input_Values`) for that tool

### Fixed

- CLI robustness and error reporting:
  - Typographical errors where a user forgets the leading dash on an option are now correctly reported instead of being silently misinterpreted. For example:
    - `./IChem --weakh d_WHb 2.4 IFP protein.mol2 ligand.mol2`
    - `d_WHb` (without `-`) is no longer consumed as a value for the previous option, it is rejected with `MoleExcept(9010106, "IChem::CONSTRUCTOR", "Unknown option or misplaced token before tool name: 'd_WHb'")`
  - Unknown IFP options or misspelled option names (e.g. `-D_WHbb`, `-D_XXX`) now produce a clear `MoleExcept(9020101, "IChem::IFP", "Unknown IFP option: ...")` instead of silently doing nothing
  - Inconsistent min/max distance overrides (e.g: `-d_Hb 4.0 -D_Hb 3.0`) are now detected and rejected early, preventing nonsensical threshold configurations from reaching the interaction engine

---

## [IChem_5.3.6] - 2025-11-27

### Added

- Refactoring IFP module:
  - IFP module refactor in `ICTools/ifp_module.cpp` and `headers/ICTools/ifp_module.hpp`:
    - New entry point `IFPModule::runIFP(const std::vector<std::string>& inputValues, const OptionMap& optValues)` called from `IChemSwitch::IFP()`
    - New configuration structures:
      - `IFPOptions` global IFP configuration
      - `InteractionOverrides` (distance / angle overrides for interactions)
    - New per ligand container `IFPEntry` storing ligand name, `Fingerprint` and bitstring (`fpString`)
    - Internal helper functions that structure the IFP pipeline:
      - `parseIFPOptions()` : parses CLI options (`Opt_Values`) into an `IFPOptions` instance
      - `configureResidueRules()` : applies solvent / cofactor rules based on `IFPOptions`
      - `computeIFPForLigand()` : runs `Interactions::detectInteractions()` and `genIFP()` on a single ligand
      - `processLigandMode()` : handles the 2 argument mode (protein + ligand file, 1 or many ligands)
      - `processReferenceMode()` : handles the 3 argument mode (protein, docked ligands, reference ligands)
      - `processFourArgMode()` : handles the 4 argument mode (two protein ligand ensembles)
      - `computeIFPsFromFiles()` : utility to compute all IFPs for a pair (protein file, ligand file)

- Extended IFP CLI: 4 argument mode
  - New usage:
    - `IChem IFP protein1.mol2 ligand1.mol2 protein2.mol2 ligand2.mol2`
      - `(protein1, ligand1)` is treated as the **docked** set
      - `(protein2, ligand2)` is treated as the **reference** set
  - Both ligand files may contain **multiple molecules**; the code computes IFPs for all ligands on both sides and prints:
    - All DOCKED bitstrings (protein1/ligand1)
    - Then all REF bitstrings (protein2/ligand2)
    - Then all pairwise Tanimoto similarities `DOCKED  REF  score`
  - Fingerprint size mismatches are detected and reported via an error instead of silently comparing incompatible fingerprints

- Unit test infrastructure:
  - New top level directory `unit_tests/` with its own `CMakeLists.txt`
  - New test target `ifp_unit_tests` in `unit_tests/ICTools/ifp_test.cpp`:
    - Exposes `chooseIFPType_forTests(bool polarOnly, bool extended, bool metalOnly)`
    - Verifies the full truth table of combinations for `--polar`, `--extended`, `--metal`
  - Tests are integrated with CTest:
    - `add_test(NAME ICTools_ifp_unit COMMAND ifp_unit_tests)` in `unit_tests/CMakeLists.txt`

### Changed

- Language and toolchain:
  - Project migrated from **C++11** to **C++20**:
    - `set(CMAKE_CXX_STANDARD 20)` in the top-level `CMakeLists.txt`
    - Updated compiler / standard requirements in the build configuration
  - Added `STATIC_ANALYSIS` CMake option to toggle GCC `-fanalyzer` and `clang-tidy`

- Build structure:
  - Core sources are now built into a reusable static library target `ichem_core`
  - The main `IChem` executable links against `ichem_core` instead of compiling every source directly
  - Unit test binaries under `unit_tests/` also link against `ichem_core`, so tests use exactly the same core code as the production binary

- IFP command line and behaviour:
  - `IChemSwitch::helpIFP()` updated:
    - Documents the 4 argument mode: `IChem [options] IFP protein1 ligand1 protein2 ligand2`
    - Document `--bitstringOFF` option added
  - The choice of IFP layout (standard / polar / extended / extended polar / metal) is centralized in:
    - `chooseIFPType(bool polarOnly, bool extended, bool metalOnly)`
    - with the result stored once in `IFPOptions::layoutIndex`

- Output ordering of IFP bitstrings and similarity lines is now consistent:
  - 3 argument and 4 argument modes both print all **docked/query** ligands first, then all **reference** ligands
  - Similarity lines always follow the pattern:
    - `DOCKED_LIGAND_NAME REF_LIGAND_NAME Tanimoto`

### Fixed

- Miscellaneous bugs (from the previous `[Unreleased]` section, now part of this release):
  - The function `getAngl_Tol_AromEF()` returning `AngT_AromFF` instead of `AngT_AromEF`
  - The function `setAngl_Tol_AromEF()` instanciating `AngT_AromFF` instead of `AngT_AromEF`
  - `Bond::~Bond()` destructor, no throwing errors and redirect the logs into the error standard
  - `Atom::getResiduName()` returns invalid reference, fixed
  - Warning implicit conversion from char to unsigned corrected
  - Return default value added in function `Molecule::getalpha()`

---

## [IChem_5.3.5] - 2025-09-05

### Added

- Refactoring:
    - Interactions::detectInteractions() replaces the old monolithic function calcInteractions() 
    - Interactions::addInteractions(): helper function to centralize creation of InterPoint entries and to update InterResults counters (N, O, Zn, CA, CZ, NZ, OD1)
    - Interactions::processAromaticInteractions(): refactored aromatic face/edge detection logic into its own function
    - Interactions::processHydrophobicInteraction(): refactored hydrophobic contact handling into its own function
    - Interactions::checkAromaticHydrophobicInteractions(): separate function that handles hydrophobic interactions between aromatic rings
    - A set of specialized helpers functions called from detectInteractions: 
        - checkMetalNitrogenSulfonamideCase(): detects metal interactions with a sulfonamide nitrogen bound to sulfur
        - checkMetalInteractions(): identifies metal ligand interactions based on the distance
        - checkHydrogenBondLigandAcceptor(): finds hydrogen bonds where the ligand acts as acceptor and protein as donor
        - checkWeakHydrogenBondLigandAcceptor(): detects weak hydrogen bonds with ligand as acceptor and protein as weak donor
        - checkWeakHydrogenBondLigandWeakAcceptor(): identifies weak hydrogen bonds when the ligand is a weak acceptor and protein is weak donor or just donor
        - checkHydrogenBondLigandDonor(): finds hydrogen bonds where the ligand acts as donor and protein as acceptor
        - checkWeakHydrogenBondLigandDonorProteinWeakAcceptor(): detects weak hydrogen bonds where the ligand is donor and protein is weak acceptor
        - checkWeakHydrogenBondLigandWeakDonorProteinAcceptor(): identifies weak hydrogen bonds with ligand as weak donor and protein as weak acceptor or just acceptor
        - checkIonicProteinInteractions(): detects ionic interactions between anionic ligand atoms and cationic protein atoms
        - checkIonicLigandInteractions(): detects ionic interactions between cationic ligand atoms and anionic protein atoms
        - checkPiCationInteraction(): identifies pi-cation interactions between ligand atoms and aromatic protein cycles
    - These helper functions break the previous long function (calcInteractions()) into focused and testable units

    - Spatial search / algorithmic data structures:
        - NeighborSearch KD-tree style index used to build a point cloud of protein atoms and perform fast neighborhood queries
        - The KD-tree can handle multiple ligand atoms in a single call returning a list of nearby protein atoms for each ligand atom efficiently and returning contact structs
        - Use of NeighborSearch::Contact (with ligand_idx, protein_idx, distance_squared) to efficiently iterate contacts
    - Minor data-structure additions:
        - resbest usage extended now used with both ligand and protein atom pointers in hydlist entries to store best hydrophobic candidate per residue
            - In the original struct with only the protein atom, we only track which protein atom was closest and its distance. But when we need to add the interaction we also need to know which ligand atom was involved. Without atmL we lose the information. it’s essential for being able to reconstruct the ligand-protein pair when flushing hydrophobic interactions
        - addInteraction increments NInter and appends to interResult.listInters using emplace_back for efficiency

- Feature IFP calculating the similarity of two IFPs (4 arguments)
- We display all the apolar interactions of type C.ar atoms between cycles
- Query optimized in neighborsearch, replace std::pair with a struct
- Nanoflann well integrated with all the interactions  

### Changed

- Refactoring:
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

- Remove the best hydrophobic interaction in aromatic interaction: we display them all now.
- Function name calcInteractions() into detectInteractions()
- nanoflann well integrated to processAromaticInteractions()
- Refactoring processAromaticInteractions()

### Removed

- Refactoring:
    - Old monolithic implementation calcInteractions() with nested grid.getAdjacency() loops and box scanning removed and replaced by detectInteractions() with neighborsearch algorithm and helpers
    - The out_lig/stdoutt logging logic has been removed
    - The old setStdRules() initialization routine and the static variables such as Interactions::Dist_H, Interactions::Dist_Hyd, Interactions::Angl_H has been replaced by a params style configuration.
    - All PLP logic in the ligand loop are removed
    - grid.getAdjacentAtoms() for hydrophobic density calculations are removed in favour of NeighborSearch

- the if block flags for pication/aromatic interactions in calcInteractionsppi() is unecessary
     + (We need to add the change for everything we added as functions and how they are named)
- CalcInteractionPpi() removed from interaction.cpp
- The function mergeSpeInts() removed as it is used within calcInteractionsppi() function

---

## [IChem_5.3.4] - 2025-06-02

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
- Sending the warnings to stderr: to prevent the pollution of warning messages while parsing mol2 files
- Change menu ints for aromatic and Pi-cation interaction length set to 5.0
- Dist_Arom and Dist_PiCation set to 5.0 from 4.0


### Removed

- ppi related files deleted (detectPPI.cpp, ppi.h), removed from CMakeLists.txt and from modifications in switch.h and switch.cpp 