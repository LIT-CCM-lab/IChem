# utils - Miscellaneous tools

`IChem utils` groups small utilities for common tasks around protein–ligand complexes.
Currently documented subcommands:

- `bsa`  : buried surface area calculation for a bound ligand
- `frag` : protein-bound ligand fragmentation (ring-core/substituent detection)

---

## Synopsis

### Buried surface area (BSA)

```bash
IChem utils bsa <protein.mol2> <ligand.mol2>
```

### Ligand fragmentation

```bash
IChem utils frag <protein.mol2> <ligand.mol2>
```

---

## Input requirements

For both subcommands:

- `protein.mol2`: protein (binding site or full protein) in MOL2 format
- `ligand.mol2`: ligand in MOL2 format

!!! important
    Protein and ligand coordinates must be in the **same coordinate frame** (same reference system) in MOL2 format.

---

## utils bsa - Buried surface area

### What it does
`IChem utils bsa` computes the **buried surface area** of the ligand upon binding to the protein.

### Command

```bash
IChem utils bsa protein.mol2 ligand.mol2 > ligand.bsa
```

### Output

The output is a **3-column table**:

```text
REF      71.4338  348.375
```

Columns:
1. protein–ligand identifier (header)
2. buried surface area (**%**)
3. ligand volume (**Å³**)

### Example (test set)

Input directory: `$ICHEM_DIR/test/BSA`  
Input files: `protein.mol2`, `ligand.mol2`  
Typical output file: `ligand.bsa`

---

## utils frag - Protein-bound ligand fragmentation

### What it does
`IChem utils frag` fragments a **protein-bound ligand** in 3D space using a method aimed at detecting substituted ring cores.

High-level workflow:
1. Detect aromatic and aliphatic rings (ring perception).
2. Classify acyclic atoms as **linkers** or **substituents** depending on whether they connect rings.
3. Keep linker atoms unchanged.
4. For substituents, cleave selected single bonds near the ring system:
   - the bond involves the closest apolar carbon (by bond distance) to any ring
   - the cleaved bond must be at least **3 bonds** away from the cyclic root atom
5. Add an anchoring atom (`Z` label) to each fragment to mark the cleavage site.

This is useful for:
- identifying core ring systems vs substituents
- building fragment libraries for bioisosteric replacement workflows

### Command

```bash
IChem utils frag protein.mol2 ligand.mol2
```

### Outputs

Several MOL2 files are written. Typical outputs (names depend on input headers and fragment index):

- `*_FRAG_<i>_MOLE.mol2`  
  Fragment extracted from the original ligand

- `*_FRAG_<i>_INTS.mol2`  
  Interaction pseudoatoms (merged mode) for the fragment

- `*.mol2`  
  Interaction pseudoatoms (merged mode) for the entire ligand

Example outputs (from test set):
- `2rh1_CAU_1_protein_FRAG_1_MOLE.mol2`
- `2rh1_CAU_1_protein_FRAG_1_INTS.mol2`
- `xxx.mol2`

!!! note
    Several fragments can be generated. Each fragment has an index:
    `header_FRAG_<index>.mol2` and the corresponding pseudoatoms
    `header_FRAG_<index>_INTS.mol2`.

### Example (test set)

Input directory: `$ICHEM_DIR/test/Frag`  
Input files: `protein.mol2`, `ligand.mol2`  
Outputs: fragment MOL2 files and pseudoatom MOL2 files as described above.

### Reference
Desaphy J, Rognan D. *J. Chem. Inf. Model.* 2014;54(7):1908–1918.  
sc-PDB-Frag: a database of protein–ligand interaction patterns for bioisosteric replacements.

---

## See also

- `ints` - interaction pseudoatoms (`commands/ints.md`)
- `IFP` - interaction fingerprints (`commands/ifp.md`)
- `volsite` - cavity detection and druggability (`commands/volsite.md`)
