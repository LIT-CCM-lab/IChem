# pdbconv - PDB processing and conversion

`IChem pdbconv` parses and processes **PDB** (or **MOL2**) structures, can automatically detect bound ligands (HET codes),
and can export ready-to-use **MOL2 representations** for downstream IChem workflows.

Typical uses:

- clean and normalize a PDB before analysis
- convert a PDB to MOL2 (molecular representation)
- split structures into protein / ligand / solvent components
- for each detected ligand, export a ligand-free protein MOL2 and cavity files (VolSite) and descriptors

---

## Synopsis

```bash
IChem pdbconv <protein.pdb|protein.mol2> <output_dir> <pdb_id>
```

- `output_dir` must exist before running the command.
- `pdb_id` is used as a prefix for output naming.

---

## Inputs

- `protein.pdb` (recommended)  
  Input PDB file to process.

- `protein.mol2` (with `--wMOL2`)  
  Use MOL2 as input; PDB-specific options are not available.

---

## Outputs

`pdbconv` writes multiple files into the given `output_dir`.

For PDB inputs, for each detected ligand, it typically produces:
- a ligand MOL2
- the corresponding ligand-free protein MOL2
- cavity MOL2 files (VolSite cavities)
- cavity descriptor file: `*_descriptor.txt`

Example (from test data):

```bash
IChem pdbconv 2RH1.pdb output 2rh1
```

Input directory: `$ICHEM_DIR/test/PDBCONV`  
Input file: `2RH1.pdb`  
Output directory: `./output`

!!! note
    The exact output filenames may depend on IChem version/build and selected options.
    The output directory must be created before executing the command.

---

## Options overview

### Top-level options

| Option | Description |
|---|---|
| `--wMOL2` | Use a MOL2 file as input. PDB options are not available. |
| `--wUnDrug` | Output undruggable cavities. |
| `--noLig` | Treat the PDB as having no ligand (disable ligand detection). |

---

## PDB options

By default, all PDB options below are enabled and **all chains are kept**.

| Option | Description |
|---|---|
| `--HARMSIZE` | Harmonize line size to 80 characters |
| `--MSEMET` | Convert `MSE` residues to `MET` |
| `--CSECYS` | Convert `CSE` residues to `CYS` |
| `--MOVHET` | Move `HETATM` records to the end of file |
| `--ALTATM` | Select alternative atoms |
| `--NUMATM` | Renumber atoms |
| `--UPDMAS` | Update the `MASTER` line |
| `--TOMOL2` | Convert to a molecular representation (instead of a flat PDB) |

!!! important
    If you use any **MOL2 options** below, you **must** also enable `--TOMOL2`
    (or use `--wMOL2` when starting from a MOL2 file).

---

## MOL2 options

These options apply to the molecular representation produced by `--TOMOL2` (or used with `--wMOL2`).

| Option | Description |
|---|---|
| `--RESTYP` | Apply residue classes (cofactor / STD_AA / MOD_AA / ligand …) |
| `--BONDSE` | Create bonds |
| `--CLNUNW` | Clean unwanted residues |
| `--MOL2TY` | Apply MOL2 atom types according to templates |
| `--SPLITM` | Split molecule into protein / ligand / solvent |
| `-SelChain <A_B_...>` | List of chains to keep, separated by underscore |
| `--SELWAT` | Select water molecules |
| `--SELLIG` | Select ligand |

---

## Examples

### Example 1 - Standard PDB processing

```bash
mkdir -p output
IChem pdbconv 2RH1.pdb output 2rh1
```

This will parse the PDB, detect ligands (unless `--noLig` is set), and export MOL2 structures and cavity-related outputs
into `output/`.

### Example 2 - Convert PDB to MOL2 representation

```bash
mkdir -p output
IChem --TOMOL2 pdbconv 2RH1.pdb output 2rh1
```

### Example 3 - Keep only specific chains (requires molecular representation)

```bash
mkdir -p output
IChem --TOMOL2 -SelChain A_B pdbconv 2RH1.pdb output 2rh1
```

### Example 4 - Use a MOL2 input (PDB options disabled)

```bash
mkdir -p output
IChem --wMOL2 pdbconv protein.mol2 output myprot
```

---

## Practical guidance

- Use `--noLig` when processing apo structures or when ligand detection is not desired.
- Use `-SelChain` to reduce the system to relevant chains before running downstream tools.
- If you plan to run `volsite` or `ints/IFP` next, prefer enabling `--TOMOL2` with appropriate MOL2 options (`--MOL2TY`, `--BONDSE`, `--SPLITM`) so the structure is in a clean molecular form.

---

## See also

- `volsite` - cavity detection and druggability prediction (`commands/volsite.md`)
- `AtomProps` - atom property annotation (`commands/atomprops.md`)
- `ints` - interaction pseudoatoms (`commands/ints.md`)
