# ints - Interaction generator

`IChem ints` detects **protein-ligand interactions** and exports them as:

- **Interaction pseudoatoms** (`*_INTS_*.mol2`) for visualization and downstream geometric analysis
- **Generic interaction fingerprints** (`.fgp`) built from **triplets of pseudoatoms** (binding-site independent)

It can also write an **interaction listing** (`-logf`) containing the same information as the interaction table produced by `IChem IFP`.

!!! important
    `ints` fingerprints are **binding-site independent** (triplets of pseudoatoms).  
    If you need **site-dependent** fingerprints whose layout depends on the residue list and ordering, use `IChem IFP`.

---

## Synopsis

### List interactions and generate pseudoatoms

```bash
IChem ints <site.mol2> <ligand.mol2>
IChem -logf <ints_table.txt> ints <site.mol2> <ligand.mol2>
IChem -type CENT ints <site.mol2> <ligand.mol2>
```

### Generate a fingerprint (`.fgp`)

```bash
IChem -fgps STD  ints <site.mol2> <ligand.mol2> <out.fgp>
IChem --small -fgps STD ints <site.mol2> <ligand.mol2> <out_small.fgp>
```

!!! note
    Options must be placed **before** the keyword `ints`, for example:  
    `IChem -logf out.txt -type CENT ints site.mol2 ligand.mol2`

---

## Inputs

- `site.mol2`  
  Binding site / active site representation (**protein MOL2**).

- `ligand.mol2`  
  Ligand MOL2 file.

---

## Outputs

Depending on the options used, `ints` can write:

### 1) Interaction listing (`-logf`)

A text file (e.g. `2rh1_ints.txt`) containing the same information as the interaction table output by `IChem IFP`.

Example:

```bash
IChem -logf 2rh1_ints.txt -type CENT ints site.mol2 ligand.mol2
```

### 2) Interaction pseudoatoms (`*_INTS_*.mol2`)

A MOL2 file describing the **properties and atomic coordinates** of interaction pseudoatoms.

Pseudoatoms are mapped according to `-type`:

- `LIG`  : pseudoatoms mapped onto **ligand-interacting atoms**
- `PROT` : pseudoatoms mapped onto **protein-interacting atoms**
- `CENT` : pseudoatoms mapped onto the **barycenter** of ligand- and protein-interacting atoms (**default**)
- `MERG` : pseudoatoms mapped onto **all** the above (LIG + PROT + CENT)

Multiple values are allowed, separated by spaces (e.g. `-type CENT PROT`).

Typical filename suffixes:

| `-type` mode | suffix |
|---|---|
| `CENT` | `_INTS_C.mol2` |
| `LIG`  | `_INTS_L.mol2` |
| `PROT` | `_INTS_P.mol2` |
| `MERG` | `_INTS_M.mol2` |

The IChem distribution also provides header templates for these modes:

| Template filename | `-type` mode |
|---|---|
| `Header_INTS_C.mol2` | `CENT` |
| `Header_INTS_L.mol2` | `LIG` |
| `Header_INTS_P.mol2` | `PROT` |
| `Header_INTS_M.mol2` | `MERG` |

### 3) Fingerprint file (`.fgp`)

If you provide an explicit output filename as the last argument, `ints` writes a **generic fingerprint** describing **triplets of interaction pseudoatoms**.

Example:

```bash
IChem -fgps STD ints site.mol2 ligand.mol2 2rh1_full.fgp
```

!!! note
    The `.fgp` file can be **appended** by running additional `ints` commands with the **same output filename**.

---

## Interaction pseudoatoms

### Pseudoatom types (atom names)

In the pseudoatom MOL2 output, the **atom name** encodes the interaction category:

| Atom name | Interaction represented |
|---|---|
| `CA`  | Hydrophobic |
| `CZ`  | Aromatic |
| `O`   | Hydrogen bond (protein atom is **acceptor**) |
| `OG`  | Hydrogen bond (protein atom is **acceptor and donor**) |
| `N`   | Hydrogen bond (protein atom is **donor**) |
| `OD1` | Ionic (protein atom is **negatively charged**) |
| `NZ`  | Ionic (protein atom is **positively charged**) |
| `ZN`  | Metal coordination |

This encoding makes it straightforward to color/select pseudoatoms in PyMOL/Chimera (e.g., select atom name `CA` for hydrophobic contacts).

### Pseudoatom positioning (`-type`)

`-type` selects where pseudoatoms are placed:

- `LIG`: on the ligand interacting atoms
- `PROT`: on the protein interacting atoms
- `CENT`: at the barycenter of the interacting ligand and protein atoms (default)
- `MERG`: include all placements above in a single output

Practical guidance:

- Use `CENT` for a **single representative point per interaction**
- Use `PROT` to annotate the **binding site pocket**
- Use `LIG` for **ligand-centric** visualization
- Use `MERG` when you need **all representations** for downstream workflows

---

## Generic fingerprints (`.fgp`)

### What the fingerprint represents

`IChem ints` can output **binding-site independent** interaction fingerprints consisting of **all possible combinations of triplets** of interaction pseudoatoms.

### Fingerprint formats (`-fgps`)

`-fgps` controls how the fingerprint vector is written:

| `-fgps` | Name | Example | Description |
|---|---|---|---|
| `STD` | Standard | `1 0 21 0 0 3` | All values written, space-separated. Simple but inefficient for long sparse vectors. |
| `SVM` | SVM format | `1:1 3:21 6:3` | Only non-zero values as `POSITION:VALUE`, space-separated. Efficient for long sparse vectors. |
| `CMP` | Compressed | `1 [1 21 [2 3` | Non-zero values are explicit; runs of zeros are encoded by `[` followed by the number of consecutive zeros. Efficient for long sparse vectors. |

Recommendations:

- For **sparse** fingerprints (many zeros), prefer `SVM` or `CMP`.
- For simple debugging or small/dense vectors, `STD` is the easiest to read.

### Small fingerprint (`--small`)

By default, the standard fingerprint can be large (e.g. **12510 integers**), many of which encode triplets describing hydrophobic contacts.  
`--small` prunes rare triplets to generate a much smaller fingerprint (e.g. **210 integers**).

Example:

```bash
IChem --small -fgps STD ints site.mol2 ligand.mol2 2rh1_small.fgp
```

---

## Options

### Output control

| Option | Default | Description |
|---|---:|---|
| `-type <modes...>` | `CENT` | Pseudoatom positioning. Multiple values allowed (space-separated): `PROT`, `LIG`, `CENT`, `MERG`. |
| `-fgps <STD|SVM|CMP>` | `STD` | Fingerprint output format (only used when writing a `.fgp` file). |
| `--small` | - | Output the pruned “small” fingerprint. |

### General options

| Option | Default | Description |
|---|---:|---|
| `-name <str>` | `prot` | Name of molecule in output file. |
| `-logf <file>` | - | Name of log file for the interaction listing. |

### Distance thresholds (Å)

| Option | Default | Description |
|---|---:|---|
| `-D_Hb <N>`  | `3.5` | Hydrogen bond maximum distance (Å). |
| `-D_Hyd <N>` | `4.5` | Hydrophobic maximum distance (Å). |
| `-D_Io <N>`  | `4.0` | Ionic maximum distance (Å). |
| `-D_Me <N>`  | `2.8` | Metal/acceptor maximum distance (Å). |
| `-D_Ar <N>`  | `5.0` | Aromatic interaction maximum distance (Å). |
| `-D_Pic <N>` | `5.0` | Pi-cation interaction maximum distance (Å). |

### Angle thresholds (radians)

| Option | Default | Description |
|---|---:|---|
| `-a_H <N>`     | `Pi`   | Hydrogen bond reference angle (rad). |
| `-at_H <N>`    | `Pi/3` | Hydrogen bond tolerance angle (rad). |
| `-a_ArFF <N>`  | `Pi`   | Aromatic face-to-face reference angle (rad). |
| `-at_ArFF <N>` | `Pi/6` | Aromatic face-to-face tolerance angle (rad). |
| `-a_ArEF <N>`  | `Pi/2` | Aromatic edge-to-face reference angle (rad). |
| `-at_ArEF <N>` | `Pi/3` | Aromatic edge-to-face tolerance angle (rad). |
| `-a_Pic <N>`   | `Pi`   | Pi-cation reference angle (rad). |
| `-at_Pic <N>`  | `Pi/6` | Pi-cation tolerance angle (rad). |

### Hydrophobic post-processing

| Option | Default | Description |
|---|---:|---|
| `--noMerge` | - | Do not merge hydrophobic interactions. |
| `--newH` | - | Use a less permissive definition of hydrophobic interaction. |

!!! note
    The authoritative defaults depend on your IChem version/build.  
    To list defaults from your installation, run: `IChem ints -h`.

---

## Examples

### Example 1 - Centered pseudoatoms and interaction listing

```bash
IChem -logf 2rh1_ints.txt -type CENT ints site.mol2 ligand.mol2
```

Expected outputs:
- `2rh1_ints.txt`
- `*_INTS_C.mol2`

### Example 2 - Pseudoatoms mapped onto protein atoms

```bash
IChem -type PROT ints site.mol2 ligand.mol2
```

Expected output:
- `*_INTS_P.mol2`

### Example 3 - Full fingerprint in SVM format

```bash
IChem -fgps SVM ints site.mol2 ligand.mol2 2rh1_full.svm.fgp
```

### Example 4 - Small fingerprint in compressed format

```bash
IChem --small -fgps CMP ints site.mol2 ligand.mol2 2rh1_small.cmp.fgp
```

---

## Troubleshooting

### No interactions detected
- Verify that `site.mol2` contains the correct binding site atoms and is in the same coordinate frame as `ligand.mol2`.
- If needed, increase relevant distance thresholds (e.g. `-D_Hb`, `-D_Hyd`) carefully and re-run.

### Too many hydrophobic contacts
- Try `--newH` (stricter hydrophobic definition).
- Consider `--noMerge` when you need to inspect individual hydrophobic contacts instead of merged ones.

---

## See also

- `IFP` - site-dependent interaction fingerprints (`commands/ifp.md`)
- Concepts / Reference - interaction definitions and geometric rules
