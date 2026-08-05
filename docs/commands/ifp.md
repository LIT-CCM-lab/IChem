# IFP - Interaction fingerprints

`IChem IFP` detects **protein–ligand interactions** between a **binding site** (protein) and a **ligand**,
then encodes them as a **site-dependent interaction fingerprint** (bit string).

`IFP` can also compare fingerprints:
- against **one or more reference ligands**, or
- between **two protein–ligand complexes**.

!!! important "Site-dependent fingerprint"
    The IFP layout depends on the **binding site residue list and ordering**.  
    Do **not** compare fingerprints generated from different binding sites (different residue sets or different ordering).

---

## Synopsis

### 1) Single complex (compute fingerprint)

```bash
IChem [options] IFP <protein.mol2> <ligand.mol2>
```

### 2) Compare to reference ligand(s)

```bash
IChem [options] IFP <protein.mol2> <ligand.mol2> <ligand_ref.mol2>
```

`<ligand_ref.mol2>` can be a **multi-MOL2** file containing multiple references.

### 3) Compare two complexes

```bash
IChem [options] IFP <protein1.mol2> <ligand1.mol2> <protein2.mol2> <ligand2.mol2>
```

!!! note
    All options must appear **before** the keyword `IFP`.

---

## Inputs

- `protein*.mol2`  
  Protein binding site / active site representation (MOL2).

- `ligand*.mol2`  
  Ligand MOL2. Ligand files can be **multi-MOL2**.

---

## What `IFP` outputs

`IFP` prints three sections to standard output:

1) **Interaction table** (one line per detected interaction)  
2) **Aromatic cycle listing** (only if aromatic interactions are detected)  
3) **Interaction fingerprint** (residue header + binary bit string)

!!! note "Angle column"
    The angle column is printed only for interaction types that have a defined angular criterion.
    For non-directional interactions (e.g., hydrophobic, ionic), the angle is typically shown as `/`.

---

## Interaction table

The interaction table lists all detected interactions with geometry information.

Columns:

- `interaction_type`: interaction category (e.g., `HBond_PROT`, `HBond_LIG`, `Hydrophobic`, `Ionic_LIG`, `Aromatic_Edge/Face`)
- `prot_atom`, `prot_atom_id`: protein atom name and atom index
- `prot_residue_index`: protein residue (name, index, chain)
- `lig_atom`, `lig_atom_id`: ligand atom name and atom index
- `lig_residue_index`: ligand residue identifier (from MOL2)
- `interaction_id`: unique identifier in the output
- `distance (Å)`: interatomic distance
- `angle (°)`: interaction angle (when applicable)

### Example (excerpt)

```bash
./IChem --basic IFP site.mol2 ligand1.mol2
```

| interaction_type   | prot_atom | prot_atom_id | prot_residue_index | lig_atom | lig_atom_id | lig_residue_index | interaction_id | distance (Å) | angle (°) |
|---|---:|---:|---|---|---:|---|---:|---:|---:|
| HBond_PROT         | ND2       | 421          | ASN 312-A          | O17      | 1           | CAU 0-XX          | 0             | 2.76961      | 159.236  |
| HBond_LIG          | OD1       | 78           | ASP 113-A          | O17      | 1           | CAU 0-XX          | 1             | 2.60795      | 134.667  |
| Hydrophobic        | CZ        | 349          | PHE 289-A          | C16      | 2           | CAU 0-XX          | 3             | 3.74597      | /        |
| Ionic_LIG          | OD1       | 78           | ASP 113-A          | N19      | 4           | CAU 0-XX          | 4             | 3.60381      | /        |
| Aromatic_Edge/Face | DuAr      | NULL         | PHE 290-A          | DuAr     | NULL        | CAU 0-XX          | 56            | 4.68033      | 105.645  |

---

## Fingerprint definition

### Interaction modes (what bits are computed)

One (and only one) interaction mode should be selected:

| Option | Meaning |
|---|---|
| `--all`   | Compute all interactions (basic + pi-cation + metal + weak H-bond) |
| `--basic` | Compute only basic interactions (other bits are set to `0`) |
| `--picat` | Compute only pi-cation bit(s) (other bits are set to `0`) |
| `--metal` | Compute only metal/acceptor bit(s) (other bits are set to `0`) |
| `--weakh` | Compute only weak H-bond bit(s) (other bits are set to `0`) |
| `--old`   | Output legacy 7-bit fingerprint format (basic only) |

Examples:

```bash
IChem --all   IFP protein.mol2 ligand.mol2
IChem --basic IFP protein.mol2 ligand.mol2
IChem --picat IFP protein.mol2 ligand.mol2
IChem --metal IFP protein.mol2 ligand.mol2
IChem --weakh IFP protein.mol2 ligand.mol2
IChem --old   IFP protein.mol2 ligand.mol2
```

!!! note
    The fingerprint structure is always kept: bits not enabled by the selected mode are set to `0`.

---

## Detection rules

Interactions are detected when the distance is within a window:

> **min ≤ distance ≤ max**, plus angle constraints when applicable.

!!! important "Precision"
    Distance thresholds accept **at most 3 decimals**.

### Maximum distance thresholds (Å)

| Option | Default | Interaction |
|---|---:|---|
| `-D_Hb <N>`  | 3.5 | H-bond |
| `-D_Hyd <N>` | 4.5 | Hydrophobic |
| `-D_Io <N>`  | 4.0 | Ionic |
| `-D_Me <N>`  | 2.8 | Metal/Acceptor |
| `-D_Ar <N>`  | 5.0 | Aromatic |
| `-D_Pic <N>` | 5.0 | Pi-cation |
| `-D_WHb <N>` | 3.5 | Weak H-bond |

### Minimum distance thresholds (Å)

| Option | Default | Interaction |
|---|---:|---|
| `-d_Hb <N>`  | 2.3 | H-bond |
| `-d_Hyd <N>` | 3.2 | Hydrophobic |
| `-d_Io <N>`  | 2.3 | Ionic |
| `-d_Me <N>`  | 1.8 | Metal/Acceptor |
| `-d_Ar <N>`  | 3.2 | Aromatic |
| `-d_Pic <N>` | 2.3 | Pi-cation |
| `-d_WHb <N>` | 2.3 | Weak H-bond |

### Angle thresholds (degrees)

| Option | Default | Meaning |
|---|---:|---|
| `-a_H <N>`     | 180 | H-bond reference angle |
| `-at_H <N>`    | 60  | H-bond tolerance |
| `-a_ArFF <N>`  | 180 | Aromatic face-to-face reference angle |
| `-at_ArFF <N>` | 30  | Aromatic face-to-face tolerance |
| `-a_ArEF <N>`  | 90  | Aromatic edge-to-face reference angle |
| `-at_ArEF <N>` | 60  | Aromatic edge-to-face tolerance |
| `-a_Pic <N>`   | 180 | Pi-cation reference angle |
| `-at_Pic <N>`  | 30  | Pi-cation tolerance |

Example (tuning hydrophobic window):

```bash
IChem -d_Hyd 3.34 -D_Hyd 4.81 --basic IFP site.mol2 ligand.mol2
```

---

## Bit layout

### Legacy 7-bit layout (`--old`)

| Bit | Interaction type |
|---:|---|
| 1 | Hydrophobic |
| 2 | Aromatic (face-to-face) |
| 3 | Aromatic (edge-to-face) |
| 4 | H-bond (protein donor) |
| 5 | H-bond (ligand donor) |
| 6 | Ionic (protein charged +) |
| 7 | Ionic (ligand charged +) |

### Default 11-bit layout

The default fingerprint expands the basic 7-bit block with additional bits:

**Bits 1–7 (basic)**

| Bit | Interaction type |
|---:|---|
| 1 | Hydrophobic |
| 2 | Aromatic (face-to-face) |
| 3 | Aromatic (edge-to-face) |
| 4 | H-bond (protein donor) |
| 5 | H-bond (ligand donor) |
| 6 | Ionic (protein charged +) |
| 7 | Ionic (ligand charged +) |

**Bit 8 (pi-cation)**

| Bit | Interaction type |
|---:|---|
| 8 | Pi-cation |

**Bit 9 (metal/acceptor)**

| Bit | Interaction type |
|---:|---|
| 9 | Metal/Acceptor |

**Bits 10–11 (weak H-bond)**

| Bit | Interaction type |
|---:|---|
| 10 | Weak H-bond (protein donor) |
| 11 | Weak H-bond (ligand donor) |

---

## Fingerprint output format

After the interaction table (and optional aromatic cycle section), `IFP` prints:

1) A **residue header row** (the residues defining the fingerprint layout)  
2) The **binary fingerprint string** (`0`/`1`) concatenating the bit blocks for each residue in order

Example (header + bitstring):

| A M82 | A V86 | A W109 | A T110 | A D113 | A V114 | A L115 | A V117 | A T118 | A C191 | A F193 | A T195 | A Y199 | A A200 | A I201 | A S203 | A S204 | A I205 | A S207 | A F208 | A W286 | A F289 | A F290 | A N293 | A Y308 | A N312 | A Y316 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|

```text
0000000000000000000000000010000000000000000000000000000100000000000000000000010000000000100000000001010000000000000000000000000000000001000000010000000000
```

!!! tip
    Each residue contributes a fixed number of bits (7 for `--old`, 11 for default).
    The final fingerprint string is the concatenation of residue blocks in the header order.

---

## Aromatic cycle section

If aromatic interactions exist, `IFP` may print an additional section listing the ligand atoms involved in aromatic cycles.

Example:

| ==== Ligand atoms involved in the aromatic cycle ==== | |
|---|---|
| Interaction 56 | C13 C12 C11 C10 C8 C9 |

---

## Usage notes

- Always use explicit `.mol2` extensions in command lines.
- Ensure MOL2 files are valid and in the same coordinate frame.
- For reproducible comparisons, keep the same binding site residue definition and ordering.
- If you need to compute fingerprint similarities across many ligands, combine `IFP` with [`sims`](sims.md).

---

## See also

- `ints` - interaction pseudoatoms and generic fingerprints (`commands/ints.md`)
- `sims` - fingerprint similarity metrics (`commands/sims.md`)
- `AtomProps` - atom property annotation (`commands/atomprops.md`)
- `grim` - interaction-based alignment (`commands/grim.md`)
