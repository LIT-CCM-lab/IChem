# IFP — Interaction FingerPrint


`IChem IFP` detects protein-ligand interactions between a **binding site** and a **ligand**, then encodes them as an **interaction fingerprint**. It can also compute similarity against reference ligand, and can compare two protein-ligand pairs.

---

## Usage

The IFP is executed from the command line using one of the following syntaxes:

```bash
IChem [options] IFP protein ligand
# Single complex analysis: computes the interaction fingerprint of one protein-ligand complex
```bash

IChem [options] IFP protein ligand ligand_ref
# Comparison with a reference ligand: computes the fingerprint and compares it to a reference ligand
```bash

IChem [options] IFP protein1 ligand1 protein2 ligand2
# Comparison of two complexes: computes and compares the interaction fingerprints of two protein-ligand complexes
```

**Parameters**

- **IChem**: Main executable of the program
- **[options]**: Command-line options described below
- **IFP**: Name of the command
- **protein**: Protein structure file in MOL2 format
- **ligand**: Ligand structure file in MOL2 format
- **ligand_ref**: Reference ligand structure file in MOL2 format
- **protein1 ligand1 protein2 ligand2**: Two protein-ligand complexes used to compute and compare their interaction fingerprints

Note:

- Ligand files can be **multi-mol2**

---

## IFP Options

In the option, there are main options in the IFP as follows:

### Main interaction modes

One of these options has to be used necessary, but restricted to only one option:

| Option    | Meaning                                                          |
| --------- | ---------------------------------------------------------------- |
| `--all`   | Compute **all** interactions (basic and pi-cation and metal and weakH) |
| `--basic` | Compute **only basic** interactions (others = `0`)               |
| `--picat` | Compute **only pi-cation** bit(s) (others = `0`)                 |
| `--metal` | Compute **only metal/acceptor** bit(s) (others = `0`)            |
| `--weakh` | Compute **only weak H-bond** bit(s) (others = `0`)               |
| `--old`   | Output **legacy 7-bit format**                                   |

Examples:
```bash
IChem --all   IFP protein.mol2 ligand.mol2
IChem --basic IFP protein.mol2 ligand.mol2
IChem --picat IFP protein.mol2 ligand.mol2
IChem --metal IFP protein.mol2 ligand.mol2
IChem --weakh IFP protein.mol2 ligand.mol2
IChem --old   IFP protein.mol2 ligand.mol2
```

## Detection rules

Interactions are detected when the distance is within the allowed window:

> **min ≤ distance ≤ max**, plus angle constraints when applicable

### Maximum thresholds distances Å

| Option     | Default | Interaction    |
| ---------- | ------: | -------------- |
| `-D_Hb N`  |     3.5 | H-bond         |
| `-D_Hyd N` |     4.5 | Hydrophobic    |
| `-D_Io N`  |     4.0 | Ionic          |
| `-D_Me N`  |     2.8 | Metal/Acceptor |
| `-D_Ar N`  |     5.0 | Aromatic       |
| `-D_Pic N` |     5.0 | Pi-cation      |
| `-D_WHb N` |     3.5 | Weak H-bond    |

### Minimum thresholds distances Å

| Option     | Default | Interaction    |
| ---------- | ------: | -------------- |
| `-d_Hb N`  |     2.3 | H-bond         |
| `-d_Hyd N` |     3.2 | Hydrophobic    |
| `-d_Io N`  |     2.3 | Ionic          |
| `-d_Me N`  |     1.8 | Metal/Acceptor |
| `-d_Ar N`  |     3.2 | Aromatic       |
| `-d_Pic N` |     2.3 | Pi-cation      |
| `-d_WHb N` |     2.3 | Weak H-bond    |

### Angles in degrees

| Option       | Default | Meaning                |
| ------------ | ------: | ---------------------- |
| `-a_H N`     |     180 | H-bond angle           |
| `-at_H N`    |      60 | H-bond tolerance       |
| `-a_ArFF N`  |     180 | Aromatic face-to-face  |
| `-at_ArFF N` |      30 | Face-to-face tolerance |
| `-a_ArEF N`  |      90 | Aromatic edge-to-face  |
| `-at_ArEF N` |      60 | Edge-to-face tolerance |
| `-a_Pic N`   |     180 | Pi-cation angle        |
| `-at_Pic N`  |      30 | Pi-cation tolerance    |


Examples:
```bash
IChem -d_Hyd 3.34 -D_Hyd 4.81 --basic IFP site.mol2 ligand.mol2
```

---

## Fingerprint formats and bits

### Fingerprint formats

The IFP has a default format of 11 bits using the main options as described above

- **Default**: **11-bit** fixed structure format: `[basic][pi-cation][metal][weakh]`

We kept also the old version that is coded in 7 bits:

- `--old` : legacy **7-bit** format: `[basic]`

Notes:

- Note that in old version, pi-cation, metal and weak H-bond interactions are not detected
- The structure is **always kept** (unused bits = `0`)

### Bit layout per residue

#### Legacy 7-bit layout

**[--old]** 

| Bit | Interaction type          |
| --: | ------------------------- |
|  1  | Hydrophobic               |
|  2  | Aromatic (face-to-face)   |
|  3  | Aromatic (edge-to-face)   |
|  4  | H-bond (protein donor)    |
|  5  | H-bond (ligand donor)     |
|  6  | Ionic (protein charged +) |
|  7  | Ionic (ligand charged +)  |

#### Default 11-bit layout


**[--basic]** 

| Bit | Interaction type          |
| --: | ------------------------- |
|  1  | Hydrophobic               |
|  2  | Aromatic (face-to-face)   |
|  3  | Aromatic (edge-to-face)   |
|  4  | H-bond (protein donor)    |
|  5  | H-bond (ligand donor)     |
|  6  | Ionic (protein charged +) |
|  7  | Ionic (ligand charged +)  |

**[--picat]**

| Bit | Interaction type          |
| --: | ------------------------- |
|  8  | Pi-cation                 |

**[--metal]**

| Bit | Interaction type          |
| --: | ------------------------- |
|  9  | Metal/Acceptor            |

**[--weakh]**

| Bit | Interaction type          |
| --: | ------------------------- |
|  10 | Weak H-bond (protein donor) |
|  11 | Weak H-bond (ligand donor)  |

---

## Example

./IChem --basic IFP site.mol2 ligand1.mol2 

| interaction_type       | prot_atom | prot_atom_id | prot_residue_index | lig_atom | lig_atom_id | lig_residue_index | interaction_id  | distance (Å) | angle (°) |
|------------------------|-----------|--------------|-------------|----------|---------------|-------|-----|--------------|-----------|
| HBond_PROT             | ND2       | 421          | ASN 312-A   | O17      | 1     | CAU 0-XX  | 0   | 2.76961      | 159.236   |
| HBond_LIG              | OD1       | 78           | ASP 113-A   | O17      | 1     | CAU 0-XX  | 1   | 2.60795      | 134.667   |
| Hydrophobic            | CZ        | 349          | PHE 289-A   | C16      | 2     | CAU 0-XX  | 3   | 3.74597      | /         |
| Ionic_LIG              | OD1       | 78           | ASP 113-A   | N19      | 4     | CAU 0-XX  | 4   | 3.60381      | /         |
| Ionic_LIG              | OD2       | 79           | ASP 113-A   | N19      | 4     | CAU 0-XX  | 5   | 2.94042      | /         |
| Hydrophobic            | CZ3       | 46           | TRP 109-A   | C21      | 6     | CAU 0-XX  | 6   | 4.43435      | /         |
| Hydrophobic            | CG2       | 64           | THR 110-A   | C21      | 6     | CAU 0-XX  | 7   | 4.42961      | /         |
| Hydrophobic            | CB        | 163          | PHE 193-A   | C21      | 6     | CAU 0-XX  | 8   | 4.07199      | /         |
| Hydrophobic            | CH2       | 47           | TRP 109-A   | C22      | 7     | CAU 0-XX  | 9   | 3.72321      | /         |
| Hydrophobic            | CZ        | 438          | TYR 316-A   | C22      | 7     | CAU 0-XX  | 10  | 4.46257      | /         |
| Hydrophobic            | CG2       | 125          | VAL 117-A   | C15      | 8     | CAU 0-XX  | 11  | 4.0788       | /         |
| Hydrophobic            | CE2       | 348          | PHE 289-A   | C15      | 8     | CAU 0-XX  | 13  | 4.0318       | /         |
| HBond_LIG              | OG        | 248          | SER 203-A   | N7       | 16    | CAU 0-XX  | 23  | 3.31671      | 127.808   |
| Hydrophobic            | CG2       | 185          | THR 195-A   | C1       | 20    | CAU 0-XX  | 24  | 4.32061      | /         |
| Hydrophobic            | CE2       | 348          | PHE 289-A   | C13      | 10    | CAU 0-XX  | 25  | 4.35445      | /         |
| Aromatic_Edge/Face     | DuAr      | NULL         | PHE 290-A   | DuAr     | NULL  | CAU 0-XX  | 56  | 4.68033      | 105.645   |
| Hydrophobic            | CZ3       | 327          | TRP 286-A   | C16      | 2     | CAU 0-XX  | 57  | 4.05558      | /         |
| Hydrophobic            | CG1       | 124          | VAL 117-A   | C12      | 11    | CAU 0-XX  | 60  | 3.98025      | /         |
| Hydrophobic            | CG1       | 89           | VAL 114-A   | C11      | 12    | CAU 0-XX  | 61  | 3.83646      | /         |
| Hydrophobic            | CB        | 288          | SER 207-A   | C10      | 13    | CAU 0-XX  | 62  | 3.68825      | /         |
| Hydrophobic            | CE2       | 168          | PHE 193-A   | C5       | 22    | CAU 0-XX  | 76  | 3.73         | /         |
| Hydrophobic            | CE2       | 348          | PHE 289-A   | C5       | 22    | CAU 0-XX  | 79  | 3.91129      | /         |
| Hydrophobic            | CZ        | 349          | PHE 289-A   | C5       | 22    | CAU 0-XX  | 80  | 3.96371      | /         |
| Hydrophobic            | CG2       | 90           | VAL 114-A   | C13      | 10    | CAU 0-XX  | 84  | 3.81147      | /         |
| Hydrophobic            | CE2       | 168          | PHE 193-A   | C4       | 18    | CAU 0-XX  | 88  | 4.16691      | /         |
| Hydrophobic            | CE2       | 202          | TYR 199-A   | C2       | 19    | CAU 0-XX  | 91  | 3.86117      | /         |
| Hydrophobic            | CE2       | 168          | PHE 193-A   | C6       | 21    | CAU 0-XX  | 96  | 3.49585      | /         |
| Hydrophobic            | CE2       | 368          | PHE 290-A   | C8       | 14    | CAU 0-XX  | 102 | 3.51967      | /         |


| ==== Ligand atoms involved in the aromatic cycle ==== |        |
|-------------------------------------------------------|--------|
| Interaction 56                                        | C13 C12 C11 C10 C8 C9 |


| A M82 | A V86 | A W109 | A T110 | A D113 | A V114 | A L115 | A V117 | A T118 | A C191 | A F193 | A T195 | A Y199 | A A200 | A I201 | A S203 | A S204 | A I205 | A S207 | A F208 | A W286 | A F289 | A F290 | A N293 | A Y308 | A N312 | A Y316 |
|-------|-------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|



```
0000000000000000000000000010000000000000000000000000000100000000000000000000010000000000100000000001010000000000000000000000000000000001000000010000000000
```

### Output format explanation


The program produces three sections:

1. **Interaction Table**  
    - Lists all detected protein-ligand interactions
    - Columns:
        - interaction_type: HBond_PROT, HBond_LIG, Hydrophobic, Ionic_LIG, Aromatic_Edge/Face
        - prot_atom: Protein atom involved
        - prot_atom_id: Protein atom index
        - prot_residue_index: Protein residue name and index
        - lig_atom: Ligand atom involved
        - lig_atom_id: Ligand atom index
        - lig_residue_index: Ligand residue name and index
        - interaction_id: Interaction id
        - distance (Å): Interatomic distance
        - Angle (°): Interaction angle

2. **Ligand Aromatic Cycle atoms** if aromatic interactions exist
    - List of ligand atoms involved in aromatic interactions

3. **Binary Interaction Fingerprint**  
    - A 0/1 string representing interaction presence for each protein residue

---

## Usage notes
- **All options must appear before** the keyword `IFP`
- Distance thresholds **(e.g D_Hyd, d_Hyd)** accept **at most 3 decimals**
- The file should explicitly mention the extension .mol2
- Ensure MOL2 files are valid


