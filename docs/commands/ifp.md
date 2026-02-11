# IFP — Interaction FingerPrint


`IChem IFP` detects protein-ligand interactions between a **binding site** and a **ligand**, then encodes them as an **interaction fingerprint**. It can also compute similarity against reference ligand, and can compare two protein-ligand pairs.

---

## Application

1. Input one (or two) protein-ligand pairs in MOL2
2. Detect interactions using **geometric rules** (optionally min/max distance and angles can be changed)
3. Encode interactions into a **fixed fingerprint structure**:
   * default: **11 bits per residue** (`IChem --basic IFP protein ligand`)
   * legacy: **7 bits per residue** using the subcommand (`IChem --old IFP protein ligand`)
4. Output:
   * an interaction table
   * residue header
   * fingerprint bit strings
   * similarity scores if using the reference the following commands:
        * (`IChem --basic IFP protein ligand reference_ligand `)
        * (`IChem --basic IFP protein1 ligand1 protein2 ligand2 `)


---

## Usage

The IFP is executed from the command line using one of the following syntaxes:

```text
IChem [options] IFP protein ligand
# Single complex analysis: Computes the interaction fingerprint of one protein-ligand complex

IChem [options] IFP protein ligand ligand_ref
# Comparison with a reference ligand: Computes the fingerprint and compares it to a reference ligand

IChem [options] IFP protein1 ligand1 protein2 ligand2
# Comparison of two complexes: Computes and compares the interaction fingerprints of two protein-ligand complexes
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