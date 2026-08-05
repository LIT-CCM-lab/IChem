# AtomProps - Atom property annotation

`IChem AtomProps` reads one or more **MOL2** files and prints, for each atom:

- the **atom identity** (name, index, residue name/id, chain)
- the **atom type** (MOL2/SYBYL/Tripos type as interpreted by IChem)
- the **IChem atom properties** used by interaction detection (e.g., acceptor, hydrophobic, aromatic, cation…)
- the **atom charge**

This command is mainly used to **inspect** how IChem interprets a molecule (ligand or protein fragment) before running other modules such as `IFP` or `ints`.

---

## Synopsis

```bash
IChem AtomProps <mol.mol2>
IChem AtomProps <mol1.mol2> <mol2.mol2> ...
```

!!! note
    The command name is case-sensitive in some environments. Use `IChem AtomProps ...` as shown above.

---

## Inputs

- `<mol*.mol2>`  
  One or more MOL2 files (ligands, cofactors, small molecules, or any MOL2 you want to inspect).

---

## Output

`AtomProps` writes its report to **standard output**.

For each input file, it prints a header identifying the molecule and then a block per atom.

### Example output (single file)

```text
./IChem AtomProps ligand1.mol2
Molecule 1 with name REF from file: ligand1.mol2

ATOM    O17  |1    |CAU   0-XX
         type : O.3
         props : Acceptor
         charge: 0

ATOM    C16  |2    |CAU   0-XX
         type : C.3
         props : Hydrophobic|Apolar
         charge: 0
...
```

### Atom block format

Each atom is reported using the following structure:

- `ATOM <atom_name> |<atom_index>|<residue> <resid>-<chain>`
- `type : <atom_type>`
- `props : <prop_1>|<prop_2>|...` (empty if no property assigned)
- `charge: <value>`

---

## IChem atom properties

The `props` field is the set of chemical labels assigned by IChem and used by interaction detection.

From typical outputs, you may see properties such as:

- `Acceptor`
- `Hydrophobic`
- `Apolar`
- `Aromatic`
- `Cation`

!!! note
    The property vocabulary depends on the IChem version and the internal typing rules.  
    Some atoms (often hydrogens or certain nitrogens) may show an empty `props` field.

---

## Using multiple input files

You can pass multiple MOL2 files in a single command. `AtomProps` prints one section per molecule, in the order provided.

Example:

```bash
IChem AtomProps ligand1.mol2 ligand2.mol2 cofactor.mol2
```

This is useful to compare how IChem assigns properties across several ligands or to validate a set of inputs before running screening workflows.

---

## Practical guidance

- If you see unexpected `props` (e.g., a missing `Acceptor`), first verify that:
  - the MOL2 atom types are correct
  - the molecule contains explicit atoms as expected (e.g., protonation state)
- If `charge` is always `0`, your MOL2 may not contain partial charges or IChem may be configured to ignore them for this report.
- An empty `props` field is not necessarily an error; it can be normal for atoms that do not match a property class.

---

## See also

- `ints` - interaction pseudoatoms and generic fingerprints (`commands/ints.md`)
- `IFP` - site-dependent interaction fingerprints (`commands/ifp.md`)
