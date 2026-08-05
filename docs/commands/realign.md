# realign - Rigid molecular alignment

`IChem realign` applies a rigid-body **rotation/translation** (rigid alignment) deduced from two molecules
to one or more additional molecules. It is commonly used to transfer the alignment computed on
interaction pseudoatoms (e.g., from `grim`) to the corresponding ligand and protein coordinates.

Typical use case:
- `grim` aligns interaction pseudoatoms and outputs an aligned pseudoatom set (`GRIM_ints.mol2`).
- `realign` transfers the same transformation to the original ligand/protein MOL2 files, producing aligned `rot_*.mol2` files.

---

## Synopsis

```bash
IChem realign <rigidM.mol2> <mobilM.mol2> <applied1.mol2> [applied2.mol2 ...]
```

Where:
- `rigidM`  : reference molecule (target orientation)
- `mobilM`  : comparison molecule (source orientation)
- `applied*`: molecule(s) to which the deduced rotation/translation is applied

---

## Inputs

- `rigidM.mol2`  
  Reference structure defining the target coordinate system.

- `mobilM.mol2`  
  Structure to be aligned onto `rigidM`. The transformation is deduced from the mapping between `rigidM` and `mobilM`.

- `applied*.mol2`  
  One or more MOL2 files that will receive the rotation/translation (e.g., original ligand and protein structures).

---

## Outputs

`realign` writes rotated/transformed MOL2 files, typically prefixed with `rot_`:

- `rot_<applied1>.mol2`
- `rot_<applied2>.mol2`
- ...

If `--wMob` is used, it also outputs the aligned `mobilM`.

Example output set (from test data):
- `rot_4amj_INTS_M.mol2`
- `rot_4amj_lig.mol2`
- `rot_4amj_prot.mol2`

---

## Options

### Matching strategy (`-gmatch`)

By default, `realign` performs an atom-by-atom match without enforcing atom typing constraints.
`-gmatch` restricts matching to specific criteria.

| Option | Default | Description |
|---|---:|---|
| `-gmatch <NAME|ATMN|MOL2|CALP>` | `NAME` | Matching criterion used to build the alignment mapping |
| `--wMob` | - | Also output the aligned `mobilM` |
| `-rule <R>` | - | Matching rule (build/version dependent) |
| `-i '<pairs>'` | - | Explicit atom index mapping string (see below) |

#### `-gmatch` modes

- `NAME` : atom name matching
- `ATMN` : atomic name matching
- `MOL2` : MOL2 type matching
- `CALP` : C-alpha atom matching (protein only)

!!! note
    The exact behavior of `NAME` vs `ATMN` may depend on how atom names are stored in your MOL2 files and the IChem build.
    If you need strict control, use an explicit index mapping with `-i`.

### Explicit atom index mapping (`-i`)

If you want to match only selected atoms (or enforce a custom pairing), use the index string.

Format:
- `'<refIndex>-<mobIndex>|<refIndex>-<mobIndex>|...'`

Example:

```bash
IChem -i '2-3|1-6|23-160' realign rigid.mol2 mobile.mol2 applied.mol2
```

Meaning:
- match atom 2 (reference) with atom 3 (mobile)
- match atom 1 (reference) with atom 6 (mobile)
- match atom 23 (reference) with atom 160 (mobile)

This is useful when:
- atom ordering differs between files
- you only want to align on a subset (e.g., common scaffold atoms)

---

## Example workflow: transfer GRIM alignment to ligand and protein

Input directory: `$ICHEM_DIR/test/REALIGN`  
Input files:
- `GRIM_ints.mol2` (aligned pseudoatoms from GRIM)
- `4amj_INTS_M.mol2` (original pseudoatoms of complex to fit)
- `4amj_lig.mol2`, `4amj_prot.mol2` (original ligand/protein coordinates)

Command:

```bash
IChem --wMob realign GRIM_ints.mol2 4amj_INTS_M.mol2 4amj_lig.mol2 4amj_prot.mol2
```

What happens:
- The transform aligning `4amj_INTS_M.mol2` onto `GRIM_ints.mol2` is computed.
- The same transform is applied to `4amj_lig.mol2` and `4amj_prot.mol2`.
- Outputs `rot_4amj_lig.mol2` and `rot_4amj_prot.mol2` contain the aligned complex.

---

## Practical guidance

- Use `--wMob` when you also want the aligned version of the `mobilM` file for debugging/visual inspection.
- If alignment fails or looks incorrect:
  - verify that `rigidM` and `mobilM` represent the *same kind* of object (e.g., both are merged interaction pseudoatoms)
  - try a stricter `-gmatch` criterion (`MOL2` or `CALP` for proteins)
  - use `-i` to enforce a controlled mapping on a known subset

---

## See also

- `grim` - interaction graph matching and alignment (`commands/grim.md`)
- `ints` - interaction pseudoatoms generation (`commands/ints.md`)
