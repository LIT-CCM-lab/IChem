# volsite - Cavity detection and druggability prediction

`IChem volsite` (VolSite) automatically detects **surface cavities** on a protein structure (MOL2) and estimates their **ligandability** (structural druggability).
Optionally, it can generate **cavity-based pharmacophores** and export descriptors used by the VolSite SVM model.

VolSite can run in two modes:

- **Unrestricted mode**: detect all cavities on the protein surface (protein only)
- **Ligand-restricted mode**: detect cavities around a bound ligand (protein + ligand)

!!! important
    The VolSite druggability model is **only valid** for the standard VolSite parameters:  
    `step=1.5`, `boxS=20`, `b=55`, `n=5`.  
    If you change these parameters, the predicted druggability value should be considered **not meaningful**.

---

## Synopsis

### Unrestricted mode (protein only)

```bash
IChem volsite <protein.mol2>
```

### Ligand-restricted mode (protein + ligand)

```bash
IChem volsite <protein.mol2> <ligand.mol2>
```

---

## Inputs

- `protein.mol2`  
  Protein binding site or full protein surface representation (MOL2).

- `ligand.mol2` (optional)  
  Bound ligand (MOL2). This activates ligand-restricted mode.

!!! note
    In ligand-restricted workflows, the ligand may be a real ligand pose or even a **single atom** with user-defined coordinates
    if you only need to constrain cavity detection around a point of interest.

---

## Outputs

VolSite produces:

- One `CAVITY_...mol2` file **per detected cavity** (ranked from largest to smallest)
- A global summary CSV: **`VolSite_stat.csv`**
- Optional additional outputs depending on flags (`--desc`, `--svm`, `--pharm`)

### Cavity MOL2 files

#### Unrestricted mode
Files are named:

- `CAVITY_Nx_ALL.mol2`

Where:
- `x` is the cavity index (1 = largest cavity)

Example output set:
- `CAVITY_N1_ALL.mol2`, `CAVITY_N2_ALL.mol2`, …

#### Ligand-restricted mode
Files are named:

- `CAVITY_Nx_y.mol2`

Where:
- `x` is the cavity index (1 = largest cavity)
- `y` is a truncation distance (Å) defining which cavity points are kept based on proximity to ligand heavy atoms

Typical `y` values:
- `4`, `6`, `8`, `12`, and `ALL`

If `y = ALL`, no truncation is applied.

---

## How it works (concepts)

VolSite discretizes space using a regular grid and detects connected sets of grid points (“cavity points”) that satisfy
buriedness and neighborhood criteria. Cavity points are then annotated with pharmacophoric categories (CA, CZ, O, …)
based on nearby protein atoms.

Key internal concepts reflected in outputs:

- **Cavity point**: a grid point belonging to the cavity volume
- **Buriedness / projection value**: how enclosed a point is by the protein surface
- **Cavity composition**: counts of cavity point types (CA, CZ, O, …)
- **Druggability**: empirical estimate derived from cavity descriptors

---

## Options

### Geometry and detection parameters

| Option | Default | Unit | Description |
|---|---:|---|---|
| `-step <N>` | `1.5` | Å | Edge length of each grid box |
| `-boxS <N>` | `20` | Å | Edge length of the main box |
| `-b <N>` | `55` | % | Minimal threshold for buriedness |
| `-n <N>` | `5` | - | Minimal neighbors for buried cavity boxes |
| `-nPTS <N>` | `70` | cubes | Minimal number of cubes to consider a cavity |
| `-NPTS <N>` | `400` | cubes | Maximal number of cubes to consider a cavity |
| `-name <str>` | - | - | PDB name used to label output cavity names |
| `-drog <N>` | - | - | Observed druggability (label/metadata, when used) |

### Input composition

| Option | Default | Description |
|---|---:|---|
| `--dna` | - | Consider DNA as part of the protein |
| `--cofactor` | - | Consider cofactors as part of the protein |
| `--solvent` | - | Consider solvent molecules as part of the protein |
| `--hydrogen` | - | Consider hydrogens |

### Descriptors and ML exports

| Option | Default | Description |
|---|---:|---|
| `--desc` | - | Write the cavity descriptor file (e.g. `*_descriptor.txt`) |
| `--svm` | - | Build an SVM property file (used by the druggability model) |

### Pharmacophore generation

| Option | Default | Description |
|---|---:|---|
| `--pharm` | - | Generate a cavity-based pharmacophore (e.g. `.chm`) |
| `--outExclu` | - | Add exclusion spheres to the pharmacophore file |

!!! note
    Pharmacophore feature filtering can be controlled with `-p*` options (e.g. `-pCA`, `-pCZ`, `-pO`, `-pN`, `-pOD1`, `-pNZ`, `-pZn`),
    which apply PLP interaction energy thresholds to keep/remove features. Use values recommended by the VolSite methodology when possible.

---

## Interpreting `VolSite_stat.csv`

`VolSite_stat.csv` summarizes detected cavities and their properties.

### Unrestricted mode (protein only)

Typical columns include:

- `Name` : protein name (from MOL2 header by default)
- `NCav` : cavity number (1..N, largest to smallest)
- `Size` : usually `ALL` (no truncation)
- `Buriedness` : average buriedness (%)
- `NPts` : number of cavity points (grid points)
- `Volume` : cavity volume (Å³)
- `CA`, `CZ`, `O`, `OG`, `OD1`, `N`, `NZ` : counts of cavity point types
- `DU` : dummy points (no protein atom within 4.5 Å)
- `Drugg` : estimated druggability (druggable if `Drugg > 0`, undruggable if `Drugg < 0`)

Cavity point types correspond to pharmacophoric categories:

- `CA` : hydrophobic cavity points
- `CZ` : aromatic cavity points
- `O`  : hydrogen-bond acceptor points
- `OG` : hydrogen-bond acceptor/donor points
- `OD1`: negatively ionizable points
- `N`  : hydrogen-bond donor points
- `NZ` : positively ionizable points
- `DU` : dummy cavity points (no nearby protein atom)

### Ligand-restricted mode (protein + ligand)

In ligand-restricted mode, VolSite outputs multiple truncation sizes (`4`, `6`, `8`, `12`, `ALL`) for each cavity.
An additional column can appear:

- `Recovery` : percentage of the binding site enclosing the ligand

!!! note
    In ligand-restricted mode, druggability is typically reported only for the **non-truncated** cavities (`Size = ALL`).

---

## Descriptor output (`--desc`)

`--desc` writes a descriptor file (e.g. `2rh1_protein_descriptor.txt`) containing, for each cavity, a vector of **89 real values**
used by the VolSite SVM druggability predictor.

The descriptor includes:
- volume
- percentage of each cavity point type (CZ, CA, O, OG, OD1, N, NZ, DU)
- distributions of **projection values** per point type

### Projection value and buriedness

The *projection value* is the number of regularly spaced 8 Å vectors emitted from each cavity point that intercept the protein surface
(maximum 120). It relates to buriedness as follows:

| Projection value | Buriedness (%) |
|---:|---:|
| `< 40` | `< 33.3` |
| `40-50` | `33.3-41.6` |
| `50-60` | `41.6-50.0` |
| `60-70` | `50.0-58.3` |
| `70-80` | `58.3-66.6` |
| `80-90` | `66.6-75.0` |
| `90-100` | `75.0-83.3` |
| `100-110` | `83.3-91.6` |
| `110-120` | `91.6-99.9` |
| `120` | `100` |

---

## Pharmacophore generation (`--pharm`)

`--pharm` generates cavity-based pharmacophores ready to be used in external tools (e.g., BIOVIA, LigandScout, Shaper2).

Typical outputs include multiple files such as:
- `Pharmacophore.chm`
- `Pharmacophore.pml`
- `Pharmacophore.mol2`
- `Pharmacophore.plp` (summarizes PLP interaction energies, kcal/mol, per feature)

`--outExclu` adds exclusion spheres to the generated pharmacophores.

---

## Examples

### Example 1 - Detect all cavities (unrestricted)

```bash
IChem volsite protein.mol2
```

Typical outputs:
- `CAVITY_N1_ALL.mol2`, `CAVITY_N2_ALL.mol2`, …
- `VolSite_stat.csv`

### Example 2 - Detect cavities around a bound ligand (ligand-restricted)

```bash
IChem volsite protein.mol2 ligand.mol2
```

Typical outputs:
- `CAVITY_N1_4.mol2`, `CAVITY_N1_6.mol2`, `CAVITY_N1_8.mol2`, `CAVITY_N1_12.mol2`, `CAVITY_N1_ALL.mol2`
- `VolSite_stat.csv`

### Example 3 - Export descriptors for druggability model

```bash
IChem --desc volsite protein.mol2
```

Outputs include:
- `VolSite_stat.csv`
- `*_descriptor.txt` (one descriptor vector per cavity)

### Example 4 - Generate cavity-based pharmacophore

```bash
IChem --pharm volsite protein.mol2 ligand.mol2
```

Add exclusion spheres:

```bash
IChem --pharm --outExclu volsite protein.mol2 ligand.mol2
```

---

## Practical guidance

- For large, accessible cavities (e.g., protein-protein interfaces), standard settings may produce cavities that look too small.
  A common approach is to decrease buriedness threshold (e.g., `-b 45`) until cavity points fill the expected interface area.
- If you want accessory molecules to participate in cavity shape, use `--dna`, `--cofactor`, `--solvent`, and/or `--hydrogen`.
- If you plan to compare druggability values across targets, keep parameters fixed to the standard model settings.

---

## See also

- `AtomProps` - atom property annotation (`commands/atomprops.md`)
- `ints` - interaction pseudoatoms and generic fingerprints (`commands/ints.md`)
- `IFP` - interaction fingerprints (`commands/ifp.md`)
