# Python API - Configuration (IFPConfig)

`IFPConfig` controls:
- the interaction **profile**
- geometric **thresholds** (min/max distances, angle)
- additional filters and legacy options

## Create a config

```python
import ichem_ifp
cfg = ichem_ifp.IFPConfig()
```

### Profiles
- basic - hydrophobic, H-bond, ionic, aromatic
- old - legacy profile (7-bit)
- metal - metal only
- weakh - weak H-bonds only
- picat - pi-cation only

#### Example:

```python
cfg = ichem_ifp.IFPConfig()
cfg.basic = True
```
- Important: Combining profiles is not allowed, only one of the profile must be selected

### Distance thresholds (Å)

- Convention:

    - d_* = minimum distance
    - D_* = maximum distance

    - d_Hb / D_Hb - H-bond
    - d_Hyd / D_Hyd - hydrophobic
    - d_Io / D_Io - ionic
    - d_Ar / D_Ar - aromatic
    - d_Me / D_Me - metal
    - d_WHb / D_WHb - weak H-bond
    - d_Pic / D_Pic - pi-cation

#### Example:

```python 
cfg = ichem_ifp.IFPConfig()
cfg.basic = True
cfg.D_Hyd = 4.5
```

- Angles (°) and tolerances (°):
    - a_H, at_H - H-bond
    - a_ArFF, at_ArFF - aromatic face-face
    - a_ArEF, at_ArEF - aromatic edge-face
    - a_Pic, at_Pic - pi-cation

```python
cfg = ichem_ifp.IFPConfig()
cfg.basic = True
cfg.a_H = 120
cfg.at_H = 30
```

### More options

The following options are enabled by default:
- `includeSolvent` - include solvent molecules in interaction analysis
- `includeCofactor` - include cofactors in interaction analysis
- `oldHydrophobic` - use the legacy hydrophobic interaction model

- If `oldHydrophobic` option is set to false, the interaction is accepted only if the protein atom is surrounded within 4.5 Å by non-hydrogen atoms of which more than 50% are hydrophobic


#### Example:
```python
cfg = ichem_ifp.IFPConfig()
cfg.basic = True
cfg.includeSolvent = False
cfg.includeCofactor = False
cfg.oldHydrophobic = False
```

## API references

This document describes the available functions for computing interaction fingerprints, interactions, and similarity scores

### `compute_ifp_interactions`

```python
compute_ifp_interactions(
    protein_file: str,
    ligand_file: str,
    cfg: IFPConfig = IFPConfig()
) -> list[LigandInteractions]
```

#### Description
Computes detailed interaction information between a protein and a ligand structure

#### Parameters

| Name | Type | Description |
|------|------|------------|
| `protein_file` | `str` | Path to the protein structure file |
| `ligand_file` | `str` | Path to the ligand structure file |
| `cfg` | `IFPConfig` | Configuration object |

#### Returns

- `list[LigandInteractions]` - A list of computed interaction objects

---

### `compute_ifp_fingerprints`

```python
compute_ifp_fingerprints(
    protein_file: str,
    ligand_file: str,
    cfg: IFPConfig = IFPConfig()
) -> list[LigandFingerprint]
```
#### Description
Generates interaction fingerprints for a protein-ligand complex

#### Parameters

| Name | Type | Description |
|------|------|------------|
| `protein_file` | `str` | Path to the protein structure file |
| `ligand_file` | `str` | Path to the ligand structure file |
| `cfg` | `IFPConfig` | Configuration object |

#### Returns

- `list[LigandFingerprint]` - A list of computed fingerprint objects

---

### `compute_ifp_tanimoto`

```python
compute_ifp_tanimoto(
    protein_file: str,
    ligand_file: str,
    reference_file: str,
    cfg: IFPConfig = IFPConfig()
) -> list[TanimotoScore]
```

#### Description
Calculates the Tanimoto similarity score between a ligand and a reference structure

#### Parameters

| Name | Type | Description |
|------|------|------------|
| `protein_file` | `str` | Path to the protein structure file |
| `ligand_file` | `str` | Path to the ligand structure file |
| `reference_file` | `str` | Path to the reference ligand structure |
| `cfg` | `IFPConfig` | Configuration object |

#### Returns

- `list[TanimotoScore]` - A list of similarity score objects

---


### `compute_ifp_tanimoto_ensembles`

```python
compute_ifp_tanimoto_ensembles(
    protein1_file: str,
    ligand1_file: str,
    protein2_file: str,
    ligand2_file: str,
    cfg: IFPConfig = IFPConfig()
) -> list[TanimotoScore]
```

#### Description
Computes the Tanimoto similarity score between two protein-ligand ensembles

#### Parameters

| Name | Type | Description |
|------|------|------------|
| `protein1_file` | `str` | First protein structure file |
| `ligand1_file` | `str` | First ligand structure file |
| `protein2_file` | `str` | Second protein structure file |
| `ligand2_file` | `str` | Second ligand structure file |
| `cfg` | `IFPConfig` | Configuration object |

#### Returns

- `list[TanimotoScore]` - A list of ensemble similarity score objects

---


### Notes


!!! note "File inputs"
    All functions accept file paths as input. Files should be in MOL2

!!! note "Configuration"
    The `IFPConfig` object allows customization of interaction detection and fingerprint parameters.
    A profile must be selected in `IFPConfig` **(basic, old, metal, weakh, picat)**. Exactly one profile is allowed. Profiles cannot be combined
