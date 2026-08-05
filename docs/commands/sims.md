# sims - Fingerprint similarity

`IChem sims` computes similarity (or distance) scores between **interaction fingerprints** using a selection of common metrics.
It supports pairwise comparison and “matrix mode” (all-against-all) depending on the input form.

Fingerprints of different textual formats are accepted, including:

- `STD` : standard dense format (e.g., `1 0 21 0 0 3`)
- `SVM` : sparse `POSITION:VALUE` format (e.g., `1:1 3:21 6:3`)
- `CMP` : compressed format with zero runs (e.g., `1 [1 21 [2 3`)

---

## Synopsis

`sims` supports three modes.

### Mode 1 - Compare two fingerprints

```bash
IChem sims <ref_fp> <comp_fp>
```

### Mode 2 - Matrix mode (all-against-all)

```bash
IChem sims <file>
```

Where `<file>` contains multiple fingerprints concatenated (e.g., one per line or one per block depending on format/export).

### Mode 3 - Compare two complexes / interaction sets

```bash
IChem sims <RefInt> <CompInt>
```

This mode can operate on fingerprints or on interaction listings, depending on options (`--wInts`, `--small`).

---

## Options

| Option | Default | Description |
|---|---:|---|
| `--wInts` | - | Use interactions instead of fingerprints (mode 3 only) |
| `--small` | - | Use small fingerprint (mode 3 only) |
| `--binary` | - | Set when the fingerprint is binary |
| `-metric <TC|HM|RT|FT|DI|SO>` | `TC` | Similarity metric selection |

### Metrics

| Code | Metric | Notes |
|---|---|---|
| `TC` | Tanimoto | Default; common for binary fingerprints |
| `HM` | Hamming distance | Distance (lower is more similar) |
| `RT` | Reference Tversky | Asymmetric similarity |
| `FT` | Fit Tversky | Asymmetric similarity |
| `DI` | Dice | Similarity |
| `SO` | Soergel | Distance-like metric |

!!! note
    Some metrics are distances (e.g., Hamming) while others are similarities.
    Always interpret scores in the context of the chosen metric.

---

## Inputs

### Mode 1 (two fingerprints)
- `<ref_fp>`: reference fingerprint file
- `<comp_fp>`: comparison fingerprint file

### Mode 2 (matrix mode)
- `<file>`: a file containing a collection of fingerprints to compare (all-against-all)

### Mode 3 (two complexes / interaction sets)
- `<RefInt>` and `<CompInt>`: reference and comparison inputs used for comparison (build/workflow dependent)
  - if `--wInts` is enabled: use **interaction listings** (mode 3 only)
  - otherwise: use fingerprints, optionally `--small`

---

## Outputs

`sims` writes a table to standard output with:

- reference identifier
- comparison identifier
- similarity score (or distance, depending on metric)

Typical pairwise output:

```text
FP1     FP2     0.809524
```

In matrix mode, it outputs all pairs (including self-comparisons), e.g.:

```text
FP1     FP1     1
FP1     FP2     0.809524
FP1     FP3     0.652174
FP2     FP1     0.809524
...
```

---

## Examples

### Example 1 - Pairwise similarity between two fingerprints

```bash
IChem --binary sims FP1.txt FP2.txt > sim.txt
```

Output (`sim.txt`):

```text
FP1     FP2     0.809524
```

### Example 2 - Full similarity matrix (all-against-all)

```bash
IChem --binary sims FP.txt > matrix.txt
```

Where `FP.txt` contains the fingerprints to compare (e.g., FP1, FP2, FP3).

### Example 3 - Mode 3 with small fingerprints

```bash
IChem --small sims RefInt CompInt
```

### Example 4 - Mode 3 using interaction listings instead of fingerprints

```bash
IChem --wInts sims RefInt CompInt
```

---

## Practical guidance

- Use `--binary` when your fingerprints represent presence/absence (0/1). This is common for IFP-like bitstrings.
- Prefer `-metric TC` (Tanimoto) for binary fingerprints and virtual screening similarity.
- Use matrix mode when you need clustering or to build pairwise similarity networks.
- If your workflow uses `ints` small fingerprints, enable `--small` to match the generation mode.

---

## See also

- `IFP` - site-dependent interaction fingerprints (`commands/ifp.md`)
- `ints` - generic interaction fingerprints (`commands/ints.md`)
