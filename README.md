## Input

It reads from two output formats of the similarity tool:

- **`scores.tsv`**: the verbose format, with columns `peptide 1`, `peptide 2`, `score`,
  and optionally `charge 1/2`, `m/z 1/2`, `iRT 1/2`.
- **`peptides.tsv` + `scores.npy`**: compact `(i, j, score)` rows that index into the
  peptide table. Add `--filtered` if `peptides.tsv` is a filtered subset whose original row
  numbers are in an `Unnamed: 0` column.

Both are loaded into the same set of columns. Each pair also gets its sequences in a
normalised form (oxidised Met as `O`), the Levenshtein distance, the peptide lengths, the
proline position and the type of amino-acid substitution.

## Command line

```bash
python plotting.py --list                                   # all plots and their options

python plotting.py --scores scores.tsv --plot score_distribution --out figures

python plotting.py --peptides peptides.tsv --npy scores.npy \
    --plot aa_swap_heatmap --opt aggfunc=count \
    --preset talk --figsize 12 9 --cmap-count magma --format pdf --out figures

python plotting.py --scores scores.tsv --threshold 0.7 --remove MO_IL \
    --plot all --opt vmin=0 --opt vmax=1 --cache scores.pkl
```

| Flag | Purpose |
|---|---|
| `--plot NAME` | plot to make; repeat it for several, or use `all` (all except `mirror`) |
| `--opt key=value` | plot-specific option, e.g. `aa=M`, `sub_type=A/S`, `xlim=0,30`, `vmin=0` |
| `--threshold`, `--remove {MO_IL,MO,IL}`, `--charge` | filter the pairs before plotting |
| `--preset {paper,talk,notebook}` | style preset |
| `--figsize W H`, `--dpi`, `--format`, `--color`, `--cmap`, `--cmap-count`, `--palette`, `--font-scale`, `--no-title`, `--transparent` | override parts of the preset |
| `--cache file.pkl` | save the loaded data so later runs skip the slow loading step |

**Plots**

| Group | Plots |
|---|---|
| Distributions | `score_distribution`, `levenshtein_distribution`, `irt_diff_histogram`, `peptide_length` |
| Score vs. features | `levenshtein_violin`, `feature_hexbin`, `binned_mean`, `score_by_combo`, `irt_ccs_hexbin` (needs CCS) |
| Heatmaps | `heatmap_mean`, `heatmap_count`, `aa_position_heatmap`, `aa_swap_heatmap`, `c_terminal_aa_swap_heatmap`, `swap_position_heatmap`, `ox_sub_distance` |
| Other | `aa_position_violin`, `swap_aa_vs_input`, `mirror` (Koina; takes `pair=PEP1,PEP2` or `top=N`) |

## Python / notebooks

Every plot function returns `(fig, ax)` and saves a file only when you pass
`output_path`. Because the figure is returned, you can keep editing it afterwards.

```python
import indistinguishability as ind

df = ind.load_scores(scores_tsv="scores.tsv")        # or peptides=..., npy=...
style = ind.get_style("paper", cmap_mean="magma")    # preset + overrides
fig, ax = ind.plot_heatmap(df[df.score >= 0.7], rows="first_p_1", vmin=0.7, vmax=1,
                           style=style, output_path="figures")
```

For a residue-by-residue swap heatmap, `plot_aa_swap_heatmap` accepts
`position="N"` (default, first two residues), `"C"` (last two residues), or
`"any"` (two swapped residues anywhere, including nonadjacent positions). A pair
must otherwise match. It accepts `aggfunc="count"` for pair counts or the default
`"mean"` for mean scores. On the command line, use for example
`--plot aa_swap_heatmap --opt position=C`.

The plot expects normalized `peptide_1_O` and `peptide_2_O` columns. If you already
have a DataFrame with `peptide 1` and `peptide 2`, prepare a copy as follows:

```python
plot_df = df.copy()
for n in (1, 2):
    plot_df[f"peptide_{n}_O"] = ind.normalise_mods(plot_df[f"peptide {n}"])
fig, ax = ind.plot_aa_swap_heatmap(plot_df, position="C")
```

By default I and L are treated as equivalent across the entire peptide; use
`collapse_il=False` to distinguish them. Normalize supported modification tags before
plotting; pairs with unrecognized bracketed modifications cannot be classified.
Swap annotations are recalculated when sequences are present. For compatibility,
N-terminal precomputed `swap_aa_row` and `swap_aa_col` can still be plotted without
sequences, but their provenance and I/L settings cannot be verified.
`plot_c_terminal_aa_swap_heatmap` and its CLI name remain as compatibility aliases.

## Layout

- `plotting.py` is the command-line tool and plot registry (`PLOTS`, `make_plot`).
- `src/indistinguishability/` contains the modules:
  - `loading.py`: loaders (`load_scores`) and filters
  - `plots.py` and `spectra.py`: plots, including the Koina mirror plots in `spectra.py`
  - `style.py`: styling (`PlotStyle`, `get_style`, presets)
  - `sequence.py`, `substitutions.py`, `rearrangements.py`: sequence and substitution annotation
  - `reporting.py`: printed summaries
