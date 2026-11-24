import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

import indistinguishability as ind
import plotting


class TerminalSwapHeatmapTests(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def test_terminal_annotation_and_existing_n_terminal_default(self):
        df = pd.DataFrame(
            {
                "peptide_1_O": ["PEPTAC", "ACPEPT", "PEPTAA", "PEPTAC", "PEPTA[mod]C"],
                "peptide_2_O": ["PEPTCA", "CAPEPT", "PEPTAA", "PEPSCA", "PEPTCA[mod]"],
            }
        )
        c = ind.annotate_terminal_swap(df, terminus="C", verbose=False)
        n = ind.annotate_terminal_swap(df, verbose=False)
        self.assertEqual(c["is_swap"].tolist(), [True, False, False, False, False])
        self.assertEqual(n["is_swap"].tolist(), [False, True, False, False, False])
        self.assertEqual((c.loc[0, "swap_aa_row"], c.loc[0, "swap_aa_col"]), ("C", "A"))
        with self.assertRaisesRegex(ValueError, "terminus must"):
            ind.annotate_terminal_swap(df, terminus="X", verbose=False)

    def test_c_terminal_heatmap_recalculates_annotations_and_means(self):
        df = pd.DataFrame(
            {
                "peptide_1_O": ["PEPTAC", "PEPTCA", "ACPEPT"],
                "peptide_2_O": ["PEPTCA", "PEPTAC", "CAPEPT"],
                "score": [0.2, 0.8, 0.99],
                "swap_aa_row": [None, None, "C"],
                "swap_aa_col": [None, None, "A"],
            }
        )
        original = df.copy(deep=True)
        with patch("indistinguishability.plots.sns.heatmap", wraps=sns.heatmap) as heatmap:
            fig, ax = ind.plot_c_terminal_aa_swap_heatmap(df, show_groups=False)
        self.assertIs(ax.figure, fig)
        self.assertAlmostEqual(heatmap.call_args.args[0].loc["C", "A"], 0.5)
        pd.testing.assert_frame_equal(df, original)

    def test_count_output_and_registry(self):
        df = pd.DataFrame(
            {"peptide_1_O": ["PEPTAC", "PEPTCA"],
             "peptide_2_O": ["PEPTCA", "PEPTAC"], "score": [0.2, 0.8]}
        )
        with tempfile.TemporaryDirectory() as output:
            with patch("indistinguishability.plots.sns.heatmap", wraps=sns.heatmap) as heatmap:
                plotting.make_plot(
                    "c_terminal_aa_swap_heatmap", df, aggfunc="count",
                    output_path=output, show_groups=False
                )
            self.assertEqual(heatmap.call_args.args[0].loc["C", "A"], 2)
            self.assertTrue(Path(output, "c_terminal_aa_swap_count_heatmap.png").is_file())
        self.assertIn("collapse_il", plotting.plot_options("c_terminal_aa_swap_heatmap"))

    def test_il_collapsing_and_missing_swaps(self):
        df = pd.DataFrame(
            {"peptide_1_O": ["IPAC", "PEPTIL"],
             "peptide_2_O": ["LPCA", "PEPTLI"], "score": [0.3, 0.9]}
        )
        collapsed = ind.annotate_terminal_swap(df, terminus="C", verbose=False)
        separate = ind.annotate_terminal_swap(
            df, terminus="C", collapse_il=False, verbose=False
        )
        self.assertEqual(collapsed["is_swap"].tolist(), [True, False])
        self.assertEqual(separate["is_swap"].tolist(), [False, True])
        with patch("indistinguishability.plots.sns.heatmap", wraps=sns.heatmap) as heatmap:
            ind.plot_c_terminal_aa_swap_heatmap(df, collapse_il=False, show_groups=False)
        self.assertAlmostEqual(heatmap.call_args.args[0].loc["L", "I"], 0.9)
        with self.assertRaisesRegex(ValueError, "no last-two-residue swap"):
            ind.plot_c_terminal_aa_swap_heatmap(df.iloc[[1]])

    def test_normalisation_before_annotation(self):
        df = pd.DataFrame(
            {"peptide 1": ["PEPTM[UNIMOD:35]C[UNIMOD:4]"],
             "peptide 2": ["PEPTC[UNIMOD:4]M[UNIMOD:35]"], "score": [0.7]}
        )
        for n in (1, 2):
            df[f"peptide_{n}_O"] = ind.normalise_mods(df[f"peptide {n}"])
        result = ind.annotate_terminal_swap(df, terminus="C", verbose=False)
        self.assertTrue(result.loc[0, "is_swap"])
        with patch("indistinguishability.plots.sns.heatmap", wraps=sns.heatmap) as heatmap:
            ind.plot_c_terminal_aa_swap_heatmap(df, show_groups=False)
        self.assertAlmostEqual(heatmap.call_args.args[0].loc["C", "O"], 0.7)


if __name__ == "__main__":
    unittest.main()
