import tempfile
import unittest
import subprocess
import sys
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

    def test_any_position_classification(self):
        df = pd.DataFrame(
            {
                "peptide_1_O": [
                    "ACPEPT", "PEPTAC", "PACGT", "APCGT", "APCGT",
                    "APCGT", "AC", "A", "AP[mod]C", "XAP", "ACG",
                ],
                "peptide_2_O": [
                    "CAPEPT", "PEPTCA", "PCAGT", "CPAGT", "CPTGA",
                    "CPCGT", "CA", "A", "CP[mod]A", "AXP", "CA",
                ],
            }
        )
        result = ind.annotate_terminal_swap(df, terminus="any", verbose=False)
        self.assertEqual(
            result["is_swap"].tolist(),
            [True, True, True, True, False, False, True, False, False, True, False],
        )
        self.assertEqual((result.loc[2, "swap_aa_row"], result.loc[2, "swap_aa_col"]), ("C", "A"))
        self.assertTrue(pd.isna(result.loc[9, "swap_aa_row"]))

    def test_position_modes_recalculate_and_aggregate(self):
        df = pd.DataFrame(
            {
                "peptide_1_O": ["ACPEPT", "PEPTAC", "PACGT", "APCGT"],
                "peptide_2_O": ["CAPEPT", "PEPTCA", "PCAGT", "CPAGT"],
                "score": [0.1, 0.3, 0.5, 0.9],
                "swap_aa_row": ["C"] * 4,
                "swap_aa_col": ["A"] * 4,
            },
            index=[3, 3, 5, 7],
        )
        original = df.copy(deep=True)
        for position, expected in (("N", 0.1), ("C", 0.3), ("any", 0.45)):
            with self.subTest(position=position):
                with patch("indistinguishability.plots.sns.heatmap", wraps=sns.heatmap) as heatmap:
                    ind.plot_aa_swap_heatmap(df, position=position, show_groups=False)
                self.assertAlmostEqual(heatmap.call_args.args[0].loc["C", "A"], expected)
        pd.testing.assert_frame_equal(df, original)

    def test_any_position_il_collapse_and_count(self):
        df = pd.DataFrame(
            {
                "peptide_1_O": ["IAPCG", "PLIQ", "PILQ"],
                "peptide_2_O": ["LCPAG", "PILQ", "PLIQ"],
                "score": [0.2, 0.4, 0.6],
            }
        )
        collapsed = ind.annotate_terminal_swap(df, terminus="any", verbose=False)
        separate = ind.annotate_terminal_swap(
            df, terminus="any", collapse_il=False, verbose=False
        )
        self.assertEqual(collapsed["is_swap"].tolist(), [True, False, False])
        self.assertEqual(separate["is_swap"].tolist(), [False, True, True])
        with patch("indistinguishability.plots.sns.heatmap", wraps=sns.heatmap) as heatmap:
            ind.plot_aa_swap_heatmap(
                df, position="any", collapse_il=False, aggfunc="count", show_groups=False
            )
        self.assertEqual(heatmap.call_args.args[0].loc["L", "I"], 2)

    def test_position_validation_and_legacy_annotations(self):
        df = pd.DataFrame({"score": [0.4], "swap_aa_row": ["C"], "swap_aa_col": ["A"]})
        with self.assertWarnsRegex(UserWarning, "cannot be verified"):
            ind.plot_aa_swap_heatmap(df, show_groups=False)
        for position in ("C", "any"):
            with self.subTest(position=position), self.assertRaisesRegex(KeyError, "peptide_1_O"):
                ind.plot_aa_swap_heatmap(df, position=position)
        with self.assertRaisesRegex(ValueError, "position must"):
            ind.plot_aa_swap_heatmap(df, position="invalid")
        self.assertIn("position", plotting.plot_options("aa_swap_heatmap"))

    def test_any_position_cli_and_output_name(self):
        df = pd.DataFrame(
            {"peptide_1_O": ["APCGT"], "peptide_2_O": ["CPAGT"], "score": [0.8]}
        )
        with tempfile.TemporaryDirectory() as output:
            plotting.make_plot(
                "aa_swap_heatmap", df, position="any", output_path=output, show_groups=False
            )
            self.assertTrue(Path(output, "any_position_aa_swap_mean_heatmap.png").is_file())

    def test_cli_position_option_end_to_end(self):
        df = pd.DataFrame(
            {"peptide 1": ["APCGT"], "peptide 2": ["CPAGT"], "score": [0.8]}
        )
        with tempfile.TemporaryDirectory() as output:
            scores = Path(output, "scores.tsv")
            df.to_csv(scores, sep="\t", index=False)
            result = subprocess.run(
                [
                    sys.executable, "plotting.py", "--scores", str(scores),
                    "--plot", "aa_swap_heatmap", "--opt", "position=any",
                    "--out", output,
                ],
                capture_output=True,
                text=True,
                check=True,
            )
            self.assertIn("1 two-residue swaps", result.stdout)
            self.assertTrue(Path(output, "any_position_aa_swap_mean_heatmap.png").is_file())


if __name__ == "__main__":
    unittest.main()
