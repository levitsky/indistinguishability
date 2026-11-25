"""Koina spectrum prediction + spectrum_utils mirror plots.

``koinapy`` and ``spectrum_utils`` are optional (install via the ``spectra`` extra) and
imported lazily, so the rest of the package works without them. The Koina model is created
once and cached; override the default via ``get_model(...)`` or pass ``model=`` explicitly.
"""

import numpy as np
import pandas as pd

from .style import PlotStyle, get_style, new_axes

DEFAULT_MODEL_NAME = "Prosit_2025_intensity_40PTM"
DEFAULT_KOINA_SERVER = "koina.wilhelmlab.org:443"

#: Fragment-ion colours for the mirror plot (``y``, ``b``, unannotated).
ION_COLORS = {"y": "#F33B16", "b": "#022873", "?": "#BEC7DA"}

_MODEL_CACHE = {}


def get_model(model_name=DEFAULT_MODEL_NAME, server=DEFAULT_KOINA_SERVER):
    """Return a cached ``koinapy.Koina`` client, creating it on first use."""
    from koinapy import Koina

    key = (model_name, server)
    if key not in _MODEL_CACHE:
        _MODEL_CACHE[key] = Koina(model_name, server)
    return _MODEL_CACHE[key]


def predict_spectrum(peptide: str,
                     charge: int = 2,
                     ce: int = 30,
                     fragmentation: str = "HCD",
                     model=None) -> dict:
    """
    Predict MS/MS spectrum using Koina (koinapy client).
    Returns dict with arrays: mz[], intensity[], annotation[], precursor_mz.
    """
    if model is None:
        model = get_model()
    inputs = pd.DataFrame({
        "peptide_sequences": [peptide],
        "precursor_charges": [charge],
        "collision_energies": [ce],
        "fragmentation_types": [fragmentation]
    })
    pred = model.predict(inputs)
    mz = pred["mz"].to_numpy()
    intensity = pred["intensities"].to_numpy()
    annotation = pred["annotation"].to_numpy()
    # Koina does not return the precursor m/z; spectrum_utils only needs a placeholder.
    precursor_mz = pred["mz"].iloc[0]
    return {
        "mz": mz,
        "intensity": intensity,
        "annotation": annotation,
        "precursor_mz": precursor_mz
    }


def koina_to_spectrum_utils(pred: dict, peptide: str, charge: int):
    from spectrum_utils.spectrum import MsmsSpectrum

    spec = MsmsSpectrum(
        identifier=peptide,
        precursor_mz=pred["precursor_mz"],
        precursor_charge=charge,
        mz=np.array(pred["mz"]),
        intensity=np.array(pred["intensity"])
    )
    spec.annotate_proforma(
        proforma_str=peptide,
        fragment_tol_mass=0.02,
        fragment_tol_mode="Da",
        ion_types="by",
        max_ion_charge=charge,
    )
    # A peak can match several fragments within the tolerance (e.g. y2 and y4^2 at
    # 272.17); the plot uses the first, so put the closest mass match first.
    for interp in spec.annotation:
        interp.fragment_annotations.sort(
            key=lambda a: abs(a.mz_delta[0]) if a.mz_delta else 0
        )
    return spec


def _ion_label(annotation):
    """b/y label with charge as '+' marks (``y2+``, ``y6++``); blank for isotopes, neutral
    losses and unknowns. spectrum_utils' default only labels charge 1, leaving
    higher-charge peaks coloured but unlabelled."""
    if (
        annotation.ion_type[0] not in "by"
        or annotation.neutral_loss is not None
        or annotation.isotope != 0
    ):
        return ""
    return f"{annotation.ion_type}{'+' * annotation.charge}"


def plot_mirror(peptide_top, peptide_bottom, charge=2, ce=30, *, style=None, ax=None,
                output_path=None, name="mirror_plot", title=None, model=None,
                ion_colors=None, figsize=None):
    """Mirror plot of the Koina-predicted spectra of two peptides (ProForma strings, e.g.
    ``M[UNIMOD:35]``). Returns ``(fig, ax)``; saved only when `output_path` is given."""
    import spectrum_utils.plot as sup

    s = style if isinstance(style, PlotStyle) else get_style(style or "notebook")
    pred1 = predict_spectrum(peptide_top, charge, ce, model=model)
    pred2 = predict_spectrum(peptide_bottom, charge, ce, model=model)
    spec1 = koina_to_spectrum_utils(pred1, peptide_top, charge)
    spec2 = koina_to_spectrum_utils(pred2, peptide_bottom, charge)
    fig, ax = new_axes(s, ax, figsize)
    sup.colors.update(ion_colors or ION_COLORS)
    sup.mirror(spec1, spec2, spectrum_kws={"annot_fmt": _ion_label}, ax=ax)
    ax.spines[["right", "top"]].set_visible(False)
    if s.show_title:
        ax.set_title(title or f"{peptide_top} vs {peptide_bottom} ({charge}+)",
                     fontsize=s.title_size)
    if output_path is not None:
        print(f"saved {s.save(fig, output_path, name)}")
    return fig, ax