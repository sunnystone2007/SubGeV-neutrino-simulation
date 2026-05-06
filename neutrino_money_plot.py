#!/usr/bin/env python3
"""Automate the ROOT money-plot peak/1σ extraction workflow.

This script mirrors the behaviour of ``peak68.C`` but runs fully in Python:

* finds ROOT files under ``1.5meter_5meter_reconstruction`` (configurable)
* loads the canvas that contains the three reconstruction curves
* determines the peak position of each curve
* applies the same binary-search interval finder used in ``peak68.C``
* optionally produces energy-vs-peak error-band graphs per flavor in a ROOT file

The output is a compact text table (and optional CSV) listing the peak and
68.27% containment interval for each histogram.
"""

from __future__ import annotations

import argparse
import fnmatch
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

import matplotlib.pyplot as plt

TARGET_FRAC = 0.682689492137  # Gaussian 1σ coverage

try:
    import ROOT  # type: ignore
except ImportError as exc:  # pragma: no cover - PyROOT is an external dep
    raise SystemExit(
        "PyROOT is required but not available. Install ROOT with Python bindings"
        " or source the corresponding environment before running this script."
    ) from exc


ROOT.gROOT.SetBatch(True)

DEFAULT_DATA_DIR = Path("1.5meter_reconstruction/antinumu_reconstruction_0326")
DEFAULT_PATTERN = "comparison_with_lepton_neutron_event*_0.5meter.root"

FLAVOR_TAGS = ("antinumu", "antinue", "nue", "numu")
GRAPH_COLORS = [
    ROOT.kBlue + 1,
    ROOT.kRed + 1,
    ROOT.kGreen + 2,
    ROOT.kMagenta + 2,
    ROOT.kOrange + 7,
]
FILL_STYLES = [3004, 3005, 3006, 3007, 3008]
YMAX_DEGREES = 150.0
MATPLOTLIB_COLORS = ('#1f77b4', '#d62728','#2ca02c','#9467bd', '#ff7f0e')


@dataclass
class IntervalResult:
    """Package the ROI information for reporting/CSV writing."""

    file_path: Path
    hist_name: str
    peak_x: float
    x_low: float
    x_high: float
    err_low: float
    err_high: float
    included_fraction: float
    flavor: str | None = None
    energy: float | None = None
    curve_key: str | None = None
    curve_label: str | None = None


def integrate_histogram(hist: "ROOT.TH1", x_min: float, x_max: float) -> float:
    """Integrate counts between ``x_min`` and ``x_max`` using bin overlaps."""

    if x_max <= x_min:
        return 0.0
    axis = hist.GetXaxis()
    nbins = hist.GetNbinsX()
    total = 0.0
    for ibin in range(1, nbins + 1):
        low = axis.GetBinLowEdge(ibin)
        high = axis.GetBinUpEdge(ibin)
        overlap = max(0.0, min(high, x_max) - max(low, x_min))
        width = high - low
        if overlap <= 0 or width <= 0:
            continue
        frac = overlap / width
        content = hist.GetBinContent(ibin)
        if content < 0.0:
            content = 0.0
        total += content * frac
    return total


def determine_flavor_tag(path: Path) -> str:
    for part in reversed(path.parts):
        lower = part.lower()
        for tag in FLAVOR_TAGS:
            if tag in lower:
                return tag
    parent = path.parent.name
    return parent if parent else path.stem


def extract_energy_from_name(path: Path) -> float | None:
    stem = path.stem
    patterns = (
        r"event([0-9]+(?:\.[0-9]+)?)",
        r"([0-9]+(?:\.[0-9]+)?)ge?v",
        r"reconstruction([0-9]+(?:\.[0-9]+)?)",
        r"_([0-9]+(?:\.[0-9]+)?)_",
    )
    for pattern in patterns:
        match = re.search(pattern, stem, re.IGNORECASE)
        if match:
            try:
                return float(match.group(1))
            except ValueError:
                continue
    # fall back to first floating-point number present
    floats = re.findall(r"[0-9]+\.[0-9]+", stem)
    if floats:
        try:
            return float(floats[0])
        except ValueError:
            pass
    return None


def sanitize_root_name(name: str) -> str:
    sanitized = re.sub(r"[^0-9A-Za-z_]+", "_", name)
    return sanitized.strip("_") or "graph"


def compute_interval_th1(
    hist: "ROOT.TH1", x_peak: float, frac: float = TARGET_FRAC
) -> tuple[float, float, float, float, float]:
    """Port of ``compute_interval_TH1`` from ``peak68.C``."""

    nbins = hist.GetNbinsX()
    if nbins <= 0:
        raise RuntimeError(f"Histogram {hist.GetName()} has no bins")

    centers: List[float] = []
    contents: List[float] = []
    for idx in range(1, nbins + 1):
        centers.append(hist.GetBinCenter(idx))
        content = hist.GetBinContent(idx)
        if content < 0.0:
            content = 0.0
        contents.append(content)
    axis = hist.GetXaxis()
    total_area = integrate_histogram(hist, axis.GetXmin(), axis.GetXmax())
    if total_area <= 0.0:
        raise RuntimeError(f"Histogram {hist.GetName()} is empty")

    target_area = frac * total_area
    peak_bin = axis.FindBin(x_peak) - 1
    peak_bin = max(0, min(nbins - 1, peak_bin))

    def find_left_cross(level: float) -> float:
        i = peak_bin
        while i > 0 and contents[i] >= level:
            i -= 1

        if contents[i] >= level:
            return axis.GetXmin()

        j = i + 1
        x1, y1 = centers[i], contents[i]
        x2, y2 = centers[j], contents[j]
        if abs(y2 - y1) < 1e-12:
            return x1
        t = (level - y1) / (y2 - y1)
        return x1 + t * (x2 - x1)

    def find_right_cross(level: float) -> float:
        i = peak_bin
        while i < nbins - 1 and contents[i] >= level:
            i += 1

        if contents[i] >= level:
            return axis.GetXmax()

        j = i - 1
        x1, y1 = centers[j], contents[j]
        x2, y2 = centers[i], contents[i]
        if abs(y2 - y1) < 1e-12:
            return x2
        t = (level - y1) / (y2 - y1)
        return x2 + t * (x1 - x2)

    max_height = max(contents)
    lo, hi = 0.0, max_height
    best_level = 0.0
    best_x_low = axis.GetXmin()
    best_x_high = axis.GetXmax()
    interval_area = 0.0

    for _ in range(48):  # binary search
        mid = 0.5 * (lo + hi)
        x_low = find_left_cross(mid)
        x_high = find_right_cross(mid)
        if x_high <= x_low:
            hi = mid
            continue

        area = integrate_histogram(hist, x_low, x_high)
        if area >= target_area:
            best_level = mid
            best_x_low = x_low
            best_x_high = x_high
            interval_area = area
            lo = mid
        else:
            hi = mid

    err_low = x_peak - best_x_low
    err_high = best_x_high - x_peak
    included_fraction = interval_area / total_area if total_area > 0 else 0.0

    return (
        float(best_x_low),
        float(best_x_high),
        float(err_low),
        float(err_high),
        float(included_fraction),
    )


def histogram_peak(hist: "ROOT.TH1") -> float:
    peak_bin = hist.GetMaximumBin()
    nbins = hist.GetNbinsX()
    if nbins <= 0:
        return hist.GetBinCenter(peak_bin)

    max_content = hist.GetBinContent(peak_bin)
    tol = max(abs(max_content) * 1e-9, 1e-12)
    tied_centers: List[float] = []
    for ibin in range(1, nbins + 1):
        content = hist.GetBinContent(ibin)
        if abs(content - max_content) <= tol:
            tied_centers.append(hist.GetBinCenter(ibin))

    if tied_centers:
        return min(tied_centers)

    return hist.GetBinCenter(peak_bin)


def extract_histograms(canvas: "ROOT.TCanvas") -> List["ROOT.TH1"]:
    primitives = canvas.GetListOfPrimitives()
    histograms: List["ROOT.TH1"] = []
    for obj in primitives:
        if not obj:
            continue
        if obj.InheritsFrom("TH1"):
            clone = obj.Clone()
            clone.SetDirectory(0)
            histograms.append(clone)
    return histograms


def pickup_histograms(
    hists: Sequence["ROOT.TH1"],
    patterns: Sequence[str] | None,
    max_count: int,
) -> List["ROOT.TH1"]:
    if patterns:
        selected = [h for h in hists if any(fnmatch.fnmatch(h.GetName(), p) for p in patterns)]
    else:
        selected = list(hists)
    if max_count > 0:
        return selected[:max_count]
    return selected


def _collect_histograms_in_directory(directory) -> List["ROOT.TH1"]:
    histograms: List["ROOT.TH1"] = []
    for key in directory.GetListOfKeys():
        obj = key.ReadObj()
        if not obj:
            continue
        if obj.InheritsFrom("TH1"):
            clone = obj.Clone()
            clone.SetDirectory(0)
            histograms.append(clone)
        elif obj.InheritsFrom("TCanvas"):
            histograms.extend(extract_histograms(obj))
        elif obj.InheritsFrom("TDirectory"):
            histograms.extend(_collect_histograms_in_directory(obj))
    return histograms


def analyze_file(
    file_path: Path,
    *,
    canvas_name: str,
    hist_patterns: Sequence[str] | None,
    max_hists: int,
) -> List[IntervalResult]:
    file = ROOT.TFile.Open(str(file_path))
    if not file or file.IsZombie():
        raise RuntimeError(f"Unable to open {file_path}")
    flavor_tag = determine_flavor_tag(file_path)
    file_energy = extract_energy_from_name(file_path)
    try:
        canvas = file.Get(canvas_name) if canvas_name else None
        if not canvas:
            for key in file.GetListOfKeys():
                obj = key.ReadObj()
                if obj.InheritsFrom("TCanvas"):
                    canvas = obj
                    break

        histograms: List["ROOT.TH1"] = []
        if canvas:
            histograms = extract_histograms(canvas)
        if not histograms:
            histograms = _collect_histograms_in_directory(file)
        if not histograms:
            raise RuntimeError(
                f"No histograms found in canvas {canvas_name} or anywhere in {file_path}"
            )

        selected = pickup_histograms(histograms, hist_patterns, max_hists)
        if not selected:
            raise RuntimeError(
                f"No histograms match the requested patterns {hist_patterns or ['<first>']} in {file_path}"
            )

        results: List[IntervalResult] = []
        for idx, hist in enumerate(selected):
            x_peak = histogram_peak(hist)
            x_low, x_high, err_low, err_high, fraction = compute_interval_th1(hist, x_peak)
            curve_key = f"curve{idx + 1}"
            label = hist.GetTitle() or hist.GetName()
            results.append(
                IntervalResult(
                    file_path=file_path,
                    hist_name=hist.GetName(),
                    peak_x=float(x_peak),
                    x_low=x_low,
                    x_high=x_high,
                    err_low=err_low,
                    err_high=err_high,
                    included_fraction=fraction,
                    flavor=flavor_tag,
                    energy=file_energy,
                    curve_key=curve_key,
                    curve_label=label,
                )
            )
        return results
    finally:
        file.Close()


def iter_root_files(paths: Sequence[Path], pattern: str) -> Iterable[Path]:
    for base in paths:
        if base.is_file():
            if fnmatch.fnmatch(base.name, pattern):
                yield base
            continue
        if not base.is_dir():
            raise FileNotFoundError(f"Input path does not exist: {base}")
        for candidate in sorted(base.rglob(pattern)):
            if candidate.is_file():
                yield candidate


def write_csv(path: Path, results: Sequence[IntervalResult]) -> None:
    header = (
        "file,flavor,curve,curve_label,energy,peak,x_low,x_high,err_low,err_high,fraction"
    )
    lines = [header]
    for res in results:
        energy = "" if res.energy is None else f"{res.energy:.6g}"
        lines.append(
            f"{res.file_path},{res.flavor or ''},{res.curve_key or ''},{res.curve_label or ''},{energy},"
            f"{res.peak_x:.6g},{res.x_low:.6g},{res.x_high:.6g},"
            f"{res.err_low:.6g},{res.err_high:.6g},{res.included_fraction:.6g}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def render_table(results: Sequence[IntervalResult]) -> str:
    headers = [
        "file",
        "flavor",
        "curve",
        "label",
        "energy(GeV)",
        "peak",
        "x_low",
        "x_high",
        "-err",
        "+err",
        "frac",
    ]
    rows = []
    for res in results:
        rel_path = res.file_path
        try:
            rel_path = res.file_path.relative_to(Path.cwd())
        except ValueError:
            pass
        energy = "" if res.energy is None else f"{res.energy:.3f}"
        rows.append(
            [
                str(rel_path),
                res.flavor or "",
                res.curve_key or "",
                res.curve_label or res.hist_name,
                energy,
                f"{res.peak_x:.4f}",
                f"{res.x_low:.4f}",
                f"{res.x_high:.4f}",
                f"{res.err_low:.4f}",
                f"{res.err_high:.4f}",
                f"{100.0 * res.included_fraction:5.2f}%",
            ]
        )
    widths: List[int] = []
    for idx, name in enumerate(headers):
        column_values = [name]
        column_values.extend(row[idx] for row in rows)
        widths.append(max(len(val) for val in column_values))

    header = " ".join(name.ljust(widths[idx]) for idx, name in enumerate(headers))
    sep = "-" * len(header)
    lines = [header, sep]
    for row in rows:
        lines.append(" ".join(row[idx].ljust(widths[idx]) for idx in range(len(headers))))
    return "\n".join(lines)


def group_results_for_bands(
    results: Sequence[IntervalResult],
) -> tuple[Dict[str, Dict[str, Dict[str, object]]], int]:
    grouped: Dict[str, Dict[str, Dict[str, object]]] = {}
    skipped = 0
    for res in results:
        if res.energy is None:
            skipped += 1
            continue
        flavor = res.flavor or determine_flavor_tag(res.file_path)
        curve_key = res.curve_key or res.hist_name
        hist_map = grouped.setdefault(flavor, {})
        entry = hist_map.setdefault(
            curve_key,
            {
                "label": res.curve_label or res.hist_name,
                "entry_map": {},
            },
        )
        entry_map: Dict[float, IntervalResult] = entry.setdefault("entry_map", {})  # type: ignore[assignment]
        energy = res.energy
        existing = entry_map.get(energy)
        if existing is None or res.peak_x < existing.peak_x:
            entry_map[energy] = res

    for hist_map in grouped.values():
        for payload in hist_map.values():
            entry_map = payload.pop("entry_map", {})
            payload["entries"] = sorted(
                entry_map.values(), key=lambda r: r.energy or 0.0
            )
    return grouped, skipped


def filter_by_flavor(
    results: Sequence[IntervalResult], flavor_filters: Sequence[str] | None
) -> List[IntervalResult]:
    if not flavor_filters:
        return list(results)
    filters = [f.lower() for f in flavor_filters]
    kept: List[IntervalResult] = []
    for res in results:
        flavor_value = (res.flavor or determine_flavor_tag(res.file_path)).lower()
        if any(f in flavor_value for f in filters):
            kept.append(res)
    return kept


def filter_by_energy(
    results: Sequence[IntervalResult], energy_min: float | None, energy_max: float | None
) -> List[IntervalResult]:
    if energy_min is None and energy_max is None:
        return list(results)
    kept: List[IntervalResult] = []
    for res in results:
        energy = res.energy
        if energy is None:
            continue
        if energy_min is not None and energy < energy_min:
            continue
        if energy_max is not None and energy > energy_max:
            continue
        kept.append(res)
    return kept


def plot_error_bands_matplotlib(
    flavor: str, hist_map: Dict[str, Dict[str, object]], output_path: Path
) -> None:
    fig, ax = plt.subplots(figsize=(8, 6))
    plotted = False
    for idx, (curve_key, payload) in enumerate(hist_map.items()):
        entries = payload.get("entries", [])
        entries = sorted(entries, key=lambda r: r.energy or 0.0)
        if not entries:
            continue
        energies = [res.energy or 0.0 for res in entries]
        peaks = [res.peak_x for res in entries]
        lower = [peak - res.err_low for peak, res in zip(peaks, entries)]
        upper = [peak + res.err_high for peak, res in zip(peaks, entries)]
        color = MATPLOTLIB_COLORS[idx % len(MATPLOTLIB_COLORS)]
        ax.plot(energies, peaks, color=color, marker='o')
        ax.fill_between(energies, lower, upper, color=color, alpha=0.3)
        plotted = True

    if not plotted:
        plt.close(fig)
        return

    ax.set_xlabel('Neutrino energy [GeV]')
    ax.set_ylabel('Peak angle [deg]')
    ax.set_title(f"{flavor} error bands")
    ax.set_ylim(0.0, YMAX_DEGREES)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
    print(f"Saved matplotlib error-band plot to {output_path}")


def write_error_band_root(path: Path, results: Sequence[IntervalResult]) -> None:
    grouped, skipped = group_results_for_bands(results)
    if not grouped:
        print("error-band: no entries with parsed energy available; skipping ROOT output")
        return

    root_file = ROOT.TFile(str(path), "RECREATE")
    for flavor, hist_map in grouped.items():
        flavor_dir = root_file.mkdir(sanitize_root_name(flavor))
        flavor_dir.cd()
        graphs = []
        for idx, (curve_key, payload) in enumerate(hist_map.items()):
            entries = sorted(payload["entries"], key=lambda r: r.energy or 0.0)
            label = payload["label"]
            graph_name = sanitize_root_name(f"{flavor}_{curve_key}_band") or f"band_{idx}"
            graph = ROOT.TGraphAsymmErrors(len(entries))
            for point_idx, res in enumerate(entries):
                energy = res.energy or 0.0
                graph.SetPoint(point_idx, energy, res.peak_x)
                graph.SetPointError(point_idx, 0.0, 0.0, res.err_low, res.err_high)

            color = GRAPH_COLORS[idx % len(GRAPH_COLORS)] if GRAPH_COLORS else ROOT.kBlue + 1
            fill_style = FILL_STYLES[idx % len(FILL_STYLES)] if FILL_STYLES else 3004
            graph.SetName(graph_name)
            graph.SetTitle(f"{flavor} {label};Neutrino energy [GeV];Peak angle [deg]")
            graph.SetLineColor(color)
            graph.SetMarkerColor(color)
            graph.SetLineWidth(2)
            graph.SetMarkerStyle(20)
            graph.SetFillColor(color)
            graph.SetFillStyle(fill_style)

            canvas = ROOT.TCanvas(f"c_{graph_name}", graph.GetTitle(), 800, 600)
            graph.GetYaxis().SetRangeUser(0.0, YMAX_DEGREES)
            graph.Draw("APE3")
            canvas.Write(f"{graph_name}_canvas")
            graph.Write(graph_name)
            graphs.append((graph, label))

        if graphs:
            combo_name = sanitize_root_name(f"{flavor}_combined_band")
            combo_canvas = ROOT.TCanvas(
                f"c_{combo_name}", f"{flavor} error bands", 900, 650
            )
            combo_canvas.cd()
            overlay_graphs = []
            for idx, (graph, label) in enumerate(graphs):
                overlay = graph.Clone(f"{graph.GetName()}_overlay_{idx}")
                overlay.SetFillColor(graph.GetFillColor())
                overlay.SetFillStyle(graph.GetFillStyle())
                overlay.GetYaxis().SetRangeUser(0.0, YMAX_DEGREES)
                draw_opt = "AP3" if idx == 0 else "3 SAME"
                overlay.Draw(draw_opt)
                overlay.Draw("P SAME")
                overlay_graphs.append(overlay)
            combo_canvas.Write(f"{combo_name}_canvas")

        png_path = path.parent / f"{sanitize_root_name(flavor)}_error_bands.png"
        plot_error_bands_matplotlib(flavor, hist_map, png_path)
        root_file.cd()
    root_file.Close()
    if skipped:
        print(f"Note: skipped {skipped} histogram entries without a parsable energy in the filename")
    print(f"Wrote error-band graphs to {path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Scan ROOT money plots, determine the peak of each curve, and compute"
            " the 68.27% containment interval using the peak68 algorithm."
        )
    )
    parser.add_argument(
        "--inputs",
        nargs="+",
        type=Path,
        default=[DEFAULT_DATA_DIR],
        help=(
            "Paths to search for ROOT files. Each path can be a directory (searched recursively)"
            f" or a specific ROOT file. Default: {DEFAULT_DATA_DIR}."
        ),
    )
    parser.add_argument(
        "--pattern",
        default=DEFAULT_PATTERN,
        help=(
            "Glob used to pick ROOT files under each --inputs entry"
            f" (default: {DEFAULT_PATTERN})."
        ),
    )
    parser.add_argument(
        "--canvas",
        default="first_reconstruction_result_on_neutrino",
        help="Canvas name that stores the comparison histograms",
    )
    parser.add_argument(
        "--hist",
        nargs="*",
        default=None,
        help="Optional fnmatch patterns for histogram selection. If omitted, take the first three curves",
    )
    parser.add_argument(
        "--max-hists",
        type=int,
        default=3,
        help="Maximum number of histograms to analyse per file (default: 3)",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=DEFAULT_DATA_DIR / "neutrino_money_plot.csv",
        help=(
            "Path to write the CSV summary. Set --no-csv to disable saving even though a default"
            f" path ({DEFAULT_DATA_DIR / 'neutrino_money_plot.csv'}) exists."
        ),
    )
    parser.add_argument(
        "--band-root",
        type=Path,
        default=DEFAULT_DATA_DIR / "neutrino_money_plot_error_bands_antinue_5meter.root",
        help=(
            "ROOT file to store error-band graphs (energy vs peak with 1σ asymmetric errors)."
            " Use --no-band-root to skip writing."
        ),
    )
    parser.add_argument(
        "--no-csv",
        action="store_true",
        help="Disable CSV output (overrides the default CSV path).",
    )
    parser.add_argument(
        "--no-band-root",
        action="store_true",
        help="Disable ROOT error-band output (overrides the default band-root path).",
    )
    parser.add_argument(
        "--flavor",
        nargs="*",
        default=None,
        help=(
            "Optional flavor filter (e.g. antinue). Only results matching these substrings"
            " will be kept and plotted."
        ),
    )
    parser.add_argument(
        "--energy-min",
        type=float,
        default=None,
        help="Minimum neutrino energy (GeV) to keep when aggregating results.",
    )
    parser.add_argument(
        "--energy-max",
        type=float,
        default=None,
        help="Maximum neutrino energy (GeV) to keep when aggregating results.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    csv_path = None if args.no_csv else args.csv
    band_path = None if args.no_band_root else args.band_root
    for path in (csv_path, band_path):
        if path:
            path.parent.mkdir(parents=True, exist_ok=True)

    files = list(iter_root_files(args.inputs, args.pattern))
    if not files:
        raise SystemExit(
            "No ROOT files found. Double-check --inputs paths and the --pattern glob."
        )

    skip_paths = {p.resolve() for p in (band_path,) if p is not None}
    filtered_files = [f for f in files if f.resolve() not in skip_paths]
    if not filtered_files:
        raise SystemExit(
            "Only output files (e.g. band-root) matched the search; adjust --pattern/--inputs."
        )

    all_results: List[IntervalResult] = []
    for file_path in filtered_files:
        results = analyze_file(
            file_path,
            canvas_name=args.canvas,
            hist_patterns=args.hist,
            max_hists=args.max_hists,
        )
        all_results.extend(results)

    flavor_results = filter_by_flavor(all_results, args.flavor)
    if not flavor_results:
        raise SystemExit("No histograms matched the specified flavor filter.")

    final_results = filter_by_energy(
        flavor_results, args.energy_min, args.energy_max
    )
    if not final_results:
        raise SystemExit("No histograms remained after applying the energy range filter.")

    print(render_table(final_results))

    if csv_path:
        write_csv(csv_path, final_results)
        print(f"Wrote CSV to {csv_path}")

    if band_path:
        write_error_band_root(band_path, final_results)


if __name__ == "__main__":
    main()
