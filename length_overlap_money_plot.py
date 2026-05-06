#!/usr/bin/env python3
"""Compare reconstruction peaks for different detector lengths.

This script mirrors the peak/1σ extraction from ``neutrino_money_plot.py`` but it
is tailored for files named ``comparison_with_lepton_neutron_event<E>_<L>meter.root``.
For each ROOT file it

* pulls the histogram named ``reco_{energy}GeV_with_noise``
* computes the peak location and 68% containment interval
* groups the results by detector length (1.0, 1.5, 5.0 m, etc.)
* produces a compact text/CSV summary and an overlay plot showing the three
  lengths on the same axes.

Example:

    ./length_overlap_money_plot.py \\
        --inputs . \\
        --pattern 'comparison_with_lepton_neutron_event*_*.root' \\
        --csv results/length_overlay.csv \\
        --plot results/length_overlay.png
"""

from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

try:  # pragma: no cover - optional dependency when --no-plot is used
    import matplotlib.pyplot as plt
except ImportError:  # pragma: no cover
    plt = None

try:
    import ROOT  # noqa: F401  # pylint: disable=unused-import
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "PyROOT is required but not available. Install/activate ROOT before running"
    ) from exc

from neutrino_money_plot import (  # type: ignore
    IntervalResult,
    analyze_file,
    iter_root_files,
)

LENGTH_PATTERN = re.compile(
    r"event(?P<energy>[0-9]+(?:\.[0-9]+)?)_+(?P<length>[0-9]+(?:\.[0-9]+)?)meter",
    re.IGNORECASE,
)
# Match neutrino_money_plot.py palette so peaks and error bands share colors
# Keep palette in sync with neutrino_money_plot.py but swap green/red ordering
# so the green series is plotted last (and vice versa for red).
MATPLOTLIB_COLORS = ('#1f77b4', '#d62728', '#2ca02c', '#9467bd', '#ff7f0e')
OVERLAP_ROUND_DIGITS = 4
MARKER_OFFSET_STEP = 0.003  # Horizontal jitter (GeV) applied when points coincide
ENERGY_GRID = [round(0.1 * idx, 10) for idx in range(1, 11)]  # 0.1 ... 1.0 GeV


@dataclass
class LengthResult:
    interval: IntervalResult
    energy_text: str
    energy_value: float
    length_text: str
    length_value: float


def parse_energy_and_length(path: Path) -> tuple[str, float, str, float]:
    match = LENGTH_PATTERN.search(path.stem)
    if not match:
        raise ValueError(
            f"Unable to parse energy/length from file name: {path.name}. Expected an 'event'"
            " token followed by energy and '<length>meter'."
        )
    energy_text = match.group('energy')
    length_text = match.group('length')
    try:
        energy_value = float(energy_text)
        length_value = float(length_text)
    except ValueError as exc:  # pragma: no cover - defensive guard
        raise ValueError(f"Invalid numeric value in file name {path.name}") from exc
    return energy_text, energy_value, length_text, length_value


def format_hist_patterns(energy_text: str, energy_value: float, template: str) -> List[str]:
    candidates = {
        energy_text,
        f"{energy_value:g}",
        f"{energy_value:.1f}",
        f"{energy_value:.2f}",
    }
    return [template.format(energy=token) for token in sorted(candidates)]


def collect_results(
    files: Sequence[Path],
    *,
    canvas_name: str,
    hist_template: str,
) -> List[LengthResult]:
    results: List[LengthResult] = []
    for file_path in files:
        try:
            energy_text, energy_value, length_text, length_value = parse_energy_and_length(
                file_path
            )
        except ValueError as exc:
            print(exc)
            continue

        patterns = format_hist_patterns(energy_text, energy_value, hist_template)
        try:
            entries = analyze_file(
                file_path,
                canvas_name=canvas_name,
                hist_patterns=patterns,
                max_hists=1,
            )
        except RuntimeError as exc:
            print(f"[warn] {file_path}: {exc}")
            continue
        if not entries:
            print(f"[warn] {file_path}: no histogram matched {patterns}")
            continue
        interval = entries[0]
        if interval.energy is None:
            interval.energy = energy_value  # type: ignore[attr-defined]
        results.append(
            LengthResult(
                interval=interval,
                energy_text=energy_text,
                energy_value=energy_value,
                length_text=length_text,
                length_value=length_value,
            )
        )
    return results


def filter_by_length(results: Sequence[LengthResult], allowed: Sequence[float] | None) -> List[LengthResult]:
    if not allowed:
        return list(results)
    kept: List[LengthResult] = []
    for res in results:
        for target in allowed:
            if math.isclose(res.length_value, target, rel_tol=0.0, abs_tol=1e-3):
                kept.append(res)
                break
    return kept


def filter_by_energy_grid(results: Sequence[LengthResult], grid: Sequence[float]) -> List[LengthResult]:
    allowed = list(grid)
    kept: List[LengthResult] = []
    for res in results:
        for target in allowed:
            if math.isclose(res.energy_value, target, rel_tol=0.0, abs_tol=1e-3):
                kept.append(res)
                break
    return kept


def group_by_length(results: Sequence[LengthResult]) -> Dict[float, List[LengthResult]]:
    grouped: Dict[float, List[LengthResult]] = {}
    for res in results:
        grouped.setdefault(res.length_value, []).append(res)
    for payload in grouped.values():
        payload.sort(key=lambda item: item.energy_value)
    return grouped


def report_missing_energies(
    grouped: Dict[float, List[LengthResult]], grid: Sequence[float]
) -> None:
    for length_value, entries in grouped.items():
        missing: List[float] = []
        for target in grid:
            found = any(
                math.isclose(item.energy_value, target, rel_tol=0.0, abs_tol=1e-3)
                for item in entries
            )
            if not found:
                missing.append(target)
        if missing:
            formatted = ", ".join(f"{val:.1f}" for val in missing)
            print(
                f"[warn] Length {length_value:.2f} m is missing energies: {formatted} GeV"
            )


def render_table(results: Sequence[LengthResult]) -> str:
    headers = [
        "file",
        "length(m)",
        "energy(GeV)",
        "peak",
        "x_low",
        "x_high",
        "-err",
        "+err",
        "frac",
    ]
    rows: List[List[str]] = []
    for res in sorted(results, key=lambda item: (item.length_value, item.energy_value)):
        interval = res.interval
        rel_path = interval.file_path
        try:
            rel_path = interval.file_path.relative_to(Path.cwd())
        except ValueError:
            pass
        rows.append(
            [
                str(rel_path),
                f"{res.length_value:.3g}",
                f"{res.energy_value:.3g}",
                f"{interval.peak_x:.4f}",
                f"{interval.x_low:.4f}",
                f"{interval.x_high:.4f}",
                f"{interval.err_low:.4f}",
                f"{interval.err_high:.4f}",
                f"{100.0 * interval.included_fraction:5.2f}%",
            ]
        )
    widths = [0] * len(headers)
    for idx, name in enumerate(headers):
        values = [name]
        values.extend(row[idx] for row in rows)
        widths[idx] = max(len(val) for val in values) if values else len(name)
    header = " ".join(name.ljust(widths[i]) for i, name in enumerate(headers))
    sep = "-" * len(header)
    lines = [header, sep]
    for row in rows:
        lines.append(" ".join(row[i].ljust(widths[i]) for i in range(len(headers))))
    return "\n".join(lines)


def write_csv(path: Path, results: Sequence[LengthResult]) -> None:
    lines = [
        "file,length_m,energy_GeV,peak,x_low,x_high,err_low,err_high,fraction"
    ]
    for res in sorted(results, key=lambda item: (item.length_value, item.energy_value)):
        interval = res.interval
        lines.append(
            f"{interval.file_path},{res.length_value:.6g},{res.energy_value:.6g},"
            f"{interval.peak_x:.6g},{interval.x_low:.6g},{interval.x_high:.6g},"
            f"{interval.err_low:.6g},{interval.err_high:.6g},{interval.included_fraction:.6g}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def compute_overlap_offsets(
    grouped: Dict[float, List[LengthResult]],
    *,
    digits: int = OVERLAP_ROUND_DIGITS,
    step: float = MARKER_OFFSET_STEP,
) -> Dict[float, List[float]]:
    """Return x-offsets so markers don't hide each other when peaks coincide."""

    offsets: Dict[float, List[float]] = {
        length_value: [0.0] * len(entries) for length_value, entries in grouped.items()
    }
    clusters: Dict[tuple[float, float], List[tuple[float, int]]] = {}
    for length_value in sorted(grouped):
        for idx, res in enumerate(grouped[length_value]):
            key = (round(res.energy_value, digits), round(res.interval.peak_x, digits))
            clusters.setdefault(key, []).append((length_value, idx))

    for members in clusters.values():
        count = len(members)
        if count <= 1:
            continue
        span = step * (count - 1)
        start = -span / 2.0
        for pos, (length_value, idx) in enumerate(members):
            offsets[length_value][idx] = start + pos * step
    return offsets


def plot_overlay(grouped: Dict[float, List[LengthResult]], output_path: Path) -> None:
    if plt is None:
        raise RuntimeError(
            "matplotlib is not available in this environment. Install it or rerun with --no-plot."
        )
    fig, ax = plt.subplots(figsize=(8, 5))
    plotted = False
    offsets = compute_overlap_offsets(grouped)

    for idx, length_value in enumerate(sorted(grouped)):
        entries = grouped[length_value]
        if not entries:
            continue

        jitter = offsets.get(length_value, [])
        energies = [res.energy_value + (jitter[i] if i < len(jitter) else 0.0) for i, res in enumerate(entries)]
        peaks = [res.interval.peak_x for res in entries]
        lower = [peak - res.interval.err_low for peak, res in zip(peaks, entries)]
        upper = [peak + res.interval.err_high for peak, res in zip(peaks, entries)]

        color = MATPLOTLIB_COLORS[idx % len(MATPLOTLIB_COLORS)]
        label = f"{length_value:.1f} m"

        # Change 'marker' here to customize the point shape:
        # 'o' = circle, 's' = square, '^' = triangle, '*' = star, 'D' = diamond
        MARKERS = ['o', 's', '^', '*', 'D', 'v']
        # Inside the loop:
        current_marker = MARKERS[idx % len(MARKERS)]
        ax.plot(energies, peaks, marker=current_marker, color=color, label=label)
        # ax.plot(energies, peaks, marker='o', color=color, label=label)

        ax.fill_between(energies, lower, upper, color=color, alpha=0.25)
        plotted = True

    if not plotted:
        plt.close(fig)
        print("[warn] Nothing to plot; skipping overlay figure")
        return

    ax.set_xlabel('Neutrino energy (GeV)')
    ax.set_ylabel('Peak angle (degrees)')
    ax.set_title('Reco peak ±1σ vs detector length')
    ax.grid(True, alpha=0.3)
    ax.legend(title='Length')
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved overlay plot to {output_path}")

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'inputs',
        nargs='*',
        type=Path,
        default=[Path('.')],
        help='Directories or ROOT files to scan (default: current directory).',
    )
    parser.add_argument(
        '--pattern',
        default='comparison_with_lepton_neutron_event*_*.root',
        help='Glob used to select ROOT files under each input path.',
    )
    parser.add_argument(
        '--canvas',
        default='first_reconstruction_result_on_neutrino',
        help='Canvas name that contains the reco histograms.',
    )
    parser.add_argument(
        '--hist-template',
        default='reco_{energy}GeV_with_noise',
        help='Template for the histogram name (energy placeholder is {energy}).',
    )
    parser.add_argument(
        '--energy-grid',
        type=float,
        nargs='*',
        default=ENERGY_GRID,
        help='List of neutrino energies (GeV) to keep. Default: 0.1 ... 1.0 in 0.1 steps.',
    )
    parser.add_argument(
        '--lengths',
        type=float,
        nargs='*',
        default=None,
        help='Optional subset of detector lengths to keep (meters).',
    )
    parser.add_argument(
        '--csv',
        type=Path,
        default=Path('length_overlap_money_plot.csv'),
        help='Path to write the CSV summary.',
    )
    parser.add_argument(
        '--no-csv',
        action='store_true',
        help='Disable writing the CSV summary.',
    )
    parser.add_argument(
        '--plot',
        type=Path,
        default=Path('length_overlap_money_plot.png'),
        help='Path to the overlay plot (PNG).',
    )
    parser.add_argument(
        '--no-plot',
        action='store_true',
        help='Skip producing the overlay figure.',
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    files = list(iter_root_files(args.inputs, args.pattern))
    if not files:
        raise SystemExit('No ROOT files matched the requested inputs/pattern.')
    results = collect_results(
        files,
        canvas_name=args.canvas,
        hist_template=args.hist_template,
    )
    if not results:
        raise SystemExit('No histograms were processed successfully.')
    energy_filtered = filter_by_energy_grid(results, args.energy_grid)
    if not energy_filtered:
        raise SystemExit('No entries matched the requested energy grid (0.1-1.0 GeV).')
    filtered = filter_by_length(energy_filtered, args.lengths)
    if not filtered:
        raise SystemExit('No entries remained after applying the length filter.')
    grouped = group_by_length(filtered)
    report_missing_energies(grouped, args.energy_grid)
    table = render_table(filtered)
    print(table)
    if not args.no_csv:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        write_csv(args.csv, filtered)
        print(f"Wrote CSV to {args.csv}")
    if not args.no_plot:
        plot_overlay(grouped, args.plot)


if __name__ == '__main__':
    main()
