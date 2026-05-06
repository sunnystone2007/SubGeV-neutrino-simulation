import os

import numpy as np

import ROOT  # type: ignore


class HistogramSpectrumSampler:
    def __init__(self, root_path: str, hist_name: str, seed: int | None = None):
        self.rng = np.random.default_rng(seed)

        # Allow shortcuts such as "~/file.root"
        full_path = os.path.expanduser(root_path)

        values, edges = self._load_hist(full_path, hist_name)

        # Ensure non-negative PDF
        pdf = np.clip(values.astype(float), 0.0, None)
        total = pdf.sum()
        if total <= 0:
            raise ValueError("Histogram has zero (or negative) total integral; cannot sample.")

        self.edges = edges
        self.cdf = np.cumsum(pdf) / total  # normalized CDF over bins

    def sample_one(self) -> float:
        u = self.rng.random()
        # pick a bin index according to CDF
        i = int(np.searchsorted(self.cdf, u, side="right"))
        i = min(max(i, 0), len(self.edges) - 2)

        # uniform within chosen bin
        lo, hi = self.edges[i], self.edges[i + 1]
        return float(self.rng.uniform(lo, hi))

    def sample(self, n: int) -> np.ndarray:
        return np.array([self.sample_one() for _ in range(n)], dtype=float)

    @staticmethod
    def _load_hist(path: str, hist_name: str) -> tuple[np.ndarray, np.ndarray]:
        f = ROOT.TFile.Open(path, "READ")
        if not f or f.IsZombie():
            raise OSError(f"Failed to open ROOT file: {path}")

        obj = f.Get(hist_name)
        if not obj:
            f.Close()
            raise KeyError(f"Object '{hist_name}' not found in {path}")

        hist = HistogramSpectrumSampler._ensure_hist(obj, hist_name)
        hist.SetDirectory(0)
        f.Close()

        nbins = hist.GetNbinsX()
        values = np.array([hist.GetBinContent(i + 1) for i in range(nbins)], dtype=float)

        edges = np.array([hist.GetBinLowEdge(i + 1) for i in range(nbins)], dtype=float)
        edges = np.append(edges, hist.GetBinLowEdge(nbins + 1))

        return values, edges

    @staticmethod
    def _ensure_hist(obj, hist_name: str):
        if hasattr(obj, "InheritsFrom") and obj.InheritsFrom("TH1"):
            return obj

        get_primitives = getattr(obj, "GetListOfPrimitives", None)
        if get_primitives is not None:
            prims = get_primitives()
            if prims is not None:
                for item in prims:
                    if hasattr(item, "InheritsFrom") and item.InheritsFrom("TH1"):
                        return item

        raise TypeError(
            f"Object '{hist_name}' of type {type(obj).__name__} is not a TH1 and contains no TH1 primitives."
        )
