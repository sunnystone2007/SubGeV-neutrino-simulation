from bdb import effective
from pathlib import Path
from typing import Any, List, Tuple, Optional, Union, Sequence
from collections import Counter
import ROOT
from IPython.core.pylabtools import figsize
import sys
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
import os
import matplotlib.patches as mpatches
from numpy.ma.core import asarray
import event
import array
from pathlib import Path

from ROOT import TH2F, TMarker, TCanvas, TLatex, TLegend
from ROOT import TF1

class Plotter:

    DEFAULT_COLOR_CYCLE = (
        '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
        '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'
    )

    PDG_COLOR_MAP = {
        2112: '#8e24aa',
        2212: '#fb8c00',
        11: '#e53935',
        -11: '#e53935',
        13: '#424242',
        -13: '#424242',
        22: '#039be5',
        211: '#43a047',
        -211: '#43a047',
        111: '#43a047',
        321: '#00897b',
        -321: '#00897b',
    }

    PDG_LABEL_MAP = {
        2112: 'neutron',
        2212: 'proton',
        11: 'electron',
        -11: 'positron',
        13: 'muon',
        -13: 'antimuon',
        22: 'photon',
        211: 'pi±',
        -211: 'pi±',
        111: 'pi0',
        321: 'K±',
        -321: 'K±',
    }

    def __init__(self, event):
        self.ME = None
        self.neutron_energy = None
        self.muon_energy = None
        self.USER_COLORS = None
        self.nearest_neutron_edep = None
        print("plotter: initialization")
        self.event = event
        #self.Jump(0, 0.5)

    def _color_for_pdg(self, pdg: int) -> str:
        color = self.PDG_COLOR_MAP.get(pdg)
        if color is not None:
            return color
        color = self.PDG_COLOR_MAP.get(-pdg)
        if color is not None:
            return color
        palette = self.USER_COLORS or list(self.DEFAULT_COLOR_CYCLE)
        return palette[abs(int(pdg)) % len(palette)]

    def _label_for_pdg(self, pdg: int) -> str:
        label = self.PDG_LABEL_MAP.get(pdg)
        if label is not None:
            return label
        label = self.PDG_LABEL_MAP.get(-pdg)
        if label is not None:
            return label
        return f'PDG {pdg}'

    def Collect(self, dQthreshold):
        ## coordinate per edep
        self.xx = np.array([])
        self.yy = np.array([])
        self.zz = np.array([])
        self.cc = np.array([])
        self.tt = np.array([])
        self.ee = np.array([])  ##energy of edep
        self.ll = np.array([])  ##length of track
        self.aa = np.array([])  ##angle relative to muon direction, does not specific which ancestor here
        self.ii = np.array([])  ##ancestor of edep
        self.pp = np.array([])  ##further distinguish  on both angle and ancestor
        self.pp_2 = np.array([])
        self.ME = np.array([])
        self.ancestor_pdgs = np.array([])
        self.USER_COLORS = list(self.DEFAULT_COLOR_CYCLE)

        ## coordinate per event
        self.muon_energy = np.array([])
        self.neutron_energy = np.array([])
        self.nearest_neutron_edep = np.array([])
        self.neutron_edep_muon_relate_angle = np.array([])

        mm2m = 0.001
        mm2cm = 0.1

        ##classification of events based on muon energy
        muon_energy = self.event.read_muon_energy()
        neutron_energy = self.event.read_neutron_energy()
        self.muon_energy = np.append(self.muon_energy, muon_energy)
        self.neutron_energy = np.append(self.neutron_energy, neutron_energy)
        nearest_neutron_id = 0
        nearest_neutron_distance = 10000


        ## give value to the coordinate at the event level
        # raw_muon,mom = self.event.read_muon_direction_and_amplitude()
        # if raw_muon is None or np.all(raw_muon == 0):
        #     raise RuntimeError("no muon direction!")
        # muon_direction = raw_muon / np.linalg.norm(raw_muon)

        ##give value to the coordinate inside every event
        for i, track in enumerate(self.event.tracks):
            ME = 0
            track_origin = track
            trkId = track.GetTrackId()
            depoList = track.association['depoList']
            ancestor = track.association['ancestor']
            ancestor = self.event.tracks[ancestor].GetPDGCode()
            pdg = track.GetPDGCode()
            while pdg != -11:
                trackid = track.GetParentId()
                track = self.event.tracks[trackid]
                pdg = track.GetPDGCode()
                if pdg == -11:
                    ME =1
                    # print("some of its parent are positrons")
                    # print(f"this track is {trkId}")
                if trackid ==0 or trackid == -1:
                    break




            for di in depoList:
                depo = self.event.depos[di]
                x = (depo.GetStart().X() + depo.GetStop().X()) / 2 * mm2m
                y = (depo.GetStart().Y() + depo.GetStop().Y()) / 2 * mm2m
                z = (depo.GetStart().Z() + depo.GetStop().Z()) / 2 * mm2m
                if x ** 2 + y ** 2 + z ** 2 <= nearest_neutron_distance ** 2:
                    nearest_neutron_id = di
                    nearest_neutron_distance = (x ** 2 + y ** 2 + z ** 2) ** 0.5
                t = (depo.GetStart().T() + depo.GetStop().T()) / 2  # ns
                e = depo.GetEnergyDeposit()  # MeV
                l = depo.GetTrackLength() * mm2cm  # most are 0.5 cm
                edep_vec = np.array([x, y, z])
                # cos = np.dot(edep_vec, muon_direction) / (np.linalg.norm(edep_vec) * np.linalg.norm(muon_direction))
                # aa = abs(np.degrees(np.arccos(np.clip(cos, -1, 1))))
                if ancestor == 2112:
                    ii = 0
                elif ancestor in (13, -13):
                    ii = 1
                else:
                    ii = 2
                # if ii == 0 and aa > 10:
                ##which means a neutron edep we are able to find
                    # pp = 0
                # elif ii == 0 and aa <= 10:
                    # pp = 1
                ##which means a neutron that is wrongly excluded
                # elif ii == 1 and aa > 10:
                #     pp = 2
                ##which means a muon wrongly included
                # elif ii == 1 and aa <= 10:
                #     pp = 3
                ##which means a muon successfully excluded
                # elif ii == 2 and aa > 10:
                #     pp = 4
                ##which means another particle wrongly included
                # elif ii == 2 and aa <= 10:
                    pp = 5
                # else:
                #     pp = 6
                #print(f"pp is {pp}")
                #
                # if ii == 0 and aa > 20:
                # ##which means a neutron edep we are able to find
                #     pp_2 = 0
                # elif ii == 0 and aa <= 20:
                #     pp_2 = 1
                # ##which means a neutron that is wrongly excluded
                # elif ii == 1 and aa > 20:
                #     pp_2 = 2
                # ##which means a muon wrongly included
                # elif ii == 1 and aa <= 20:
                #     pp_2 = 3
                # ##which means a muon successfully excluded
                # elif ii == 2 and aa > 20:
                #     pp_2 = 4
                # ##which means another particle wrongly included
                # elif ii == 2 and aa <= 20:
                #     pp_2 = 5
                # else:
                #     pp_2 = 6
                #print(f"pp_2 is {pp}")
                if self.event.ChargeBirksLaw(e, l) > dQthreshold:  # detector threshold
                    self.xx = np.append(self.xx, x)
                    self.yy = np.append(self.yy, y)
                    self.zz = np.append(self.zz, z)
                    self.tt = np.append(self.tt, t)
                    self.ee = np.append(self.ee, e)
                    self.ll = np.append(self.ll, l)
                    self.cc = np.append(self.cc, self._color_for_pdg(ancestor))
                    # self.aa = np.append(self.aa, aa)
                    self.ii = np.append(self.ii, ii)
                    # self.pp = np.append(self.pp, pp)
                    # self.pp_2 = np.append(self.pp_2, pp_2)
                    self.ME = np.append(self.ME, ME)
                    self.ancestor_pdgs = np.append(self.ancestor_pdgs, ancestor)
        depo = self.event.depos[nearest_neutron_id]
        x = (depo.GetStart().X() + depo.GetStop().X()) / 2
        y = (depo.GetStart().Y() + depo.GetStop().Y()) / 2
        z = (depo.GetStart().Z() + depo.GetStop().Z()) / 2
        neutron_direction = np.array([x, y, z]) / np.linalg.norm([x, y, z])
        #print("the neutron direction is :", neutron_direction)
        # self.neutron_edep_muon_relate_angle = np.append(self.neutron_edep_muon_relate_angle, abs(np.arctan(np.dot(neutron_direction, muon_direction))))
                #print("the relative angle between the muon direction and the shortest distance muon edep is:",
                      #abs(np.arctan(np.dot(neutron_direction, muon_direction))))

    def Jump(self, entryNo, dQthreshold):
        self.event.Jump(entryNo)
        self.Collect(dQthreshold)

    def Next(self, dQthreshold):
        self.event.Next()
        self.Collect(dQthreshold)

    def Prev(self, dQthreshold):
        self.event.Prev()
        self.Collect(dQthreshold)
    def neutron_E_P_plot(self):
        """here we want to make a distribution of the neutron energy at the first vertex versus its momentum magnitude there.
        we expect to see it first varies quadratically and then bend to a linear line because as energy increase, there shall be
        a smooth transfer from classic neutron to relativistic neutron. for example, for 5mev neutron, its speed is 0.1c and for 500 mev
        neutron its speed is 0.75c, highly relavilistic.

        ps: for now only antinumus with one neutron one muon are considered
        """
        magnitude_of_momentum_list,Energy_list,edep_list = self.event.find_Energy_and_momentum()
        print(f"magnitude_of_momentum_list: {magnitude_of_momentum_list}")
        print(f"Energy_list: {Energy_list}")
        total_events= len(Energy_list)
        canvas = ROOT.TCanvas("Canvas", "Canvas", 800, 600)
        energy_array = array.array('d', Energy_list)
        momentum_array = array.array('d', magnitude_of_momentum_list)
        graph = ROOT.TGraph(total_events,energy_array,momentum_array)
        graph.GetXaxis().SetTitle("Neutron Initial Kinetic Energy at Vertex[MeV]")
        graph.GetYaxis().SetTitle("Neutron Initial Momentum Magnitude at Vertex[MeV/c]")
        graph.SetMarkerSize(10)
        graph.SetMarkerSize(1)
        graph.SetMarkerColor(ROOT.kBlue + 2)
        graph.Draw("AP")
        output_filename = f"NIKEV_NIMMV_scatter.root"
        canvas.SaveAs(output_filename)
        print("file saved at :", output_filename)
    def check_final_state(self):
        max_events = 1000
        available_events = getattr(self.event, 'nEntry', 0)
        n_events = min(max_events, available_events)
        if n_events <= 0:
            print('Momentum conservation check skipped: no events available.')
            return

        angle_differences = []
        residual_norms = []

        for i in range(n_events):
            self.event.Jump(i)
            if i % 100 == 0:
                print(f'running {i}/{n_events} events')

            neutrino_momentum = np.array([
                self.event.info.get('Px_nu', 0.0),
                self.event.info.get('Py_nu', 0.0),
                self.event.info.get('Pz_nu', 0.0),
            ], dtype=float)
        nucleus_momentum = self.event.read_final_state_nucleus()
        fSI_particle_momentum = self.event.read_final_state_momentum()

    def check_momentum_conservation_2(self):
        """a plot showing TPMM and TPMM+nucleus loss"""

        max_events = 1000
        available_events = getattr(self.event, 'nEntry', 0)
        n_events = min(max_events, available_events)
        if n_events <= 0:
            print('Momentum conservation check skipped: no events available.')
            return

        cos_with_nucleus = ROOT.TH1D(
            'angle_with_nucleus',
            'angle(incoming,outgoing) with nucleus momentum',
            50,
            0.0,
            180.0,
        )
        cos_without_nucleus = ROOT.TH1D(
            'angle_without_nucleus',
            'angle(incoming,outgoing) without nucleus momentum',
            50,
            0.0,
            180.0,
        )
        for i in range(n_events):
            self.event.Jump(i)
            if i % 100 == 0:
                print(f'running {i}/{n_events} events')
            nucleus_momentum = self.event.read_final_state_nucleus()
            nucleus_momentum = np.array(nucleus_momentum*1000, dtype=float)
            neutrino_momentum = np.array([
                self.event.info.get('Px_nu', 0.0),
                self.event.info.get('Py_nu', 0.0),
                self.event.info.get('Pz_nu', 0.0),
            ], dtype=float)
            leptonic_hadronic = self.event.get_momentum_direction_and_amplitude(i)
            if leptonic_hadronic is None:
                continue

            muon_mom, neutron_mom = leptonic_hadronic
            muon_mom = np.array(muon_mom, dtype=float)
            neutron_mom = np.array(neutron_mom, dtype=float)

            total_proton_momentum = np.zeros(3)
            for is_track, _, momentum, *_ in self.event.proton_withoutthreshold():
                if is_track:
                    total_proton_momentum += np.array(momentum, dtype=float)

            total_pion_momentum = np.zeros(3)
            for is_track, _, momentum in self.event.pion_threshold():
                if is_track:
                    total_pion_momentum += np.array(momentum, dtype=float)

            incoming_total = neutrino_momentum
            outgoing_with_nucleus = (
                muon_mom
                + neutron_mom
                + total_proton_momentum
                # + total_pion_momentum
                + nucleus_momentum
            )
            outgoing_without_nucleus = (
                muon_mom
                + neutron_mom
                + total_proton_momentum
                # + total_pion_momentum
            )

            incoming_norm = np.linalg.norm(incoming_total)
            if incoming_norm == 0.0:
                continue

            with_norm = np.linalg.norm(outgoing_with_nucleus)
            without_norm = np.linalg.norm(outgoing_without_nucleus)
            if with_norm != 0.0:
                cos_angle = np.dot(incoming_total, outgoing_with_nucleus) / (incoming_norm * with_norm)
                cos_angle = np.clip(cos_angle, -1.0, 1.0)
                cos_with_nucleus.Fill(np.degrees(np.arccos(cos_angle)))

            if without_norm != 0.0:
                cos_angle_wo = np.dot(incoming_total, outgoing_without_nucleus) / (incoming_norm * without_norm)
                cos_angle_wo = np.clip(cos_angle_wo, -1.0, 1.0)
                cos_without_nucleus.Fill(np.degrees(np.arccos(cos_angle_wo)))

        if cos_with_nucleus.GetEntries() == 0 and cos_without_nucleus.GetEntries() == 0:
            print('No qualifying events for cos(angle) comparison.')
            return

        canvas = ROOT.TCanvas('Canvas', 'Momentum Conservation', 800, 600)
        cos_with_nucleus.SetLineColor(ROOT.kBlue + 1)
        cos_with_nucleus.SetLineWidth(2)
        cos_without_nucleus.SetLineColor(ROOT.kRed + 1)
        cos_without_nucleus.SetLineStyle(2)
        cos_without_nucleus.SetLineWidth(2)

        cos_with_nucleus.Draw('HIST')
        cos_without_nucleus.Draw('HIST SAME')
        cos_with_nucleus.GetXaxis().SetTitle('angle difference [degrees]')
        cos_with_nucleus.GetYaxis().SetTitle('counts')

        legend = ROOT.TLegend(0.15, 0.75, 0.55, 0.88)
        legend.AddEntry(cos_with_nucleus, 'with nucleus momentum', 'l')
        legend.AddEntry(cos_without_nucleus, 'without nucleus momentum', 'l')
        legend.Draw()

        output_filename = 'momentum_cosine_comparison_numuevent.root'
        canvas.SaveAs(output_filename)
        print('Saved comparison histogram to', output_filename)





















    def check_momentum_conservation(self):
        """Compare total incoming and outgoing momenta for each event.

        The incoming state combines the neutrino and the struck nucleus
        (if GENIE provides the latter). The outgoing state includes the
        reconstructed muon, neutron, and any protons/pions above the energy
        thresholds defined in :meth:`Event.proton_threshold` and
        :meth:`Event.pion_threshold`. The routine records both the angular
        difference between the incoming and outgoing momentum vectors and the
        magnitude of the residual momentum (|p_in - p_out|).
        """

        max_events = 1000
        available_events = getattr(self.event, 'nEntry', 0)
        n_events = min(max_events, available_events)
        if n_events <= 0:
            print('Momentum conservation check skipped: no events available.')
            return

        angle_differences = []
        residual_norms = []

        for i in range(n_events):
            self.event.Jump(i)
            if i % 100 == 0:
                print(f'running {i}/{n_events} events')

            neutrino_momentum = np.array([
                self.event.info.get('Px_nu', 0.0),
                self.event.info.get('Py_nu', 0.0),
                self.event.info.get('Pz_nu', 0.0),
            ], dtype=float)

            nucleus_momentum = self.event.read_final_state_nucleus()
            nucleus_momentum = np.array(nucleus_momentum*1000, dtype=float)

            leptonic_hadronic = self.event.get_momentum_direction_and_amplitude(i)
            if leptonic_hadronic is None:
                continue

            muon_mom, neutron_mom = leptonic_hadronic
            muon_mom = np.array(muon_mom, dtype=float)
            neutron_mom = np.array(neutron_mom, dtype=float)

            total_proton_momentum = np.zeros(3)
            for is_track, _, momentum, *_ in self.event.proton_withoutthreshold():
                if is_track:
                    total_proton_momentum += np.array(momentum, dtype=float)

            total_pion_momentum = np.zeros(3)
            for is_track, _, momentum in self.event.pion_threshold():
                if is_track:
                    total_pion_momentum += np.array(momentum, dtype=float)

            incoming_total = neutrino_momentum
            outgoing_total = (
                muon_mom
                + neutron_mom
                + total_proton_momentum
                # + total_pion_momentum
                + nucleus_momentum
            )
            # print(f'incoming total: {incoming_total}')
            # print(f'outgoing total: {outgoing_total}=')
            # print(f"{muon_mom}(muon_mom)+{neutron_mom}"
            #       f"(neutron_mom)+{nucleus_momentum}(nucleus_momentum)+...")
            # print('\n')
            incoming_norm = np.linalg.norm(incoming_total)
            outgoing_norm = np.linalg.norm(outgoing_total)

            if incoming_norm == 0.0 or outgoing_norm == 0.0:
                continue

            cos_angle = np.dot(incoming_total, outgoing_total) / (incoming_norm * outgoing_norm)
            cos_angle = np.clip(cos_angle, -1.0, 1.0)
            angle_differences.append(np.degrees(np.arccos(cos_angle)))

            residual_norm = np.linalg.norm(incoming_total - outgoing_total)
            residual_norms.append(residual_norm)

        if not angle_differences:
            print('Momentum conservation check found no qualifying events.')
            return

        canvas = ROOT.TCanvas('Canvas', 'Momentum Conservation', 800, 600)
        canvas.Divide(1, 2)

        angle_hist = ROOT.TH1D('angle_diff', 'Incoming vs outgoing momentum angle(0.2GeV AntiNuMu)', 50, 0, 180)
        angle_hist.GetXaxis().SetTitle('angle difference [degrees]')
        angle_hist.GetYaxis().SetTitle('counts')
        for value in angle_differences:
            angle_hist.Fill(value)

        residual_max = max(residual_norms) if residual_norms else 0.0
        upper_bound = residual_max * 1.1 if residual_max > 0.0 else 1.0
        residual_hist = ROOT.TH1D('residual_norm', 'Momentum residual magnitude(0.2GeV AntiNuMu)', 50, 0.0, upper_bound)
        residual_hist.GetXaxis().SetTitle('|p_{in} - p_{out}| [MeV/c]')
        residual_hist.GetYaxis().SetTitle('counts')
        for value in residual_norms:
            residual_hist.Fill(value)

        canvas.cd(1)
        angle_hist.Draw('HIST')
        canvas.cd(2)
        residual_hist.Draw('HIST')

        filepath = '/home/frunshi/frunshiworking/neutrino_reconstruction/momentum_conservation_check/results_by_simulation0.2GeVAntiNuMu.root'
        canvas.SaveAs(filepath)
        print('results_by_simulation saved at :', filepath)
        

    def hadronic_momentum_fraction(self, max_events: int = 9000,
                                   energy_min_GeV: float = 0.2,
                                   energy_max_GeV: float = 1.0,
                                   bin_width_GeV: float = 0.1) -> None:
        total_entries = getattr(self.event, 'nEntry', 0)
        n_events = min(max_events, total_entries)

        min_mev = energy_min_GeV * 1000.0
        max_mev = energy_max_GeV * 1000.0
        step_mev = bin_width_GeV * 1000.0

        # Construct bins such as 200-250, 250-300, ..., 950-1000 MeV.
        bin_edges = np.arange(min_mev, max_mev + step_mev, step_mev)
        bin_summaries = [
            {
                'muon': 0.0,
                'neutron': 0.0,
                'proton': 0.0,
                'pion': 0.0,
                'nucleus': 0.0,
            }
            for _ in range(len(bin_edges) - 1)
        ]
        bin_counts = np.zeros(len(bin_edges) - 1, dtype=int)

        # Order chosen so the pion/gamma/others component is stacked last (drawn on top).
        species = ('muon', 'neutron', 'proton', 'nucleus', 'pion')

        for idx in range(n_events):
            self.event.Jump(idx)
            if idx % 100 == 0:
                print(f'running {idx}/{n_events} events')

            e_nu = self.event.info.get('E_nu', 0.0)
            if e_nu < min_mev or e_nu >= max_mev:
                continue

            bin_index = int((e_nu - min_mev) // step_mev)
            if not (0 <= bin_index < len(bin_summaries)):
                continue

            leptonic_hadronic = self.event.get_momentum_direction_and_amplitude(idx)
            if leptonic_hadronic is None:
                continue

            muon_mom, neutron_mom = leptonic_hadronic
            muon_mag = float(np.linalg.norm(muon_mom))
            neutron_mag = float(np.linalg.norm(neutron_mom))

            total_proton_momentum = np.zeros(3)
            for is_track, _, momentum, *_ in self.event.proton_withoutthreshold():
                if is_track:
                    total_proton_momentum += np.array(momentum, dtype=float)
            proton_mag = float(np.linalg.norm(total_proton_momentum))

            nucleus_momentum = np.array(self.event.read_final_state_nucleus(), dtype=float) * 1000.0
            nucleus_mag = float(np.linalg.norm(nucleus_momentum))

            neutrino_momentum = np.array(
                [
                    self.event.info.get('Px_nu', 0.0),
                    self.event.info.get('Py_nu', 0.0),
                    self.event.info.get('Pz_nu', 0.0),
                ],
                dtype=float,
            )
            outgoing_sum = (
                np.array(muon_mom, dtype=float)
                + np.array(neutron_mom, dtype=float)
                + total_proton_momentum
                + nucleus_momentum
            )
            pion_momentum = neutrino_momentum - outgoing_sum
            pion_mag = float(np.linalg.norm(pion_momentum))

            total_mag = muon_mag + neutron_mag + proton_mag + pion_mag + nucleus_mag
            if total_mag <= 0.0:
                continue

            bin_counts[bin_index] += 1
            summary = bin_summaries[bin_index]
            summary['muon'] += muon_mag
            summary['neutron'] += neutron_mag
            summary['proton'] += proton_mag
            summary['pion'] += pion_mag
            summary['nucleus'] += nucleus_mag

        if not bin_counts.any():
            print('hadronic_momentum_fraction: no events passed the selection.')
            return

        bin_centers_GeV = (bin_edges[:-1] + bin_edges[1:]) / 2 / 1000.0
        stacked_values = {key: [] for key in species}
        for summary in bin_summaries:
            total = sum(summary.values())
            if total <= 0.0:
                for key in species:
                    stacked_values[key].append(0.0)
                continue
            for key in species:
                stacked_values[key].append(summary[key] / total * 100.0)

        edges_GeV = array.array('d', (bin_edges / 1000.0).tolist())
        histograms = {}
        color_map = {
            'muon': ROOT.kViolet + 1,
            'neutron': ROOT.kBlue + 1,
            'proton': ROOT.kOrange + 1,
            'pion': ROOT.kGreen + 2,
            'nucleus': ROOT.kRed + 1,
        }

        for key in species:
            hist = ROOT.TH1F(f'h_{key}_fraction', f'{key} fraction', len(bin_edges) - 1, edges_GeV)
            histograms[key] = hist
            values = stacked_values[key]
            for bin_idx, value in enumerate(values, start=1):
                hist.SetBinContent(bin_idx, value)
            hist.SetFillColor(color_map.get(key, ROOT.kGray + 1))
            hist.SetLineColor(ROOT.kBlack)

        stack = ROOT.THStack('hadronic_fraction_stack', 'Momentum composition vs E_{#nu}')
        legend = ROOT.TLegend(0.7, 0.55, 0.9, 0.88)
        for key in species:
            hist = histograms[key]
            stack.Add(hist)
            if key == 'pion':
                legend_label = 'pion, gamma and other particles'
            else:
                legend_label = f'{key.capitalize()} share'
            legend.AddEntry(hist, legend_label, 'f')

        canvas = ROOT.TCanvas('c_hadronic_fraction', 'Hadronic momentum fraction', 900, 600)
        stack.Draw('HIST')
        stack.SetMinimum(0.0)
        stack.SetMaximum(100.0)
        stack.GetXaxis().SetTitle('Neutrino energy E_{#nu} [GeV]')
        stack.GetYaxis().SetTitle('Momentum contribution [%]')
        legend.Draw()

        output_path = '~/frunshiworking/reference_plot/hadronic_momentum_fraction.root'
        canvas.SaveAs(output_path)
        print(f'Saved hadronic momentum fraction plot to {output_path}')





    def antinumu_event_muon_spectrum(self):
        n_events = 8000
        muon_energy_array = []
        for i in range(n_events):
            self.event.Jump(i)
            if i%256==0:
                print(f'event: {i}')
            muon_energy=self.event.read_muon_energy()
            muon_energy_array.append(muon_energy)
        canvas = ROOT.TCanvas("muon_energy_spectrum_numu", "muon_energy_spectrum_numu", 800, 600)
        h1 = ROOT.TH1F("muon_energy_spectrum_numu", "muon_energy_spectrum", 100, 0, 1000)
        for value in muon_energy_array:
            h1.Fill(value)
        h1.SetTitle("muon energy spectrum")
        h1.Draw("HIST")
        save_path = "/home/frunshi/frunshiworking/neutrino_reconstruction/muon_energy_spectrum_numu.root"
        canvas.SaveAs(save_path)

    def momentum_fraction(self):
        neutron_momentum = self.plot_particle_momentum_distribution(2112)
        proton_momentum = self.plot_particle_momentum_distribution(2212)
        pion_momentum = self.plot_particle_momentum_distribution(111,211,-211)
        muon_momentum = self.plot_particle_momentum_distribution(13,-13)
        sum = neutron_momentum + proton_momentum+ pion_momentum+ muon_momentum




        neutron_fraction = neutron_momentum / sum
        proton_fraction = proton_momentum / sum
        pion_fraction = pion_momentum / sum
        muon_fraction = muon_momentum / sum
        stack = ROOT.THStack('hadronic_fraction_stack', 'Momentum composition vs E_{#nu}')
        legend = ROOT.TLegend(0.7, 0.55, 0.9, 0.88)
        for key in species:
            hist = histograms[key]
            stack.Add(hist)
            if key == 'pion':
                legend_label = 'pion, gamma and other particles'
            else:
                legend_label = f'{key.capitalize()} share'
            legend.AddEntry(hist, legend_label, 'f')

        canvas = ROOT.TCanvas('c_hadronic_fraction', 'Hadronic momentum fraction', 900, 600)
        stack.Draw('HIST')
        stack.SetMinimum(0.0)
        stack.SetMaximum(100.0)
        stack.GetXaxis().SetTitle('Neutrino energy E_{#nu} [GeV]')
        stack.GetYaxis().SetTitle('Momentum contribution [%]')
        legend.Draw()

        output_path = '~/frunshiworking/reference_plot/hadronic_momentum_fraction.root'
        canvas.SaveAs(output_path)
        print(f'Saved hadronic momentum fraction plot to {output_path}')







    def hadronic_energy_fraction(self, max_events: int = 9000,
                                 energy_min_GeV: float = 0.2,
                                 energy_max_GeV: float = 1.0,
                                 bin_width_GeV: float = 0.1) -> None:
        """Plot stacked energy fractions of final-state particles vs neutrino energy."""

        total_entries = getattr(self.event, 'nEntry', 0)
        n_events = min(max_events, total_entries)

        min_mev = energy_min_GeV * 1000.0
        max_mev = energy_max_GeV * 1000.0
        step_mev = bin_width_GeV * 1000.0

        bin_edges = np.arange(min_mev, max_mev + step_mev, step_mev)
        bin_summaries = [
            {
                'muon': 0.0,
                'neutron': 0.0,
                'proton': 0.0,
                'nucleus': 0.0,
                'pion': 0.0,
            }
            for _ in range(len(bin_edges) - 1)
        ]
        bin_counts = np.zeros(len(bin_edges) - 1, dtype=int)
        species = ('muon', 'neutron', 'proton', 'nucleus', 'pion')

        def _final_state_nucleus_ke() -> float:
            genie_tree = getattr(self.event, 'genieTree', None)
            if genie_tree is None:
                return 0.0

            genie_tree.GetEntry(self.event.currentEntry)
            status = getattr(genie_tree, 'StdHepStatus', None)
            p4 = getattr(genie_tree, 'StdHepP4', None)
            count = int(getattr(genie_tree, 'StdHepN', 0))
            if status is None or p4 is None:
                print("no energy in nucleus found")
                return 0.0

            for part_idx in range(count):
                try:
                    if status[part_idx] != 15:
                        continue
                except (IndexError, TypeError):
                    break

                base = part_idx * 4
                try:
                    px = float(p4[base]) * 1000.0
                    py = float(p4[base + 1]) * 1000.0
                    pz = float(p4[base + 2]) * 1000.0
                    energy = float(p4[base + 3]) * 1000.0
                except (TypeError, IndexError):
                    print("no energy in nucleus found")
                    return 0.0

                momentum_sq = px * px + py * py + pz * pz
                rest_mass_sq = max(energy * energy - momentum_sq, 0.0)
                rest_mass = math.sqrt(rest_mass_sq)
                # print(f"the returned energy is {energy-rest_mass} GeV")
                return max(energy - rest_mass, 0.0)

            return 0.0

        for idx in range(n_events):
            self.event.Jump(idx)
            if idx % 100 == 0:
                print(f'running {idx}/{n_events} events')

            e_nu = self.event.info.get('E_nu', 0.0)
            if e_nu < min_mev or e_nu >= max_mev:
                continue

            bin_index = int((e_nu - min_mev) // step_mev)
            if not (0 <= bin_index < len(bin_summaries)):
                continue

            muon_energy = 0.0
            neutron_energy = 0.0
            proton_energy = 0.0
            for particle in getattr(self.event.vertex, 'Particles', []):
                pdg = particle.GetPDGCode()
                mom = particle.GetMomentum()
                total_energy = float(mom.E())
                mass = float(mom.M())
                kinetic_energy = max(total_energy - mass, 0.0)

                if pdg in (-13, 13):
                    muon_energy += total_energy
                elif pdg == 2112:
                    neutron_energy += kinetic_energy
                elif pdg == 2212:
                    proton_energy += kinetic_energy

            nucleus_energy = _final_state_nucleus_ke()*1000
            assigned_energy = muon_energy + neutron_energy + proton_energy + nucleus_energy
            pion_energy = max(e_nu - assigned_energy, 0.0)
            total_energy = assigned_energy + pion_energy
            if total_energy <= 0.0:
                continue

            bin_counts[bin_index] += 1
            summary = bin_summaries[bin_index]
            summary['muon'] += muon_energy
            summary['neutron'] += neutron_energy
            summary['proton'] += proton_energy
            summary['nucleus'] += nucleus_energy
            summary['pion'] += pion_energy

        if not bin_counts.any():
            print('hadronic_energy_fraction: no events passed the selection.')
            return

        bin_centers_GeV = (bin_edges[:-1] + bin_edges[1:]) / 2 / 1000.0
        stacked_values = {key: [] for key in species}
        for summary in bin_summaries:
            total = sum(summary.values())
            if total <= 0.0:
                for key in species:
                    stacked_values[key].append(0.0)
                continue
            for key in species:
                stacked_values[key].append(summary[key] / total * 100.0)

        edges_GeV = array.array('d', (bin_edges / 1000.0).tolist())
        histograms = {}
        color_map = {
            'muon': ROOT.kViolet + 1,
            'neutron': ROOT.kBlue + 1,
            'proton': ROOT.kOrange + 1,
            'pion': ROOT.kGreen + 2,
            'nucleus': ROOT.kRed + 1,
        }

        for key in species:
            hist = ROOT.TH1F(f'h_{key}_energy_fraction', f'{key} energy fraction', len(bin_edges) - 1, edges_GeV)
            histograms[key] = hist
            values = stacked_values[key]
            for bin_idx, value in enumerate(values, start=1):
                hist.SetBinContent(bin_idx, value)
            hist.SetFillColor(color_map.get(key, ROOT.kGray + 1))
            hist.SetLineColor(ROOT.kBlack)

        stack = ROOT.THStack('hadronic_energy_fraction_stack', 'Energy composition vs E_{#nu}')
        legend = ROOT.TLegend(0.7, 0.55, 0.9, 0.88)
        for key in species:
            hist = histograms[key]
            stack.Add(hist)
            if key == 'pion':
                legend_label = 'pion, gamma and other particles'
            else:
                legend_label = f'{key.capitalize()} share'
            legend.AddEntry(hist, legend_label, 'f')

        canvas = ROOT.TCanvas('c_hadronic_energy_fraction', 'Hadronic energy fraction', 900, 600)
        stack.Draw('HIST')
        stack.SetMinimum(0.0)
        stack.SetMaximum(100.0)
        stack.GetXaxis().SetTitle('Neutrino energy E_{#nu} [GeV]')
        stack.GetYaxis().SetTitle('Energy contribution [%]')
        legend.Draw()

        output_path = '~/frunshiworking/reference_plot/hadronic_energy_fraction.root'
        canvas.SaveAs(output_path)
        print(f'Saved hadronic energy fraction plot to {output_path}')

    def plot_particle_energy_distribution(
        self,
        *particle_pdgs: int,
        max_events: int = 9000,
        energy_min_GeV: float = 0.2,
        energy_max_GeV: float = 1.0,
        bin_width_GeV: float = 0.1,
        output_path: Optional[Union[str, Path]] = None,
    ) -> Optional[Path]:
        """Histogram average particle energy per neutrino-energy bin.

        Parameters mirror ``hadronic_*_fraction`` so the scan covers the same
        neutrino-energy window (default 0.2--1.0 GeV) while looping over up to
        9000 events. Each entry contributes the summed kinetic energy carried by
        the requested PDG codes (single or multiple) inside the event, averaged
        within the neutrino-energy spectrum bins.
        """

        if not particle_pdgs:
            particle_pdgs = (2112,)
        particle_pdgs = tuple(int(pdg) for pdg in particle_pdgs)

        total_entries = getattr(self.event, 'nEntry', 0)
        if total_entries <= 0:
            print('plot_particle_energy_distribution: no events available.')
            return None

        n_events = min(max_events, total_entries)
        if n_events <= 0:
            print('plot_particle_energy_distribution: no events selected.')
            return None

        if energy_max_GeV <= energy_min_GeV or bin_width_GeV <= 0:
            raise ValueError('Invalid neutrino energy bounds or bin width.')

        bin_edges_GeV = np.arange(energy_min_GeV, energy_max_GeV + bin_width_GeV, bin_width_GeV)
        if len(bin_edges_GeV) < 2:
            raise ValueError('Neutrino energy range is too small for histogram bins.')

        energy_sums = np.zeros(len(bin_edges_GeV) - 1, dtype=float)
        bin_counts = np.zeros(len(bin_edges_GeV) - 1, dtype=int)
        
        print(f"looking at energy spectrum of particle{particle_pdgs}...")

        for idx in range(n_events):
            self.event.Jump(idx)
            if idx % 100 == 0:
                print(f'running {idx}/{n_events} events')

            e_nu_GeV = self.event.info.get('E_nu', 0.0) / 1000.0
            if e_nu_GeV < energy_min_GeV or e_nu_GeV >= energy_max_GeV + 0.1:
                continue

            bin_index = np.searchsorted(bin_edges_GeV, e_nu_GeV, side='right') - 1
            bin_index = max(0, min(bin_index, len(energy_sums) - 1))
            if not (0 <= bin_index < len(energy_sums)):
                continue

            total_energy = 0.0
            for pdg in particle_pdgs:
                total_energy += float(self.event.particle_energy_per_event(pdg))
            energy_sums[bin_index] += total_energy
            bin_counts[bin_index] += 1

        if not bin_counts.any():
            print('plot_particle_energy_distribution: no events passed the selection.')
            return None

        averages = np.divide(energy_sums, bin_counts, out=np.zeros_like(energy_sums), where=bin_counts > 0)
        edges_array = array.array('d', bin_edges_GeV.tolist())
        pdg_label = '_'.join(str(pdg) for pdg in particle_pdgs)
        pdg_title = ','.join(str(pdg) for pdg in particle_pdgs)
        hist_name = f'h_particle_energy_pdg{pdg_label}'
        hist_title = f'Average particle energy (pdg={pdg_title}) vs E_{{#nu}} numu events'
        hist = ROOT.TH1F(hist_name, hist_title, len(edges_array) - 1, edges_array)
        for bin_idx, value in enumerate(averages, start=1):
            hist.SetBinContent(bin_idx, value)
        hist.SetMaximum(500)
        hist.GetXaxis().SetTitle('Neutrino energy E_{#nu} [GeV]')
        hist.GetYaxis().SetTitle('Average particle energy [MeV]')
        hist.SetLineColor(ROOT.kAzure + 2)
        hist.SetFillColor(ROOT.kAzure - 9)
        hist.SetFillStyle(3001)

        canvas = ROOT.TCanvas('c_particle_energy_distribution', 'Particle energy vs neutrino energy', 900, 600)
        hist.Draw('HIST')

        if output_path is None:
            output_path = f'~/frunshiworking/reference_plot/particle_energy_distribution_{pdg_label}.root'
        output_path = Path(output_path).expanduser()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        canvas.SaveAs(str(output_path))
        print(f'Saved particle energy distribution to {output_path}')
        return output_path

    def plot_particle_momentum_distribution(
        self,
        *particle_pdgs: int,
        max_events: int = 9000,
        energy_min_GeV: float = 0.2,
        energy_max_GeV: float = 1.0,
        bin_width_GeV: float = 0.1,
        output_path: Optional[Union[str, Path]] = None,
    ) -> Optional[Path]:
        """Histogram average particle momentum amplitude per neutrino-energy bin."""

        if not particle_pdgs:
            particle_pdgs = (2112,)
        particle_pdgs = tuple(int(pdg) for pdg in particle_pdgs)

        total_entries = getattr(self.event, 'nEntry', 0)
        if total_entries <= 0:
            print('plot_particle_momentum_distribution: no events available.')
            return None

        n_events = min(max_events, total_entries)
        if n_events <= 0:
            print('plot_particle_momentum_distribution: no events selected.')
            return None

        if energy_max_GeV <= energy_min_GeV or bin_width_GeV <= 0:
            raise ValueError('Invalid neutrino energy bounds or bin width.')

        bin_edges_GeV = np.arange(energy_min_GeV, energy_max_GeV + bin_width_GeV, bin_width_GeV)
        if len(bin_edges_GeV) < 2:
            raise ValueError('Neutrino energy range is too small for histogram bins.')

        momentum_sums = np.zeros(len(bin_edges_GeV) - 1, dtype=float)
        bin_counts = np.zeros(len(bin_edges_GeV) - 1, dtype=int)

        for idx in range(n_events):
            self.event.Jump(idx)
            if idx % 100 == 0:
                print(f'running {idx}/{n_events} events')

            e_nu_GeV = self.event.info.get('E_nu', 0.0) / 1000.0
            if e_nu_GeV < energy_min_GeV or e_nu_GeV >= energy_max_GeV + 0.1:
                continue

            bin_index = np.searchsorted(bin_edges_GeV, e_nu_GeV, side='right') - 1
            bin_index = max(0, min(bin_index, len(momentum_sums) - 1))
            if not (0 <= bin_index < len(momentum_sums)):
                continue

            total_momentum = 0.0
            for pdg in particle_pdgs:
                total_momentum += float(self.event.particle_momentum_amplitude_per_event(pdg))
            momentum_sums[bin_index] += total_momentum
            bin_counts[bin_index] += 1

        if not bin_counts.any():
            print('plot_particle_momentum_distribution: no events passed the selection.')
            return None

        averages = np.divide(momentum_sums, bin_counts, out=np.zeros_like(momentum_sums), where=bin_counts > 0)
        edges_array = array.array('d', bin_edges_GeV.tolist())
        pdg_label = '_'.join(str(pdg) for pdg in particle_pdgs)
        pdg_title = ','.join(str(pdg) for pdg in particle_pdgs)
        hist_name = f'h_particle_momentum_pdg{pdg_label}'
        hist_title = f'Average particle momentum amplitude (pdg={pdg_title}) vs E_{{#nu}}antinumu events'
        hist = ROOT.TH1F(hist_name, hist_title, len(edges_array) - 1, edges_array)
        for bin_idx, value in enumerate(averages, start=1):
            hist.SetBinContent(bin_idx, value)
        hist.GetXaxis().SetTitle('Neutrino energy E_{#nu} [GeV]')
        hist.GetYaxis().SetTitle('Average particle momentum amplitude [MeV/c]')
        hist.SetLineColor(ROOT.kViolet + 6)
        hist.SetFillColor(ROOT.kViolet - 5)
        hist.SetFillStyle(3002)

        canvas = ROOT.TCanvas('c_particle_momentum_distribution', 'Particle momentum vs neutrino energy', 900, 600)
        hist.Draw('HIST')

        if output_path is None:
            output_path = f'~/frunshiworking/reference_plot/particle_momentum_distribution_{pdg_label}.root'
        output_path = Path(output_path).expanduser()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        canvas.SaveAs(str(output_path))
        print(f'Saved particle momentum distribution to {output_path}')
        return averages









    def identify_neutrino_process(self, max_events=None, quiet=False, return_details=False):
        """Classify the underlying interaction type for a batch of events.

        Parameters
        ----------
        max_events : Optional[int]
            Cap the number of entries scanned. When ``None`` the full tree is
            processed.
        quiet : bool
            Suppress the default console printout when ``True``.

        Returns
        -------
        Counter
            Counts for each interaction label encountered. Keys are the
            uppercase interaction tags used by GENIE ("QES", "RES", etc.).
        other_details : Counter, optional
            When ``return_details`` is ``True`` a second ``Counter`` with the
            raw reaction names that were grouped into ``OTHER`` is returned.
        """

        total_entries = getattr(self.event, 'nEntry', 0)
        if total_entries <= 0:
            if not quiet:
                print('identify_neutrino_process: no events available to scan.')
            return Counter()

        if max_events is None:
            n_events = total_entries
        else:
            n_events = min(max_events, total_entries)

        interaction_counts = Counter({
            'QES': 0,
            'RES': 0,
            'MEC': 0,
            'DIS': 0,
            'COH': 0,
            'OTHER': 0,
        })
        other_breakdown = Counter()

        for i in range(n_events):
            self.event.Jump(i)
            reaction_name = self.event.vertex.GetReaction() or ''
            reaction_upper = reaction_name.upper()
            if 'QES' in reaction_upper:
                interaction_counts['QES'] += 1
            elif 'RES' in reaction_upper:
                interaction_counts['RES'] += 1
            elif 'MEC' in reaction_upper:
                interaction_counts['MEC'] += 1
            elif 'DIS' in reaction_upper:
                interaction_counts['DIS'] += 1
            elif 'COH' in reaction_upper:
                interaction_counts['COH'] += 1
            else:
                interaction_counts['OTHER'] += 1
                cleaned_name = reaction_name.strip() or 'UNKNOWN'
                other_breakdown[cleaned_name] += 1

        if not quiet:
            print(f'Processed {n_events} events (available: {total_entries}).')
            print(f"  Quasi Elastic Scattering (QES): {interaction_counts['QES']}")
            print(f"  Resonance (RES): {interaction_counts['RES']}")
            print(f"  Meson Exchange Current (MEC): {interaction_counts['MEC']}")
            print(f"  Deep Inelastic Scattering (DIS): {interaction_counts['DIS']}")
            print(f"  Coherent (COH): {interaction_counts['COH']}")
            if interaction_counts['OTHER']:
                print(f"  Other/Unclassified: {interaction_counts['OTHER']}")
                for name, count in other_breakdown.most_common(5):
                    print(f"    - {name}: {count}")

        if return_details:
            return interaction_counts, other_breakdown

        return interaction_counts
















    def neutron_scattering_process(self):
        elastic = []
        inelastic = []
        for i in range(9000):
            self.Jump(i,0)
            if i %100 ==0:
                print(i)
            neutron_energy = self.event.read_neutron_energy()
            subprocess = self.event.identify_neutron_scattering_type()
            # print(subprocess)
            if subprocess is not None:
                if subprocess == 121: ##inelastic scattering
                    inelastic.append(neutron_energy)
                elif subprocess == 111:
                    elastic.append(neutron_energy)

        all_energies = np.array(elastic + inelastic,dtype = float)
        E_min = 0
        E_max = np.max(all_energies)
        nbins = 60
        c =ROOT.TCanvas("Kinetic Energy of neutron", "Elastic vs Inelastic", 900, 650)

        hEl = ROOT.TH1D(
            "Elastic",
            "Elastic(sbuprocess code 111)",
            nbins,
            E_min,
            E_max,
        )
        hIn = ROOT.TH1D(
            "Inelastic",
            "Inelastic(subprocess code 121)",
            nbins,
            E_min,
            E_max,
        )
        print("number of all elastic events:",len(elastic))
        for E in elastic:
            hEl.Fill(E)
        for E in inelastic:
            hIn.Fill(E)

        # build stacked histogram: elastic on bottom, inelastic on top
        hEl.SetFillColor(ROOT.kBlue - 9)
        hEl.SetLineColor(ROOT.kBlue + 1)
        hIn.SetFillColor(ROOT.kRed - 4)
        hIn.SetLineColor(ROOT.kRed + 1)
        hEl.SetLineWidth(2)
        hIn.SetLineWidth(2)

        stack = ROOT.THStack("neutron_scatter_stack", "")
        stack.Add(hEl, "HIST")
        stack.Add(hIn, "HIST")
        stack.Draw("HIST")
        stack.GetXaxis().SetTitle("Neutron Kinetic energy[MeV]")
        stack.GetYaxis().SetTitle("Number of events")
        stack.SetTitle("Stack plot Elastic vs Inelastic Neutron Scattering Antinue events")

        leg = ROOT.TLegend(0.65, 0.75, 0.88, 0.88)
        leg.AddEntry(hEl, "Elastic scattering", "f")
        leg.AddEntry(hIn, "Inelastic scattering", "f")
        leg.Draw()
        filepath = f"~/frunshiworking/neutrino_reconstruction/antinue_neutron_process.root"
        c.SaveAs(filepath)
        print(f"Plot saved to {filepath}")






    def two_reconstruction_method_compare(self, neutrino_energy, energy_threshold=None ,require_michel_neutron=False):


        ## notice that the full event version only differ with two line in event: 3237
        n_events =1000
        neutron_reconstruction_list = []
        lepton_reconstruction_list = []
        without_neutron_reconstruction_list = []
        noise_free_reconstruction_list = []
        reconstruction_direction = []
        counter_neutrino_direction_reconstruction = 0
        counter_get_momentum_direction_and_amplitude = 0
        counter_without_neutron_reconstruction = 0
        counter_neutrino_direction_reconstruction_without_noise = 0
        counter_feature_requirement = 0
        filter_breakdown = Counter()
        threshold = self.event.energy_threshold if energy_threshold is None else energy_threshold
        cone_method = self.event.true_cone_angle_list_antinumu_brutal_95

        for i in range(n_events):
            if i%100 == 0:
                print(f"running {i}/1000 events")
            self.event.Jump(i)
            if require_michel_neutron:
                has_michel = self.event.has_michel_electron()
                has_neutron = self.event.has_primary_vertex_neutron()
                if not (has_michel and has_neutron):
                    counter_feature_requirement += 1
                    filter_breakdown['missing_michel_or_neutron'] += 1
                    continue
            result= self.event.neutrino_direction_reconstruction(
                i,
                cone_method=cone_method,
                energy_threshold=threshold,
            )
            if result is None or result == (None, None):
                counter_neutrino_direction_reconstruction += 1
                reason = getattr(self.event, 'last_reco_failure_reason', 'unknown')
                filter_breakdown[f"neutrino_reco:{reason}"] += 1
                continue
            neutron_reconstruction,lepton_reconstruction = result
            neutron_reconstruction_list.append(neutron_reconstruction)
            lepton_reconstruction_list.append(lepton_reconstruction)

            noise_free_result = self.event.neutrino_direction_reconstruction_without_noises_antinumu(
                i,
                cone_method=cone_method,
                energy_threshold=threshold,
            )
            if noise_free_result is None or noise_free_result == (None, None):
                counter_neutrino_direction_reconstruction_without_noise += 1
                reason_nf = getattr(self.event, 'last_reco_failure_reason', 'unknown')
                filter_breakdown[f"neutrino_reco_noise_free:{reason_nf}"] += 1
            else:
                noise_free_angle, _ = noise_free_result
                noise_free_reconstruction_list.append(noise_free_angle)

            result_without_neutron = self.event.neutrino_direction_reconstruction_without_neutron(i)
            if result_without_neutron is None or result_without_neutron == (None, None):
                counter_without_neutron_reconstruction += 1
                reason_wo = getattr(self.event, 'last_reco_failure_reason', 'unknown')
                filter_breakdown[f"neutrino_reco_wo_neutron:{reason_wo}"] += 1
            else:
                without_neutron_reconstruction, _ = result_without_neutron
                without_neutron_reconstruction_list.append(without_neutron_reconstruction)

            #read out the momentum of all particles from simulation
            simulation_direction = self.event.simulation_readout()
            angle_difference = np.degrees(np.arccos(simulation_direction[2]))
            reconstruction_direction.append(angle_difference)
        # print(f"neutron_reconstruction_list: {neutron_reconstruction_list}")

        print("neutrino_direction_reconstruction filter out",counter_neutrino_direction_reconstruction)
        print("neutrino_direction_reconstruction_without_noise filter out",counter_neutrino_direction_reconstruction_without_noise)
        print("neutrino_direction_reconstruction_without_neutron filter out",counter_without_neutron_reconstruction)
        print("get_momentum_direction_and_amplitude filter out",counter_get_momentum_direction_and_amplitude)
        if require_michel_neutron:
            print("Michel+neutron requirement filter out", counter_feature_requirement)
        if filter_breakdown:
            print("Filter breakdown (reason: count):")
            for reason, count in filter_breakdown.most_common():
                print(f"  {reason}: {count}")
        neutron_muon_resolution = 0
        muon_only_resolution = 0
        without_neutron_resolution = 0
        noise_free_resolution = 0
        for value in neutron_reconstruction_list:
            neutron_muon_resolution += value
        if neutron_reconstruction_list:
            neutron_muon_resolution /= len(neutron_reconstruction_list)
        for value in lepton_reconstruction_list:
            muon_only_resolution += value
        if lepton_reconstruction_list:
            muon_only_resolution /= len(lepton_reconstruction_list)
        for value in without_neutron_reconstruction_list:
            without_neutron_resolution += value
        if without_neutron_reconstruction_list:
            without_neutron_resolution /= len(without_neutron_reconstruction_list)
        for value in noise_free_reconstruction_list:
            noise_free_resolution += value
        if noise_free_reconstruction_list:
            noise_free_resolution /= len(noise_free_reconstruction_list)

        canvas = ROOT.TCanvas("first_reconstruction_result_on_neutrino","two method of neutrino reconstruction comparison",800,600)
        # hist1 = ROOT.TH1D("Lepton_only_reco","Lepton_only_reco",40,0,180)
        title = (f"{neutrino_energy} GeV anti muon neutrino direction reconstruction ")
        name = f"reco_{neutrino_energy}GeV_with_noise"

        hist2 = ROOT.TH1D(name, title, 40, 0, 180)
        hist4 = ROOT.TH1D("TRPM",title ,40,0,180)
        hist3 = ROOT.TH1D(f"h_lepton_proton_reco_{neutrino_energy}GeV", title, 40, 0, 180)
        # hist5 = ROOT.TH1D(f"h_noise_free_reco_{neutrino_energy}GeV", title, 40, 0, 180)

        print("the length of reconstruction with noise:",len(neutron_reconstruction_list))
        for value in neutron_reconstruction_list:
            hist2.Fill(value)
        # for value in lepton_reconstruction_list:
        #     hist1.Fill(value)
        for value in reconstruction_direction:
            hist3.Fill(value)
        for value in without_neutron_reconstruction_list:
            hist4.Fill(value)
        # for value in noise_free_reconstruction_list:
        #     hist5.Fill(value)

        legend = ROOT.TLegend(0.60, 0.75, 0.88, 0.88)
        # legend.AddEntry(hist1, "Lepton(Muon) Only Method", "l")
        legend.AddEntry(hist2, "R1 ", "l")
        legend.AddEntry(hist3, "R2 ", "l")
        legend.AddEntry(hist4, "R3 ", "l")
        # legend.AddEntry(hist5, "R1*", "l")
        # hist1.SetLineColor(ROOT.kBlue + 1)
        hist4.SetLineColor(ROOT.kGreen + 1)
        hist3.SetLineColor(ROOT.kRed + 1)
        hist2.SetLineColor(ROOT.kBlue + 2)
        # hist5.SetLineColor(ROOT.kAzure + 2)
        hist2.SetLineWidth(6)
        hist3.SetLineWidth(6)
        hist4.SetLineWidth(6)
        # hist5.SetLineWidth(6)
        # hist1.SetTitle("Lepton Neutron Reconstruction Comparison")
        # ROOT freezes axis limits after the first Draw, so align the maximum
        # before overlaying the different reconstruction strategies (R1/R2/R3).
        max_bin_content = max(
            hist2.GetMaximum(),
            hist3.GetMaximum(),
            hist4.GetMaximum(),
            # hist5.GetMaximum(),
        )
        hist2.SetMaximum(max_bin_content * 1.05 if max_bin_content else 1.0)

        hist2.Draw("HIST")
        hist3.Draw("HIST SAME")
        hist4.Draw("HIST SAME")
        # hist5.Draw("HIST SAME")
        legend.Draw()
        filepath = f"~/frunshiworking/1.5meter_reconstruction/antinumu_reconstruction_0326/comparison_with_lepton_neutron_event{round(neutrino_energy,2)}_{self.event.inspection_length_energy_side}meter.root"
        canvas.SaveAs(filepath)
        print(f"Plot saved to {filepath}")
        return muon_only_resolution, neutron_muon_resolution

    def _collect_reconstruction_distributions(
        self,
        *,
        n_events: int = 1000,
        require_michel_neutron: bool = False,
        verbose: bool = False,
    ) -> dict:
        """Gather reconstruction angles shared by comparison plots."""

        neutron_reconstruction_list: List[float] = []
        lepton_reconstruction_list: List[float] = []
        without_neutron_reconstruction_list: List[float] = []
        reconstruction_direction: List[float] = []

        counter_neutrino_direction_reconstruction = 0
        counter_without_neutron_reconstruction = 0
        counter_feature_requirement = 0
        filter_breakdown = Counter()

        for i in range(n_events):
            if verbose and i % 100 == 0:
                print(f'collecting event {i}/{n_events}')
            self.event.Jump(i)

            if require_michel_neutron:
                has_michel = self.event.has_michel_electron()
                has_neutron = self.event.has_primary_vertex_neutron()
                if not (has_michel and has_neutron):
                    counter_feature_requirement += 1
                    filter_breakdown['missing_michel_or_neutron'] += 1
                    continue

            result = self.event.neutrino_direction_reconstruction(i)
            if result is None or result == (None, None):
                counter_neutrino_direction_reconstruction += 1
                reason = getattr(self.event, 'last_reco_failure_reason', 'unknown')
                filter_breakdown[f'neutrino_reco:{reason}'] += 1
                continue

            neutron_reconstruction, lepton_reconstruction = result
            neutron_reconstruction_list.append(neutron_reconstruction)
            lepton_reconstruction_list.append(lepton_reconstruction)

            result_wo = self.event.neutrino_direction_reconstruction_without_neutron(i)
            if result_wo is None or result_wo == (None, None):
                counter_without_neutron_reconstruction += 1
                reason_wo = getattr(self.event, 'last_reco_failure_reason', 'unknown')
                filter_breakdown[f'neutrino_reco_wo_neutron:{reason_wo}'] += 1
            else:
                without_neutron_reconstruction, _ = result_wo
                without_neutron_reconstruction_list.append(without_neutron_reconstruction)

            simulation_direction = self.event.simulation_readout()
            angle_difference = np.degrees(np.arccos(simulation_direction[2]))
            reconstruction_direction.append(angle_difference)

        muon_only_resolution = (
            float(np.mean(lepton_reconstruction_list)) if lepton_reconstruction_list else None
        )
        neutron_muon_resolution = (
            float(np.mean(neutron_reconstruction_list)) if neutron_reconstruction_list else None
        )
        without_neutron_resolution = (
            float(np.mean(without_neutron_reconstruction_list))
            if without_neutron_reconstruction_list
            else None
        )

        return {
            'lepton_only': lepton_reconstruction_list,
            'lepton_neutron': neutron_reconstruction_list,
            'lepton_proton_pion': without_neutron_reconstruction_list,
            'truth': reconstruction_direction,
            'means': {
                'lepton_only': muon_only_resolution,
                'lepton_neutron': neutron_muon_resolution,
                'lepton_proton_pion': without_neutron_resolution,
            },
            'counters': {
                'neutrino_reco_fail': counter_neutrino_direction_reconstruction,
                'without_neutron_fail': counter_without_neutron_reconstruction,
                'feature_requirement_fail': counter_feature_requirement,
            },
            'filter_breakdown': filter_breakdown,
        }

    def _prepare_output_path(self, output_path: Optional[Union[str, Path]]) -> Path:
        """Return an expanded output path, creating directories when needed."""

        resolved = Path(output_path).expanduser()
        resolved.parent.mkdir(parents=True, exist_ok=True)
        return resolved

    def lepton_only_reco_hist_antinue(self,neutrino_energy):
        angle_deviation_list =[]
        bins = 90
        for i in range(1000):

            if i%100==0:
                print(i)
            self.Jump(i,0.5)
            angle_deviation = self.event.lepton_only_reconstruction_antinue_per_event()
            # print("angle deviation is :",angle_deviation)
            angle_deviation_list.append(angle_deviation)
        save_path = self._prepare_output_path(
            f'~/frunshiworking/antinumu_seperateplot/lepton_only_reco_{neutrino_energy:.2f}GeV_antinue.root'
        )
        canvas = ROOT.TCanvas('c_lepton_only_reco', 'Lepton-only reconstruction', 800, 600)
        hist = ROOT.TH1D(
            f'h_lepton_only_reco_{neutrino_energy:.2f}GeV',
            f'Lepton-only reconstruction;Angle difference [deg];Counts ({bins} bins)',
            bins,
            0,
            180,
        )
        for value in angle_deviation_list:
            hist.Fill(value)
        hist.SetLineColor(ROOT.kBlue + 1)
        hist.SetLineWidth(2)
        hist.Draw('HIST')
        canvas.SaveAs(str(save_path))
        print(f'Lepton-only histogram saved to {save_path}')


    def neutrino_direction_leptononly_vs_neutroninclude_antinue(
        self,neutrino_energy,
        *,
        n_events: int = 1000,
        bins: int = 40,
        output_path: Optional[Union[str, Path]] = None,
    ):
        """Overlay lepton-only, hadrons w/o n, and lepton+neutron angles for anti-nu_e CC events."""

        max_events = min(n_events, self.event.nEntry)
        lepton_only_angles: List[float] = []
        TRPM_angles:List[float] = []
        neutron_inclusion_angles: List[float] = []
        neutron_exclusion_angles: List[float] = []
        # neutron_exclusion_angles2:List[float] = []
        failures = 0

        for idx in range(max_events):
            if idx % 100 == 0:
                print(f'processing  nu_e event {idx}/{max_events}')

            result = self.event.neutrino_direction_reconstruction_antinue(idx)
            if result is None:
                failures += 1
                continue

            TRPM, _, neutron_exclusion_angle,neutron_inclusion_angle, _, _ = result
            # lepton_only_angles.append(lepton_only_angle)
            TRPM_angles.append(TRPM)
            neutron_exclusion_angles.append(neutron_exclusion_angle)
            neutron_inclusion_angles.append(neutron_inclusion_angle)
            # neutron_exclusion_angles.append(neutron_exclusion_angle)

        # if not lepton_only_angles or not neutron_inclusion_angles:
        #     print('No valid antinu_e reconstruction angles collected.')
        #     return None, None, None

        canvas = ROOT.TCanvas(
            'c_antinue_lepton_vs_neutron',
            ' nu_e neutron_exclusion trackable proton vs all proton(both cases neutron absent)',
            900,
            700,
        )
        # hist_lepton = ROOT.TH1D(
        #     'h_antinue_trackable proton',
        #     f'Anti-#nu_{{e}} reconstruction ({neutrino_energy:.2f} GeV);Angle difference [deg];Counts',
        #     int(bins),
        #     0.0,
        #     180.0,
        # )
        hist_TRPM = ROOT.TH1D(
            'hist_TRPM',
            'TRPM',
            int(bins),
            0.0,
            180.0,)


        hist_neutron = ROOT.TH1D(
            'h_antinue_all proton',
            f'# anti nu_{{e}} reconstruction ({neutrino_energy:.2f} GeV);Angle difference [deg];Counts',
            int(bins),
            0.0,
            180.0,
        )

        hist_exclusion = ROOT.TH1D(
            'h_anti nue_neutron_exclude',
            f'# anti nu_{{e}} reconstruction ({neutrino_energy:.2f} GeV);Angle difference [deg];Counts',
            int(bins),
            0.0,
            180.0,
        )

        # for value in neutron_exclusion_angles:
        #     hist_lepton.Fill(value)
        for value in TRPM_angles:
            hist_TRPM.Fill(value)
        # for value in neutron_exclusion_angles2:
        #     hist.Fill(value)
        for value in neutron_exclusion_angles:
            hist_exclusion.Fill(value)
        for value in neutron_inclusion_angles:
            hist_neutron.Fill(value)


        # ROOT freezes the axis range when the first histogram is drawn, so
        # explicitly set the maximum to fit all overlaid histograms.
        max_bin_content = max(
            # hist_lepton.GetMaximum(),
            hist_neutron.GetMaximum(),
            hist_TRPM.GetMaximum()
            # hist_exclusion.GetMaximum(),
        )
        hist_neutron.SetMaximum(max_bin_content * 1.05 if max_bin_content else 1.0)

        # hist_lepton.SetLineColor(ROOT.kBlue + 1)
        # hist_lepton.SetLineWidth(2)
        hist_TRPM.SetMaximum(max_bin_content * 1.05 if max_bin_content else 1.0)
        hist_TRPM.SetLineColor(ROOT.kGreen)
        hist_TRPM.SetLineWidth(5)
        hist_neutron.SetLineColor(ROOT.kRed + 1)
        hist_neutron.SetLineWidth(5)
        hist_exclusion.SetLineColor(ROOT.kBlue + 1)
        hist_exclusion.SetLineWidth(5)

        # hist_exclusion.SetLineColor(ROOT.kBlack + 2)
        # hist_exclusion.SetLineWidth(2)

        hist_neutron.Draw('HIST')
        hist_exclusion.Draw('HIST SAME')
        # hist_lepton.Draw('HIST SAME')
        hist_TRPM.Draw('HIST SAME')

        # mean_lepton = float(np.mean(lepton_only_angles))
        mean_neutron = float(np.mean(neutron_inclusion_angles))
        mean_exclusion = float(np.mean(neutron_exclusion_angles))
        mean_TRPM = float(np.mean(TRPM_angles))
        legend = ROOT.TLegend(0.55, 0.68, 0.88, 0.88)
        # legend.AddEntry(hist_lepton, f'other track-like particle+trackable proton (mean={mean_lepton:.1f} deg)', 'l')
        legend.AddEntry(hist_TRPM, f'TRPM (mean={mean_TRPM:.2f} deg)','l',)
        legend.AddEntry(
            hist_neutron,
            f'other track-like particle+ neutron(mean={mean_neutron:.1f} deg)',
            'l',
        )

        legend.AddEntry(
            hist_exclusion,
            f'Lepton+hadron-neutron (mean={mean_exclusion:.1f} deg)',
            'l',
        )
        legend.Draw()

        if failures:
            print(f'Skipped {failures} events due to reconstruction failure.')

        default_path = f'~/frunshiworking/1.5meter_reconstruction/antinue_reco_0326/antinue_neutrino_reconstruction{neutrino_energy}_2_0326.root'
        save_path = self._prepare_output_path(output_path or default_path)
        canvas.SaveAs(str(save_path))
        print(f'-nu_e proton verification histogram saved to {save_path}')

        return lepton_only_angles, neutron_inclusion_angles, save_path


    def neutrino_direction_lepton_vs_neutron_only(
        self, neutrino_energy,
        *,
        n_events: int = 1000,
        bins: int = 40,
        output_path: Optional[Union[str, Path]] = None):

        max_events = min(n_events, self.event.nEntry)
        TRPM_angles: list[float] = []
        lepton_only_angles: List[float] = []
        neutron_angles:List[float] = []
        failures = 0

        for idx in range(max_events):
            if idx % 100 == 0:
                print(f'processing nu_e event {idx}/{max_events}')

            result = self.event.neutrino_direction_reconstruction_antinue(idx)
            neutron_angle_only = self.event.reconstruct_neutron_only(idx)
            if result is None:
                failures += 1
                continue
            TRPM, _, neutron_exclusion_angle, _, _, _ = result
            # lepton_only_angles.append(lepton_only_angle)
            TRPM_angles.append(TRPM)
            lepton_only_angles.append(neutron_exclusion_angle)
            neutron_angles.append(neutron_angle_only)





    def lepton_only_reco_hist(
        self,
        neutrino_energy: float,
        *,
        n_events: int = 1000,
        require_michel_neutron: bool = False,
        bins: int = 40,
        output_path: Optional[Union[str, Path]] = None,
        distributions: Optional[dict] = None,
        verbose: bool = False,
    ) -> Tuple[Optional[dict], Optional[Path]]:
        """Plot the lepton-only reconstruction angle distribution individually."""

        data = distributions or self._collect_reconstruction_distributions(
            n_events=n_events,
            require_michel_neutron=require_michel_neutron,
            verbose=verbose,
        )
        lepton_only = data['lepton_only']
        if not lepton_only:
            print('No lepton-only reconstruction entries; skipping plot.')
            return data, None

        save_path = self._prepare_output_path(
            f'~/frunshiworking/antinumu_seperateplot/lepton_only_reco_{neutrino_energy:.2f}GeV.root',
        )
        canvas = ROOT.TCanvas('c_lepton_only_reco', 'Lepton-only reconstruction', 800, 600)
        hist = ROOT.TH1D(
            f'h_lepton_only_reco_{neutrino_energy:.2f}GeV',
            f'Lepton-only reconstruction;Angle difference [deg];Counts ({bins} bins)',
            bins,
            0,
            180,
        )
        for value in lepton_only:
            hist.Fill(value)
        hist.SetLineColor(ROOT.kBlue + 1)
        hist.SetLineWidth(2)
        hist.Draw('HIST')

        mean_value = data['means'].get('lepton_only')
        if mean_value is not None:
            latex = ROOT.TLatex(0.6, 0.85, f'#LT#Delta#theta#GT = {mean_value:.2f} deg')
            latex.SetNDC(True)
            latex.Draw()

        canvas.SaveAs(str(save_path))
        print(f'Lepton-only histogram saved to {save_path}')
        return data, save_path

    def count_radiative_influence(self,number_of_noises_point,radius_of_inspection,energy_threshold=None):
        threshold = self.event.energy_threshold if energy_threshold is None else energy_threshold
        count = 0
        log_path = f"/home/frunshi/frunshiworking/plots/radiative_influence_log_{energy_threshold}.txt"
        with open(log_path, "w", encoding="utf-8") as f:
            for i in range(9000):
                if i%100==0:
                    print(f"processing {i}/9000")

                self.Jump(i,threshold)
                indication = self.event.identify_radiactive_noise_influence(
                    number_of_noises_point,
                    radius_of_inspection,
                    threshold,
                )
                count += indication
            f.write(f"applying an energy threshold of {threshold}MeV\n")
            f.write(f"in 9000 events of antinumu, {count / 9000} events picked noises into reconstruction\n")





    def lepton_neutron_reco_hist(
        self,
        neutrino_energy: float,
        *,
        n_events: int = 1000,
        require_michel_neutron: bool = False,
        bins: int = 40,
        output_path: Optional[Union[str, Path]] = None,
        distributions: Optional[dict] = None,
        verbose: bool = False,
    ) -> Tuple[Optional[dict], Optional[Path]]:
        """Plot the combined lepton+neutron reconstruction histogram."""

        data = distributions or self._collect_reconstruction_distributions(
            n_events=n_events,
            require_michel_neutron=require_michel_neutron,
            verbose=verbose,
        )
        lepton_neutron = data['lepton_neutron']
        if not lepton_neutron:
            print('No lepton+neutron reconstruction entries; skipping plot.')
            return data, None

        save_path = self._prepare_output_path(

            f'~/frunshiworking/antinumu_seperateplot/lepton_neutron_reco_{neutrino_energy:.2f}GeV.root',
        )

        canvas = ROOT.TCanvas('c_lepton_neutron_reco', 'Lepton+neutron reconstruction', 800, 600)
        hist = ROOT.TH1D(
            f'h_lepton_neutron_reco_{neutrino_energy:.2f}GeV',
            'Lepton + Neutron reconstruction;Angle difference [deg];Counts',
            bins,
            0,
            180,
        )
        for value in lepton_neutron:
            hist.Fill(value)
        hist.SetLineColor(ROOT.kRed + 1)
        hist.SetLineWidth(2)
        hist.Draw('HIST')

        mean_value = data['means'].get('lepton_neutron')
        if mean_value is not None:
            latex = ROOT.TLatex(0.6, 0.85, f'#LT#Delta#theta#GT = {mean_value:.2f} deg')
            latex.SetNDC(True)
            latex.Draw()

        canvas.SaveAs(str(save_path))
        print(f'Lepton+neutron histogram saved to {save_path}')
        return data, save_path

    def truth_reco_hist(
        self,
        neutrino_energy: float,
        *,
        n_events: int = 1000,
        require_michel_neutron: bool = False,
        bins: int = 40,
        output_path: Optional[Union[str, Path]] = None,
        distributions: Optional[dict] = None,
        verbose: bool = False,
    ) -> Tuple[Optional[dict], Optional[Path]]:
        """Plot the TRPM/simulation reference angle distribution on its own canvas."""

        data = distributions or self._collect_reconstruction_distributions(
            n_events=n_events,
            require_michel_neutron=require_michel_neutron,
            verbose=verbose,
        )
        truth_angles = data['truth']
        if not truth_angles:
            print('No TRPM entries collected; skipping truth plot.')
            return data, None

        save_path = self._prepare_output_path(
            f'~/frunshiworking/antinumu_seperateplot/trpm_truth_{neutrino_energy:.2f}GeV.root',
        )
        canvas = ROOT.TCanvas('c_truth_reco', 'Truth/TRPM distribution', 800, 600)
        hist = ROOT.TH1D(
            f'h_truth_reco_{neutrino_energy:.2f}GeV',
            'True direction (TRPM);Angle difference [deg];Counts',
            bins,
            0,
            180,
        )
        for value in truth_angles:
            hist.Fill(value)
        hist.SetLineColor(ROOT.kGreen + 2)
        hist.SetLineWidth(2)
        hist.Draw('HIST')
        canvas.SaveAs(str(save_path))
        print(f'TRPM histogram saved to {save_path}')
        return data, save_path

    def lepton_proton_pion_reco_hist(
        self,
        neutrino_energy: float,
        *,
        n_events: int = 1000,
        require_michel_neutron: bool = False,
        bins: int = 40,
        output_path: Optional[Union[str, Path]] = None,
        distributions: Optional[dict] = None,
        verbose: bool = False,
    ) -> Tuple[Optional[dict], Optional[Path]]:
        """Plot the lepton+proton(+pion) reconstruction histogram individually."""

        data = distributions or self._collect_reconstruction_distributions(
            n_events=n_events,
            require_michel_neutron=require_michel_neutron,
            verbose=verbose,
        )
        without_neutron = data['lepton_proton_pion']
        if not without_neutron:
            print('No lepton+proton+pion reconstruction entries; skipping plot.')
            return data, None

        save_path = self._prepare_output_path(
            f'~/frunshiworking/antinumu_seperateplot/lepton_proton_pion_reco_{neutrino_energy:.2f}GeV.root',
        )
        canvas = ROOT.TCanvas('c_lepton_proton_pion_reco', 'Lepton + proton(+pion) reconstruction', 800, 600)
        hist = ROOT.TH1D(
            f'h_lepton_proton_pion_reco_{neutrino_energy:.2f}GeV',
            'Lepton + Proton(+Pion);Angle difference [deg];Counts',
            bins,
            0,
            180,
        )
        for value in without_neutron:
            hist.Fill(value)
        hist.SetLineColor(ROOT.kMagenta + 1)
        hist.SetLineWidth(2)
        hist.Draw('HIST')

        mean_value = data['means'].get('lepton_proton_pion')
        if mean_value is not None:
            latex = ROOT.TLatex(0.6, 0.85, f'#LT#Delta#theta#GT = {mean_value:.2f} deg')
            latex.SetNDC(True)
            latex.Draw()

        canvas.SaveAs(str(save_path))
        print(f'Lepton+proton(+pion) histogram saved to {save_path}')
        return data, save_path

    def TRPM_antinue(self):
        """TRPM stands for True Reconstructable Particle Momentum"""
        max_events = min(1000,self.event.nEntry)
        angles = []
        for entry in range(max_events):
            if entry % 100 == 0:
                print(f"processing event {entry}/{max_events}")





    def target_nucleus_direction_histogram(self, neutrino_energy: float = 0.2, n_events: int = 1000, bins: int = 40):
        """Draw a histogram of the target nucleus zenith angle for the sample."""

        max_events = min(n_events, self.event.nEntry)
        angles = []
        missing = 0

        for entry in range(max_events):
            if entry % 100 == 0:
                print(f"processing event {entry}/{max_events}")
            self.event.Jump(entry)
            angle = self.event.info.get('target_nucleus_zenith_angle_deg')
            if angle is None:
                missing += 1
                continue
            angles.append(angle)

        if not angles:
            print("No target nucleus direction information available for histogram.")
            return None

        canvas = ROOT.TCanvas("target_nucleus_direction", "Target nucleus direction", 800, 600)
        hist = ROOT.TH1D(
            f"h_target_nucleus_direction_{neutrino_energy}GeV",
            f"Target nucleus zenith angle ({neutrino_energy} GeV sample)",
            bins,
            0,
            180,
        )

        for value in angles:
            hist.Fill(value)

        hist.GetXaxis().SetTitle("θ(target nucleus, +z) [deg]")
        hist.GetYaxis().SetTitle("Counts")
        hist.Draw("HIST")

        output_dir = Path.home() / "frunshiworking" / "neutrino_reconstruction" / "target_nucleus"
        output_dir.mkdir(parents=True, exist_ok=True)
        filepath = output_dir / f"target_nucleus_direction_{neutrino_energy:.1f}GeV.root"
        canvas.SaveAs(str(filepath))
        mean_angle = sum(angles) / len(angles)

        print(
            f"Histogram saved to {filepath}. Used {len(angles)} events, "
            f"skipped {missing} with undefined direction. Mean angle: {mean_angle:.2f} degrees."
        )

        return hist
    def muon_minus_decay_rate(self):
        n_events =1000
        muon_decay_number = 0
        for i in range(n_events):
            # muon decay indicator will be zero if capture by nucleus and
            # one if decay.https://github.com/ECCE-EIC/geant4/blob/master/source/physics_lists/builders/OrderingParameterTable
            self.event.Jump(i)
            muon_decay_indicator = self.event.muon_decay_indicator()
            muon_decay_number+= muon_decay_indicator
        print("number of decay", muon_decay_number)
        rate_of_decay = muon_decay_number / n_events
        print("muon decay rate:", rate_of_decay)


    def draw_process_histogram(self,neutrino_energy):
        """here we want to see how many events are elastic and how many are not"""
        n_events = 1000
        compton = 0  #SubprocessCode 12
        Scintillation = 0   #22
        neutron_inelastic = 0   #121
        neutron_capture = 0   #403
        proton_inelastic = 0   #131
        hadron_elastic = 0  #111
        transportation = 0  #91
        ionization = 0  #2
        rayleigh = 0    #11
        compton_scattering = 0
        Bremsstrahlung = 0 # 3




        process_list = []
        sub_process_list = []
        other_process_list = []
        for i in range(n_events):
            self.event.Jump(i)
            if i%100 ==0:
                print(f"processing event: {i}/1000")
            process_list_event,sub_process_list_event= self.event.print_neutron_process()
            for process in process_list_event:
                process_list.append(process)
            for sub_process in sub_process_list_event:
                sub_process_list.append(sub_process)
        for value in sub_process_list:
            if value == 12:
                compton += 1
            elif value == 22 :
                Scintillation += 1
            elif value ==121 :
                neutron_inelastic += 1
            elif value == 403 :
                neutron_capture += 1
            elif value ==131 :
                proton_inelastic += 1
            elif value ==111:
                hadron_elastic +=1
            elif value == 3:
                Bremsstrahlung += 1
            elif value ==2:
                ionization += 1
            elif value == 91:
                transportation += 1
            elif value ==11:
                rayleigh += 1
            elif value == 13:
                compton_scattering += 1

            else :
                other_process_list.append(value)
        canvas = ROOT.TCanvas("Neutron_interaction_process", f"Neutron_interaction_process({neutrino_energy}GeV Neutrino)", 800, 600)
        hist1 = ROOT.TH1D("subprocess of Neutron interaction","subprocess of Neutron interaction",500,0,500)
        for value in sub_process_list:
            hist1.Fill(value)
        hist1.SetLineColor(ROOT.kBlue + 2)
        hist1.SetTitle("subprocess of Neutron interaction")
        hist1.Draw("HIST")
        filepath = f"~/frunshiworking/neutrino_reconstruction_result/Neutron_process_{neutrino_energy}GeV.root"
        canvas.SaveAs(filepath)
        print(f"Plot saved to {filepath}")
        print(f"in 1000 events there are {compton} compton and {Scintillation} gamma conversions"
              f"\n {neutron_inelastic} neutron inelastic and {proton_inelastic} proton inelastic"
              f"\n {ionization} ionization and {rayleigh} rayleigh"
              f"\n {Bremsstrahlung} Bremsstrahlung and {hadron_elastic} hadron elastic"
              f"\n {transportation} transportation "
              f"\n other processes include {other_process_list}")


    def neutron_energy_spectrum(self):
        neutron_energy_list=[]
        n_events = 9000
        for i in range(n_events):
            self.event.Jump(i)
            neutron_energy = self.event.read_neutron_energy()
            neutron_energy_list.append(neutron_energy)
        print(f"Neutron energy spectrum: {neutron_energy_list}")
        canvas = ROOT.TCanvas("neutron_energy_spectrum","neutron_energy_spectrum",800,600)
        h1= ROOT.TH1F("neutron_energy_spectrum","neutron_energy_spectrum",100,0,1000)
        for value in neutron_energy_list:
            h1.Fill(value)
        h1.SetTitle("neutron energy spectrum")
        h1.Draw("HIST")
        save_path = "/home/frunshi/frunshiworking/neutrino_reconstruction/neutron_energy_spectrum.root"
        canvas.SaveAs(save_path)

    def antinumu_neutron_process_histogram(
        self,
        *,
        max_events: int = 9000,
        save_path: Optional[Union[str, Path]] = None,
    ) -> Tuple[Optional[float], Optional[Path]]:
        """Scan anti-\nu_\mu events, summarize neutron processes, and plot energies.

        Returns
        -------
        (ratio_percent, output_path)
            ``ratio_percent`` is the percentage of elastic (subprocess 111)
            interactions w.r.t. the number of inelastic (subprocess 121)
            interactions. ``output_path`` records where the histogram was
            written (``None`` if the scan produced no qualifying neutrons).
        """

        total_entries = getattr(self.event, 'nEntry', 0)
        if total_entries <= 0:
            print('No events loaded; neutron process scan skipped.')
            return None, None

        scan_limit = min(max_events, total_entries)
        elastic_count = 0
        inelastic_count = 0
        neutron_energies_gev: List[float] = []
        processed_events = 0

        for entry in range(scan_limit):
            self.event.Jump(entry)
            if entry%100==0:
                print(f"{entry} events processed")
            if self.event.info.get('nu_pdg') != -14:
                continue

            processed_events += 1
            _, subprocesses = self.event.print_neutron_process()
            for code in subprocesses:
                if code == 111:
                    elastic_count += 1
                elif code == 121:
                    inelastic_count += 1

            neutron_energy_mev = self.event.read_neutron_energy()
            if neutron_energy_mev is not None:
                neutron_energies_gev.append(neutron_energy_mev / 1000.0)

        if processed_events == 0:
            print('No anti-nu_mu events encountered in the requested range.')
            return None, None

        if inelastic_count == 0:
            ratio_percent = None
            print('No inelastic (121) interactions observed; ratio undefined.')
        else:
            ratio_percent = (elastic_count / inelastic_count) * 100.0
            print(
                f'Elastic (111) interactions are {ratio_percent:.2f}% of '
                f'inelastic (121) interactions.'
            )

        energies_in_range = [e for e in neutron_energies_gev if 0.1 <= e <= 1.0]
        if not energies_in_range:
            print('No neutron energies fell within 0.1-1.0 GeV; nothing to plot.')
            return ratio_percent, None

        if save_path is None:
            output_dir = Path(self.event.plotpath)
            output_dir.mkdir(parents=True, exist_ok=True)
            save_path = output_dir / 'antinumu_neutron_energy_hist.root'
        else:
            save_path = Path(save_path).expanduser()
            save_path.parent.mkdir(parents=True, exist_ok=True)

        canvas = ROOT.TCanvas('c_antinumu_neutron_energy', 'Anti-nu_mu Neutron Spectrum', 900, 900)
        canvas.cd()

        pad_top = ROOT.TPad('pad_energy', '', 0.0, 0.38, 1.0, 1.0)
        pad_top.SetBottomMargin(0.02)
        pad_top.Draw()
        pad_top.cd()
        ROOT.gStyle.SetOptStat(0)
        bin_edges = array.array('d', [0.1 * i for i in range(1, 11)])
        energy_hist = ROOT.TH1F(
            'h_antinumu_neutron_energy',
            'Anti-#nu_{#mu} Neutron Energy Spectrum;Neutron kinetic energy [GeV];Counts per 0.1 GeV',
            len(bin_edges) - 1,
            bin_edges,
        )
        for energy in energies_in_range:
            energy_hist.Fill(energy)
        energy_hist.SetLineColor(ROOT.kAzure + 4)
        energy_hist.SetFillColor(ROOT.kAzure - 4)
        energy_hist.SetLineWidth(2)
        energy_hist.Draw('HIST')

        canvas.cd()
        pad_bottom = ROOT.TPad('pad_summary', '', 0.0, 0.0, 1.0, 0.38)
        pad_bottom.SetTopMargin(0.02)
        pad_bottom.SetBottomMargin(0.3)
        pad_bottom.SetGridy(True)
        pad_bottom.Draw()
        pad_bottom.cd()

        max_count = max(inelastic_count, elastic_count, 1)
        y_lim = max_count * 1.3
        frame_hist = ROOT.TH1F('h_interaction_frame', '', 2, 0.0, 2.0)
        frame_hist.SetMaximum(y_lim)
        frame_hist.SetMinimum(-y_lim)
        frame_hist.SetStats(False)
        xaxis = frame_hist.GetXaxis()
        xaxis.SetBinLabel(1, 'Inelastic (121)')
        xaxis.SetBinLabel(2, 'Elastic (111)')
        xaxis.CenterLabels(True)
        frame_hist.GetYaxis().SetTitle('Total interactions')
        frame_hist.Draw('AXIS')

        h_inelastic = ROOT.TH1F('h_inelastic_counts', '', 2, 0.0, 2.0)
        h_inelastic.SetBinContent(1, inelastic_count)
        h_inelastic.SetFillColor(ROOT.kRed - 4)
        h_inelastic.SetLineColor(ROOT.kRed + 1)
        h_inelastic.SetBarWidth(0.4)
        h_inelastic.SetBarOffset(-0.1)
        h_inelastic.Draw('BAR SAME')

        h_elastic = ROOT.TH1F('h_elastic_counts', '', 2, 0.0, 2.0)
        h_elastic.SetBinContent(2, -elastic_count)
        h_elastic.SetFillColor(ROOT.kBlue - 4)
        h_elastic.SetLineColor(ROOT.kBlue + 1)
        h_elastic.SetBarWidth(0.4)
        h_elastic.SetBarOffset(0.1)
        h_elastic.Draw('BAR SAME')

        legend = ROOT.TLegend(0.15, 0.75, 0.5, 0.92)
        legend.AddEntry(h_inelastic, 'Inelastic (121)', 'f')
        legend.AddEntry(h_elastic, 'Elastic (111)', 'f')
        legend.Draw()

        canvas.SaveAs(str(save_path))
        print(f'Saved neutron energy histogram to {save_path}')

        return ratio_percent, save_path


    def _track_span_length_cm(self, event_obj, track):
        """Return the summed path length between the first and last energy deposits for a track."""
        depo_indices = track.association.get('depoList', [])
        if not depo_indices:
            return None
        mm2cm = 0.1
        length_cm = 0.0
        for di in depo_indices:
            depo = event_obj.depos[di]
            length_cm += max(depo.GetTrackLength(), 0.0) * mm2cm
        if length_cm <= 0.0:
            return None
        return length_cm

    def _primary_particle_energy_length_histograms(
        self,
        *,
        pdg_codes: Tuple[int, ...],
        label: str,
        file_path: str,
        max_events: int,
        bins: int,
        output_path=None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """General routine to build primary-particle energy/length ROOT histograms.

        Parameters
        ----------
        pdg_codes: tuple of int
            PDG codes identifying the particle species of interest.
        label: str
            Short label inserted into histogram names/titles.
        file_path: str
            Input ROOT file path.
        max_events: int
            Maximum number of events to scan.
        bins: int
            Number of histogram bins.
        output_path: Optional[str | Path]
            Target ROOT file for the histograms. When omitted, defaults to
            ``self.event.plotpath`` with a filename derived from the label.
        """
        resolved_path = Path(file_path).expanduser()
        if not resolved_path.exists():
            raise FileNotFoundError(f"Requested file does not exist: {resolved_path}")

        event_obj = self.event
        using_external_event = False
        try:
            current_path = Path(event_obj.fileName).expanduser()
            same_file = current_path.resolve(strict=False) == resolved_path.resolve(strict=False)
        except Exception:
            same_file = False

        if not same_file:
            event_obj = event.Event(str(resolved_path), evgen=self.event.evgen)
            using_external_event = True

        energies = []
        lengths = []
        total_entries = min(max_events, event_obj.nEntry)
        if total_entries <= 0:
            print(f"No events available for {label} scan.")
            return np.asarray(energies), np.asarray(lengths)

        original_entry = event_obj.currentEntry if not using_external_event else None

        pdg_set = set(pdg_codes)
        for entry in range(total_entries):
            event_obj.Jump(entry)
            if not getattr(event_obj, 'vertex', None):
                continue

            track_lookup = {track.GetTrackId(): track for track in event_obj.tracks}

            for particle in event_obj.vertex.Particles:
                pdg = particle.GetPDGCode()
                if pdg not in pdg_set:
                    continue

                track_id = particle.GetTrackId()
                if track_id < 0:
                    continue

                track = track_lookup.get(track_id)
                if track is None or track.GetParentId() != -1:
                    continue

                if not track.association.get('depoList'):
                    continue

                mom = particle.GetMomentum()
                if not mom:
                    continue

                length_cm = self._track_span_length_cm(event_obj, track)
                if length_cm is None:
                    continue

                kinetic_energy = mom.E() - mom.M()
                energies.append(float(kinetic_energy))
                lengths.append(length_cm)

            if (entry + 1) % 500 == 0 or entry + 1 == total_entries:
                print(f"Processed {entry + 1}/{total_entries} events for primary {label} scan.")

        if original_entry is not None and 0 <= original_entry < event_obj.nEntry:
            event_obj.Jump(original_entry)

        if not energies:
            print(f"No primary {label}s with energy deposits were found in the scanned events.")
            return np.asarray(energies), np.asarray(lengths)

        energy_array = np.asarray(energies, dtype=float)
        length_array = np.asarray(lengths, dtype=float)
        ratio_array = np.divide(
            energy_array,
            length_array,
            out=np.zeros_like(energy_array, dtype=float),
            where=length_array > 0.0
        )

        energy_min, energy_max = float(energy_array.min()), float(energy_array.max())
        length_min, length_max = float(length_array.min()), float(length_array.max())
        ratio_min, ratio_max = float(ratio_array.min()), float(ratio_array.max())

        def _range_with_padding(vmin, vmax):
            if vmax <= vmin:
                vmax = vmin + 1.0
            span = vmax - vmin
            return vmin - 0.05 * span, vmax + 0.05 * span

        e_low, e_high = _range_with_padding(energy_min, energy_max)
        l_low, l_high = _range_with_padding(length_min, length_max)
        r_low, r_high = _range_with_padding(ratio_min, ratio_max)

        label_tag = label.replace(' ', '_')
        hist_ke_name = f'h_{label_tag}_primary_ke'
        hist_len_name = f'h_{label_tag}_primary_span'
        hist_ratio_name = f'h_{label_tag}_primary_ke_over_span'

        energy_hist = ROOT.TH1F(
            hist_ke_name,
            f"Primary {label} kinetic energy;Kinetic energy [MeV];Counts",
            bins,
            e_low,
            e_high
        )
        length_hist = ROOT.TH1F(
            hist_len_name,
            f"Primary {label} energy-deposition span;Span length [cm];Counts",
            bins,
            l_low,
            l_high
        )
        ratio_hist = ROOT.TH1F(
            hist_ratio_name,
            f"Primary {label} energy/length;Energy/length [MeV/cm];Counts",
            bins,
            r_low,
            r_high
        )

        for energy, length, ratio in zip(energy_array, length_array, ratio_array):
            energy_hist.Fill(float(energy))
            length_hist.Fill(float(length))
            ratio_hist.Fill(float(ratio))

        if output_path is None:
            output_dir = Path(self.event.plotpath)
            output_dir.mkdir(parents=True, exist_ok=True)
            output_file = output_dir / f'primary_{label_tag}_energy_length.root'
        else:
            output_file = Path(output_path).expanduser()
            output_file.parent.mkdir(parents=True, exist_ok=True)

        canvas = ROOT.TCanvas(
            f'c_{label_tag}_energy_length',
            f'Primary {label} energy and length',
            1500,
            500
        )
        canvas.Divide(3, 1)

        canvas.cd(1)
        energy_hist.SetLineColor(ROOT.kAzure + 2)
        energy_hist.Draw('HIST')

        canvas.cd(2)
        length_hist.SetLineColor(ROOT.kOrange + 7)
        length_hist.Draw('HIST')

        canvas.cd(3)
        ratio_hist.SetLineColor(ROOT.kGreen + 2)
        ratio_hist.Draw('HIST')

        canvas.Update()

        root_file = ROOT.TFile(str(output_file), 'RECREATE')
        energy_hist.Write()
        length_hist.Write()
        ratio_hist.Write()
        canvas.Write()
        root_file.Close()

        print(f"Saved {label} energy/length histograms to {output_file}")
        return energy_array, length_array

    def proton_energy_length_histograms(
        self,
        file_path="~/frunshiworking/data_antinumu/antinumu/edep_antinumu_merged.root",
        max_events=9000,
        bins=300,
        output_path=None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Build ROOT histograms for primary proton kinetic energy, span length, and energy/length."""
        return self._primary_particle_energy_length_histograms(
            pdg_codes=(2212,),
            label='proton',
            file_path=file_path,
            max_events=max_events,
            bins=bins,
            output_path=output_path,
        )


    def pion_energy_length_histograms(
        self,
        file_path="~/frunshiworking/data_antinumu/antinumu/edep_antinumu_merged.root",
        max_events=9000,
        bins=60,
        output_path=None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Build ROOT histograms for primary charged-pion kinetic energy, span length, and energy/length.

        Both pi+ and pi- emerging from the neutrino vertex are included."""
        return self._primary_particle_energy_length_histograms(
            pdg_codes=(211, 111, -211),
            label='pion',
            file_path=file_path,
            max_events=max_events,
            bins=bins,
            output_path=output_path,
        )


    def electron_energy_length_histograms(
        self,
        file_path="~/frunshiworking/data_nue/newevent_nue/edep_nue_all.root",
        max_events=9000,
        bins=60,
        output_path=None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Build ROOT histograms for primary electron/positron kinetic energy, span length, and energy/length."""
        return self._primary_particle_energy_length_histograms(
            pdg_codes=(11, -11),
            label='electron',
            file_path=file_path,
            max_events=max_events,
            bins=bins,
            output_path=output_path,
        )

    def draw_electron_shower(self, entry, energy_threshold=0.0):
        """Plot Y-Z positions of electron deposits, colored by subprocess."""

        if entry is not None:
            self.event.Jump(entry)

        electron_track, track_index, descendant_indices = self.event.get_primary_electron_track()
        if electron_track is None or not descendant_indices:
            print('draw_electron_shower: primary electron track not found for this event.')
            return

        mm2cm = 0.1
        subprocess_labels = {
            '2': 'Ionization',
            '3': 'Bremsstrahlung',
            '22': 'Scintillation',
            '401': 'StepLimiter',
            '111': 'HadronElastic',
        }

        # Cache points per descendant track
        track_point_cache = {}
        for track_idx in descendant_indices:
            track = self.event.tracks[track_idx]
            points = []
            coords = []
            labels = []
            for point in track.Points:
                pos = point.GetPosition()
                vec = np.array([pos.X() * mm2cm, pos.Y() * mm2cm, pos.Z() * mm2cm])
                subprocess = point.GetSubprocess() or 'UNKNOWN'
                if isinstance(subprocess, bytes):
                    subprocess = subprocess.decode('utf-8', errors='ignore')
                if not isinstance(subprocess, str):
                    subprocess_code = str(subprocess)
                else:
                    subprocess_code = subprocess
                label = subprocess_labels.get(subprocess_code, 'Other')
                coords.append(vec)
                labels.append(label)
                points.append({'coords': vec, 'label': label})
            if coords:
                track_point_cache[track_idx] = {
                    'coords': np.vstack(coords),
                    'labels': labels,
                }

        deposits_by_label = {}
        max_energy = 0.0
        all_y = []
        all_z = []

        has_deposits = False
        for track_idx in descendant_indices:
            track = self.event.tracks[track_idx]
            depo_list = track.association.get('depoList', [])
            if not depo_list:
                continue
            cache = track_point_cache.get(track_idx)
            for di in depo_list:
                if di < 0 or di >= len(self.event.depos):
                    continue
                depo = self.event.depos[di]
                energy = float(depo.GetEnergyDeposit())
                if energy <= energy_threshold:
                    continue
                start, stop = depo.GetStart(), depo.GetStop()
                midpoint = np.array([
                    (start.X() + stop.X()) / 2 * mm2cm,
                    (start.Y() + stop.Y()) / 2 * mm2cm,
                    (start.Z() + stop.Z()) / 2 * mm2cm,
                ])

                label = 'Other'
                if cache:
                    distances = np.linalg.norm(cache['coords'] - midpoint, axis=1)
                    if distances.size:
                        idx_min = int(np.argmin(distances))
                        label = cache['labels'][idx_min]

                bucket = deposits_by_label.setdefault(label, {'z': [], 'y': [], 'e': []})
                bucket['z'].append(midpoint[2])
                bucket['y'].append(midpoint[1])
                bucket['e'].append(energy)
                all_y.append(midpoint[1])
                all_z.append(midpoint[2])
                max_energy = max(max_energy, energy)
                has_deposits = True

        if not has_deposits:
            print('draw_electron_shower: no electron deposits passed the selection.')
            return

        all_y = np.asarray(all_y)
        all_z = np.asarray(all_z)
        y_margin = max((all_y.max() - all_y.min()) * 0.1, 5.0)
        z_margin = max((all_z.max() - all_z.min()) * 0.1, 5.0)

        canvas = ROOT.TCanvas(
            f"electron_shower_{self.event.currentEntry}",
            f"Electron shower event {self.event.currentEntry}",
            900,
            700,
        )

        frame = ROOT.TH2F(
            'electron_shower_frame',
            f"Electron deposits (entry {self.event.currentEntry})",
            100,
            all_z.min() - z_margin,
            all_z.max() + z_margin,
            100,
            all_y.min() - y_margin,
            all_y.max() + y_margin,
        )
        frame.GetXaxis().SetTitle('Z [cm]')
        frame.GetYaxis().SetTitle('Y [cm]')
        frame.Draw()

        color_cycle = {
            'Ionization': ROOT.kRed + 1,
            'Bremsstrahlung': ROOT.kBlue + 1,
            'Scintillation': ROOT.kGreen + 2,
            'StepLimiter': ROOT.kMagenta + 2,
            'HadronElastic': ROOT.kOrange + 7,
            'Other': ROOT.kGray + 2,
        }

        graphs = []
        for idx, (label, data) in enumerate(deposits_by_label.items()):
            graph = ROOT.TGraph(
                len(data['z']),
                array.array('d', data['z']),
                array.array('d', data['y']),
            )
            color = color_cycle.get(label, color_cycle['Other'])
            graph.SetMarkerStyle(20)
            graph.SetMarkerColor(color)
            graph.SetMarkerSize(0.5)
            graph.Draw('P SAME')
            graphs.append((graph, label, data))

        legend = ROOT.TLegend(0.12, 0.72, 0.45, 0.9)
        legend_markers = []
        for graph, label, data in graphs:
            legend_label = label or 'Other'
            proxy_marker = ROOT.TMarker(0, 0, 20)
            proxy_marker.SetMarkerColor(graph.GetMarkerColor())
            proxy_marker.SetMarkerSize(max(graph.GetMarkerSize() * 1.5, 0.6))
            legend_markers.append(proxy_marker)
            legend.AddEntry(
                proxy_marker,
                f"{legend_label} ({len(data['z'])} hits)",
                'p'
            )
        legend.Draw()

        OUTPUT_DIR = Path('root_plot')
        OUTPUT_DIR.mkdir(exist_ok=True)
        output_path = f"~/frunshiworking/root_plot/electron_shower_0.8GeVantinue/electron_shower_event_{self.event.currentEntry}.root"
        canvas.SaveAs(str(output_path))
        print(f"draw_electron_shower: saved to {output_path}")



    def neutron_E_E_plot(self, max_events=None):
        """Scatter plot of neutron kinetic energy at the vertex vs. summed deposits.

        Parameters
        ----------
        max_events : int, optional
            Maximum number of events to scan. Defaults to every entry in the
            underlying `Event` instance.
        """

        if not hasattr(self.event, 'nEntry'):
            print('neutron_E_E_plot: event object has no entries to iterate.')
            return

        total_entries = getattr(self.event, 'nEntry', 0)
        if total_entries <= 0:
            print('neutron_E_E_plot: no events available.')
            return

        limit = min(total_entries, max_events) if max_events is not None else total_entries
        # limit =1000
        print(f"the number of events we are dealing with: {limit}")
        neutron_ke = []
        summed_edep = []

        for entry in range(limit):
            self.event.Jump(entry)
            if entry%100==0:
                print(f"event: {entry}")
            energy_sum, _ = self.event.sum_neutron_energy_deposition()
            if energy_sum > 150:
                continue  # exclude smeared events with unrealistically large summed deposits
            # print(f"sum of the recognized neutron energies: {energy_sum}")

            ke_candidates = []
            vertex = getattr(self.event, 'vertex', None)
            particles = getattr(vertex, 'Particles', []) if vertex is not None else []
            for particle in particles:
                pdg = particle.GetPDGCode()
                momentum = particle.GetMomentum()
                if momentum is None:
                    continue

                try:
                    kinetic_energy = float(momentum.E() - momentum.M())
                except TypeError:
                    continue

                if pdg == 2112:
                    ke_candidates.append(kinetic_energy)

            ke_at_vertex = sum(ke_candidates) if ke_candidates else 0.0

            neutron_ke.append(ke_at_vertex)
            summed_edep.append(energy_sum)

        if not neutron_ke:
            print('neutron_E_E_plot: no neutron entries collected; nothing to draw.')
            return

        graph = ROOT.TGraph(
            len(neutron_ke),
            array.array('d', neutron_ke),
            array.array('d', summed_edep),
        )
        graph.SetName('g_neutron_ke_vs_edep')
        graph.SetTitle(f'Neutron KE vs summed deposits {limit}events(numu Flatten spectrum);Neutron initial kinetic energy (MeV);Summed neutron energy deposits (MeV)')
        graph.SetMarkerStyle(20)
        graph.SetMarkerSize(0.2)
        graph.SetMarkerColor(ROOT.kBlue + 2)

        canvas = ROOT.TCanvas('c_neutron_ke_vs_edep', 'Neutron KE scatter', 900, 700)
        graph.Draw('AP')

        fit_range_min = min(neutron_ke)
        fit_range_max = max(neutron_ke)
        if fit_range_max > fit_range_min:
            fit_func = TF1('f_lin_neutron_ke_vs_edep', '[0] + [1]*x', fit_range_min, fit_range_max)
            fit_func.SetLineColor(ROOT.kRed + 1)
            fit_func.SetLineWidth(2)
            fit_func.FixParameter(0, 0)
            graph.Fit(fit_func, 'Q')
            intercept = fit_func.GetParameter(0)
            slope = fit_func.GetParameter(1)
            label = ROOT.TLatex(0.15, 0.85, f'y = {slope:.2f}x')
            label.SetNDC(True)
            label.SetTextSize(0.03)
            label.SetTextColor(ROOT.kRed + 1)
            label.Draw()
            print(f'neutron_E_E_plot: linear fit slope={slope:.4f} intercept={intercept:.4f}')
        else:
            print('neutron_E_E_plot: insufficient spread in KE values for linear fit.')

        output_dir = Path('1.5meter_reconstruction')
        output_dir.mkdir(exist_ok=True)
        output_path = output_dir / 'neutron_ke_vs_edep_1.5m_4cm_composite_neutron_0.6thresh_antinumu.root'
        root_file = ROOT.TFile(str(output_path), 'RECREATE')
        graph.Write()
        canvas.Write()
        root_file.Close()

        print(f'neutron_E_E_plot: saved scatter plot to {output_path}')


    def neutron_E_E_plot_2d(self, max_events=10000):
            """Single TH2Poly: 10x10 MeV bins for Y<50, 25x25 MeV bins for Y>=50.
            Curve: peak below 50, weighted mean above 50. Red line at Y=50.
            """

            if not hasattr(self.event, 'nEntry'):
                print('event object has no entries.')
                return

            total_entries = getattr(self.event, 'nEntry', 0)
            limit = min(total_entries, max_events) if max_events is not None else total_entries

            # TH2Poly for display: variable bin sizes per region
            poly = ROOT.TH2Poly(
                'h2_poly',
                'Neutron Energy Correlation(#bar{#nu}_{#mu});'
                'Initial Kinetic Energy (MeV);Summed Deposits (MeV)',
                0, 600, 0, 675,
            )
            # Low region: 60 x 5 = 300 bins of 10x10 MeV  (Y: 0-50)
            for yi in range(5):
                y1, y2 = yi * 10, (yi + 1) * 10
                for xi in range(60):
                    poly.AddBin(xi * 10, y1, (xi + 1) * 10, y2)
            # High region: 24 x 25 = 600 bins of 25x25 MeV (Y: 50-675)
            for yi in range(25):
                y1, y2 = 50 + yi * 25, 50 + (yi + 1) * 25
                for xi in range(24):
                    poly.AddBin(xi * 25, y1, (xi + 1) * 25, y2)

            # TH2D helpers kept for curve computation (need GetBinContent(xbin,ybin))
            hist_low = ROOT.TH2D(
                'h2_low_curve', '', 60, 0, 600, 5, 0, 50,
            )
            x_edges_high = array.array('d', [i * 25 for i in range(25)])
            y_edges_high = array.array('d', [50 + i * 25 for i in range(26)])
            hist_high = ROOT.TH2D(
                'h2_high_curve', '',
                len(x_edges_high) - 1, x_edges_high,
                len(y_edges_high) - 1, y_edges_high,
            )

            training_pairs = []
            for entry in range(limit):
                self.event.Jump(entry)
                if entry % 100 == 0:
                    print("processing entry", entry)
                energy_sum, _ = self.event.sum_neutron_energy_deposition()

                ke_candidates = []
                vertex = getattr(self.event, 'vertex', None)
                particles = getattr(vertex, 'Particles', []) if vertex is not None else []
                for particle in particles:
                    if particle.GetPDGCode() == 2112:
                        momentum = particle.GetMomentum()
                        if momentum:
                            ke_candidates.append(momentum.E() - momentum.M())
                ke_at_vertex = sum(ke_candidates) if ke_candidates else 0.0

                poly.Fill(ke_at_vertex, energy_sum)
                if energy_sum < 50.0:
                    hist_low.Fill(ke_at_vertex, energy_sum)
                else:
                    hist_high.Fill(ke_at_vertex, energy_sum)
                training_pairs.append((energy_sum, ke_at_vertex))

            ROOT.gStyle.SetOptStat(0)
            ROOT.gStyle.SetPalette(ROOT.kBird)

            output_dir = Path('paper_work')
            output_dir.mkdir(exist_ok=True)

            # ── Draw ──
            canvas = ROOT.TCanvas('c_poly', 'c_poly', 900, 800)
            canvas.SetRightMargin(0.15)
            poly.SetMinimum(0.5)
            poly.GetZaxis().SetTitle('Number of Events per Bin')
            poly.GetZaxis().SetTitleOffset(1.2)
            poly.Draw('COLZ')
            canvas.SetLogz(1)
            canvas.Update()

            # Red line at Y = 50
            red_line = ROOT.TLine(0, 50.0, 600, 50.0)
            red_line.SetLineColor(ROOT.kRed)
            red_line.SetLineWidth(2)
            red_line.Draw('SAME')

            # Curve helpers
            def _make_curve(hist, compute_peak):
                nx = hist.GetNbinsX()
                ny = hist.GetNbinsY()
                pts = []
                for ybin in range(1, ny + 1):
                    yc = hist.GetYaxis().GetBinCenter(ybin)
                    if compute_peak:
                        best_c, best_x = 0.0, None
                        for xbin in range(1, nx + 1):
                            xc = hist.GetXaxis().GetBinCenter(xbin)
                            if xc <= 10:
                                continue
                            c = hist.GetBinContent(xbin, ybin)
                            if c > best_c:
                                best_c, best_x = c, xc
                        if best_x is not None:
                            pts.append((best_x, yc))
                    else:
                        w_sum, w_tot = 0.0, 0.0
                        for xbin in range(1, nx + 1):
                            xc = hist.GetXaxis().GetBinCenter(xbin)
                            if xc <= yc:
                                continue
                            c = hist.GetBinContent(xbin, ybin)
                            if c <= 0:
                                continue
                            w_sum += c * xc
                            w_tot += c
                        if w_tot > 0:
                            pts.append((w_sum / w_tot, yc))
                return pts

            def _draw_graph(pts, name):
                if not pts:
                    return None
                g = ROOT.TGraph(len(pts))
                for i, (x, y) in enumerate(pts):
                    g.SetPoint(i, x, y)
                g.SetName(name)
                g.SetLineColor(ROOT.kOrange + 7)
                g.SetLineWidth(3)
                g.SetLineStyle(2)
                g.SetMarkerStyle(21)
                g.SetMarkerColor(ROOT.kOrange + 7)
                g.SetMarkerSize(0.7)
                g.Draw('LP SAME')
                return g

            curve_low  = _make_curve(hist_low,  compute_peak=True)
            curve_high = _make_curve(hist_high, compute_peak=False)
            graph_low  = _draw_graph(curve_low,  'g_curve_low')
            graph_high = _draw_graph(curve_high, 'g_curve_high')

            canvas.Update()
            out_stem  = 'neutron_2D_poly_antinumu_0.5'
            png_path  = output_dir / f'{out_stem}.png'
            root_path = output_dir / f'{out_stem}.root'
            canvas.SaveAs(str(png_path))
            rf = ROOT.TFile(str(root_path), 'RECREATE')
            poly.Write()
            if graph_low  is not None: graph_low.Write()
            if graph_high is not None: graph_high.Write()
            rf.Close()
            print(f'Saved {png_path} and {root_path}')

            refs_low  = (curve_low,  red_line, graph_low)
            refs_high = (curve_high, red_line, graph_high)

            # Save peak curve points to text
            for refs, stem in ((refs_low, 'low'), (refs_high, 'high')):
                pts = refs[0]
                if pts:
                    txt_path = output_dir / f'neutron_peak_curve_{stem}_antinumu_1.5_meter.txt'
                    with txt_path.open('w', encoding='ascii') as fh:
                        fh.write('# Initial Kinetic Energy (MeV)\tSummed Deposits (MeV)\n')
                        for x_val, y_val in pts:
                            fh.write(f'{x_val:.6f}\t{y_val:.6f}\n')
                    print(f'Saved curve points to {txt_path}')

            # Save training data
            training_dir = Path('training_data')
            training_dir.mkdir(exist_ok=True)
            training_path = training_dir / 'neutron_summed_vs_ke_antinumu.csv'
            with training_path.open('w', encoding='ascii') as fh:
                fh.write('summed_energy_MeV,initial_ke_MeV\n')
                for e_sum, ke in training_pairs:
                    fh.write(f'{e_sum:.6f},{ke:.6f}\n')
            print(f'Saved training data to {training_path}')

    def neutron_blip_ancestry_breakdown(self, max_events=10000):
        """Apply the same spatial/cone/track filter as sum_neutron_energy_deposition.

        Reports:
          1. Blip count and deposited energy broken down by the track's own PDG type
             (muon, proton, pion, neutron, gamma, electron).
          2. Fraction of blips (count and energy) whose ultimate ancestor or direct
             parent is a neutron.
        """

        if not hasattr(self.event, 'nEntry'):
            print('neutron_blip_ancestry_breakdown: event object has no entries.')
            return

        total_entries = getattr(self.event, 'nEntry', 0)
        limit = min(total_entries, max_events) if max_events is not None else total_entries

        def _category(pdg):
            if pdg in (13, -13):       return 'muon'
            if pdg == 2212:            return 'proton'
            if pdg in (211, -211, 111): return 'pion'
            if pdg == 2112:            return 'neutron'
            if pdg == 22:              return 'gamma'
            if pdg in (11, -11):       return 'electron'
            return 'other'

        categories = ['muon', 'proton', 'pion', 'neutron', 'gamma', 'electron', 'other']
        blip_count       = {c: 0   for c in categories}
        blip_energy      = {c: 0.0 for c in categories}
        anc_blip_count   = {c: 0   for c in categories}
        anc_blip_energy  = {c: 0.0 for c in categories}
        other_anc_pdg_count  = {}   # raw PDG → blip count  for "other" ancestors
        other_anc_pdg_energy = {}   # raw PDG → energy sum  for "other" ancestors

        total_blips = 0
        total_energy = 0.0
        neutron_ancestor_blips  = 0
        neutron_ancestor_energy = 0.0
        neutron_parent_blips    = 0
        neutron_parent_energy   = 0.0
        mm2cm = 0.1

        for entry in range(limit):
            self.event.Jump(entry)
            if entry % 100 == 0:
                print(f'processing entry {entry}')

            ev = self.event
            muon_energy = ev.read_muon_energy()
            muon_dir, _ = ev.read_muon_direction_and_amplitude()
            cone = ev.pairing_muon_energy_with_muon_cone(muon_energy)
            threshold = ev.energy_threshold

            track_proton_pion = set()
            for particle in ev.vertex.Particles:
                pdg = particle.GetPDGCode()
                if pdg == 2212:
                    mom = particle.GetMomentum()
                    if mom and (mom.E() - mom.M()) > ev.proton_trackable_length:
                        track_proton_pion.add(particle.GetTrackId())
                elif pdg in (111, 211, -211):
                    track_proton_pion.add(particle.GetTrackId())

            for track in ev.tracks:
                trackid  = track.GetTrackId()
                pdg_code = track.GetPDGCode()
                if pdg_code not in [22, 2112, 2212, 211, -211, 111, 11, -11, 13, -13]:
                    continue
                if trackid in track_proton_pion:
                    continue

                parent_id = track.GetParentId()
                if parent_id != -1 and ev.tracks[parent_id].GetTrackId() in track_proton_pion:
                    continue

                cat = _category(pdg_code)
                ancestor_idx = track.association.get('ancestor', None)
                ancestor_pdg = ev.tracks[ancestor_idx].GetPDGCode() if ancestor_idx is not None else pdg_code
                parent_pdg   = ev.tracks[parent_id].GetPDGCode() if parent_id != -1 else None

                for di in track.association.get('depoList', []):
                    if not (0 <= di < len(ev.depos)):
                        continue
                    depo = ev.depos[di]
                    edep = depo.GetEnergyDeposit()
                    if edep < threshold:
                        continue
                    start, stop = depo.GetStart(), depo.GetStop()
                    midpoint = np.array([
                        (start.X() + stop.X()) / 2 * mm2cm,
                        (start.Y() + stop.Y()) / 2 * mm2cm,
                        (start.Z() + stop.Z()) / 2 * mm2cm,
                    ])
                    midpointlength = np.linalg.norm(midpoint)
                    angle = np.degrees(ev.compute_angle_between(muon_dir, midpoint))
                    if (angle < cone
                            or midpointlength > ev.inspection_length_energy_side * 100
                            or midpointlength < ev.minimum_length):
                        continue

                    total_blips  += 1
                    total_energy += edep
                    blip_count[cat]  += 1
                    blip_energy[cat] += edep
                    anc_cat = _category(ancestor_pdg)
                    anc_blip_count[anc_cat]  += 1
                    anc_blip_energy[anc_cat] += edep
                    if anc_cat == 'other':
                        other_anc_pdg_count[ancestor_pdg]  = other_anc_pdg_count.get(ancestor_pdg, 0) + 1
                        other_anc_pdg_energy[ancestor_pdg] = other_anc_pdg_energy.get(ancestor_pdg, 0.0) + edep

                    if ancestor_pdg == 2112:
                        neutron_ancestor_blips  += 1
                        neutron_ancestor_energy += edep
                    if parent_pdg == 2112:
                        neutron_parent_blips  += 1
                        neutron_parent_energy += edep

        print(f'\n--- Blip ancestry breakdown ({limit} events) ---')
        print(f'Total blips passing filter : {total_blips}   total energy: {total_energy:.2f} MeV')
        if total_blips == 0:
            return

        header = f'{"Type":<10} {"Blips":>8} {"Blip %":>8} {"Energy(MeV)":>13} {"Energy %":>9}'
        sep    = '-' * 52

        print(f'\nBy track particle type:')
        print(header)
        print(sep)
        for cat in categories:
            n = blip_count[cat]
            e = blip_energy[cat]
            print(f'{cat:<10} {n:>8d} {100*n/total_blips:>7.1f}% {e:>13.2f} {100*e/total_energy:>8.1f}%')

        print(f'\nBy ultimate ancestor particle type:')
        print(header)
        print(sep)
        for cat in categories:
            n = anc_blip_count[cat]
            e = anc_blip_energy[cat]
            print(f'{cat:<10} {n:>8d} {100*n/total_blips:>7.1f}% {e:>13.2f} {100*e/total_energy:>8.1f}%')

        if other_anc_pdg_count:
            print(f'\n  "other" ancestor breakdown by raw PDG code:')
            print(f'  {"PDG":>10} {"Blips":>8} {"Blip %":>8} {"Energy(MeV)":>13} {"Energy %":>9}')
            print('  ' + '-' * 52)
            for pdg, n in sorted(other_anc_pdg_count.items(), key=lambda x: -x[1]):
                e = other_anc_pdg_energy[pdg]
                print(f'  {pdg:>10} {n:>8d} {100*n/total_blips:>7.1f}% {e:>13.2f} {100*e/total_energy:>8.1f}%')

        print(f'\n{"Ancestry cut":<26} {"Blips":>8} {"Blip %":>8} {"Energy(MeV)":>13} {"Energy %":>9}')
        print('-' * 68)
        print(f'{"Neutron as ancestor":<26} {neutron_ancestor_blips:>8d} {100*neutron_ancestor_blips/total_blips:>7.1f}% {neutron_ancestor_energy:>13.2f} {100*neutron_ancestor_energy/total_energy:>8.1f}%')
        print(f'{"Neutron as direct parent":<26} {neutron_parent_blips:>8d} {100*neutron_parent_blips/total_blips:>7.1f}% {neutron_parent_energy:>13.2f} {100*neutron_parent_energy/total_energy:>8.1f}%')

    def neutron_E_E_plot_2d_uniform(self, max_events=10000):
        """Single 10x10 MeV 2D histogram of neutron KE vs summed deposits,
        full Y range, with pion/proton ancestor exclusion applied via
        sum_neutron_energy_deposition."""

        if not hasattr(self.event, 'nEntry'):
            print('event object has no entries.')
            return

        total_entries = getattr(self.event, 'nEntry', 0)
        limit = min(total_entries, max_events) if max_events is not None else total_entries

        # X: initial neutron KE  0-600 MeV, 60 bins x 10 MeV
        # Y: summed deposits     0-700 MeV, 70 bins x 10 MeV
        hist = ROOT.TH2D(
            'h2_uniform',
            'Neutron Energy Correlation(#bar{#nu}_{#mu});'
            'Initial Kinetic Energy (MeV);Summed Deposits (MeV);Counts',
            60, 0, 600,
            70, 0, 700,
        )

        training_pairs = []
        for entry in range(limit):
            self.event.Jump(entry)
            if entry % 100 == 0:
                print(f'processing entry {entry}')
            energy_sum, _ = self.event.sum_neutron_energy_deposition()

            ke_candidates = []
            vertex = getattr(self.event, 'vertex', None)
            particles = getattr(vertex, 'Particles', []) if vertex is not None else []
            for particle in particles:
                if particle.GetPDGCode() == 2112:
                    momentum = particle.GetMomentum()
                    if momentum:
                        ke_candidates.append(momentum.E() - momentum.M())
            ke_at_vertex = sum(ke_candidates) if ke_candidates else 0.0

            hist.Fill(ke_at_vertex, energy_sum)
            training_pairs.append((energy_sum, ke_at_vertex))

        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetPalette(ROOT.kBird)

        canvas = ROOT.TCanvas('c_uniform', 'c_uniform', 900, 800)
        canvas.SetRightMargin(0.15)
        hist.GetZaxis().SetTitle('Number of Events per Bin')
        hist.GetZaxis().SetTitleOffset(1.2)
        hist.SetMinimum(0.5)
        hist.Draw('COLZ')
        canvas.SetLogz(1)
        canvas.Update()

        # red line at Y = 50 MeV marking the low/high boundary
        x_min = hist.GetXaxis().GetXmin()
        x_max = hist.GetXaxis().GetXmax()
        red_line = ROOT.TLine(x_min, 50.0, x_max, 50.0)
        red_line.SetLineColor(ROOT.kRed)
        red_line.SetLineWidth(2)
        red_line.Draw('SAME')

        # below Y=50: peak (mode) per Y bin; at/above Y=50: weighted mean per Y bin
        nbins_x = hist.GetNbinsX()
        nbins_y = hist.GetNbinsY()
        curve_points = []
        for ybin in range(1, nbins_y + 1):
            y_center = hist.GetYaxis().GetBinCenter(ybin)
            if y_center < 50.0:
                max_content, peak_x = 0.0, None
                for xbin in range(1, nbins_x + 1):
                    x_center = hist.GetXaxis().GetBinCenter(xbin)
                    if x_center <= 10:
                        continue
                    content = hist.GetBinContent(xbin, ybin)
                    if content > max_content:
                        max_content = content
                        peak_x = x_center
                if peak_x is not None:
                    curve_points.append((peak_x, y_center))
            else:
                w_sum, w_tot = 0.0, 0.0
                for xbin in range(1, nbins_x + 1):
                    x_center = hist.GetXaxis().GetBinCenter(xbin)
                    if x_center <= y_center:
                        continue
                    content = hist.GetBinContent(xbin, ybin)
                    if content <= 0:
                        continue
                    w_sum += content * x_center
                    w_tot += content
                if w_tot > 0:
                    curve_points.append((w_sum / w_tot, y_center))

        graph = None
        if curve_points:
            graph = ROOT.TGraph(len(curve_points))
            for idx, (x_val, y_val) in enumerate(curve_points):
                graph.SetPoint(idx, x_val, y_val)
            graph.SetName('g_curve_uniform')
            graph.SetLineColor(ROOT.kOrange + 7)
            graph.SetLineWidth(3)
            graph.SetLineStyle(2)
            graph.SetMarkerStyle(21)
            graph.SetMarkerColor(ROOT.kOrange + 7)
            graph.SetMarkerSize(0.7)
            graph.Draw('LP SAME')

        canvas.Update()

        output_dir = Path('paper_work')
        output_dir.mkdir(exist_ok=True)
        out_stem = 'neutron_2D_uniform_10x10_antinumu'

        png_path  = output_dir / f'{out_stem}.png'
        root_path = output_dir / f'{out_stem}.root'
        canvas.SaveAs(str(png_path))
        rf = ROOT.TFile(str(root_path), 'RECREATE')
        hist.Write()
        if graph is not None:
            graph.Write()
        rf.Close()
        print(f'Saved {png_path} and {root_path}')

        # curve points to text
        if curve_points:
            txt_path = output_dir / f'neutron_peak_curve_uniform_antinumu.txt'
            with txt_path.open('w', encoding='ascii') as fh:
                fh.write('# Initial Kinetic Energy (MeV)\tSummed Deposits (MeV)\n')
                for x_val, y_val in curve_points:
                    fh.write(f'{x_val:.6f}\t{y_val:.6f}\n')
            print(f'Saved curve points to {txt_path}')

        # training data
        training_dir = Path('training_data')
        training_dir.mkdir(exist_ok=True)
        training_path = training_dir / 'neutron_summed_vs_ke_uniform_antinumu.csv'
        with training_path.open('w', encoding='ascii') as fh:
            fh.write('summed_energy_MeV,initial_ke_MeV\n')
            for e_sum, ke in training_pairs:
                fh.write(f'{e_sum:.6f},{ke:.6f}\n')
        print(f'Saved training data to {training_path}')

    def neutron_e_e_one_plus_one_dimensional(self, low_sum=20.0, high_sum=30.0, max_events=10000):
            """Histogram neutron KE for events whose summed deposits fall into [low_sum, high_sum]."""

            if not hasattr(self.event, 'nEntry'):
                print('event object has no entries.')
                return

            total_entries = getattr(self.event, 'nEntry', 0)
            limit = min(total_entries, max_events) if max_events is not None else total_entries

            hist = ROOT.TH1D(
                'h_neutron_ke_filtered',
                f'Neutron KE for {low_sum}-{high_sum} MeV deposits;Initial Kinetic Energy (MeV);Events',
                100, 0, 600
            )

            filled_events = 0
            for entry in range(limit):
                self.event.Jump(entry)
                if entry % 100 == 0:
                    print('processing entry', entry)

                energy_sum, _ = self.event.sum_neutron_energy_deposition()
                if energy_sum < low_sum or energy_sum > high_sum:
                    continue

                ke_candidates = []
                vertex = getattr(self.event, 'vertex', None)
                particles = getattr(vertex, 'Particles', []) if vertex is not None else []
                for particle in particles:
                    if particle.GetPDGCode() == 2112:
                        momentum = particle.GetMomentum()
                        if momentum:
                            ke_candidates.append(momentum.E() - momentum.M())

                ke_at_vertex = sum(ke_candidates) if ke_candidates else 0.0
                hist.Fill(ke_at_vertex)
                filled_events += 1

            if filled_events == 0:
                print(f'neutron_e_e_one_plus_one_dimensional: no events with summed deposits in {low_sum}-{high_sum} MeV.')
                return

            canvas = ROOT.TCanvas('c_neutron_ke_filtered', 'Neutron KE Histogram', 900, 700)
            ROOT.gStyle.SetOptStat(1110)
            hist.SetLineColor(ROOT.kAzure + 2)
            hist.SetLineWidth(2)
            hist.Draw('HIST')
            canvas.Update()

            output_dir = Path('1.5meter_reconstruction_energy_reco')
            output_dir.mkdir(exist_ok=True)
            basename = f'neutron_ke_hist_{int(low_sum)}to{int(high_sum)}_numu'
            output_path = output_dir / f'{basename}.root'
            root_file = ROOT.TFile(str(output_path), 'RECREATE')
            hist.Write()
            canvas.Write()
            canvas.SaveAs(str(output_path).replace('.root', '.png'))
            root_file.Close()

            print(f'Saved neutron KE histogram (filtered {low_sum}-{high_sum} MeV deposits) to {output_path}')

    def neutron_E_E_plot_2d_EV(self,max_events = 10000):
            """2D Histogram of neutron kinetic energy vs. summed deposits."""

            if not hasattr(self.event, 'nEntry'):
                print('event object has no entries.')
                return

            total_entries = getattr(self.event, 'nEntry', 0)
            limit = min(total_entries, max_events) if max_events is not None else total_entries

            # 1. Initialize the 2D Histogram (TH2D)
            # Bins, X-min, X-max, Bins, Y-min, Y-max
            hist2d = ROOT.TH2D(
                'h2_neutron_ke_vs_edep',
                'Neutron Energy Correlation(nue);Initial Kinetic Energy (MeV);Summed Deposits (MeV);Counts',
                60, 0, 600,  # 100 bins for KE from 0 to 150
                60, 0, 600  # 100 bins for Edep from 0 to 150
            )

            for entry in range(limit):
                self.event.Jump(entry)
                if entry%100==0:
                    print("processing entry",entry)
                energy_sum, _ = self.event.sum_neutron_energy_deposition_EV()

                ke_candidates = []
                vertex = getattr(self.event, 'vertex', None)
                particles = getattr(vertex, 'Particles', []) if vertex is not None else []
                for particle in particles:
                    if particle.GetPDGCode() == 2112:
                        momentum = particle.GetMomentum()
                        if momentum:
                            ke_candidates.append(momentum.E() - momentum.M())

                ke_at_vertex = sum(ke_candidates) if ke_candidates else 0.0
                # 2. Fill the 2D histogram
                hist2d.Fill(ke_at_vertex, energy_sum)



            # 3. Drawing and Styling
            canvas = ROOT.TCanvas('c_2d', '2D Histogram', 900, 800)
            # 1. Adjust the margins to make room for the color bar
            canvas.SetRightMargin(0.15)

            # 2. Add a title to the color bar (Z-axis)
            hist2d.GetZaxis().SetTitle("Number of Events per Bin")
            hist2d.GetZaxis().SetTitleOffset(1.2)

            # 4. Draw with the 'Z' option
            ROOT.gStyle.SetOptStat(0)  # Hide the stats box for a cleaner look
            ROOT.gStyle.SetPalette(ROOT.kBird)  # Modern color palette

            # 'COLZ' draws the 2D histogram with a color-coded frequency scale
            hist2d.SetMinimum(0.5)
            hist2d.Draw('COLZ')
            canvas.SetLogz(1)
            canvas.Update()

            # Determine peak (low sums) or mean (high sums) KE per summed deposit bin
            peak_points = []
            nbins_y = hist2d.GetNbinsY()
            nbins_x = hist2d.GetNbinsX()
            for ybin in range(1, nbins_y + 1):
                y_center = hist2d.GetYaxis().GetBinCenter(ybin)
                if y_center < 50.0:
                    max_content = 0.0
                    peak_x_center = None
                    for xbin in range(1, nbins_x + 1):
                        x_center = hist2d.GetXaxis().GetBinCenter(xbin)
                        if x_center <= y_center:
                            continue
                        content = hist2d.GetBinContent(xbin, ybin)
                        if content <= 0:
                            continue
                        if content > max_content:
                            max_content = content
                            peak_x_center = x_center
                    if peak_x_center is not None:
                        peak_points.append((peak_x_center, y_center))
                else:
                    weighted_sum = 0.0
                    weight_total = 0.0
                    for xbin in range(1, nbins_x + 1):
                        x_center = hist2d.GetXaxis().GetBinCenter(xbin)
                        if  x_center <= y_center:
                            continue
                        content = hist2d.GetBinContent(xbin, ybin)
                        if content <= 0:
                            continue
                        weighted_sum += content * x_center
                        weight_total += content
                    if weight_total > 0:
                        peak_points.append((weighted_sum / weight_total, y_center))

            # Draw peak curve
            peak_graph = None
            if peak_points:
                peak_graph = ROOT.TGraph(len(peak_points))
                for idx, (x_val, y_val) in enumerate(peak_points):
                    peak_graph.SetPoint(idx, x_val, y_val)
                peak_graph.SetName('g_neutron_peak_curve_nue')
                peak_graph.SetLineColor(ROOT.kOrange + 7)
                peak_graph.SetLineWidth(3)
                peak_graph.SetLineStyle(2)
                peak_graph.SetMarkerStyle(21)
                peak_graph.SetMarkerColor(ROOT.kOrange + 7)
                peak_graph.SetMarkerSize(0.7)
                peak_graph.Draw('LP SAME')
            else:
                print('neutron_E_E_plot_2d_EV: no peak points found; skipping peak curve.')

            # Save output
            output_dir = Path('paper_work')
            output_dir.mkdir(exist_ok=True)
            output_path = output_dir / 'neutron_2D_hist_nue_peak_0.5meter.root'

            root_file = ROOT.TFile(str(output_path), 'RECREATE')
            hist2d.Write()
            if peak_graph is not None:
                peak_graph.Write()
            canvas.SaveAs(str(output_path).replace('.root', '.png'))
            root_file.Close()

            if peak_points:
                peak_txt_path = output_dir / 'neutron_peak_curve_points_nue_peak_0.5meter.txt'
                with peak_txt_path.open('w', encoding='ascii') as fh:
                    fh.write('# Peak curve points: Initial Kinetic Energy (MeV)\tSummed Deposits (MeV)\n')
                    for x_val, y_val in peak_points:
                        fh.write(f'{x_val:.6f}\t{y_val:.6f}\n')
                print(f'Saved peak curve points to {peak_txt_path}')

            print(f'Saved 2D Histogram to {output_path}')










    def neutron_E_E_plot_EV(self,max_events = 10000):

        if not hasattr(self.event, 'nEntry'):
            print('neutron_E_E_plot_EV: event object has no entries to iterate.')
            return

        total_entries = getattr(self.event, 'nEntry', 0)
        if total_entries <= 0:
            print('neutron_E_E_plot_EV: no events available.')
            return

        limit = min(total_entries, max_events) if max_events is not None else total_entries
        # limit =1000
        print(f"the number of events we are dealing with: {limit}")
        neutron_ke = []
        summed_edep = []

        for entry in range(limit):
            self.event.Jump(entry)
            if entry%100==0:
                print(f"event: {entry}")
            energy_sum, _ = self.event.sum_neutron_energy_deposition_EV()
            if energy_sum > 150:
                continue  # exclude smeared events with unrealistically large summed deposits
            # print(f"sum of the recognized neutron energies: {energy_sum}")

            ke_candidates = []
            vertex = getattr(self.event, 'vertex', None)
            particles = getattr(vertex, 'Particles', []) if vertex is not None else []
            for particle in particles:
                pdg = particle.GetPDGCode()
                momentum = particle.GetMomentum()
                if momentum is None:
                    continue
                kinetic_energy = float(momentum.E() - momentum.M())
                if pdg == 2112:
                    ke_candidates.append(kinetic_energy)

            ke_total = sum(ke_candidates) if ke_candidates else 0.0


            neutron_ke.append(ke_total)
            summed_edep.append(energy_sum)

        if not neutron_ke:
            print('neutron_E_E_plot: no neutron entries collected; nothing to draw.')
            return

        graph = ROOT.TGraph(
            len(neutron_ke),
            array.array('d', neutron_ke),
            array.array('d', summed_edep),
        )
        graph.SetName('g_neutron_ke_vs_edep')
        graph.SetTitle(f'Neutron KE vs summed deposits {limit}events(nue Flatten spectrum) 2 meter radius;Neutron initial kinetic energy (MeV);Summed neutron energy deposits (MeV)')
        graph.SetMarkerStyle(20)
        graph.SetMarkerSize(0.2)
        graph.SetMarkerColor(ROOT.kBlue + 2)

        canvas = ROOT.TCanvas('c_neutron_ke_vs_edep', 'Neutron KE scatter', 900, 700)
        graph.Draw('AP')

        fit_range_min = min(neutron_ke)
        fit_range_max = max(neutron_ke)
        if fit_range_max > fit_range_min:
            fit_func = TF1('f_lin_neutron_ke_vs_edep', '[0] + [1]*x', fit_range_min, fit_range_max)
            fit_func.SetLineColor(ROOT.kRed + 1)
            fit_func.SetLineWidth(2)
            fit_func.FixParameter(0, 0)
            graph.Fit(fit_func, 'Q')
            intercept = fit_func.GetParameter(0)
            slope = fit_func.GetParameter(1)
            label = ROOT.TLatex(0.15, 0.85, f'y = {slope:.2f}x')
            label.SetNDC(True)
            label.SetTextSize(0.03)
            label.SetTextColor(ROOT.kRed + 1)
            label.Draw()
            print(f'neutron_E_E_plot: linear fit slope={slope:.4f} intercept={intercept:.4f}')
        else:
            print('neutron_E_E_plot: insufficient spread in KE values for linear fit.')

        output_dir = Path('1.5meter_reconstruction')
        output_dir.mkdir(exist_ok=True)
        output_path = output_dir / 'neutron_ke_vs_edep_nue_2m_4cm.root'
        root_file = ROOT.TFile(str(output_path), 'RECREATE')
        graph.Write()
        canvas.Write()
        root_file.Close()

        print(f'neutron_E_E_plot: saved scatter plot to {output_path}')


    def Draw_timing_AfterFilter(self, neutrino_energy,nth_event,axis='yz', value='time', energy='GeV',
             markerSize=3, cmap='jet', vmax=2000):
        """we want to know how exactly the neutron bounces around, so we plot the third
        plot in Draw function, then try to show the timing information of neutron energy
        deposition"""
        self.event.Jump(self.event.currentEntry)
        if self.event.SelectProtonEvent() == 0:
            print(f"Event {self.event.currentEntry}: no primary proton; skipping Draw.")
            return
        self.Collect(0.5)
        mapping = {'x':self.xx, 'y':self.yy, 'z':self.zz,'t':self.tt}
        mask = self.ii == 0
        x_coordinates = self.zz[mask]
        y_coordinates = self.yy[mask]
        times =self.tt[mask]
        times=np.insert(times, 0, 0)
        x_coordinates=np.insert(x_coordinates, 0, 0)
        y_coordinates=np.insert(y_coordinates, 0, 0)
        sort_indices = np.argsort(times)
        x_sorted = x_coordinates[sort_indices]
        y_sorted = y_coordinates[sort_indices]
        canvas = ROOT.TCanvas(f"C{self.event.currentEntry}", "Neutron Timing", 800, 600)
        x_min = x_sorted.min()
        y_min = y_sorted.min()
        x_max = x_sorted.max()
        y_max = y_sorted.max()
        x_margin = (x_max - x_min)*0.1
        y_margin = (y_max - y_min)*0.1
        color = self.tt[mask]
        hist_title = f"{self.event.currentEntry} event {neutrino_energy}GeV antinumu Neutron Timing"
        hist = ROOT.TH2F('n_scatter',hist_title, 100, x_min-x_margin, x_max+x_margin,
                         100,y_min-y_margin, y_max+y_margin)
        for x,y,c in zip(x_sorted,y_sorted,color):
            hist.Fill(x,y,c)
        ROOT.gStyle.SetPalette(ROOT.kRainBow)
        hist.GetZaxis().SetRangeUser(0, vmax)
        hist.Draw("COLZ")  # "COL" for color, "Z" to draw the color palette legend

        graph_trajectory = ROOT.TGraph(len(x_sorted),x_sorted,y_sorted)
        graph_trajectory.SetLineStyle(2)  # 2 is a dashed line in ROOT
        graph_trajectory.SetLineWidth(2)
        graph_trajectory.Draw("L SAME")

        marker_vertex = ROOT.TMarker(0, 0, 20)  # 20 is a solid circle
        marker_vertex.SetMarkerColor(ROOT.kBlue)
        marker_vertex.SetMarkerSize(1.5)
        marker_vertex.Draw("SAME")
        dir_vec = self.event.read_neutron_direction()
        if dir_vec is not None:
            comp = {'x': 0, 'y': 1, 'z': 2}
            dir_vec_2d = np.array([dir_vec[0][comp[axis[1]]], dir_vec[0][comp[axis[0]]]])

            norm = np.linalg.norm(dir_vec_2d)
            if norm > 0:
                vec = dir_vec_2d / norm
                # Scale the arrow to be 20% of the plot's smaller dimension for visibility
                arrow_scale = min(x_max - x_min, y_max - y_min) * 0.2

                arrow = ROOT.TArrow(0, 0, vec[0] * arrow_scale, vec[1] * arrow_scale, 0.02, "|>")
                arrow.SetLineWidth(2)
                arrow.SetLineColor(ROOT.kBlue)
                arrow.SetFillColor(ROOT.kBlue)
                arrow.Draw("SAME")

        # 4. --- Legend and Saving ---
        legend = ROOT.TLegend(0.1, 0.91, 0.45, 0.99)  # (x1, y1, x2, y2) in canvas coordinates
        legend.SetNColumns(2)
        legend.AddEntry(marker_vertex, "Vertex (0,0)", "p")
        legend.AddEntry(graph_trajectory, "Time-ordered path", "l")
        legend.Draw("SAME")

        # Create the output path and save the canvas
        filepath = f"/home/frunshi/frunshiworking/neutron_timing_withthreshold/{neutrino_energy}GeVantinumu{nth_event}thEvent.root"
        Path(filepath).parent.mkdir(parents=True, exist_ok=True)
        canvas.SaveAs(filepath)
        print(f"Plot saved to {filepath}")
        # print(f"x_sorted{x_sorted}, y_sorted{y_sorted},times{times}")
        # plt.plot(x_sorted, y_sorted, linestyle='--', color='gray', zorder=1)
        #
        # plt.scatter(0, 0, c="blue", s=15, label="neutron direction at the vertex")
        # plt.scatter(mapping[axis[1]][mask],mapping[axis[0]][ma
        # sk],c=self.cc[mask],s=markerSize)
        # dir_vec = self.event.read_neutron_direction()
        # if dir_vec is None:
        #     print(f"Event {self.event.currentEntry}: no neutron direction, skipping arrow.")
        # else:
        #     # print(f"--  DEBUG  -- the dir_vec is {dir_vec}")
        #     comp = {'x': 0, 'y': 1, 'z': 2}
        #     # print(f"dir_vec: {dir_vec}")
        #     dir_vec = dir_vec[0]
        #     dx = dir_vec[comp[axis[1]]]
        #     dy = dir_vec[comp[axis[0]]]
        #     vec = np.array([dx, dy])
        #     if np.linalg.norm(vec) == 0 or np.linalg.norm(vec) is None:
        #         return None
        #     vec = vec / np.linalg.norm(vec)
        #     plt.arrow(0, 0,
        #     vec[0], vec[1],
        #     color='blue',
        #     length_includes_head=True,
        #     head_width=0.1,
        #     head_length=0.1
        #     )
        # filepath = f"/home/frunshi/frunshiworking/neutron_timing/{neutrino_energy}GeVantinumu{nth_event}thEvent.pdf"
        # Path(filepath).parent.mkdir(parents=True, exist_ok=True)
        # plt.savefig(filepath)



    def Draw(self, axis='yz', value='time', energy='GeV',
             markerSize=3, cmap='jet', vmax=2000, entry=None,
             thresholds=None, output_dir=None, focus=None):
        """Render 2D views of a single event with configurable energy thresholds."""
        # focus argument retained for backward compatibility

        target_entry = self.event.currentEntry if entry is None else entry
        self.event.Jump(target_entry)
        self.Collect(0)
        mapping = {'x': self.xx, 'y': self.yy, 'z': self.zz, 'e': self.ee, 't': self.tt}

        thresholds = tuple(thresholds or (0.3, 0.4, 0.5, 0.6, 0.7))
        zoom_threshold = thresholds[0] if thresholds else 0.0
        zoom_mask = self.ee >= zoom_threshold if thresholds else np.isfinite(self.ee)

        positive_times = self.tt[self.tt > 0]
        if positive_times.size:
            vmin = positive_times.min()
            vmax = positive_times.max()
            if vmax <= vmin:
                vmax = vmin * 1.01
            time_norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)
        else:
            time_norm = None

        if self.event.evgen == 'Genie':
            plane_ylim = (-4, 4)
            plane_xlim = (-2, 8)
        elif self.event.evgen == 'Marley':
            plane_ylim = (-1, 1)
            plane_xlim = (-1, 0.5)
        else:
            plane_ylim = plane_xlim = None

        zoom_xlim = zoom_ylim = None
        zoom_x = mapping[axis[1]][zoom_mask]
        zoom_y = mapping[axis[0]][zoom_mask]
        if zoom_x.size > 0 and zoom_y.size > 0:
            finite_mask = np.isfinite(zoom_x) & np.isfinite(zoom_y)
            if np.any(finite_mask):
                zoom_x = zoom_x[finite_mask]
                zoom_y = zoom_y[finite_mask]

                def _pad(values):
                    vmin = values.min()
                    vmax = values.max()
                    span = max(vmax - vmin, 1e-3)
                    pad = max(span * 0.1, 0.02)
                    return vmin - pad, vmax + pad

                zoom_xlim = _pad(zoom_x)
                zoom_ylim = _pad(zoom_y)

        fig, axes = plt.subplots(3, 2, figsize=(10, 12), dpi=100)
        axes = axes.flatten()

        try:
            fig.suptitle(self.event.vertex.GetReaction())
            scatter_kwargs_time = dict(
                cmap=cmap,
                s=markerSize,
                edgecolors='none',
                linewidths=1,
                zorder=1,
            )
            if time_norm is not None:
                scatter_kwargs_time['norm'] = time_norm
            scatter_kwargs_color = dict(
                s=markerSize,
                edgecolors='none',
                linewidths=1,
                zorder=1,
            )

            main_axis = axes[0]
            collected_graphs = {}
            threshold_axes = axes[1:5]
            direction_axis = axes[5]
            main_mask = np.isfinite(mapping[axis[1]]) & np.isfinite(mapping[axis[0]])
            plot_time = main_axis.scatter(
                mapping[axis[1]],
                mapping[axis[0]],
                c=self.tt,
                **scatter_kwargs_time,
            )
            main_axis.set_title('All deposits (time)')
            fig.colorbar(plot_time, orientation='vertical', label='time [ns]', ax=main_axis, pad=0.02, shrink=1.0)
            if np.any(main_mask):
                collected_graphs['all_deposits'] = (
                    mapping[axis[1]][main_mask],
                    mapping[axis[0]][main_mask],
                )

            used_thresholds = list(thresholds[:len(threshold_axes)])
            for axis_obj, threshold in zip(threshold_axes, used_thresholds):
                mask = self.ee >= threshold
                x_vals = mapping[axis[1]][mask]
                y_vals = mapping[axis[0]][mask]
                if x_vals.size == 0:
                    axis_obj.set_title(f'No deposits ≥ {threshold:.1f} MeV')
                    axis_obj.text(
                        0.5, 0.5, 'No data',
                        transform=axis_obj.transAxes,
                        ha='center', va='center', fontsize=8,
                    )
                    continue

                collected_graphs[f'threshold_{threshold:.1f}'] = (x_vals, y_vals)
                colors = self.cc[mask] if self.cc.size else None
                scatter = axis_obj.scatter(
                    x_vals,
                    y_vals,
                    c=colors,
                    **scatter_kwargs_color,
                )
                axis_obj.set_title(f'threshold ≥ {threshold:.1f} MeV (ancestor)')

                pdgs = self.ancestor_pdgs[mask] if self.ancestor_pdgs.size else []
                unique_pdgs = []
                if len(pdgs):
                    seen = set()
                    for pdg in pdgs:
                        pdg_int = int(pdg)
                        if pdg_int in seen:
                            continue
                        seen.add(pdg_int)
                        unique_pdgs.append(pdg_int)
                if unique_pdgs:
                    legend_handles = []
                    for pdg in unique_pdgs:
                        legend_handles.append(
                            mpatches.Patch(
                                color=self._color_for_pdg(pdg),
                                label=self._label_for_pdg(pdg),
                            )
                        )
                    axis_obj.legend(handles=legend_handles, loc='upper right', fontsize=8, frameon=False)

            for axis_obj in threshold_axes[len(used_thresholds):]:
                axis_obj.axis('off')

            effective_xlim = zoom_xlim if zoom_xlim is not None else plane_xlim
            effective_ylim = zoom_ylim if zoom_ylim is not None else plane_ylim

            def _format_axes(target_axes, limits):
                xlim, ylim = limits
                for axis_obj in target_axes:
                    if axis_obj not in axes:
                        continue
                    if xlim is not None and ylim is not None:
                        axis_obj.set_xlim(*xlim)
                        axis_obj.set_ylim(*ylim)
                    axis_obj.tick_params(axis='both', direction='in', length=2)
                    axis_obj.set_xlabel(f'{axis[1]} [m]')

            _format_axes((main_axis,), (plane_xlim, plane_ylim))
            _format_axes(tuple(threshold_axes), (effective_xlim, effective_ylim))

            if len(threshold_axes) > 1:
                zoom_axis = threshold_axes[1]
                zoom_extent = 0.1
                zoom_axis.set_xlim(-zoom_extent, zoom_extent)
                zoom_axis.set_ylim(-zoom_extent, zoom_extent)
                zoom_axis.set_title(f"{zoom_axis.get_title()} (±0.1 m)")

            _format_axes((direction_axis,), (plane_xlim, plane_ylim))
            direction_axis.set_title('Initial neutron/proton directions')
            direction_axis.axhline(0, color='lightgray', linewidth=0.8)
            direction_axis.axvline(0, color='lightgray', linewidth=0.8)
            comp = {'x': 0, 'y': 1, 'z': 2}
            arrow_scale = None
            if effective_xlim and effective_ylim:
                arrow_scale = min(effective_xlim[1] - effective_xlim[0], effective_ylim[1] - effective_ylim[0]) * 0.25
            else:
                arrow_scale = 0.5
            legend_handles_dir = []
            legend_seen = set()
            for particle in getattr(self.event.vertex, 'Particles', []):
                pdg = particle.GetPDGCode()
                if pdg not in (2112, 2212):
                    continue
                momentum = particle.GetMomentum()
                if momentum is None:
                    continue
                vec = np.array([momentum.Px(), momentum.Py(), momentum.Pz()])
                if not np.any(vec):
                    continue
                vec = vec / np.linalg.norm(vec)
                plane_vec = np.array([
                    vec[comp[axis[1]]],
                    vec[comp[axis[0]]],
                ])
                if not np.any(plane_vec):
                    continue
                norm = np.linalg.norm(plane_vec)
                if norm == 0:
                    continue
                plane_vec = plane_vec / norm * arrow_scale
                color = self._color_for_pdg(pdg)
                direction_axis.arrow(
                    0,
                    0,
                    plane_vec[0],
                    plane_vec[1],
                    color=color,
                    head_width=arrow_scale * 0.15,
                    head_length=arrow_scale * 0.15,
                    length_includes_head=True,
                )
                if pdg not in legend_seen:
                    legend_seen.add(pdg)
                    legend_handles_dir.append(mpatches.Patch(color=color, label=self._label_for_pdg(pdg)))
            if legend_handles_dir:
                direction_axis.legend(handles=legend_handles_dir, loc='upper right', fontsize=8, frameon=False)

            for idx, axis_obj in enumerate(axes):
                if idx % 2 == 0:
                    axis_obj.set_ylabel(f'{axis[0]} [m]')

            label_axis = threshold_axes[1] if len(threshold_axes) > 1 else (threshold_axes[0] if threshold_axes else main_axis)
            entries = []
            for particle in getattr(self.event.vertex, 'Particles', []):
                trkId = particle.GetTrackId()
                if trkId < 0:
                    continue
                momentum = particle.GetMomentum()
                if momentum is None:
                    continue
                kinetic_energy = momentum.E() - momentum.M()
                text = f'{particle.GetName()}: {kinetic_energy:.1f} MeV'
                entries.append((text, self._color_for_pdg(particle.GetPDGCode())))

            for idx, (text, color) in enumerate(entries):
                label_axis.text(
                    0.02,
                    0.98 - idx * 0.05,
                    text,
                    transform=label_axis.transAxes,
                    va='top',
                    color=color,
                    fontsize=8,
                )

            save_dir = Path(output_dir) if output_dir else Path(self.event.plotpath) / 'individual_display'
            save_dir.mkdir(parents=True, exist_ok=True)
            filename_stem = f'particle_{value}_{axis}_{energy}_evt_{self.event.currentEntry:05d}'
            output_path = save_dir / f'{filename_stem}.png'
            fig.savefig(output_path, dpi=200)

            root_dir = save_dir / 'root'
            root_dir.mkdir(exist_ok=True)
            root_path = root_dir / f'{filename_stem}.root'
            root_file = ROOT.TFile(str(root_path), 'RECREATE')
            try:
                for graph_name, (x_vals, y_vals) in collected_graphs.items():
                    if x_vals.size == 0:
                        continue
                    graph = ROOT.TGraph(
                        len(x_vals),
                        array.array('d', x_vals.astype(float)),
                        array.array('d', y_vals.astype(float)),
                    )
                    graph.SetName(f'{graph_name}_{axis}')
                    graph.SetTitle(f'{graph_name} ({axis})')
                    graph.Write()
            finally:
                root_file.Close()

            print(
                f'Draw: saved figure for entry {self.event.currentEntry} to {output_path} '
                f'and ROOT graphs to {root_path}'
            )
        finally:
            plt.close(fig)




    def Draw_pion(self, axis='yz', value='time', energy='GeV',
                  markerSize=3, cmap='jet', vmax=2000):
        """Wrapper that draws views focused on pion energy depositions."""
        return self.Draw(axis=axis, value=value, energy=energy,
                         markerSize=markerSize, cmap=cmap, vmax=vmax, focus='pion')




    def electron_energy_spectrum(self):
            n_events = 9000
            electron_energy_array = []
            for i in range(n_events):
                self.event.Jump(i)
                if i % 256 == 0:
                    print(f'event: {i}')
                electron_energy = self.event.read_electron_energy()
                electron_energy_array.append(electron_energy)
            canvas = ROOT.TCanvas("Anti-electron neutrino", "Anti-electron neutrino", 800, 600)
            h1 = ROOT.TH1F("electron_energy_spectrum_antinue", "electron_energy_spectrum", 100, 0, 1000)
            for value in electron_energy_array:
                h1.Fill(value)
            h1.SetTitle("electron energy spectrum")
            h1.Set_xlabel("electron kinetic energy [MeV]")
            h1.Set_ylabel("Counts")
            h1.Draw("HIST")
            save_path = "/home/frunshi/frunshiworking/neutrino_reconstruction/muon_energy_spectrum_antinue.root"
            canvas.SaveAs(save_path)







    def Draw_N_M(self, axis='yz', value='time', energy='GeV',
             markerSize=3, cmap='jet', vmax=2000):
        # load event & collect deposit info
        self.event.Jump(self.event.currentEntry,0)
        self.Collect(0)
        mapping = {'x': self.xx, 'y': self.yy, 'z': self.zz, 'e': self.ee, 't': self.tt}
        mask = self.ee >= 0.5
        mask_ME = (self.ee >= 0.5) & (self.ME ==1)

        fig, axes = plt.subplots(3, 2, figsize=(10, 12), dpi=100)
        axes = axes.flatten()

        ax1, ax2, ax3, ax4, ax5, ax6 = axes[:6]

        fig.suptitle(self.event.vertex.GetReaction())

        # — Panel 1: time scatter —
        plot_time = ax1.scatter(
            mapping[axis[1]], mapping[axis[0]],
            c=self.tt,
            cmap='viridis',
            norm=matplotlib.colors.LogNorm(),
            s=15,
            edgecolors='none',
            linewidths=1,
            zorder=1
        )

        # — Panel 2: energy scatter —
        plot_energy = ax2.scatter(
            mapping[axis[1]], mapping[axis[0]],
            c=self.ee * 2,
            cmap='viridis',
            vmax=2.5,
            s=15,
            edgecolors='none',
            linewidths=1,
            zorder=1
        )
        x_start, x_end = ax2.get_xlim()
        y_start, y_end = ax2.get_ylim()

        ax2.set_xlabel('Energy [GeV]')
        ax2.set_title('Energy [MeV/cm]')
        fig.colorbar(plot_energy, orientation='vertical', label="MeV/cm", ax=ax2, pad=0.02, shrink=1.0)
        fig.colorbar(plot_time, orientation ='vertical', label = "time[ns]",ax=ax1, pad=0.02, shrink=1.0)

        # — Panel 3: neutron direction + filtered scatter —
        dir_vec = self.event.read_neutron_direction()
        midpoint, electron_direction, counter, id, midpoint_end = self.event.find_Michel_Electron()
        y_michel = midpoint[1] / 100
        z_michel = midpoint[2] / 100
        y_michel_end = midpoint_end[1] / 100
        z_michel_end = midpoint_end[2] / 100
        if dir_vec is None:
            print(f"Event {self.event.currentEntry}: no neutron direction, skipping arrow.")
        else:
            # print(f"--  DEBUG  -- the dir_vec is {dir_vec}")
            comp = {'x': 0, 'y': 1, 'z': 2}
            # print(f"--- DEBUG --- At Event {self.event.currentEntry}, the value of 'axis' is: '{axis}'")
            dx = dir_vec[comp[axis[1]]]
            dy = dir_vec[comp[axis[0]]]
            vec = np.array([dx, dy])
            vec = vec / np.linalg.norm(vec)

            ax3.arrow(
                0, 0,
                vec[0], vec[1],
                color='blue',
                length_includes_head=True,
                head_width=0.1,
                head_length=0.1,
            )
            ax3.scatter(z_michel, y_michel, c="green", label="michel electron created")
            ax3.scatter(z_michel_end, y_michel_end, c="red", label="michel electron stopped")
            ax3.scatter(0, 0, c="blue", s=15, label="neutron direction")
            ax3.set_xlabel("z direction [m]")
            ax3.set_ylabel("y direction [m]")
            ax3.legend()
            ax3.scatter(
                mapping[axis[1]][mask],
                mapping[axis[0]][mask],
                c=self.cc[mask],
                s=markerSize
            )
            ax3.set_xlabel('z-y plane,color-coated with different ancestor, [m])')

        # ax4.scatter(
        #     mapping[axis[1]][mask],
        #     mapping[axis[0]][mask],
        #     c="black",
        #     s=markerSize
        # )
        #
        x_ax3_start, x_ax3_end = ax3.get_xlim()
        y_ax3_start, y_ax3_end = ax3.get_ylim()
        #
        # ax4.set_xlabel('Ancestor ID')
        # ax4.set_xlim(x_ax3_start, x_ax3_end)
        # ax4.set_ylim(y_ax3_start, y_ax3_end)
        # ax4.set_xlabel("zoom out version of selected edep(z-y plane, [m])")

        ax4.scatter(mapping[axis[1]][mask], mapping[axis[0]][mask],c=self.tt[mask], s=markerSize)
        ax4.scatter(z_michel, y_michel, c="green", label="michel electron created")
        ax4.scatter(z_michel_end, y_michel_end, c="red", label="michel electron stopped")
        all_z = [z_michel, z_michel_end]
        all_y = [y_michel, y_michel_end]
        min_y = min(all_y)
        max_y = max(all_y)
        min_z = min(all_z)
        max_z = max(all_z)
        padding =0.5
        ax4.set_xlim(min_z - abs(min_z-max_z)*padding, max_z + abs(max_z-min_z)*padding)
        ax4.set_ylim(min_y - abs(min_y-max_y)*padding, max_y + abs(max_y-min_y)*padding)
        ax4.set_xlabel('z-y plane,zoom into Michel Electron [m])')
        ax4.set_ylabel('y direction [m]')
        ax4.legend()



        data = np.asarray(self.event.neutron_edep_information())
        if data is None or len(data) == 0:
            return
        arr = np.asarray(data)
        if any(len(row) < 3 for row in data):
            return

        all_length = data[:, 0]
        all_energy = data[:, 1]
        all_angle = data[:, 2]


        # plot 5 section


        # plot_energy_angle = ax5.scatter(
        #     all_energy, all_angle,
        #     c=all_length,
        #     cmap='viridis',
        #     s=15,
        #     edgecolors='none',
        #     linewidths=1,
        #     zorder=1
        # )
        # ax5.set_xlim(0, 3)
        # ax5.set_ylim(0, 3.1416)
        # fig.colorbar(plot_energy_angle, label="length[cm]", ax=ax5, pad=0.02, shrink=1.0)
        # ax5.set_xlabel("energy[MeV]")
        # ax5.set_ylabel("angle[rad]")

        ax5.scatter(z_michel, y_michel, c="green", label="michel electron created")
        ax5.scatter(z_michel_end, y_michel_end, c="red", label="michel electron stopped")
        ax5.scatter(0, 0, c="blue", s=15, label="neutron direction")
        ax5.set_xlabel("z direction [m]")
        ax5.set_ylabel("y direction [m]")
        ax5.legend()

        ax5.scatter(
            mapping[axis[1]][mask],
            mapping[axis[0]][mask],
            c=self.tt[mask],
            s=markerSize
        )

        ax5.set_xlim(x_ax3_start, x_ax3_end)
        ax5.set_ylim(y_ax3_start, y_ax3_end)
        ax5.set_xlabel('z-y plane,timing of edep above 0.5MeV threshold)')
        ax5.set_ylabel("y direction [m]")
        ax5.legend()


        # plot 6 section

        # plot_length_angle = ax6.scatter(
        #     all_length, all_angle,
        #     c=all_energy,
        #     cmap='viridis',
        #     s=15,
        #     edgecolors='none',
        #     linewidths=1,
        #     zorder=1
        # )
        # ax6.set_xlim(0, 400)
        # ax6.set_ylim(0, 3.1416)
        # fig.colorbar(plot_length_angle, label="energy[MeV]", ax=ax6, pad=0.02, shrink=1.0)
        # ax6.set_xlabel("Length[cm]")
        # ax6.set_ylabel("angle[rad]")

        plot6 = ax6.scatter(mapping[axis[1]][mask_ME], mapping[axis[0]][mask_ME], c=self.tt[mask_ME], s=markerSize)
        ax6.set_xlim(x_ax3_start, x_ax3_end)
        ax6.set_ylim(y_ax3_start, y_ax3_end)
        ax6.set_xlabel("z direction [m]")
        ax6.set_ylabel("y direction [m]")
        fig.colorbar(plot6, orientation='vertical', label="time[ns]", ax=ax6, pad=0.02, shrink=1.0)
        ax6.legend()

        ax6.set_xlabel('z-y plane,timing of Michel Electron Children edep')
        ax6.legend()

        for ax in (ax1, ax2):
            if self.event.evgen == 'Genie':
                ax.set_ylim(-4, 4)
                ax.set_xlim(-2, 8)
                xpos, ypos, interval = -1.8, 3.6, 0.35
            elif self.event.evgen == 'Marley':
                ax.set_ylim(-1, 1)
                ax.set_xlim(-1, 0.5)
                xpos, ypos, interval = -0.95, 0.9, 0.1
            else:
                print("Unknown event generator!")
                sys.exit()
            ax.tick_params(axis='both', direction='in', length=2)
            ax.set_xlabel(f'{axis[1]} [m]')
            if ax is ax1:
                ax.set_ylabel(f'{axis[0]} [m]')

        entries = []
        for particle in self.event.vertex.Particles:
            trkId = particle.GetTrackId()
            if trkId < 0:
                continue
            mom = particle.GetMomentum()
            KE = mom.E() - mom.M()
            name = f'{particle.GetName()}: {KE :.1f} MeV'
            color = self._color_for_pdg(particle.GetPDGCode())
            entries.append((name, color))

        for idx, (text, color) in enumerate(entries):
                ax3.text(
                    0.02,
                    0.98 - idx * 0.05,
                    text,
                    transform=ax3.transAxes,
                    va='top',
                    color=color,
                    fontsize=8
                )

        filepath_1= f"{self.event.plotpath}/developing/BeforeAfterThreshold/particle_{value}_{axis}_{energy}_evt_{self.event.currentEntry}.pdf"
        Path(filepath_1).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(filepath_1)
        plt.clf()
        plt.close()


    # For GeV events
    def hist_dEdx(self):
        plt.hist(self.ee * 2, range=(0, 6), bins=100)
        plt.xlabel('dE/dx [MeV/cm]')
        plt.draw()
        plt.savefig(self.event.plotpath + '/dE_dx_evt_%d.pdf' % self.event.currentEntry)
        plt.clf()  # important to clear figure
        plt.close()
        #plt.show()

    # For low E MeV events, edep length is smaller than 0.5cm (5mm step limit in gdml)
    def hist_LowE_dEdx(self):
        plt.hist(np.divide(self.ee, self.ll), range=(0, 30), bins=500)
        plt.xlabel('dE/dx [MeV/cm]')
        plt.draw()
        plt.savefig(self.event.plotpath + '/LowE_dE_dx_evt_%d.pdf' % self.event.currentEntry)
        plt.clf()  # important to clear figure
        plt.close()

    def print_muon_direction(self):
        """a function inputs the order of the event and return the direction of muon energy deposition relative to
         muon direction"""
        n_events = 1000
        muon_angle = []
        energy_low = 250
        energy_high = 350
        energy_mean = (energy_low + energy_high) / 2000.0

        for idx in range(n_events):
            self.event.Jump(idx)
            self.Collect(0.5)
            CurrentEnergy = self.muon_energy[0]
            print(f"the muon energy is {CurrentEnergy} for the {idx}th event")
            if energy_low < CurrentEnergy < energy_high:
                raw_angle =self.aa
                angle =np.array(raw_angle)
                print(f"in this event, there is in total {len(raw_angle)}energy depositions")
                trig =0
                for j in range(len(angle)):
                    if np.array(self.ii[j] == 1):
                        #this judges if it is a muon energy deposition
                        muon_angle.append(angle[j])
        muon_angle = np.array(muon_angle)
        print(f"the muon angle in {n_events} events in the range{energy_low}~{energy_high} is {muon_angle}")
        file_name = f"{energy_mean}GeV_muon_angle_distribution.txt"
        np.savetxt(file_name, muon_angle)
    # def print_muon_per_event(self):
    #     """this time we want to have a right definition of purity and efficiency, to do that instead fo printing out all muon
    #     energy depositions we decide to print the muon energy deposition with the largest angle to the muon direction per event
    #     in this case we select one muon energy deposition and one neutron energy deposition per event"""
    #     n_events = 1000
    #     muon_angle = []
    #     energy_low = 250
    #     energy_high = 350
    #     energy_mean = (energy_low + energy_high) / 2000
    #     for idx in range(n_events):
    #         self.event.Jump(idx)
    #         self.Collect(0.5)
    #         currentenergy = self.muon_energy[0]
    #         muon_direction= np.array(self.event.read_muon_direction())
    #         angle_max = 0
    #         if energy_low < currentenergy < energy_high:
    #             for





    def print_nearest_muon(self,energy_for_muon):
        """a function inputs the order of the event and return the direction of muon energy deposition relative to the muon track,
        check event.py for detail function on how to find the first discontinuous m-edep and the angle in a certain function"""
        count = 0
        n_events = 9000
        energy_low = energy_for_muon*1000-25
        energy_high = energy_for_muon*1000+25
        energy_mean = round((energy_low + energy_high) / 2000,4)
        muon_total_position = []
        muon_total_direction = []
        relative_angle = []
        for idx in range(n_events):
            self.event.Jump(idx)
            self.Collect(0.5)
            trig = self.event.check_N_M_only_event()
            currentenergy = self.muon_energy[0]
            angle_Medep = 0
            if energy_low < currentenergy < energy_high and trig == 1:
                print("this is a legal event for muon edep")
                count+=1
                result = self.event.find_nearest_discontinuous_medep()
                if result is None:
                    # This block runs if the function returned None.
                    # You can't find the medep, so set default values.
                    print("Warning: find_nearest_discontinuous_medep() returned None. Using default values.")
                    angle_Medep = 0  # A safe default
                    position = np.array([100, 10000, 10000])  # A safe default
                else:
                    # 3. If the result is NOT None, it's safe to unpack it
                    angle_Medep, position = result
                relative_angle.append(angle_Medep)
                muon_total_position.append(position)
        filename = f"/home/frunshi/frunshiworking/P_E_plot/{energy_mean}GeV_medep_angle_distribution.txt"
        np.savetxt(filename, relative_angle)
        print(f'{count}event recorded on Medep')
        return count


    def print_nearest_neutron(self,energy_for_muon):
        """a function inputs the order of the event and return the direction of nearest neutron energy deposition relative to the muon track"""
        n_events = 9000
        nearest_neutron = []
        mm2m =0.001
        ##-----------------------------------------------------------------##
        # notice that the muon energy is usually 200MeV less than the energy of primary neutrinos
        energy_low = energy_for_muon*1000-25
        energy_high = energy_for_muon*1000+25
        energy_mean = round((energy_low + energy_high) / 2000,5)
        neutron_total_position = []
        neutron_total_direction = []
        relative_angle = []
        count = 0
        for idx in range(n_events):
            self.event.Jump(idx)
            self.Collect(0.5)
            currentenergy = self.muon_energy[0]
            muon_direction,mom = self.event.read_muon_direction_and_amplitude()
            trig = self.event.check_N_M_only_event()
            # print("the muon energy is :", currentene
            # rgy)
            # print("the muon direction is :", muon_direction)
            nearest_neutron_distance = 10000
            nearest_x,nearest_y,nearest_z =0,0,0
            if energy_low < currentenergy < energy_high and trig == 1:
                # print("the energy qualified")
                count+=1
                for i, track in enumerate(self.event.tracks):
                    depoList = track.association['depoList']
                    ancestor = track.association['ancestor']
                    pdg_check =self.event.tracks[ancestor].GetPDGCode()
                    if pdg_check == 2112 :
                        # print("exist neutron track")
                        for di in depoList:
                            depo = self.event.depos[di]
                            x = (depo.GetStart().X() + depo.GetStop().X()) / 2 * mm2m
                            y = (depo.GetStart().Y() + depo.GetStop().Y()) / 2 * mm2m
                            z = (depo.GetStart().Z() + depo.GetStop().Z()) / 2 * mm2m
                            if x ** 2 + y ** 2 + z ** 2 <= nearest_neutron_distance ** 2:
                                nearest_neutron_distance = (x ** 2 + y ** 2 + z ** 2) ** 0.5
                                nearest_x, nearest_y, nearest_z = x, y, z
            if nearest_x == 0 and nearest_y == 0 and nearest_z == 0:
                # print("there is a event where we aren't able to find the nearest neutron position")
                continue
            neutron_position = np.array([nearest_x, nearest_y, nearest_z])
            neutron_total_position.append(neutron_position)

            neutron_direction = neutron_position / np.linalg.norm(neutron_position)
            neutron_total_direction.append(neutron_direction)
            # print("the neutron array is ",neutron_direction)
            cos = neutron_direction[0]*muon_direction[0]+neutron_direction[1]*muon_direction[1]+neutron_direction[2]*muon_direction[2]
            # print ("cosine is :", cos)
            relative_angle =abs(np.degrees(np.arccos(cos)))
            # relative_angle.append(abs(np.degrees(np.arccos(cos))))
            # print("the neutron angle in this event is :",relative_angle)
            nearest_neutron.append(relative_angle)

        file_name = f"/home/frunshi/frunshiworking/P_E_plot/{energy_mean}GeV_nedep_distribution.txt"
        np.savetxt(file_name, nearest_neutron)
        print(f'{count}event recorded on Nedep')
        return neutron_total_position, neutron_total_direction,relative_angle

    def percentage_of_NEAREARLY(self, n_events: int = 9000, energy_step: int = 100,
                                 max_energy: int = 900, time_tolerance_ns: float = 0.8):
        """Loop over events and report how often the nearest neutron deposition is also the earliest.

        Args:
            n_events: Maximum number of events to scan.
            energy_step: Width of each kinetic-energy bin in MeV.
            max_energy: Upper edge of the final bin (MeV).
            time_tolerance_ns: Absolute tolerance (ns) used when comparing the
                nearest and earliest interaction times.

        Returns:
            List of dictionaries, each describing a single energy bin with the
            total number of valid events, how many had matching nearest/earliest
            times, and the resulting percentage.
        """

        if energy_step <= 0:
            raise ValueError("energy_step must be positive")

        n_entries = getattr(self.event, 'nEntry', n_events)
        limit = min(n_events, n_entries)

        # Build bins such as (0, 100), (100, 200), ..., (800, 900)
        energy_bins: List[Tuple[int, int]] = []
        start = 0
        while start < max_energy:
            energy_bins.append((start, min(start + energy_step, max_energy)))
            start += energy_step

        stats = {energy_range: {"total": 0, "matches": 0} for energy_range in energy_bins}

        for idx in range(limit):
            self.event.Jump(idx)
            if idx%100==0:
                print(f'{idx}events scanned')

            if self.event.SelectNeutronEvent() != 1:
                continue  # Restrict to single-neutron events for a well-defined earliest hit

            neutron_energy = self.event.read_neutron_energy()
            if neutron_energy is None:
                continue

            bin_key = None
            for low, high in energy_bins:
                # Treat the upper edge of the final bin as inclusive to catch 900 MeV entries
                if low <= neutron_energy < high or (high == max_energy and math.isclose(neutron_energy, high, abs_tol=1e-9)):
                    bin_key = (low, high)
                    break

            if bin_key is None:
                continue

            earliest_info = self.event.earliest_interaction_time()
            if earliest_info is None:
                continue
            earliest_time, _ = earliest_info
            if earliest_time is None or earliest_time >= 1.0e11:
                # Skip events where no neutron energy deposition was found
                continue

            muon_energy = self.event.read_muon_energy()
            # Ensure the downstream helper sees the muon energy we just read
            self.muon_energy = np.array([muon_energy])
            nearest_info = self.event.find_nearest_discontinuous_neutron(self)
            # print("nearest_info", nearest_info)
            if nearest_info is None:
                continue
            nearest_position, _, nearest_time = nearest_info
            if nearest_position is None or nearest_time is None:
                continue
            if np.linalg.norm(np.asarray(nearest_position)) < 1e-9:
                continue

            stats[bin_key]["total"] += 1
            if math.isclose(earliest_time, nearest_time, abs_tol=time_tolerance_ns):
                stats[bin_key]["matches"] += 1

        summary = []
        total_number_of_events = 0
        print("percentage of Nearest N-Edep(Neutron Energy deposition) is actually the earliest")
        print("Only considering one Neutron (from Neutrino vertex) Events")
        for low,high in energy_bins:
            counts = stats[(low, high)]
            total_number_of_events += counts["total"]
        print(f"there are in total {total_number_of_events} Antinumu events in consideration")
        for low, high in energy_bins:
            counts = stats[(low, high)]
            total = counts["total"]
            matches = counts["matches"]
            percentage = (matches / total * 100.0) if total else 0.0
            summary.append({
                "range": (low, high),
                "total_events": total,
                "matching_events": matches,
                "percentage": percentage,
            })
            print(f"Neutron Kinetic Energy {low:3d}-{high:3d} MeV: {percentage:6.2f}% ({matches}/{total})")

        return summary
    def find_nearest_particles(self, muon_angle_cut):
        """
        Efficiently finds the nearest muon- and neutron-induced depositions in a single pass.

        For each event that passes the muon energy filter, this function identifies the
        nearest valid deposition from a muon ancestor and a neutron ancestor. A
        deposition is considered "valid" if it lies outside the specified angle cut.

        It calculates:
        - Purity: Fraction of events where the nearest valid deposition is a neutron.
        - Efficiency: Fraction of events where at least one valid neutron deposition is found.

        Args:
            muon_angle_cut (float): The angle in degrees. Depositions inside
                                    this cone are ignored for all particles.

        Returns:
            tuple: A tuple containing the purity and efficiency values.
        """
        n_events = 200
        energy_low = 250
        energy_high = 350
        mm2m = 0.001

        # --- Lists to store final results (optional, for debugging) ---
        muon_total_position, muon_total_direction, muon_total_angle = [], [], []
        neutron_total_position, neutron_total_direction, neutron_total_angle = [], [], []

        # --- Counters for metrics ---
        event_counter = 0  # Total events passing the energy filter
        effective_events = 0  # Events with at least one valid neutron
        pure_events = 0  # Events where the nearest valid deposition is a neutron

        # --- 1. SINGLE EVENT LOOP ---
        for idx in range(n_events):
            self.event.Jump(idx)
            self.Collect(0.5)

            if not self.muon_energy:
                continue

            current_energy = self.muon_energy[0]

            # Filter for events with primary muons in the desired energy range
            if energy_low < current_energy < energy_high:
                event_counter += 1
                muon_direction ,mom= np.array(self.event.read_muon_direction_and_amplitude())
                muon_direction /= np.linalg.norm(muon_direction)

                # --- Initialize "nearest" variables for this event ---
                nearest_muon_dist_sq = float('inf')
                nearest_muon_pos = None
                nearest_neutron_dist_sq = float('inf')
                nearest_neutron_pos = None

                # --- 2. UNIFIED TRACK PROCESSING ---
                nearest_neutron_angle =0
                for track in self.event.tracks:
                    if not track.association['depoList']:
                        continue

                    ancestor = track.association['ancestor']
                    pdg_check = self.event.tracks[ancestor].GetPDGCode()
                    depo = self.event.depos[track.association['depoList'][0]]

                    x = (depo.GetStart().X() + depo.GetStop().X()) / 2 * mm2m
                    y = (depo.GetStart().Y() + depo.GetStop().Y()) / 2 * mm2m
                    z = (depo.GetStart().Z() + depo.GetStop().Z()) / 2 * mm2m
                    dist_sq = x ** 2 + y ** 2 + z ** 2

                    # Calculate angle relative to the primary muon direction
                    direction_of_edep = np.array([x, y, z]) / np.sqrt(dist_sq)
                    dot_prod = np.clip(np.dot(direction_of_edep, muon_direction), -1.0, 1.0)
                    angle = np.degrees(np.arccos(dot_prod))

                    # To be considered, any deposition MUST be outside the cone
                    if abs(pdg_check) == 2112:
                        if dist_sq < nearest_neutron_dist_sq:
                            nearest_neutron_angle = angle

                    if angle >= muon_angle_cut:
                        # Check for MUON ancestor
                        if abs(pdg_check) == 13:
                            if dist_sq < nearest_muon_dist_sq:
                                nearest_muon_dist_sq = dist_sq
                                nearest_muon_pos = [x, y, z]

                        # Check for NEUTRON ancestor
                        elif pdg_check == 2112:
                            # print("there is a neutron detected")
                            if dist_sq < nearest_neutron_dist_sq:
                                nearest_neutron_dist_sq = dist_sq
                                nearest_neutron_pos = [x, y, z]

                # --- 4. METRIC CALCULATION FOR THIS EVENT ---

                # -- Purity: Is the nearest valid deposition a neutron? --
                if nearest_neutron_dist_sq < nearest_muon_dist_sq:
                    pure_events += 1
                # print(f"nearest_neutron_dist_sq: {nearest_neutron_dist_sq}")
                # -- Efficiency: Was at least one valid neutron found? --
                if nearest_neutron_angle > muon_angle_cut:
                    effective_events += 1

        # --- 5. CALCULATE FINAL METRICS ---
        if event_counter == 0:
            return 0.0, 0.0

        # Divide by the total number of events that passed the energy filter
        purity = pure_events / event_counter
        efficiency = effective_events / event_counter

        print(f"Total events analyzed: {event_counter}")
        print(f"Pure events (nearest valid depo is neutron): {pure_events}")
        print(f"Effective events (at least one valid neutron found): {effective_events}")
        print(f"Purity: {purity:.4f}")
        print(f"Efficiency: {efficiency:.4f}")

        return purity, efficiency


    # def print_nearest_muon(self,angle_cut):
    #
    #
    #     """this function inputs the angle cut, and under the cut first exclude all the energy depositions
    #     within the cut and then figure out which energy deposition is the nearest, like the function print_nearest_neutron,
    #     eventually we want it to print out three lists"""
    #     n_events = 1000
    #     energy_low = 250
    #     energy_high = 350
    #     mm2m =0.001
    #     muon_total_position = []
    #     muon_total_distance = []
    #     for idx in range(n_events):
    #         self.event.Jump(idx)
    #         self.Collect(0.5)
    #         currentenergy = self.muon_energy[0]
    #         muon_direction= np.array(self.event.read_muon_direction())
    #         if energy_low < currentenergy < energy_high:
    #             print("the energy qualified")
    #             nearest_muon_distance = 10000
    #             nearest_x, nearest_y, nearest_z = 0, 0, 0
    #             for i, track in enumerate(self.event.tracks):
    #                 depoList = track.association['depoList']
    #                 ancestor = track.association['ancestor']
    #                 pdg_check =self.event.tracks[ancestor].GetPDGCode()
    #                 if abs(pdg_check) == 13 :
    #                         depo = self.event.depos[depoList[0]]
    #                         x = (depo.GetStart().X() + depo.GetStop().X()) / 2 * mm2m
    #                         y = (depo.GetStart().Y() + depo.GetStop().Y()) / 2 * mm2m
    #                         z = (depo.GetStart().Z() + depo.GetStop().Z()) / 2 * mm2m
    #                         direction_of_edep = np.array([x,y,z])/np.linalg.norm([x,y,z])
    #                         angle_difference = abs(np.degrees(np.arccos(np.dot(direction_of_edep,muon_direction))))
    #                         if angle_difference > angle_cut:
    #                             if x ** 2 + y ** 2 + z ** 2 <= nearest_muon_distance ** 2:
    #                                 nearest_muon_distance = (x ** 2 + y ** 2 + z ** 2) ** 0.5
    #                                 nearest_x, nearest_y, nearest_z = x, y, z
    #             if nearest_x == 0 and nearest_y == 0 and nearest_z == 0:
    #                 print("there is a event where we arent able to find the nearest neutron position")
    #                 continue
    #             else:
    #                 muon_total_distance.append(nearest_muon_distance)
    #                 muon_total_position.append([nearest_x, nearest_y, nearest_z])
    #
    #     return muon_total_position, muon_total_distance
    def Draw_angle(self, axis='yz', value='charge', energy='GeV', markerSize=3, cmap='jet', vmax=2000):
        # particle, timing, dE/dx
        mapping = {'x': self.xx, 'y': self.yy, 'z': self.zz}

        fig, axes = plt.subplots(2, 3, figsize=(15, 8), dpi=600)
        ax1, ax2, ax3, ax4, ax5, ax6 = axes.flatten()
        fig.suptitle(self.event.vertex.GetReaction())
        cb_ax = fig.add_axes([.94, .124, .02, .754])
        midpoint, electron_direction, counter, id,midpoint_end = event.find_Michel_Electron()
        y_michel = midpoint[1]
        z_michel = midpoint[2]
        y_michel_end = midpoint_end[1]
        z_michel_end = midpoint_end[2]
        ax1.scatter(z_michel,y_michel,c= "green",label="michel electron")
        ax1.scatter(z_michel_end,y_michel_end,c= "red",label="michel electron(at the end )")
        # particle plot
        ax1.scatter(mapping[axis[1]], mapping[axis[0]], c=self.cc, s=markerSize)

        if value == 'time':
            # timing plot
            norm = matplotlib.colors.LogNorm(vmin=1, vmax=vmax)
            plot_12 = ax2.scatter(mapping[axis[1]], mapping[axis[0]], c=self.tt, cmap=cmap, vmin=1, vmax=vmax,
                                  norm=matplotlib.colors.LogNorm(), s=markerSize)
            cb_ax.set_xlabel('ns')

        elif value == 'charge':
            # charge plot
            norm = matplotlib.colors.LogNorm(vmin=1, vmax=vmax)
            if energy == 'GeV':
                plot_12 = ax2.scatter(mapping[axis[1]], mapping[axis[0]], c=self.ee * 2, cmap=cmap, vmax=4,
                                      s=markerSize)
            elif energy == 'MeV':
                # plot_12 = ax2.scatter(mapping[axis[1]], mapping[axis[0]], c=np.divide(self.ee, self.ll), cmap=cmap, vmax=4, s=markerSize)
                plot_12 = ax2.scatter(mapping[axis[1]], mapping[axis[0]], c=self.ee, cmap=cmap, vmax=2.5, norm=norm,
                                      s=markerSize)
            cb_ax.set_xlabel('MeV')

        fig.colorbar(plot_12, orientation='vertical', cax=cb_ax)
        # fig.tight_layout()
        for ax in (ax1, ax2, ax3, ax4, ax5, ax6):
            if self.event.evgen == 'Genie':
                ax.set_ylim(-0.5, 0.5)
                ax.set_xlim(-0.3, 0.5)
                xpos = -1.8
                ypos = 3.6
                interval = 0.35
            else:
                print("Unknown event generator!")
                sys.exit()
            ax.tick_params(axis='y', direction='in', length=2)
            ax.tick_params(axis='x', direction='in', length=2)
            ax.set_xlabel(f'{axis[1]} [m]')
            if ax == ax1:
                ax.set_ylabel(f'{axis[0]} [m]')
        for particle in self.event.vertex.Particles:
            trkId = particle.GetTrackId()
            if trkId < 0:
                continue
            name = particle.GetName()
            color = self._color_for_pdg(particle.GetPDGCode())
            mom = particle.GetMomentum()
            KE = mom.E() - mom.M()
            name = '%s: %.1f MeV' % (name, KE)

            ax1.text(xpos, ypos, name, color=color)
            ypos -= interval
        #####
        #draw the arrow
        comp = {'x': 0, 'y': 1, 'z': 2}

        dn = self.event.read_neutron_direction()
        dm ,mom= self.event.read_muon_direction_and_amplitude()
        if dn is None or dm is None:
            return
        for vec, color, label in ((dn, 'blue', 'neutron'), (dm, 'red', 'muon')):
            v = np.array([vec[comp['z']], vec[comp['y']]])
            v = v / np.linalg.norm(v)
            ax2.arrow(0, 0, v[0], v[1],
                      color=color, length_includes_head=True,
                      head_width=0.1, head_length=0.1)
            # legend via patches
        handles = [
            mpatches.Patch(color='blue', label='neutron(arrow)'),
            mpatches.Patch(color='red', label='muon(arrow)'),
        ]
        ax2.legend(handles=handles, loc='upper right', fontsize=8)
        # ax2.plot()

        dathreshold = 10
        mask = self.aa > dathreshold
        all_colors = np.where(self.ii == 0, 'red', np.where(self.ii == 1, 'black', 'grey'))
        print("self.ii is ", self.ii)
        ax3.scatter(
            mapping[axis[1]][mask],
            mapping[axis[0]][mask],
            s=markerSize,
            c=all_colors[mask],  # pick only the masked ones
        )
        handles_2 = [
            mpatches.Patch(color='red', label='neutron'),
            mpatches.Patch(color='black', label='misidentify muon'),
            mpatches.Patch(color='grey', label='other primary particles'),
        ]
        ax3.legend(handles=handles_2, loc='upper right', fontsize=8)
        ax3.set_title("particle survived under the filter(10 degrees)")

        conditions = [
            (self.pp == 0),
            (self.pp == 1),
            (self.pp == 2),

        ]
        choices = [
            'red',  # for 0
            'black',  # for 1
            'green',  # for 3

        ]
        conditions_6 = [
            (self.pp_2 == 0),
            (self.pp_2 == 1),
            (self.pp_2 == 2),
        ]

        # any element not matching one of the above gets 'grey'
        all_colors_3 = np.select(conditions, choices, default='grey')
        all_colors_6 = np.select(conditions_6, choices, default='grey')

        dathreshold_2 = 20
        mask2 = self.aa > dathreshold_2
        ax4.scatter(mapping[axis[1]][mask2],
                    mapping[axis[0]][mask2],
                    s=markerSize,
                    c=all_colors[mask2])
        ax4.legend(handles=handles_2, loc='upper right', fontsize=8)
        ax4.set_title("particle survived by the filter(20 degrees)")
        mask3 = ~np.isin(self.pp, [3, 4, 5])
        mask6 = ~np.isin(self.pp_2, [3, 4, 5])

        ax5.scatter(mapping[axis[1]][mask3],
                    mapping[axis[0]][mask3],
                    s=markerSize,
                    c=all_colors_3[mask3]
                    )
        handles_3 = [
            mpatches.Patch(color='red', label='neutron edep able to find'),
            mpatches.Patch(color='black', label='neutron edep wrongly exclude'),
            mpatches.Patch(color='green', label='muon wrongly included'),
            mpatches.Patch(color='grey', label='other primary particle edep')
        ]
        ax5.legend(handles=handles_3, loc='upper right', fontsize=8)
        ax5.set_title("particle purity and efficiency(10 degrees)")

        ax6.scatter(mapping[axis[1]][mask6],
                    mapping[axis[0]][mask6],
                    s=markerSize,
                    c=all_colors_6[mask6]
                    )

        ax6.legend(handles=handles_3, loc='upper right', fontsize=8)
        ax6.set_title("particle purity and efficiency(20 degrees)")

        fig.savefig(
            self.event.plotpath + '/particle_%s_%s_%s_evt_%d.pdf' % (value, axis, energy, self.event.currentEntry))
        plt.clf()  # important to clear figure
        plt.close()
        # plt.show()









    def hist_dx(self):
        plt.hist(self.ll, range=(0, 4), bins=100)
        plt.xlabel('dx [cm]')
        plt.draw()
        plt.savefig(self.event.plotpath + '/single_edep_dx_evt_%d.pdf' % self.event.currentEntry)
        plt.clf()  # important to clear figure
        plt.close()
        return self.ll  #cm

    def hist_trklength(self):
        all_tracks_length = np.array([])
        for i, track in enumerate(self.event.tracks):
            all_tracks_length = np.append(all_tracks_length, track.length['selfDepo'])

        plt.hist(all_tracks_length, range=(0, 6), bins=100)
        plt.xlabel('track length [cm]')
        plt.draw()
        plt.savefig(self.event.plotpath + '/primary_track_length_evt_%d.pdf' % self.event.currentEntry)
        plt.clf()  # important to clear figure
        plt.close()
        return all_tracks_length  #cm

    def evt_maxdtdr(self):  # per evt has one max dt
        nedeps = len(self.tt)
        max_edep_dr, max_edep_dt = 0, 0
        for iedep in range(nedeps):
            for jedep in range(iedep + 1, nedeps):
                temp_dt = abs(np.float128(self.tt[iedep] - self.tt[jedep]))  # ns
                # unit m
                temp_dr = math.sqrt(pow(np.float128(self.xx[iedep] - self.xx[jedep]), 2) + pow(
                    np.float128(self.yy[iedep] - self.yy[jedep]), 2) + pow(np.float128(self.zz[iedep] - self.zz[jedep]),
                                                                           2))
                if temp_dt > max_edep_dt:
                    max_edep_dt = temp_dt
                if temp_dr > max_edep_dr:
                    max_edep_dr = temp_dr

        return [max_edep_dr * 100, max_edep_dt]  #cm, ns


    def neutron_money_plot(self, neutrino_energy, cone_method, bins: int = 90,
                           thresholds: Optional[Sequence[float]] = None):
        """Plot angular differences for one or multiple energy thresholds."""

        threshold_scan = [None] if thresholds is None else list(thresholds)
        ROOT.gROOT.SetBatch(True)
        canvas = ROOT.TCanvas("c_angle_diff", "Angle Difference", 800, 600)
        legend = ROOT.TLegend(0.58, 0.7, 0.88, 0.88)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        colors = [ROOT.kBlue + 1, ROOT.kRed + 1, ROOT.kGreen + 2,
                  ROOT.kMagenta + 1, ROOT.kOrange + 7]

        histograms = []
        for idx, thr in enumerate(threshold_scan):
            cos_list, counter, angles = self.event.collect_cos(
                cone_method=cone_method ,
                energy_threshold= thr
            )
            if not angles:
                continue
            hist_name = f"h_angle_diff_{idx}"
            hist_title = (
                f"Neutron money plot {neutrino_energy}GeV"
                f" (threshold {self.event._resolve_energy_threshold(thr):.2f} MeV)"
            )
            hist = ROOT.TH1F(hist_name, hist_title, bins, 0, 180)
            for angle in angles:
                hist.Fill(angle)
            color = colors[idx % len(colors)]
            hist.SetLineColor(color)
            hist.SetLineWidth(2)
            draw_option = "HIST" if not histograms else "HIST SAME"
            hist.Draw(draw_option)
            legend.AddEntry(
                hist,
                f"thr {self.event._resolve_energy_threshold(thr):.2f} MeV ({counter} events)",
                "l"
            )
            histograms.append(hist)

        if not histograms:
            print("neutron money plot: no angle data available")
            return None

        histograms[0].GetXaxis().SetTitle("Angle difference between real and reco [deg]")
        histograms[0].GetYaxis().SetTitle("Counts")
        if len(histograms) > 1:
            legend.Draw()

        output_dir = Path('~/frunshiworking/plots').expanduser()
        output_dir.mkdir(parents=True, exist_ok=True)
        if thresholds is None:
            root_path = output_dir / f'neutron_money_plot_{neutrino_energy}GeVneutrino_2_0.2.root'
        else:
            thr_label = "_".join(f"{self.event._resolve_energy_threshold(t):.2f}" for t in threshold_scan)
            thr_label = thr_label.replace('.', 'p')
            root_path = output_dir / f'neutron_money_plot_{neutrino_energy}GeVneutrino_thresh_{thr_label}.root'
        root_file = ROOT.TFile(str(root_path), "RECREATE")
        for hist in histograms:
            hist.Write()
        if len(histograms) > 1:
            legend.Write("legend")
        canvas.Write("angle_difference_canvas")
        root_file.Close()


    def neutron_money_plot_u_overlap_cone_angle(self, neutrino_energy,
                                                cone_method1,
                                                cone_method2,
                                                bins: int = 90):
        """Overlay two muon-neutrino money plots with different cone methods."""

        result1 = self.event.collect_cos(cone_method=cone_method1)
        result2 = self.event.collect_cos(cone_method=cone_method2)

        angles1 = result1[2] if result1 else []
        count1 = result1[1] if result1 else 0
        angles2 = result2[2] if result2 else []
        count2 = result2[1] if result2 else 0

        if not angles1 and not angles2:
            print("neutron money plot overlap: no angle data available for either cone method")
            return None

        ROOT.gROOT.SetBatch(True)
        canvas = ROOT.TCanvas("c_angle_cone_overlap", "Angle Difference Overlap", 800, 600)
        hist_title = f"Muon overlap {neutrino_energy}GeV (cone methods)( OneMuon-OneNeutron Event)"
        hist_1 = ROOT.TH1F("h_angle_cone1", hist_title, bins, 0, 180)
        hist_2 = ROOT.TH1F("h_angle_cone2", hist_title, bins, 0, 180)

        for angle in angles1:
            hist_1.Fill(angle)
        for angle in angles2:
            hist_2.Fill(angle)

        hist_1.SetLineColor(ROOT.kBlue + 1)
        hist_1.SetLineWidth(2)
        hist_2.SetLineColor(ROOT.kGreen + 2)
        hist_2.SetLineWidth(2)

        hist_1.GetXaxis().SetTitle("Angle difference between real and reco [deg]")
        hist_1.GetYaxis().SetTitle("Counts")
        hist_1.Draw("HIST")
        if angles2:
            hist_2.Draw("HIST SAME")

        legend = ROOT.TLegend(0.58, 0.7, 0.88, 0.88)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        if angles1:
            legend.AddEntry(hist_1, f"Cone method : Purity and Efficiency Plot({count1} events)", "l")
        if angles2:
            legend.AddEntry(hist_2, f"Cone method : Brutal Selection ({count2} events)", "l")
        legend.Draw()

        output_dir = Path('~/frunshiworking/plots').expanduser()
        output_dir.mkdir(parents=True, exist_ok=True)
        root_path = output_dir / f'neutron_money_plot_u_cone_overlap_{neutrino_energy}GeV.root'
        root_file = ROOT.TFile(str(root_path), "RECREATE")
        hist_1.Write("cone_method_purity_efficiency")
        hist_2.Write("cone_method_brutal_selection")
        canvas.Write("cone_overlap_canvas")
        root_file.Close()




    # def count_neutron_antinue(self):
    #     events_with
    #     for i in range(1000):
    #         self.Jump(i,0.5)



    def neutron_money_plot_EV(self, neutrino_energy, bins: int = 45):

        angles = self.event.collect_cos_EV()
        if not angles:
            print("neutron money plot: no angle data available")
            return None
        counter = len(angles)
        ROOT.gROOT.SetBatch(True)
        canvas = ROOT.TCanvas("c_angle_diff", "Angle Difference", 800, 600)
        hist_title = f"Neutron Money plot of{neutrino_energy}GeV neutrino ({counter} events)"
        hist = ROOT.TH1F("h_angle_diff", hist_title, bins, 0, 180)
        for angle in angles:
            hist.Fill(angle)
        hist.SetLineColor(ROOT.kBlue + 1)
        hist.SetLineWidth(2)
        hist.GetXaxis().SetTitle("Angle difference between real and reco [deg]")
        hist.GetYaxis().SetTitle("Counts")
        hist.Draw("HIST")
        root_path = f'~/frunshiworking/plots/electron_neutrino/neutron_money_plot_{neutrino_energy}GeVneutrino.root'
        root_file = ROOT.TFile(str(root_path), "RECREATE")
        hist.Write()
        root_file.Close()


    def neutron_money_plot_EV_pureblood(self,neutrino_energy , bins:int = 45):
        # pure blood is defined as use a neutron energy deposition for reconstruction,
        # not necessarily a nearest neutron
        angles = self.event.collect_cos_EV_pureblood()
        if not angles:
            print("neutron money plot: no angle data available")
            return None
        counter = len(angles)
        ROOT.gROOT.SetBatch(True)
        canvas = ROOT.TCanvas("c_angle_diff", "Angle Difference", 800, 600)
        hist_title = f"Neutron Money plot of{neutrino_energy}GeV neutrino ({counter} events)(selecting a neutron deposition for reconstruction)"
        hist = ROOT.TH1F("h_angle_diff", hist_title, bins, 0, 180)
        for angle in angles:
            hist.Fill(angle)
        hist.SetLineColor(ROOT.kBlue + 1)
        hist.SetLineWidth(2)
        hist.GetXaxis().SetTitle("Angle difference between real and reco [deg]")
        hist.GetYaxis().SetTitle("Counts")
        hist.Draw("HIST")
        root_path = f'~/frunshiworking/plots/electron_neutrino/neutron_money_plot_{neutrino_energy}GeVneutrino.root'
        root_file = ROOT.TFile(str(root_path), "RECREATE")
        hist.Write()
        root_file.Close()

    def neutron_money_plot_e_u_overleap(self, muon_plotter, electron_plotter,
                                        neutrino_energy,muon_cone, bins: int = 90):
        """Overlay muon and electron neutrino neutron money plots from separate plotters."""

        muon_data = muon_plotter.event.collect_cos(muon_cone)
        muon_angles = []
        muon_count = 0
        if muon_data:
            _, muon_count, muon_angles = muon_data

        electron_angles = electron_plotter.event.collect_cos_EV()
        electron_count = len(electron_angles)

        if not muon_angles and not electron_angles:
            print("neutron money plot overlap: no angle data available for either sample")
            return None

        ROOT.gROOT.SetBatch(True)
        canvas = ROOT.TCanvas("c_angle_overlap", "Angle Difference Overlap", 800, 600)
        hist_title = f"Neutron Money plot overlap {neutrino_energy}GeV neutrino"
        hist_mu = ROOT.TH1F("h_angle_mu", hist_title, bins, 0, 180)
        hist_e = ROOT.TH1F("h_angle_e", hist_title, bins, 0, 180)

        for angle in muon_angles:
            hist_mu.Fill(angle)
        for angle in electron_angles:
            hist_e.Fill(angle)

        hist_mu.SetLineColor(ROOT.kBlue + 1)
        hist_mu.SetLineWidth(2)
        hist_e.SetLineColor(ROOT.kRed + 1)
        hist_e.SetLineWidth(2)

        hist_mu.GetXaxis().SetTitle("Angle difference between real and reco [deg]")
        hist_mu.GetYaxis().SetTitle("Counts")
        hist_mu.Draw("HIST")
        if electron_angles:
            hist_e.Draw("HIST SAME")

        legend = ROOT.TLegend(0.58, 0.7, 0.88, 0.88)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        if muon_angles:
            legend.AddEntry(hist_mu, f"Muon neutrino ({muon_count} events)", "l")
        if electron_angles:
            legend.AddEntry(hist_e, f"Electron neutrino ({electron_count} events)", "l")
        legend.Draw()

        output_dir = Path('~/frunshiworking/plots').expanduser()
        output_dir.mkdir(parents=True, exist_ok=True)
        root_path = output_dir / f'neutron_money_plot_overlap_{neutrino_energy}GeVneutrino.root'
        root_file = ROOT.TFile(str(root_path), "RECREATE")
        hist_mu.Write("muon_angles")
        hist_e.Write("electron_angles")
        canvas.Write("overlap_canvas")
        root_file.Close()


    def count_oneelectrononeneutronevent_antinue(self):
        right_event_count = 0
        for i in range(1000):
            self.Jump(i, 0.5)
            result = self.event.selectelectronneutronevent(i)
            right_event_count += result
        print(f"there are {right_event_count} electrons neutron events ")

    def plot_neutron_earliest_depo_energy_spectrum(self):
        n_events = 9000
        earliest_depo_list = []
        for i in range(n_events):
            if i % 100 == 0:
                print(f"processing events {i}/9000")
            self.Jump(i, 0)
            earliest_depo = self.event.find_energy_of_earliest_neutron_energy_deposition()
            if earliest_depo is None:
                continue
            earliest_depo_list.append(float(earliest_depo))

        emin, emax = 0.0,6.0
        nbins = 60
        h = ROOT.TH1D(
            "hEarliestNeutronDepoE",
            "Earliest neutron energy deposition spectrum;Energy;Counts",
            nbins, emin, emax
        )
        h.SetLineWidth(2)

        for e in earliest_depo_list:
            h.Fill(e)

        c = ROOT.TCanvas("cEarliestNeutronDepoE", "Earliest neutron depo energy", 900, 650)
        h.Draw("HIST")
        c.Update()
        outfile = os.path.expanduser("~/frunshiworking/plots/earliest_depo_energy_spectrum.root")

        # Write histogram + canvas explicitly (more reliable for ROOT files)
        f = ROOT.TFile(outfile, "RECREATE")
        h.Write()
        c.Write()
        f.Close()

        print(f"[plot] Saved histogram+canvas to: {outfile}")

    def plot_earliest_depo_length_spectra(
            self,
            max_events: int = 10000,
            energy_threshold: Optional[float] = None,
            nbins: int = 100,
    ):
        """Overlay earliest neutron and low-energy proton deposition distances with linear binning."""
        import os
        import array
        import math
        import ROOT

        available_events = getattr(self.event, 'nEntry', 0)
        n_events = min(max_events, available_events)
        if n_events <= 0:
            print("earliest length spectra: no events available")
            return None

        resolved_threshold = self.event.energy_threshold if energy_threshold is None else energy_threshold
        neutron_distances = []
        low_proton_distances = []

        # --- Data Collection ---
        for i in range(n_events):
            if i % 500 == 0:
                print(f"Processing event {i}/{n_events}...")
            self.Jump(i, 0)

            earliest_distance = self.event.find_distance_of_earliest_neutron_energy_deposition(
                energy_threshold=energy_threshold
            )
            if earliest_distance is not None:
                neutron_distances.append(float(earliest_distance))

            distances = self.event.collect_low_energy_proton_deposition_distances(
                energy_threshold=energy_threshold
            )
            if distances:
                low_proton_distances.extend(float(dist) for dist in distances)

        if not neutron_distances and not low_proton_distances:
            print("earliest length spectra: no events passed selection.")
            return None

        # --- 1. Purely Linear Bin Edges ---
        # Define the range (0 to 10.0 meters as in your original log_end)
        max_range = 10.0
        bin_edges = []
        for i in range(nbins + 1):
            edge = (max_range * i) / nbins
            bin_edges.append(edge)

        edges_array = array.array('d', bin_edges)

        # --- 2. Plotting Setup ---
        canvas = ROOT.TCanvas("cEarliestDepo", "Deposition Distances", 900, 650)
        canvas.SetLogy(1)
        # Note: Y-axis remains log for dynamic range, but X-axis is now linearly binned

        legend = ROOT.TLegend(0.55, 0.7, 0.88, 0.88)
        legend.SetBorderSize(0)

        hist_specs = [
            ("hNeutron", "Neutron", neutron_distances, ROOT.kAzure + 2),
            ("hProton", "Proton", low_proton_distances, ROOT.kOrange + 7),
        ]

        drawn_hists = []
        current_max = 0.0

        for name, title, data, color in hist_specs:
            if not data: continue

            # Use the linear edges_array
            h = ROOT.TH1D(name, f";Distance to vertex [m];Normalized Density [1/m]", len(edges_array) - 1, edges_array)
            for d in data: h.Fill(d)

            # Normalize to probability density
            integral = h.Integral()
            if integral > 0:
                h.Scale(1.0 / integral, "width")

            h.SetLineColor(color)
            h.SetLineWidth(2)
            h.SetStats(0)

            current_max = max(current_max, h.GetMaximum())
            drawn_hists.append(h)
            legend.AddEntry(h, title, "l")

        # --- 3. Final Draw & Save ---
        for i, h in enumerate(drawn_hists):
            h.SetMaximum(current_max * 2.0)
            # You can adjust RangeUser here to see more or less of the linear scale
            h.GetXaxis().SetRangeUser(0.0, 5.0)
            h.Draw("HIST" if i == 0 else "HIST SAME")

        legend.Draw()
        canvas.Update()

        out_path = os.path.expanduser("~/frunshiworking/plots/earliest_depo_length_overlay_linear")
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        canvas.SaveAs(f"{out_path}.pdf")

        f_out = ROOT.TFile(f"{out_path}.root", "RECREATE")
        for h in drawn_hists: h.Write()
        canvas.Write()
        f_out.Close()
    def nearest_earliest_angle_distribution(
        self,
        neutrino_energy,
        p,
        thresholds=None,
        bins: int = 50):
        """Overlay the angle distribution for multiple neutron energy thresholds."""

        ROOT.gROOT.SetBatch(True)
        n_events = min(1000, self.event.nEntry)
        original_threshold = self.event.energy_threshold

        if not thresholds:
            threshold_scan = [None]
        else:
            # Preserve user order while removing duplicates
            seen = set()
            threshold_scan = []
            for thr in thresholds:
                if thr in seen:
                    continue
                seen.add(thr)
                threshold_scan.append(thr)

        canvas = ROOT.TCanvas("c_nearest_earliest", "Angle Distribution Canvas", 800, 600)
        legend = ROOT.TLegend(0.58, 0.7, 0.88, 0.88)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        colors = [
            ROOT.kBlack,
            ROOT.kBlue + 1,
            ROOT.kRed + 1,
            ROOT.kGreen + 2,
            ROOT.kOrange + 7,
            ROOT.kViolet + 5,
        ]

        histograms = []
        for idx, thr in enumerate(threshold_scan):
            resolved_threshold = self.event._resolve_energy_threshold(thr)
            self.event.energy_threshold = resolved_threshold
            angle_list = []

            for i in range(n_events):
                self.Jump(i, resolved_threshold)
                if i%100==0:
                    print(f"processing events {i}/1000")
                result = self.event.nearest_vs_earliest(p)
                if not result:
                    continue
                angle_distribution = result.get("angle_distribution")
                if angle_distribution is None:
                    continue
                angle_list.append(angle_distribution)

            if not angle_list:
                continue

            hist_name = f"angle_hist_thr_{idx}"
            hist_title = (
                f"Angle Distribution (nearest vs earliest) {neutrino_energy:.1f}GeV"
                f" with different energy threshold"
            )
            hist = ROOT.TH1F(hist_name, hist_title, bins, 0, 180)
            for angle in angle_list:
                hist.Fill(angle)
            color = colors[idx % len(colors)]
            hist.SetLineColor(color)
            hist.SetLineWidth(2)
            hist.SetMaximum(5000)
            draw_option = "HIST" if not histograms else "HIST SAME"
            hist.Draw(draw_option)
            legend.AddEntry(
                hist,
                f"{resolved_threshold:.2f} MeV ",
                "l"
            )
            histograms.append((hist, resolved_threshold))

        self.event.energy_threshold = original_threshold

        if not histograms:
            print("nearest/earliest angle distribution: no events passed selection")
            return None

        histograms[0][0].GetXaxis().SetTitle("angle [degrees]")
        histograms[0][0].GetYaxis().SetTitle("Counts")
        # Log axis cannot include zero, so clip the visible minimum to a small positive value
        histograms[0][0].GetXaxis().SetRangeUser(1e-2, histograms[0][0].GetXaxis().GetXmax())
        # canvas.SetLogx()
        canvas.SetLogy()
        if len(histograms) > 1:
            legend.Draw()

        thr_label = f"{histograms[0][1]:.2f}"
        if len(histograms) > 1:
            thr_label = "_".join(
                f"{thr_val:.2f}" for _, thr_val in histograms
            )
        thr_label = thr_label.replace('.', 'p')

        output_filename = (
            f"angle_between_earliest_nearest_{neutrino_energy:.1f}GeV_thr_{thr_label}.root"
        )
        root_file = ROOT.TFile(output_filename, "RECREATE")
        for hist, _ in histograms:
            hist.Write()
        if len(histograms) > 1:
            legend.Write("legend")
        canvas.Write("angle_distribution_canvas")
        root_file.Close()

        output_filename_pdf = (
            f"angle_between_earliest_nearest_{neutrino_energy:.1f}GeV_thr_{thr_label}.pdf"
        )
        canvas.SaveAs(output_filename_pdf)
        print(f"Visual plot saved to {output_filename_pdf}")







    def plot_angle_distribution(self):
        """
        this code explore the potential correlation between the muon angle and neutron angle
        """
        muon_distribution = []
        neutron_distribution = []
        relative_distribution = []
        n_events = min(10, self.event.nEntry)
        for idx in range(n_events):
            self.event.Jump(idx)
            neutron_direction = self.event.read_neutron_direction()
            muon_direction,mom = self.event.read_muon_direction_and_amplitude()
            if neutron_direction is None or muon_direction is None or np.all(neutron_direction) == 0 or np.all(
                    muon_direction) == 0:
                print(f"Event {self.event.currentEntry}: no neutron direction, skipping arrow!!!.")
                continue
            neutron_theta = np.arctan(
                (neutron_direction[0] ** 2 + neutron_direction[1] ** 2) ** 0.5 / neutron_direction[2])
            if neutron_theta <= 0:
                neutron_theta = 3.1415 + neutron_theta
            if neutron_direction[0] > 0:
                if neutron_direction[1] > 0:
                    neutron_phi = np.arctan(neutron_direction[1] / neutron_direction[0])
                else:
                    neutron_phi = 6.2832 + np.arctan(neutron_direction[1] / neutron_direction[0])
            elif neutron_direction[0] < 0:
                if neutron_direction[1] > 0:
                    neutron_phi = 3.14159 + np.arctan(neutron_direction[1] / neutron_direction[0])
                else:
                    neutron_phi = 3.14159 + np.arctan(neutron_direction[1] / neutron_direction[0])
            neutron_distribution.append([neutron_theta, neutron_phi])

            muon_theta = np.arctan((muon_direction[0] ** 2 + muon_direction[1] ** 2) ** 0.5 / muon_direction[2])
            muon_phi = 0
            if muon_theta <= 0:
                muon_theta = 3.1415 + muon_theta
            if muon_direction[0] > 0:
                if muon_direction[1] > 0:
                    muon_phi = np.arctan(muon_direction[1] / muon_direction[0])
                else:
                    muon_phi = 6.2832 + np.arctan(muon_direction[1] / muon_direction[0])
            elif muon_direction[0] < 0:
                if muon_direction[1] > 0:
                    muon_phi = 3.14159 + np.arctan(muon_direction[1] / muon_direction[0])
                else:
                    muon_phi = 3.14159 + np.arctan(muon_direction[1] / muon_direction[0])

            muon_distribution.append([muon_theta, muon_phi])
            relative_theta = muon_theta - neutron_theta
            relative_phi = muon_phi - neutron_phi
            if relative_theta < 0:
                relative_theta = relative_theta + 2 * np.pi
            if relative_phi < 0:
                relative_phi = relative_phi + 2 * np.pi
            colors = 'red' if self.event.selectmuonneutronevent(idx) == 1 else 'blue'
            relative_distribution.append((relative_theta, relative_phi, colors))

        # Convert to arrays and degrees
        neutron_thetas, neutron_phis = np.array(neutron_distribution).T
        muon_thetas, muon_phis = np.array(muon_distribution).T

        rel_thetas = [row[0] for row in relative_distribution]
        rel_phis = [row[1] for row in relative_distribution]
        colors_rel = [row[2] for row in relative_distribution]

        neutron_thetas = np.degrees(neutron_thetas)
        neutron_phis = np.degrees(neutron_phis)
        muon_thetas = np.degrees(muon_thetas)
        muon_phis = np.degrees(muon_phis)
        rel_thetas = np.degrees(rel_thetas)
        rel_phis = np.degrees(rel_phis)

        # Compute true separation angles
        separations = []
        for idx in range(n_events):
            self.event.Jump(idx)
            n_dir = self.event.read_neutron_direction()
            m_dir ,mom= self.event.read_muon_direction_and_amplitude()
            if n_dir is None or m_dir is None or np.all(n_dir) == 0 or np.all(m_dir) == 0:
                print(f"Event {self.event.currentEntry}: no neutron direction, skipping arrow.")
                continue
            cosang = np.dot(n_dir, m_dir) / (np.linalg.norm(n_dir) * np.linalg.norm(m_dir))
            cosang = np.clip(cosang, -1, 1)
            separations.append(np.degrees(np.arccos(cosang)))

        # Plotting
        fig, axs = plt.subplots(2, 2, figsize=(12, 10))
        axs[0, 0].scatter(muon_thetas, muon_phis, s=5)
        axs[0, 0].set_xlabel('Muon θ (°)')
        axs[0, 0].set_ylabel('Muon φ (°)')
        axs[0, 0].set_title('Muon Angular Distribution')

        axs[0, 1].scatter(neutron_thetas, neutron_phis, s=5)
        axs[0, 1].set_xlabel('Neutron θ (°)')
        axs[0, 1].set_ylabel('Neutron φ (°)')
        axs[0, 1].set_title('Neutron Angular Distribution')

        axs[1, 0].scatter(rel_thetas, rel_phis, c=colors_rel, s=5)
        axs[1, 0].set_xlabel('Δθ (°)')
        axs[1, 0].set_ylabel('Δφ (°)')
        axs[1, 0].set_title('Relative Angular Difference')
        red_patch = mpatches.Patch(color='red', label='Other Event')
        blue_patch = mpatches.Patch(color='blue', label='Muon-Neutron Event')
        #axs[1, 0].legend(handle)
        s = [red_patch, blue_patch], loc = ('upper right')  # or adjust `loc` as needed

        bins = np.arange(0, max(separations) + 5, 5)
        axs[1, 1].hist(separations, bins=bins)
        axs[1, 1].set_xlabel('Separation Angle (°)')
        axs[1, 1].set_ylabel('Number of Events')
        axs[1, 1].set_title('Histogram of Separation Angles (5° bins)')

        plt.tight_layout()
        plt.savefig('../angle_distribution.pdf')



    def plot_length_inverse_hist_with_angle_curve(self):
        """
        Combine the first 10 events, collect all [length, energy, angle] rows,
        then:
          • Draw an energy histogram (bin = 0.1 MeV) with bar heights = counts
            on the left y-axis.
          • Compute the mean angle in each bin, and plot ⟨angle⟩ vs. bin_center
            on a secondary right y-axis as a line.
        Saves figure as: ./plots/energy_hist_angle_curve.png
        """

        # 1) Collect (energy, angle) from the first up to 10 events
        n_events = min(1000, self.event.nEntry)
        all_energy = []
        all_angle = []
        all_length = []
        all_length_log = []
        all_length_inversesquare = []

        for idx in range(n_events):
            self.event.Jump(idx)
            data = self.event.neutron_edep_information()
            arr = np.asarray(data)

            # If not N×3, skip
            if arr.ndim != 2 or arr.shape[1] != 3:
                print(
                    f"Warning: neutron_edep_information() for event {idx} "
                    f"returned shape {arr.shape}. Skipping."
                )
                continue

            # arr[:,1] = energy, arr[:,2] = angle
            all_length.extend(arr[:, 0].tolist())
            all_energy.extend(arr[:, 1].tolist())
            all_angle.extend(arr[:, 2].tolist())
            all_length_log.extend(np.log(arr[:, 0]).tolist())
            all_length_inversesquare.extend((arr[:, 0] ** (-2)).tolist())

        # 2) If no data was collected, bail out
        if len(all_energy) == 0:
            print("No valid (energy, angle) points found in the first 10 events; no histogram made.")
            return

        all_energy = np.asarray(all_energy)
        all_angle = np.asarray(all_angle)
        all_length = np.asarray(all_length)
        all_length_inversesquare = np.asarray(all_length_inversesquare)
        all_length_log = np.asarray(all_length_log)
        #
        # # 3) Build bins of width 0.1 MeV
        # emin = all_length.min()
        # emax = all_energy.max()
        # # emin=0.5
        # emax=1.5

        #inverse case
        emin = 0
        emax = 0.075

        bin_start = emin
        bin_stop = emax
        # bin_start = np.floor(emin * 10) / 10
        # bin_stop  = np.ceil(emax * 10) / 10
        bins = np.arange(bin_start, bin_stop, 0.005)
        counts, bin_edges = np.histogram(all_length_inversesquare, bins=bins)
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

        # 4) For each bin, compute angle‐mean
        means = []
        for i in range(len(bin_centers)):
            low = bin_edges[i]
            high = bin_edges[i + 1]
            mask = (all_length_inversesquare >= low) & (all_length_inversesquare < high)

            if np.any(mask):
                angs = all_angle[mask]
                means.append(np.mean(angs))
            else:
                means.append(np.nan)

        # 5) Plot histogram on left y-axis
        fig, ax = plt.subplots(figsize=(8, 5), dpi=100)
        bars = ax.bar(
            bin_centers,
            counts,
            width=0.01,
            align="center",
            edgecolor="black",
            label="Counts"
        )
        ax.set_xlabel("inverse square length")
        ax.set_ylabel("number of points")
        ax.set_title(f"inverse square length histogram & ⟨angle⟩ curve (first {n_events} events)")

        # 6) Create a secondary y-axis for ⟨angle⟩
        ax2 = ax.twinx()
        ax2.plot(
            bin_centers,
            means,
            color="red",
            marker="o",
            linestyle="-",
            label="⟨angle⟩"
        )
        ax2.set_ylabel("mean angle")
        ax2.grid(False)  # turn off grid on the second axis

        # 7) Put a legend for the angle curve on the right
        ax2.legend(loc="upper right")

        # 8) Save to ./plots
        outdir = self.event.plotpath
        os.makedirs(outdir, exist_ok=True)
        outpath = os.path.join(outdir, "inversesquarelength_hist_angle_curve.png")
        fig.tight_layout()
        fig.savefig(outpath)
        plt.close(fig)
        print(f"  → saved combined histogram + ⟨angle⟩ curve to {outpath}")



    def stacked_energy_distribution_histogram(self, *, return_per_event: bool = False):
        n_events = self.event.nEntry
        proton_total_energy = 0
        neutron_total_energy = 0
        lepton_total_energy = 0
        pion_total_energy = 0
        gamma_total_energy = 0
        other_total_energy = 0
        per_event_breakdown = None
        if return_per_event:
            per_event_breakdown = {
                'proton': [],
                'neutron': [],
                'lepton': [],
                'pion': [],
                'gamma': [],
                'other': [],
            }
        for i in range(n_events):
            # Use the same entryNo/threshold ordering as every other Jump call.
            self.Jump(i, 0.5)
            proton_energy, neutron_energy, lepton_energy, pion_energy, gamma_energy, other_energy = self.event.read_particle_energy_fraction()
            proton_total_energy += proton_energy
            neutron_total_energy += neutron_energy
            lepton_total_energy += lepton_energy
            pion_total_energy += pion_energy
            gamma_total_energy += gamma_energy
            other_total_energy += other_energy
            if per_event_breakdown is not None:
                per_event_breakdown['proton'].append(proton_energy)
                per_event_breakdown['neutron'].append(neutron_energy)
                per_event_breakdown['lepton'].append(lepton_energy)
                per_event_breakdown['pion'].append(pion_energy)
                per_event_breakdown['gamma'].append(gamma_energy)
                per_event_breakdown['other'].append(other_energy)

        totals = (
            proton_total_energy,
            neutron_total_energy,
            lepton_total_energy,
            pion_total_energy,
            gamma_total_energy,
            other_total_energy,
        )
        if per_event_breakdown is not None:
            return totals, per_event_breakdown
        return totals





    def antinumu_neutron_distance_spectrum(self,low_neutron_energy, high_neutron_energy):
        self.Jump(0,0.5)
        distance_per_edep = []
        mm2cm =0.1
        for i in range(9000):
            self.event.Jump(i)
            if i % 500==0:
                print("processing event", i)
            neutron_energy = self.event.read_neutron_energy()
            if neutron_energy > high_neutron_energy or neutron_energy < low_neutron_energy:
                continue
            for j in range(self.event.tracks.size):
                ancestor_pdg = self.event.tracks[self.event.tracks[j].association["ancestor"]].GetPDGCode()
                if ancestor_pdg ==2112:
                    for depo_index in self.event.tracks[j].association.get('depoList', []):
                        if depo_index < 0 or depo_index >= len(self.event.depos):
                            continue
                        if  self.event.depos[depo_index].GetEnergyDeposit()<0.5:
                            continue
                        depo = self.event.depos[depo_index]
                        start, stop = depo.GetStart(), depo.GetStop()
                        midpoint = np.array([
                            (start.X() + stop.X()) / 2 * mm2cm,
                            (start.Y() + stop.Y()) / 2 * mm2cm,
                            (start.Z() + stop.Z()) / 2 * mm2cm,
                        ])
                        distance = np.linalg.norm(midpoint)
                        distance_per_edep.append(distance)
        print(f"the total number of energy deposition recorded are{len(distance_per_edep)}")
        return distance_per_edep







    # def electron_neutrino_direction_reconstruction:
    #     """loop over 1000 events of one certain energy of anti electron neutrino
    #     in each event, read electron neutrinno direction, positron direction, neutron and other particle direction
    #     then use fafhalkfjaj
    #
    #

        # """


    def plot_neutron_distance_spectrum(self,low, high):
        """Fill a ROOT histogram from the antinumu neutron distance spectrum and persist it."""
        distances = self.antinumu_neutron_distance_spectrum(low,high)
        if not distances:
            print("No neutron energy deposit distances found; nothing to plot.")
            return

        min_distance = min(distances)
        # max_distance = max(distances)
        max_distance = 1800
        nbins = 100
        hist = ROOT.TH1F(
            "h_neutron_distance",
            f"Antinumu neutron distance spectrum({low}~{high}MeV Neutron);distance from origin (cm);entries",
            nbins,
            min_distance,
            max_distance,
        )
        print(f"each bin has width{(max_distance-min_distance)/nbins}")
        for value in distances:
            hist.Fill(value)

        total_hits = len(distances)
        sorted_distances = np.sort(np.array(distances, dtype=float))
        percentages = np.arange(1, total_hits + 1, dtype=float) / total_hits * 100

        cdf_graph = ROOT.TGraph(
            total_hits,
            array.array('d', sorted_distances.tolist()),
            array.array('d', percentages.tolist())
        )
        cdf_graph.SetName("g_neutron_distance_cdf")
        cdf_graph.SetTitle("Fraction of events shorter than distance;distance from origin (cm);fraction of events (%)")
        cdf_graph.SetLineColor(ROOT.kRed)
        cdf_graph.SetLineWidth(2)

        quantiles_to_report = [50, 80, 90, 95, 99]
        for target in quantiles_to_report:
            idx = np.searchsorted(percentages, target, side="left")
            idx = min(idx, total_hits - 1)
            distance_at_percentile = sorted_distances[idx]
            print(f"{target}% of neutron deposits are within {distance_at_percentile:.2f} cm")

        output_dir = Path("~/frunshiworking/plots").expanduser()
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"antinumu_continuous_spectrum_neutron_distance_spectrum_{low}to{high}mevneutron.root"
        root_file = ROOT.TFile(str(output_path), "RECREATE")

        canvas = ROOT.TCanvas("c_neutron_distance", "Neutron distance y-logx", 900, 800)

        pad_hist = ROOT.TPad("pad_hist", "", 0, 0.35, 1, 1)
        pad_hist.SetBottomMargin(0.15)
        pad_hist.Draw()
        pad_hist.cd()
        hist.GetXaxis().SetRangeUser(0, 1800)
        hist.Draw()
        hist_xaxis = hist.GetXaxis()
        hist_xaxis.SetTitle("distance from origin (cm)")
        hist_xaxis.SetLabelSize(0.045)
        hist_xaxis.SetTitleSize(0.05)
        hist_xaxis.SetTitleOffset(1.1)
        hist_xaxis.SetNdivisions(510)

        canvas.cd()
        pad_cdf = ROOT.TPad("pad_cdf", "", 0, 0, 1, 0.35)
        pad_cdf.SetTopMargin(0.1)
        pad_cdf.SetBottomMargin(0.25)
        pad_cdf.Draw()
        pad_cdf.cd()
        cdf_graph.Draw("AL")
        cdf_graph.GetXaxis().SetLimits(min_distance, max_distance)
        cdf_graph.GetYaxis().SetRangeUser(0, 100)
        cdf_graph.GetYaxis().SetNdivisions(505)
        cdf_graph.GetYaxis().SetTitleOffset(0.7)
        cdf_graph.GetYaxis().SetTitleSize(0.09)
        cdf_graph.GetYaxis().SetLabelSize(0.08)
        cdf_xaxis = cdf_graph.GetXaxis()
        cdf_xaxis.SetTitle("distance from vertex (cm)")
        cdf_xaxis.SetTitleSize(0.12)
        cdf_xaxis.SetLabelSize(0.1)
        cdf_xaxis.SetTitleOffset(1.0)

        canvas.Update()

        hist.Write()
        cdf_graph.Write()
        canvas.Write()
        root_file.Close()
        print(f"Saved neutron distance histogram to {output_path}")






    def plot_energy_vs_angle(self):

        # How many events to combine:
        n_events = min(10, self.event.nEntry)

        # Lists to accumulate every event’s points
        all_length = []
        all_energy = []
        all_angle = []

        for idx in range(n_events):
            # 1) Jump into the idx-th event
            self.event.Jump(idx)

            # 2) Retrieve the raw data array (shape should be (N,3))
            data = self.event.neutron_edep_information()
            arr = np.asarray(data)

            # 3) Verify it’s a 2D array with 3 columns; skip if not
            if arr.ndim != 2 or arr.shape[1] != 3:
                print(
                    f"Warning: neutron_edep_information() returned shape {arr.shape} "
                    f"for event {idx}. Skipping."
                )
                continue

            # 4) Unpack length/energy/angle columns
            length = arr[:, 0]
            energy = arr[:, 1]
            angle = arr[:, 2]

            # 5) Append into our “big” lists
            all_length.extend(length.tolist())
            all_energy.extend(energy.tolist())
            all_angle.extend(angle.tolist())

        # If nothing was collected (all events skipped), warn & return
        if len(all_length) == 0:
            print("No valid points found in the first 10 events; no plot made.")
            return

        # Convert to numpy arrays for plotting
        all_length = np.asarray(all_length)
        all_energy = np.asarray(all_energy)
        all_angle = np.asarray(all_angle)

        # 6) Single scatter: angle on X, energy on Y, color by length
        plt.figure(figsize=(7, 5), dpi=100)
        sc = plt.scatter(
            all_angle,
            all_energy,
            c=all_length,
            cmap="viridis",
            s=20,
            edgecolors="",
        )
        cbar = plt.colorbar(sc)
        cbar.set_label("length")
        plt.xlabel("angle")
        plt.ylabel("energy")
        plt.title(f"Combined energy vs. angle (first {n_events} events)")

        # 7) Ensure output directory exists and save
        outdir = self.event.plotpath
        os.makedirs(outdir, exist_ok=True)
        outpath = os.path.join(outdir, "energy_vs_angle_combined.png")
        plt.tight_layout()
        plt.savefig(outpath)
        plt.close()
        print(f"  → saved combined plot to {outpath}")

    def plot_length_vs_angle(self):

        # How many events to combine:
        n_events = min(10, self.event.nEntry)

        # Lists to accumulate every event’s points
        all_length = []
        all_energy = []
        all_angle = []
        all_length_log = []
        all_length_squareinverse = []

        for idx in range(n_events):
            # 1) Jump into the idx-th event
            self.event.Jump(idx)

            # 2) Retrieve the raw data array (shape should be (N,3))
            data = self.event.neutron_edep_information()
            arr = np.asarray(data)

            # 3) Verify it’s a 2D array with 3 columns; skip if not
            if arr.ndim != 2 or arr.shape[1] != 3:
                print(
                    f"Warning: neutron_edep_information() returned shape {arr.shape} "
                    f"for event {idx}. Skipping."
                )
                continue

            # 4) Unpack length/energy/angle columns
            length = arr[:, 0]
            energy = arr[:, 1]
            angle = arr[:, 2]

            # 5) Append into our “big” lists
            all_length.extend(length.tolist())
            all_energy.extend(energy.tolist())
            all_angle.extend(angle.tolist())
            all_length_log.extend(log(length).tolist())
            all_length_squareinverse.extend((length ** (-2).tolist))

        # If nothing was collected (all events skipped), warn & return
        if len(all_length) == 0:
            print("No valid points found in the first 10 events; no plot made.")
            return

        # Convert to numpy arrays for plotting
        all_length = np.asarray(all_length)
        all_energy = np.asarray(all_energy)
        all_angle = np.asarray(all_angle)

        # 6) Single scatter: angle on X, energy on Y, color by length
        plt.figure(figsize=(7, 5), dpi=100)
        sc = plt.scatter(
            all_angle,
            all_energy,
            c=all_length,
            cmap="viridis",
            s=20,
            edgecolors="",
        )
        cbar = plt.colorbar(sc)
        cbar.set_label("length")
        plt.xlabel("angle")
        plt.ylabel("energy")
        plt.title(f"Combined energy vs. angle (first {n_events} events)")

        # 7) Ensure output directory exists and save
        outdir = self.event.plotpath
        os.makedirs(outdir, exist_ok=True)
        outpath = os.path.join(outdir, "energy_vs_angle_combined.png")
        plt.tight_layout()
        plt.savefig(outpath)
        plt.close()
        print(f"  → saved combined plot to {outpath}")

    #---------------------------------------------
    def DrawROOT(self, dim2d='yz', markerSize=0.2):
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetMarkerStyle(24)
        ROOT.gStyle.SetMarkerSize(markerSize)
        colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kMagenta]
        nColor = len(colors)

        c1 = TCanvas("c1", self.event.vertex.GetReaction(), 800, 800)

        # canvas of 5m x 5m
        dummy = TH2F("dummy", "", 100, -5, 5, 100, -5, 5)
        dummy.GetXaxis().SetTitle('[m]')
        dummy.GetYaxis().SetTitle('[m]')
        dummy.Draw()

        mm2m = 0.001
        txts = []
        markers = []
        txtX = 0.15
        txtY = 0.85
        txtSize = 0.03
        # nDepo = 0
        for i, track in enumerate(self.event.tracks):
            depoList = track.association['depoList']
            ancestor = track.association['ancestor']
            color = colors[ancestor % nColor]

            if ancestor == i:
                txt = TLatex()
                txt = txt.DrawLatexNDC(txtX, txtY, track.GetName())
                txt.SetTextColor(color)
                txt.SetTextSize(txtSize)
                txts.append(txt)
                txtY -= txtSize

            for j, di in enumerate(depoList):
                depo = self.event.depos[di]
                x = (depo.GetStart().X() + depo.GetStop().X()) / 2 * mm2m
                y = (depo.GetStart().Y() + depo.GetStop().Y()) / 2 * mm2m
                z = (depo.GetStart().Z() + depo.GetStop().Z()) / 2 * mm2m
                mapping = {'x': x, 'y': y, 'z': z}

                m = TMarker(mapping[dim2d[1]], mapping[dim2d[0]], 24)
                m.SetMarkerColor(color)

                m.Draw()
                markers.append(m)
                # nDepo += 1

        # print('depo points drawn: ', nDepo, '| total depo: ', self.event.depos.size)

        ROOT.gPad.Update()
        return c1

    def event_type(self):
        """a function to understand what types of events we are dealing with, classified by the number of hadrons"""
        n_events =1000
        neutron_include = 0
        muon_include = 0
        proton_include = 0
        multi_neutron = 0
        multi_muon = 0
        multi_proton = 0
        multi_proton_neutron = 0
        for idx in range(n_events):

            self.event.Jump(idx)
            trig_proton = self.event.SelectProtonEvent()
            trig_neutron = self.event.SelectNeutronEvent()
            trig_muon = self.event.SelectMuonEvent()
            if trig_proton > 0:
                proton_include += 1
                if trig_proton > 1 :
                    multi_proton += 1

            if trig_muon > 0:
                muon_include += 1
                if trig_muon > 1 :
                    multi_muon += 1

            if trig_neutron > 0 :
                neutron_include += 1
                if trig_neutron > 1 :
                    multi_neutron += 1

            if trig_neutron > 1 and trig_proton > 1 :
                multi_proton_neutron += 1

        print(f"among 1000 events, there are {muon_include}events with muon, {neutron_include}events with neutron, {proton_include}events with proton ")
        print(f"among 1000 events, there are {multi_muon}events with more than 1 muon")
        print(f"among 1000 events, there are {multi_neutron}events with more than 1 neutron")
        print(f"among 1000 events, there are {multi_proton}events with more than 1 proton")
        print(f"among 1000 events, there are {multi_proton_neutron}events with more than 1 neutron and 1 proton")

if __name__ == "__main__":
    event = Event(sys.argv[1])
    p = Plotter(event)
    #p.Next()
    #p.Draw('yz')
    # c1 = p.DrawROOT('xz', 0.2)
    # input('press a key to continue ...')
    for entry in range(20):
        p.event.currentEntry = entry
        #    p.Draw('yz')
        p.effeciency_plot()
    # p.plot_length_inverse_hist_with_angle_curve()
