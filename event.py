# The standard library imports below reflect historical debugging helpers that
# may not all be required by the current analysis routines. Keep them grouped
# so unused ones can be trimmed as workflows evolve.
from functools import total_ordering
from subprocess import check_output
from traceback import print_stack
from xml.sax import parse
import ROOT
from ROOT import TG4Event, TFile, TChain
from collections import Counter
from typing import Tuple, Optional
import sys, os
import numpy as np
import random
import math
from pathlib import Path
import matplotlib.pyplot as plt
from typing import Tuple, List
from numpy import floating
from numpy.f2py.crackfortran import true_intent_list
from scipy.optimize import curve_fit
from scipy.stats import norm
import noise_generator


class Event:

    # neutrino_type = 'antinue'
    neutrino_type = 'antinumu'
    # neutrino_type = 'nue'
    # neutrino_type = 'numu'

    true_cone_angle_list_antinumu_brutal_95= [58,22,12,13,14,15,17,14,14,6 ,6]
    true_cone_angle_list_antinue           = [50,50,50,43,42,55,46,49,57,57,57] # brutal selection angle
    true_cone_angle_list_nue               = [25,25,26,30,33,39,22,35,27,45,45]
    true_cone_angle_list_numu_95percent    = [13,9 ,10,10,9 ,9 ,9 ,9 ,8 ,3 ,5]

    cone_applying_for_electron = true_cone_angle_list_antinue
    cone_applying_for_muon = true_cone_angle_list_numu_95percent
    time_interval = 0 # macro second
    noise_number = 0.75 * time_interval
    noise_points = 0

    energy_threshold = 0.6 # the energy that a deposition is visible in the detector
    proton_trackable_length = 70 # 70 mev for current assumption on whether proton is identifiable or not



    # we here try different sphere for different reconstruction part
    # namely 5 meter radius for energy reconstruction and 2 meter radius for direction reconstruction
    inspection_length_energy_side = 0.5
    inspection_length_direction_side = 1.5
    minimum_length = 4 # centimeters



    def __init__(self, fileName, evgen='Genie'):
        print("event: initialization")
        self.fileName = fileName
        self.evgen = evgen
        self.ReadTree()
        self.currentEntry = 0
        # create a folder to store plots
        self.plotpath = "./plots"
        if not os.path.exists(self.plotpath):
            os.makedirs(self.plotpath)
            print("plotpath '" + self.plotpath + "' did not exist. It has been created!")

        # Track metadata about the most recent energy-deposition candidate used
        # during reconstruction so downstream diagnostics can inspect the origin.
        self._last_first_edep_source = None
        self._last_depo_outside_source = None
        self.last_reco_failure_reason = None

        # Provide access to helper modules used by optional utilities.
        self.noise_generator = noise_generator
        self._radiative_noise_sampler = None  # cache to avoid reopening ROOT histograms
        # Cache the last reconstruction results so downstream comparisons can
        # inspect noise-on/noise-off performance without re-running the logic.
        self.neutron_direction_reconstruction_result_with_noise = None
        self.neutron_direction_reconstruction_result_without_noise = None
        self.neutrino_direction_reconstruction_result_with_noise = None
        self.neutrino_direction_reconstruction_result_without_noise = None

    def _resolve_energy_threshold(self, override: Optional[float]) -> float:
        """Return the override value or fall back to the class-wide threshold."""
        return self.energy_threshold if override is None else override

    # ------------------------
    def ReadTree(self):
        # self.rootFile = TFile(self.fileName)
        self.simTree = TChain("EDepSimEvents")
        self.simTree.Add(self.fileName)
        self.nEntry = self.simTree.GetEntries()

        # Only Genie has gRooTracker, Marley doesn't
        if self.evgen == 'Genie':
            self.genieTree = TChain("DetSimPassThru/gRooTracker")
            self.genieTree.Add(self.fileName)
            # self.genieTree = self.rootFile.Get("DetSimPassThru/gRooTracker")
            if self.nEntry != self.genieTree.GetEntries():
                print("Edep-sim tree and GENIE tree number of entries do not match!")
                sys.exit()
        elif self.evgen == 'Marley':
            print("Marley events: Skip looking for gRooTracker directory.")
        else:
            print("Unknown event generator!")
            sys.exit()

        self.event = TG4Event()
        self.simTree.SetBranchAddress("Event", self.event)

        self.event = TG4Event()
        self.simTree.SetBranchAddress("Event", self.event)

    # ------------------------
    def Jump(self, entryNo):
        # if entryNo % 100 == 0:
            # print(f'reading event {entryNo}/{self.nEntry}')

        self.currentEntry = entryNo
        self.simTree.GetEntry(entryNo)
        self.ReadVertex()
        self.ReadTracks()
        self.ReadEnergyDepo('SimEnergyDeposit')
        # self.PrintVertex
        # self.PrintTracksEnergy()
        # self.SelectNeutronEvent()
        # self.PrintTracks
        # self.read_neutron_momentum()

        self.info = {}
        self.info['Event_ID'] = self.currentEntry
        self.info['E_nu'] = 0
        self.info['Px_nu'] = 0
        self.info['Py_nu'] = 0
        self.info['Pz_nu'] = 0
        self.info['E_avail'] = 0
        self.info['E_availList'] = np.zeros(8)  # lepton, proton, neutron, pi+-, pi0, gamma, alpha, others.
        self.info['E_depoTotal'] = 0
        self.info['E_depoTotal_track'] = 0
        self.info['Q_depoTotal_dots_th_75keV'] = 0
        self.info['Q_depoTotal'] = 0
        self.info['Q_depoTotal_th_75keV'] = 0
        self.info['Q_depoTotal_th_500keV'] = 0
        self.info['Q_depoTotal_MBox'] = 0
        self.info['Q_depoTotal_MBox_th_75keV'] = 0
        self.info['Q_depoTotal_MBox_th_500keV'] = 0
        self.info['L_depoTotal_avg_220PEpMeV'] = 0
        self.info['L_depoTotal_avg_180PEpMeV'] = 0
        self.info['L_depoTotal_avg_140PEpMeV'] = 0
        self.info['L_depoTotal_avg_100PEpMeV'] = 0
        self.info['L_depoTotal_avg_35PEpMeV'] = 0
        self.info['L_depoTotal_MBox_avg_220PEpMeV'] = 0
        self.info['L_depoTotal_MBox_avg_180PEpMeV'] = 0
        self.info['L_depoTotal_MBox_avg_140PEpMeV'] = 0
        self.info['L_depoTotal_MBox_avg_100PEpMeV'] = 0
        self.info['L_depoTotal_MBox_avg_35PEpMeV'] = 0
        self.info['E_depoList'] = np.zeros(8)  # lepton, proton, neutron, pi+-, pi0, gamma, alpha, others.
        self.info['E_depoList_track'] = np.zeros(8)
        self.info['Q_depoList'] = np.zeros(8)
        gamma_count = self.CountGammasFromFirstVertex()
        self.info['primary_gamma_count'] = gamma_count
        self.info['has_primary_gamma'] = gamma_count > 0
        # self.info['Q_depoList_th_75keV'] = np.zeros(8)
        # self.info['Q_depoList_th_500keV'] = np.zeros(8)
        # self.info['Q_depoList_dots_th_75keV'] = np.zeros(8)
        # self.info['Q_depoList_MBox'] = np.zeros(8)
        # self.info['Q_depoList_MBox_th_75keV'] = np.zeros(8)
        # self.info['Q_depoList_MBox_th_500keV'] = np.zeros(8)
        # self.info['L_depoList_avg_220PEpMeV'] = np.zeros(8)
        # self.info['L_depoList_avg_180PEpMeV'] = np.zeros(8)
        # self.info['L_depoList_avg_140PEpMeV'] = np.zeros(8)
        # self.info['L_depoList_avg_100PEpMeV'] = np.zeros(8)
        # self.info['L_depoList_avg_35PEpMeV'] = np.zeros(8)
        # self.info['L_depoList_MBox_avg_220PEpMeV'] = np.zeros(8)
        # self.info['L_depoList_MBox_avg_180PEpMeV'] = np.zeros(8)
        # self.info['L_depoList_MBox_avg_140PEpMeV'] = np.zeros(8)
        # self.info['L_depoList_MBox_avg_100PEpMeV'] = np.zeros(8)
        # self.info['L_depoList_MBox_avg_35PEpMeV'] = np.zeros(8)
        # self.info['N_parList'] = np.zeros(8)  # lepton, proton, neutron, pi+-, pi0, gamma, alpha, others.
        # self.info['nu_pdg'] = 0
        # self.info['nu_xs'] = self.vertex.GetCrossSection()
        # self.info['nu_proc'], self.info['nu_nucl'] = self.GetReaction()
        # self.FillEnergyInfo()
        if self.evgen == 'Genie':
            self.ReadGenie()
        elif self.evgen == 'Marley':
            self.ReadMarley()
        else:
            print("Unknown event generator!")
            sys.exit()

        self.info['gamma_momentum_fraction'] = self.GammaMomentumFraction()

    # ------------------------
    def ReadGenie(self):
        # gRooTracker info
        # sa
        # StdHepStatus: 0: initial state; 1: final state particles; others: intermediate transport
        # following assumes the first particle is always the neutrino.
        self.genieTree.GetEntry(self.currentEntry)
        self.info['nu_pdg'] = self.genieTree.StdHepPdg[0]
        self.info['E_nu'] = self.genieTree.StdHepP4[3] * 1000
        # Convert neutrino 3-momentum to MeV/c to match Geant4 units.
        self.info['Px_nu'] = self.genieTree.StdHepP4[0] * 1000
        self.info['Py_nu'] = self.genieTree.StdHepP4[1] * 1000
        self.info['Pz_nu'] = self.genieTree.StdHepP4[2] * 1000



        nucleus_p, nucleus_dir = self.get_primary_nucleus_kinematics()
        self.info['target_nucleus_momentum'] = nucleus_p
        self.info['target_nucleus_direction'] = nucleus_dir
        self.info['target_nucleus_zenith_angle_deg'] = self._compute_zenith_angle(nucleus_dir)
    def PrintGenie(self):
        self.genieTree.GetEntry(self.currentEntry)
        for i in range(self.genieTree.StdHepN):
            pdg = self.genieTree.StdHepPdg[i]
            status = self.genieTree.StdHepStatus[i]
            momentum= np.array([self.genieTree.StdHepP4[4*i],self.genieTree.StdHepP4[4*i+1],self.genieTree.StdHepP4[4*i+2]])
            print(f"pdg: {pdg}, momentum: {momentum},status: {status}")
        print("\n")
    def read_final_state_momentum(self):
        self.genieTree.GetEntry(self.currentEntry)
        momentum = [0, 0, 0]
        for i in range(self.genieTree.StdHepN):
            status = self.genieTree.StdHepStatus[i]
            if status == 15:
                momentum += np.array([self.genieTree.StdHepP4[4 * i], self.genieTree.StdHepP4[4 * i + 1],
                                     self.genieTree.StdHepP4[4 * i + 2]])
        return momentum

    def read_final_state_nucleus(self):
        self.genieTree.GetEntry(self.currentEntry)
        momentum = [0,0,0]
        for i in range(self.genieTree.StdHepN):
            status = self.genieTree.StdHepStatus[i]
            if status == 15:
                momentum = np.array([self.genieTree.StdHepP4[4 * i], self.genieTree.StdHepP4[4 * i + 1],
                                 self.genieTree.StdHepP4[4 * i + 2]])

        return momentum



    def read_electron_direction(self):
            """ here we just read the muon direction AND AMPLITUDE from
            intrisic function """
            direction_vector = []
            i = 0
            track = self.tracks[i]
            parId = track.GetParentId()
            pdg = track.GetPDGCode()
            electron_mom_amplitude = 0
            while parId != -1 or pdg not in (-11, 11):
                    i += 1
                    track = self.tracks[i]
                    parId = track.GetParentId()
                    pdg = track.GetPDGCode()

            px = track.GetInitialMomentum().X()
            py = track.GetInitialMomentum().Y()
            pz = track.GetInitialMomentum().Z()
            p_square = (px ** 2 + py ** 2 + pz ** 2) ** 0.5
            if p_square == 0.0:
                    direction_vector = np.array([0.0, 0.0, 0.0])
            else:
                    direction_vector.append(px / p_square)
                    direction_vector.append(py / p_square)
                    direction_vector.append(pz / p_square)
                    electron_mom_amplitude = p_square

            print(f"direction_vector: {direction_vector}","\n",f"electron_mom_amplitude: {electron_mom_amplitude}")
            return direction_vector, electron_mom_amplitude
    def pairing_electron_energy_with_electron_cone(self,electron_energy):
        angle = self.cone_applying_for_electron[round(electron_energy / 100)]
        return angle



    def sum_neutron_energy_deposition_EV(self,energy_threshold: Optional[float] = None):

        electron_energy = self.read_electron_energy()
        cone = self.pairing_electron_energy_with_electron_cone(electron_energy)
        electron_dir,_ = self.read_electron_direction()
        track_proton_pion = []
        threshold = self._resolve_energy_threshold(energy_threshold)
        mm2cm = 0.1
        energy_sum = 0
        position = []
        for i, particle in enumerate(self.vertex.Particles):
            pdg = particle.GetPDGCode()
            if pdg ==2212:
                mom=particle.GetMomenum()
                KE= mom.E()-mom.M()
                if KE> self.proton_trackable_length:
                    track_proton_pion.append(particle.GetTrackId())
                elif pdg in (111, 211, -211):
                    track_proton_pion.append(particle.GetTrackId())
            for track in self.tracks:
                trackid = track.GetTrackId()
                pdg_code= track.GetPDGCode()
                if pdg_code not in[22,111,211,-211,2112,2212,11,-11,13,-13]:
                    continue
                if trackid in track_proton_pion:
                    continue
                parent_id = track.GetParentId()
                if parent_id != -1 and self.tracks[parent_id].GetTrackId() in track_proton_pion:
                    continue
                for di in track.association.get('depoList', []):
                        if not (0 <= di < len(self.depos)):
                            continue
                        depo = self.depos[di]
                        edep = depo.GetEnergyDeposit()
                        # print(f"edep: {edep}")
                        if edep < threshold:
                            continue
                        start, stop = depo.GetStart(), depo.GetStop()
                        midpoint = np.array([
                            (start.X() + stop.X()) / 2 * mm2cm,
                            (start.Y() + stop.Y()) / 2 * mm2cm,
                            (start.Z() + stop.Z()) / 2 * mm2cm
                        ])
                        # print(f"midpoint: {midpoinidpoint}")
                        midpointlength = np.linalg.norm(midpoint)
                        angle = self.compute_angle_between(electron_dir, midpoint)
                        angle = np.degrees(angle)
                        # Skip depositions in the muon cone or farther than 5 m from the vertex
                        # print(f"angle: {angle},midpointlength: {midpointlength}")
                        if angle < cone or midpointlength > self.inspection_length_energy_side*100 or midpointlength < self.minimum_length:
                            continue
                        energy_sum += edep
                        position.append(midpoint)
            return energy_sum, position

    def sum_neutron_energy_deposition(self):
        """The next step is to redo the neutron momentum magnitude reconstruction in a more realistic way,
        e,g, first exclude all the energy deposition above 5 meters, then exclude all the energy deposition inside the muon cone,
        then decide for proton and pion if their energy deposition forms a track, if so, remove the track."""

        ## this code is working for one event. do not forget to jump to this event, every time call this function
        muon_energy = self.read_muon_energy()

        muon_dir,mom = self.read_muon_direction_and_amplitude()
        cone = self.pairing_muon_energy_with_muon_cone(muon_energy)
        track_proton_pion = []
        threshold = self.energy_threshold
        mm2cm = 0.1
        energy_sum = 0
        position = []
        for i, particle in enumerate(self.vertex.Particles):
            pdg = particle.GetPDGCode()
            if pdg ==2212:
                mom = particle.GetMomentum()
                KE = mom.E() - mom.M()
                if KE > self.proton_trackable_length:
                    track_proton_pion.append(particle.GetTrackId())
            elif pdg in (111,211,-211):
                track_proton_pion.append(particle.GetTrackId())
        for track in self.tracks:
            trackid = track.GetTrackId()
            pdg_code = track.GetPDGCode()
            if pdg_code not in [22,2112,2212,211,-211,111,11,-11,13,-13]:
                continue
            if trackid in track_proton_pion:
                continue
            parent_id = track.GetParentId()
            if parent_id != -1 and self.tracks[parent_id].GetTrackId() in track_proton_pion:
                continue
            for di in track.association.get('depoList', []):
                    if not (0 <= di < len(self.depos)):
                        continue
                    depo = self.depos[di]
                    edep = depo.GetEnergyDeposit()
                    # print(f"edep: {edep}")
                    if edep < threshold:
                        continue
                    start, stop = depo.GetStart(), depo.GetStop()
                    midpoint = np.array([
                        (start.X() + stop.X()) / 2 * mm2cm,
                        (start.Y() + stop.Y()) / 2 * mm2cm,
                        (start.Z() + stop.Z()) / 2 * mm2cm
                    ])
                    # print(f"midpoint: {midpoinidpoint}")
                    midpointlength = np.linalg.norm(midpoint)
                    angle = self.compute_angle_between(muon_dir, midpoint)
                    angle = np.degrees(angle)
                    # Skip depositions in the muon cone or farther than 5 m from the vertex
                    # print(f"angle: {angle},midpointlength: {midpointlength}")
                    if angle < cone or midpointlength > self.inspection_length_energy_side*100 or midpointlength < self.minimum_length:
                        continue
                    energy_sum += edep
                    position.append(midpoint)
        return energy_sum, position






    def get_primary_electron_track(self):
        """Return the primary electron ancestor track, its index, and descendants."""

        if not hasattr(self, 'tracks') or getattr(self.tracks, 'size', 0) == 0:
            return None, None, []

        for idx in range(self.tracks.size):
            track = self.tracks[idx]
            if not hasattr(track, 'association'):
                continue
            ancestor_idx = track.association.get('ancestor')
            if ancestor_idx is None or ancestor_idx < 0 or ancestor_idx >= self.tracks.size:
                continue

            ancestor_track = self.tracks[ancestor_idx]
            if ancestor_track.GetPDGCode() in (-11, 11) and ancestor_track.GetParentId() ==-1:
                descendants = [i for i in range(self.tracks.size)
                               if hasattr(self.tracks[i], 'association')
                               and self.tracks[i].association.get('ancestor') == ancestor_idx]
                return ancestor_track, ancestor_idx, descendants

        return None, None, []


    def print_electron_track_processes(self):
        """Print the Geant4 process/subprocess tags for the primary electron track."""

        electron_track, electron_index, _ = self.get_primary_electron_track()
        if electron_track is None:
            print('print_electron_track_processes: primary electron track not found.')
            return

        print(
            f"Electron track #{electron_index} (pdg={electron_track.GetPDGCode()}):"
        )

        for step, point in enumerate(electron_track.Points, start=1):
            process = point.GetProcess() or 'UNKNOWN'
            subprocess = point.GetSubprocess() or 'UNKNOWN'
            pos = point.GetPosition()
            print(
                f"  step {step:03d}: process={process}, subprocess={subprocess}, "
                f"x={pos.X():.2f} mm, y={pos.Y():.2f} mm, z={pos.Z():.2f} mm, t={pos.T():.2f} ns"
            )



    def get_primary_nucleus_kinematics(self):
        """Return initial-state nucleus momentum (MeV/c) and direction unit vector.

        GENIE encodes the target nucleus at index 1 of the flattened StdHepP4
        array. When the kinematics are unavailable (non-GENIE sample or missing
        entry) the method returns zero vectors."""

        if self.evgen != 'Genie' or not hasattr(self, 'genieTree'):
            return (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)

        # Ensure the current GENIE entry is loaded so the StdHep arrays are valid.
        self.genieTree.GetEntry(self.currentEntry)

        components = self._extract_stdhep_four_vector(3)
        if components is None:
            return (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)

        px, py, pz, _ = components
        momentum = np.array([px, py, pz])
        p_mag = math.sqrt(px ** 2 + py ** 2 + pz ** 2)
        if p_mag == 0.0:
            direction = (0.0, 0.0, 0.0)
        else:
            inv_mag = 1.0 / p_mag
            direction = (px * inv_mag, py * inv_mag, pz * inv_mag)

        return momentum, direction

    @staticmethod
    def _compute_zenith_angle(direction_vector):
        """Return the angle (degrees) between a unit vector and the +z axis."""

        if direction_vector is None:
            return None

        dx, dy, dz = direction_vector
        norm = math.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
        if norm == 0.0:
            return None

        cos_theta = max(min(dz / norm, 1.0), -1.0)
        return math.degrees(math.acos(cos_theta))

    def _extract_stdhep_four_vector(self, particle_index):
        """Return (px, py, pz, E) in MeV for a StdHep particle or None if unavailable."""

        if self.evgen != 'Genie' or not hasattr(self, 'genieTree'):
            return None

        p4 = getattr(self.genieTree, 'StdHepP4', None)
        if p4 is None:
            return None

        try:
            base = particle_index * 4
            px = float(p4[base]) * 1000.0
            py = float(p4[base + 1]) * 1000.0
            pz = float(p4[base + 2]) * 1000.0
            energy = float(p4[base + 3]) * 1000.0
            return px, py, pz, energy
        except (TypeError, IndexError):
            pass

        try:
            lorentz = p4[particle_index]
        except Exception:
            lorentz = None

        if lorentz is None and hasattr(p4, 'At'):
            try:
                lorentz = p4.At(particle_index)
            except Exception:
                lorentz = None

        if lorentz is None:
            return None

        for getter in ('Px', 'Py', 'Pz', 'E'):
            if not hasattr(lorentz, getter):
                return None

        return (
            float(lorentz.Px()) * 1000.0,
            float(lorentz.Py()) * 1000.0,
            float(lorentz.Pz()) * 1000.0,
            float(lorentz.E()) * 1000.0,
        )

    def PrintGenieTree(self, event_index: int, include_kinematics: bool = True):
        """Print the GENIE StdHep record for the provided event index."""

        if self.evgen != 'Genie':
            print(f"GENIE tree printing is unavailable for '{self.evgen}' samples.")
            return

        if not hasattr(self, 'genieTree'):
            print('GENIE gRooTracker tree has not been loaded yet.')
            return

        total_entries = self.genieTree.GetEntries()
        if total_entries <= 0:
            print('GENIE tree is empty.')
            return

        if event_index < 0 or event_index >= total_entries:
            print(f'Event index {event_index} is outside the available range [0, {total_entries - 1}].')
            return

        previous_entry = getattr(self, 'currentEntry', None)
        self.genieTree.GetEntry(event_index)

        particle_count = int(getattr(self.genieTree, 'StdHepN', 0))
        if particle_count <= 0:
            print(f'GENIE event {event_index} does not contain StdHep particles.')
            if previous_entry is not None and 0 <= previous_entry < total_entries and previous_entry != event_index:
                self.genieTree.GetEntry(previous_entry)
            return

        def branch_values(branch, cast=None):
            if branch is None:
                return []

            try:
                iterable = list(branch)
            except TypeError:
                try:
                    size = branch.size() if hasattr(branch, 'size') else particle_count
                    iterable = [branch[i] for i in range(int(size))]
                except Exception:
                    return []

            values = []
            for item in iterable[:particle_count]:
                if cast is None:
                    values.append(item)
                else:
                    try:
                        values.append(cast(item))
                    except Exception:
                        values.append(item)
            return values

        pdgs = branch_values(getattr(self.genieTree, 'StdHepPdg', None), int)
        statuses = branch_values(getattr(self.genieTree, 'StdHepStatus', None), int)
        first_mothers = branch_values(getattr(self.genieTree, 'StdHepFm', None), int)
        last_mothers = branch_values(getattr(self.genieTree, 'StdHepLm', None), int)
        first_daughters = branch_values(getattr(self.genieTree, 'StdHepFd', None), int)
        last_daughters = branch_values(getattr(self.genieTree, 'StdHepLd', None), int)

        pdg_database = ROOT.TDatabasePDG.Instance() if hasattr(ROOT, 'TDatabasePDG') else None
        status_names = {0: 'init', 1: 'final', 2: 'inter', 3: 'doca', 11: 'decay', 12: 'doca2'}

        def format_range(first, last):
            try:
                first = int(first)
                last = int(last)
            except (TypeError, ValueError):
                return '-'

            if first <= 0 and last <= 0:
                return '-'

            if last < first:
                first, last = last, first

            if first == last:
                return str(first)

            return f'{first}-{last}'

        header = f"{'idx':>3} {'status':>8} {'pdg':>8} {'name':<18} {'mothers':>10} {'daughters':>11}"
        if include_kinematics:
            header += f" {'E [MeV]':>10} {'px':>10} {'py':>10} {'pz':>10}"

        print(f'GENIE gRooTracker for event {event_index} ({particle_count} particles)')
        print(header)
        print('-' * len(header))

        for idx in range(particle_count):
            pdg_code = pdgs[idx] if idx < len(pdgs) else None
            status_code = statuses[idx] if idx < len(statuses) else None
            mother_display = format_range(
                first_mothers[idx] if idx < len(first_mothers) else None,
                last_mothers[idx] if idx < len(last_mothers) else None,
            )
            daughter_display = format_range(
                first_daughters[idx] if idx < len(first_daughters) else None,
                last_daughters[idx] if idx < len(last_daughters) else None,
            )

            if pdg_database is not None and pdg_code is not None:
                pdg_entry = pdg_database.GetParticle(int(pdg_code))
                particle_name = pdg_entry.GetName() if pdg_entry and pdg_entry.GetName() else ''
            else:
                particle_name = ''

            if not particle_name:
                particle_name = 'unknown'

            status_label = status_names.get(status_code)
            if status_label:
                status_display = f'{status_code}:{status_label}'
            elif status_code is None:
                status_display = '?'
            else:
                status_display = str(status_code)

            pdg_display = str(pdg_code) if pdg_code is not None else '?'

            row = f"{idx:>3} {status_display:>8} {pdg_display:>8} {particle_name:<18} {mother_display:>10} {daughter_display:>11}"

            if include_kinematics:
                components = self._extract_stdhep_four_vector(idx)
                if components is not None:
                    px, py, pz, energy = components
                    row += f" {energy:>10.1f} {px:>10.1f} {py:>10.1f} {pz:>10.1f}"
                else:
                    row += f" {'-':>10} {'-':>10} {'-':>10} {'-':>10}"

            print(row)

        if previous_entry is not None and 0 <= previous_entry < total_entries and previous_entry != event_index:
            self.genieTree.GetEntry(previous_entry)

    # ------------------------
    def ReadMarley(self):
        # we are mostly looking at nue anyway
        self.info['nu_pdg'] = 12
        if self.GetnuPDGFromFileName() == 'nue': self.info['nu_pdg'] = 12
        if self.GetnuPDGFromFileName() == 'numu': self.info['nu_pdg'] = 14
        if self.GetnuPDGFromFileName() == 'anue': self.info['nu_pdg'] = -12
        if self.GetnuPDGFromFileName() == 'anumu': self.info['nu_pdg'] = -14
        self.info['E_nu'] = self.GetEnuFromFileName()
        if self.currentEntry == 0:  # only print once
            print("Marley events: assert info from file name")

    # ------------------------
    def ReadVertex(self):
        primaries = np.array(self.event.Primaries)
        # print(f"Primaries: {primaries}")
        self.vertex = None
        if (primaries.size != 1 and self.evgen == 'Genie'):
            print("Number of primaries not equal to 1 (not neutrino vertex)!")
            return

        self.vertex = primaries[0]

    def CountGammasFromFirstVertex(self):
        """Return the number of gamma (pdg 22) particles at the primary vertex."""
        if not hasattr(self, 'vertex') or self.vertex is None:
            return 0

        gamma_count = 0
        for particle in self.vertex.Particles:
            if particle.GetTrackId() < 0:
                continue
            if particle.GetPDGCode() == 22:
                gamma_count += 1
        return gamma_count

    def HasGammaFromFirstVertex(self):
        """Return True if any gamma originates from the primary vertex."""
        return self.CountGammasFromFirstVertex() > 0

    def PrimaryFinalStatePDGs(self):
        """Return PDG codes for tracked particles emitted at the primary vertex."""
        if not hasattr(self, 'vertex') or self.vertex is None:
            return []

        pdgs = []
        for particle in self.vertex.Particles:
            if particle.GetTrackId() < 0:
                continue
            pdgs.append(particle.GetPDGCode())
        return pdgs

    def print_other_vertex_particles(
        self, exclude_pdgs: Tuple[int, ...] = (2212, 211, -211,111,13,-13,2112,22)
    ) -> Tuple[str, ...]:
        """Print the other tracked primaries at the vertex after excluding certain PDGs.

        Parameters
        ----------
        exclude_pdgs: tuple of int, optional
            PDG codes to skip when reporting. Defaults to removing protons and charged pions.

        Returns
        -------
        tuple of str
            Human-readable summary lines describing the remaining particles.
        """
        if not hasattr(self, 'vertex') or self.vertex is None:
            print('No primary vertex available in the current event.')
            return tuple()

        pdg_exclusions = set(exclude_pdgs or [])
        pdg_db = ROOT.TDatabasePDG.Instance() if hasattr(ROOT, 'TDatabasePDG') else None
        counts = Counter()
        kinetic_energy_map = {}
        baryon_with_ke = {3122, 3112, 3212}

        for particle in self.vertex.Particles:
            if particle.GetTrackId() < 0:
                continue

            pdg = particle.GetPDGCode()
            if pdg in pdg_exclusions:
                continue

            name = particle.GetName() if hasattr(particle, 'GetName') else ''
            if not name and pdg_db is not None:
                pdg_entry = pdg_db.GetParticle(pdg)
                if pdg_entry is not None:
                    fetched = pdg_entry.GetName()
                    if fetched:
                        name = fetched
            if not name:
                name = f'PDG {pdg}'
            key = (pdg, name)
            counts[key] += 1

            if pdg in baryon_with_ke and hasattr(particle, 'GetMomentum'):
                momentum = particle.GetMomentum()
                kinetic_energy = None
                if momentum is not None and hasattr(momentum, 'E') and hasattr(momentum, 'M'):
                    try:
                        kinetic_energy = float(momentum.E() - momentum.M())
                    except TypeError:
                        kinetic_energy = None

                if kinetic_energy is not None:
                    kinetic_energy_map.setdefault(key, []).append(kinetic_energy)

        if not counts:
            print('No additional tracked particles at the primary vertex after exclusions.')
            return tuple()

        print('Primary-vertex particles (excluding requested PDG codes):')
        lines = []
        for (pdg, name), amount in sorted(counts.items(), key=lambda item: item[0][0]):
            line = f'  {name} (PDG {pdg}): {amount}'
            print(line)
            lines.append(line)

            energies = kinetic_energy_map.get((pdg, name))
            if energies:
                formatted = ', '.join(f'{ke:.3f} GeV' for ke in energies)
                energy_line = f'    kinetic energy at vertex: {formatted}'
                print(energy_line)
                lines.append(energy_line)
        return tuple(lines)
    def GammaMomentumFraction(self):
        """Fraction of the incoming neutrino momentum carried by primary gammas.

        Gammas produced at the primary vertex are typically emitted in pairs from
        a neutral pion decay. Summing their three-momenta vectorially can cancel
        large, back-to-back contributions, driving the fraction to zero even
        though energetic gammas are present. Measuring the total scalar momentum
        (equivalent to energy for massless gammas) avoids that cancellation and
        better reflects how much of the neutrino's momentum leaves in gammas."""

        if not hasattr(self, 'vertex') or self.vertex is None:
            return 0.0

        neutrino_px = self.info.get('Px_nu', 0.0)
        neutrino_py = self.info.get('Py_nu', 0.0)
        neutrino_pz = self.info.get('Pz_nu', 0.0)
        neutrino_p_mag = math.sqrt(neutrino_px ** 2 + neutrino_py ** 2 + neutrino_pz ** 2)
        if neutrino_p_mag <= 0.0:
            return 0.0

        gamma_p_total = 0.0
        for particle in self.vertex.Particles:
            if particle.GetTrackId() < 0 or particle.GetPDGCode() != 22:
                continue

            momentum = particle.GetMomentum()
            gamma_p_total += momentum.P()

        if gamma_p_total <= 0.0:
            return 0.0

        return gamma_p_total / neutrino_p_mag

    def read_neutrino_direction(self):
        for particle in self.event.Particles:
            pdg = particle.GetPDG()
            if pdg in (12,14,16):
                print("see the neutrino")
    # --------------------------
    def GetReaction(self):
        txt_list = self.vertex.GetReaction().split(';')
        proc = ''
        nucl = 0
        for x in txt_list[2:]:
            if 'proc:' in x:
                proc = x.replace('proc:', '').replace('Weak[', '').replace('],', '')
            elif 'N:' in x:
                nucl = int(x.replace('N:', ''))
        # print(proc, nucl)
        proc_num = 0
        if (proc.startswith('CC')):
            proc_num = 10
            proc = proc.replace('CC', '')
        elif (proc.startswith('NC')):
            proc_num = 20
            proc = proc.replace('NC', '')
        else:
            proc_num = 30
        proc_dict = {
            'QES': 1,
            'RES': 2,
            'DIS': 3,
            'COH': 4,
            'MEC': 5,
        }
        proc_num += proc_dict.get(proc, 0)
        return proc_num, nucl

    # ------------------------
    def ReadTracks(self):
        self.tracks = np.array(self.event.Trajectories)
        for i in range(self.tracks.size):
            self.tracks[i].energy = {}
            self.tracks[i].association = {}
            self.tracks[i].length = {}
            self.tracks[i].energy['depoTotal'] = 0
            self.tracks[i].energy['depoTotal_charge'] = 0
            self.tracks[i].energy['depoTotal_charge_MBox'] = 0
            self.tracks[i].energy['depoTotal_charge_th_75keV'] = 0
            self.tracks[i].energy['depoTotal_charge_MBox_th_75keV'] = 0
            self.tracks[i].energy['depoTotal_charge_th_500keV'] = 0
            self.tracks[i].energy['depoTotal_charge_MBox_th_500keV'] = 0
            self.tracks[i].energy['depoTotal_light_avg_220PEpMeV'] = 0
            self.tracks[i].energy['depoTotal_light_avg_220PEpMeV_MBox'] = 0
            self.tracks[i].energy['depoTotal_light_avg_180PEpMeV'] = 0
            self.tracks[i].energy['depoTotal_light_avg_180PEpMeV_MBox'] = 0
            self.tracks[i].energy['depoTotal_light_avg_140PEpMeV'] = 0
            self.tracks[i].energy['depoTotal_light_avg_140PEpMeV_MBox'] = 0
            self.tracks[i].energy['depoTotal_light_avg_100PEpMeV'] = 0
            self.tracks[i].energy['depoTotal_light_avg_100PEpMeV_MBox'] = 0
            self.tracks[i].energy['depoTotal_light_avg_35PEpMeV'] = 0
            self.tracks[i].energy['depoTotal_light_avg_35PEpMeV_MBox'] = 0
            self.tracks[i].length['selfDepo'] = 0
            self.tracks[i].association['depoList'] = []
            self.tracks[i].association['children'] = []
            self.tracks[i].association['ancestor'] = i

        for i in range(self.tracks.size):
            track = self.tracks[i]
            parId = track.GetParentId()
            if parId == -1: continue
            self.tracks[parId].association['children'].append(i)
            while parId != -1:
                track.association['ancestor'] = parId
                parId = self.tracks[parId].GetParentId()

    # ------------------------
    def ReadEnergyDepo(self, detName):
        self.depos = np.array(self.event.SegmentDetectors[detName])

        # add depo info to tracks
        # depoList = self.FindDepoListFromTrack(1)
        # depoEnergy = np.sum([depo.GetEnergyDeposit() for depo in self.depos[depoList]])
        # print('debug: muon deposit energy:', depoEnergy)
        mm2cm = 0.1

        for i, depo in enumerate(self.depos):
            trkId = depo.Contrib[0]
            edep = depo.GetEnergyDeposit()
            trkLength = depo.GetTrackLength() * mm2cm
            Qdep = self.ChargeBirksLaw(edep, trkLength)
            Qdep_MBox = self.ChargeModifiedBoxModel(edep, trkLength)  # Modified Box Model
            track = self.tracks[trkId]
            track.association['depoList'].append(i)
            track.energy['depoTotal'] += edep
            track.energy['depoTotal_charge'] += Qdep
            track.energy['depoTotal_charge_MBox'] += Qdep_MBox
            track.length['selfDepo'] += trkLength
            if Qdep > 0.075: track.energy['depoTotal_charge_th_75keV'] += Qdep
            if Qdep > 0.5: track.energy['depoTotal_charge_th_500keV'] += Qdep
            if Qdep_MBox > 0.075: track.energy['depoTotal_charge_MBox_th_75keV'] += Qdep_MBox
            if Qdep_MBox > 0.5: track.energy['depoTotal_charge_MBox_th_500keV'] += Qdep_MBox

            # Light deposit
            # edepSecond = depo.GetSecondaryDeposit()
            # Ideally this should be the light dep
            # how the light yield is actually derived: 25k photons / MeV are simualted and 220 on average are detected
            # the overall photon collection efficiency (PCE) is
            #    220/25000=0.88% PCE
            #    180/25000=0.72% PCE
            #    140/25000=0.56% PCE
            #    100/25000=0.4% PCE
            #    35/25000=0.14% PCE

            #    220/21622=1.0% PCE
            #    180/21622=0.83% PCE
            #    140/21622=0.65% PCE
            #    100/21622=0.46% PCE
            #    35/21622=0.16% PCE
            # (dL/W_ph)*PCE = (dL/19.5eV)*PCE = number of photons (taking into account Birks model)
            Ldep = edep - Qdep
            Ldep_MBox = edep - Qdep_MBox
            nPE_220 = (Ldep * 1000000 / 19.5) * 0.01
            nPE_180 = (Ldep * 1000000 / 19.5) * 0.0083
            nPE_140 = (Ldep * 1000000 / 19.5) * 0.0065
            nPE_100 = (Ldep * 1000000 / 19.5) * 0.0046
            nPE_35 = (Ldep * 1000000 / 19.5) * 0.0016
            nPE_220_MBox = Ldep_MBox * 220
            nPE_180_MBox = Ldep_MBox * 180
            nPE_140_MBox = Ldep_MBox * 140
            nPE_100_MBox = Ldep_MBox * 100
            nPE_35_MBox = Ldep_MBox * 35
            nPE_220_detected = random.gauss(nPE_220, math.sqrt(nPE_220))
            nPE_180_detected = random.gauss(nPE_180, math.sqrt(nPE_180))
            nPE_140_detected = random.gauss(nPE_140, math.sqrt(nPE_140))
            nPE_100_detected = random.gauss(nPE_100, math.sqrt(nPE_100))
            nPE_35_detected = random.gauss(nPE_35, math.sqrt(nPE_35))
            nPE_220_MBox_detected = random.gauss(nPE_220_MBox, math.sqrt(nPE_220_MBox))
            nPE_180_MBox_detected = random.gauss(nPE_180_MBox, math.sqrt(nPE_180_MBox))
            nPE_140_MBox_detected = random.gauss(nPE_140_MBox, math.sqrt(nPE_140_MBox))
            nPE_100_MBox_detected = random.gauss(nPE_100_MBox, math.sqrt(nPE_100_MBox))
            nPE_35_MBox_detected = random.gauss(nPE_35_MBox, math.sqrt(nPE_35_MBox))
            Ldep_avg_220PEpMeV_detected = nPE_220_detected * 19.5 / 1000000 / 0.01  # MeV
            Ldep_avg_180PEpMeV_detected = nPE_180_detected * 19.5 / 1000000 / 0.0083
            Ldep_avg_140PEpMeV_detected = nPE_140_detected * 19.5 / 1000000 / 0.0065
            Ldep_avg_100PEpMeV_detected = nPE_100_detected * 19.5 / 1000000 / 0.0046
            Ldep_avg_35PEpMeV_detected = nPE_35_detected * 19.5 / 1000000 / 0.0016
            Ldep_avg_220PEpMeV_MBox_detected = nPE_220_MBox_detected / 220
            Ldep_avg_180PEpMeV_MBox_detected = nPE_180_MBox_detected / 180
            Ldep_avg_140PEpMeV_MBox_detected = nPE_140_MBox_detected / 140
            Ldep_avg_100PEpMeV_MBox_detected = nPE_100_MBox_detected / 100
            Ldep_avg_35PEpMeV_MBox_detected = nPE_35_MBox_detected / 35

            track.energy['depoTotal_light_avg_220PEpMeV'] += Ldep_avg_220PEpMeV_detected
            track.energy['depoTotal_light_avg_180PEpMeV'] += Ldep_avg_180PEpMeV_detected
            track.energy['depoTotal_light_avg_140PEpMeV'] += Ldep_avg_140PEpMeV_detected
            track.energy['depoTotal_light_avg_100PEpMeV'] += Ldep_avg_100PEpMeV_detected
            track.energy['depoTotal_light_avg_35PEpMeV'] += Ldep_avg_35PEpMeV_detected

            track.energy['depoTotal_light_avg_220PEpMeV_MBox'] += Ldep_avg_220PEpMeV_MBox_detected
            track.energy['depoTotal_light_avg_180PEpMeV_MBox'] += Ldep_avg_180PEpMeV_MBox_detected
            track.energy['depoTotal_light_avg_140PEpMeV_MBox'] += Ldep_avg_140PEpMeV_MBox_detected
            track.energy['depoTotal_light_avg_100PEpMeV_MBox'] += Ldep_avg_100PEpMeV_MBox_detected
            track.energy['depoTotal_light_avg_35PEpMeV_MBox'] += Ldep_avg_35PEpMeV_MBox_detected

        # E_tot = np.sum([depo.GetEnergyDeposit() for depo in self.depos])
        # print('total deposit energy: ', E_tot)

    # ------------------------
    def ChargeBirksLaw(self, edep, trkLength):
        # PHYS. REV. D 99, 036009 (2019)
        # return R = (dQ/dx)/(dE/dx) = dQ/dE = A/(1+kQ*dE/dx)
        # A = 0.8, kQ = 0.0972 g/(MeV*cm2) @500V/cm
        # LAr density 1.4g/cm^3
        # when considering light, need to take account excitation and ionization ratio (0.83 vs 0.17)
        return 0.83 * edep * 0.8 / (1 + 0.0972 * edep / trkLength / 1.4)

    # ------------------------
    def ChargeModifiedBoxModel(self, edep, trkLength):
        # https://lar.bnl.gov/properties/pass.html
        # return R = (dQ/dx)/(dE/dx) = dQ/dE = ln(A+B*dE/dx)/(B*dE/dx)
        # A = 0.930, B = 0.424 g/(MeV*cm2) @500V/cm
        # LAr density 1.4g/cm^3
        # when considering light, need to take account excitation and ionization ratio (0.83 vs 0.17)
        return 0.83 * edep * np.log(0.930 + 0.424 * edep / trkLength / 1.4) / (0.424 * edep / trkLength / 1.4)

    # ------------------------
    def FindDepoListFromTrack(self, trkId):
        x = [i for (i, depo)
             in enumerate(self.depos)
             if (trkId in depo.Contrib)]
        # if (trkId == depo.GetPrimaryId())]
        return x

    # ------------------------
    def PrintVertex(self):
        if self.evgen == 'Genie':
            self.info['nu_pdg'] = self.genieTree.StdHepPdg[0]
            self.info['E_nu'] = self.genieTree.StdHepP4[3] * 1000
            print(f"neutrino {self.info['nu_pdg']}: {self.info['E_nu']} MeV")
        elif self.evgen == 'Marley':
            self.ReadMarley()
        else:
            print("Unknown event generator!")
            sys.exit()

        posx = self.vertex.GetPosition().X()
        posy = self.vertex.GetPosition().Y()
        posz = self.vertex.GetPosition().Z()
        print(f"vertex @: ({posx}, {posy}, {posz}) [mm]")
        print(f"reaction: {self.vertex.GetReaction()}")
        print(f"interaction #: {self.vertex.GetInteractionNumber()}")

        print(f"{self.vertex.Particles.size()} particles at the vertex", )
        print(f'{"pdg":>8}{"name":>8}{"trkId":>6}{"mass":>10}{"KE":>10}')
        print(f'{"":>8}{"":>8}{"":>6}{"[MeV]":>10}{"[MeV]":>10}')
        print('-'*(8+8+6+10+10))
        for particle in self.vertex.Particles:
            trkId = particle.GetTrackId()
            # Skip negative trk id: in the case of Marley events,
            # this usually is the final nucleus before deexcitation that G4 doesn't track
            # the kinematics are not correct either
            if trkId < 0: continue
            pdg = particle.GetPDGCode()
            name = particle.GetName()
            # mom_direction=particle.GetMomentum.direction()
            # print("this is momentum direction:",mom_direction)
            mom = particle.GetMomentum()
            # print("this is momentum direction:",mom)
            mass = mom.M()
            KE = mom.E() - mass
        print(f'{pdg:>8d}{name:>8s}{trkId:>6d}{mass:>10.2f}{KE:>10.2f}')
        print('-' * (8 + 8 + 6 + 10 + 10))

        print(f'{self.info}')

    # ------------------------
    def FillEnergyInfo(self):

        for particle in self.vertex.Particles:
            trkId = particle.GetTrackId()
            # Skip negative trk id: in the case of Marley events,
            # this usually is the final nucleus before deexcitation that G4 doesn't track
            # the kinematics are not correct either
            if trkId < 0: continue
            pdg = particle.GetPDGCode()
            depoE = self.GetEnergyDepoWithDesendents(trkId)
            track = self.tracks[trkId]
            # if longer than 2cm, assume can reconstruct dE
            if track.length['selfDepo'] > 2:
                depoE_track = track.energy['depoTotal']
                depoQ_dots = 0
            else:
                # dots/blips
                depoE_track = 0
                depoQ_dots = track.energy['depoTotal_charge_th_75keV']
            depoQ = self.GetChargeDepoWithDesendents(trkId)[0]
            depoQ_75keV = self.GetChargeDepoWithDesendents(trkId)[1]
            depoQ_500keV = self.GetChargeDepoWithDesendents(trkId)[2]
            depoQ_MBox = self.GetChargeMBoxDepoWithDesendents(trkId)[0]
            depoQ_MBox_75keV = self.GetChargeMBoxDepoWithDesendents(trkId)[1]
            depoQ_MBox_500keV = self.GetChargeMBoxDepoWithDesendents(trkId)[2]
            depoL_220PEpMeV = self.GetLightDepoWithDesendents(trkId)[0]
            depoL_180PEpMeV = self.GetLightDepoWithDesendents(trkId)[1]
            depoL_140PEpMeV = self.GetLightDepoWithDesendents(trkId)[2]
            depoL_100PEpMeV = self.GetLightDepoWithDesendents(trkId)[3]
            depoL_35PEpMeV = self.GetLightDepoWithDesendents(trkId)[4]
            depoL_MBox_220PEpMeV = self.GetLightMBoxDepoWithDesendents(trkId)[0]
            depoL_MBox_180PEpMeV = self.GetLightMBoxDepoWithDesendents(trkId)[1]
            depoL_MBox_140PEpMeV = self.GetLightMBoxDepoWithDesendents(trkId)[2]
            depoL_MBox_100PEpMeV = self.GetLightMBoxDepoWithDesendents(trkId)[3]
            depoL_MBox_35PEpMeV = self.GetLightMBoxDepoWithDesendents(trkId)[4]
            mom = particle.GetMomentum()
            mass = mom.M()
            KE = mom.E() - mass
            self.info['E_depoTotal'] += depoE
            self.info['E_depoTotal_track'] += depoE_track
            self.info['Q_depoTotal'] += depoQ
            self.info['Q_depoTotal_th_75keV'] += depoQ_75keV
            self.info['Q_depoTotal_th_500keV'] += depoQ_500keV
            self.info['Q_depoTotal_dots_th_75keV'] += depoQ_dots
            self.info['Q_depoTotal_MBox'] += depoQ_MBox
            self.info['Q_depoTotal_MBox_th_75keV'] += depoQ_MBox_75keV
            self.info['Q_depoTotal_MBox_th_500keV'] += depoQ_MBox_500keV
            self.info['L_depoTotal_avg_220PEpMeV'] += depoL_220PEpMeV
            self.info['L_depoTotal_avg_180PEpMeV'] += depoL_180PEpMeV
            self.info['L_depoTotal_avg_140PEpMeV'] += depoL_140PEpMeV
            self.info['L_depoTotal_avg_100PEpMeV'] += depoL_100PEpMeV
            self.info['L_depoTotal_avg_35PEpMeV'] += depoL_35PEpMeV
            self.info['L_depoTotal_MBox_avg_220PEpMeV'] += depoL_MBox_220PEpMeV
            self.info['L_depoTotal_MBox_avg_180PEpMeV'] += depoL_MBox_180PEpMeV
            self.info['L_depoTotal_MBox_avg_140PEpMeV'] += depoL_MBox_140PEpMeV
            self.info['L_depoTotal_MBox_avg_100PEpMeV'] += depoL_MBox_100PEpMeV
            self.info['L_depoTotal_MBox_avg_35PEpMeV'] += depoL_MBox_35PEpMeV
            # fill E_availList: lepton, proton, neutron, pi+-, pi0, gamma, alpha, others.
            if (pdg in [13, -13, 11, -11]):
                self.info['E_avail'] += (KE + mass)
                self.info['E_availList'][0] += (KE + mass)
                self.info['E_depoList'][0] += depoE
                self.info['E_depoList_track'][0] += depoE_track
                self.info['Q_depoList'][0] += depoQ
                self.info['Q_depoList_th_75keV'][0] += depoQ_75keV
                self.info['Q_depoList_th_500keV'][0] += depoQ_500keV
                self.info['Q_depoList_dots_th_75keV'][0] += depoQ_dots
                self.info['Q_depoList_MBox'][0] += depoQ_MBox
                self.info['Q_depoList_MBox_th_75keV'][0] += depoQ_MBox_75keV
                self.info['Q_depoList_MBox_th_500keV'][0] += depoQ_MBox_500keV
                self.info['L_depoList_avg_220PEpMeV'][0] += depoL_220PEpMeV
                self.info['L_depoList_avg_180PEpMeV'][0] += depoL_180PEpMeV
                self.info['L_depoList_avg_140PEpMeV'][0] += depoL_140PEpMeV
                self.info['L_depoList_avg_100PEpMeV'][0] += depoL_100PEpMeV
                self.info['L_depoList_avg_35PEpMeV'][0] += depoL_35PEpMeV
                self.info['L_depoList_MBox_avg_220PEpMeV'][0] += depoL_MBox_220PEpMeV
                self.info['L_depoList_MBox_avg_180PEpMeV'][0] += depoL_MBox_180PEpMeV
                self.info['L_depoList_MBox_avg_140PEpMeV'][0] += depoL_MBox_140PEpMeV
                self.info['L_depoList_MBox_avg_100PEpMeV'][0] += depoL_MBox_100PEpMeV
                self.info['L_depoList_MBox_avg_35PEpMeV'][0] += depoL_MBox_35PEpMeV
                self.info['N_parList'][0] += 1
            elif (pdg == 2212):  # proton
                self.info['E_avail'] += KE
                self.info['E_availList'][1] += KE
                self.info['E_depoList'][1] += depoE
                self.info['E_depoList_track'][1] += depoE_track
                self.info['Q_depoList'][1] += depoQ
                self.info['Q_depoList_th_75keV'][1] += depoQ_75keV
                self.info['Q_depoList_th_500keV'][1] += depoQ_500keV
                self.info['Q_depoList_dots_th_75keV'][1] += depoQ_dots
                self.info['Q_depoList_MBox'][1] += depoQ_MBox
                self.info['Q_depoList_MBox_th_75keV'][1] += depoQ_MBox_75keV
                self.info['Q_depoList_MBox_th_500keV'][1] += depoQ_MBox_500keV
                self.info['L_depoList_avg_220PEpMeV'][1] += depoL_220PEpMeV
                self.info['L_depoList_avg_180PEpMeV'][1] += depoL_180PEpMeV
                self.info['L_depoList_avg_140PEpMeV'][1] += depoL_140PEpMeV
                self.info['L_depoList_avg_100PEpMeV'][1] += depoL_100PEpMeV
                self.info['L_depoList_avg_35PEpMeV'][1] += depoL_35PEpMeV
                self.info['L_depoList_MBox_avg_220PEpMeV'][1] += depoL_MBox_220PEpMeV
                self.info['L_depoList_MBox_avg_180PEpMeV'][1] += depoL_MBox_180PEpMeV
                self.info['L_depoList_MBox_avg_140PEpMeV'][1] += depoL_MBox_140PEpMeV
                self.info['L_depoList_MBox_avg_100PEpMeV'][1] += depoL_MBox_100PEpMeV
                self.info['L_depoList_MBox_avg_35PEpMeV'][1] += depoL_MBox_35PEpMeV
                self.info['N_parList'][1] += 1
            elif (pdg == 2112):  # neutron
                self.info['E_avail'] += KE
                self.info['E_availList'][2] += KE
                self.info['E_depoList'][2] += depoE
                self.info['E_depoList_track'][2] += depoE_track
                self.info['Q_depoList'][2] += depoQ
                self.info['Q_depoList_th_75keV'][2] += depoQ_75keV
                self.info['Q_depoList_th_500keV'][2] += depoQ_500keV
                self.info['Q_depoList_dots_th_75keV'][2] += depoQ_dots
                self.info['Q_depoList_MBox'][2] += depoQ_MBox
                self.info['Q_depoList_MBox_th_75keV'][2] += depoQ_MBox_75keV
                self.info['Q_depoList_MBox_th_500keV'][2] += depoQ_MBox_500keV
                self.info['L_depoList_avg_220PEpMeV'][2] += depoL_220PEpMeV
                self.info['L_depoList_avg_180PEpMeV'][2] += depoL_180PEpMeV
                self.info['L_depoList_avg_140PEpMeV'][2] += depoL_140PEpMeV
                self.info['L_depoList_avg_100PEpMeV'][2] += depoL_100PEpMeV
                self.info['L_depoList_avg_35PEpMeV'][2] += depoL_35PEpMeV
                self.info['L_depoList_MBox_avg_220PEpMeV'][2] += depoL_MBox_220PEpMeV
                self.info['L_depoList_MBox_avg_180PEpMeV'][2] += depoL_MBox_180PEpMeV
                self.info['L_depoList_MBox_avg_140PEpMeV'][2] += depoL_MBox_140PEpMeV
                self.info['L_depoList_MBox_avg_100PEpMeV'][2] += depoL_MBox_100PEpMeV
                self.info['L_depoList_MBox_avg_35PEpMeV'][2] += depoL_MBox_35PEpMeV
                self.info['N_parList'][2] += 1
            elif (pdg in [211, -211]):  # pi+/-
                self.info['E_avail'] += (KE + mass)
                self.info['E_availList'][3] += (KE + mass)
                self.info['E_depoList'][3] += depoE
                self.info['E_depoList_track'][3] += depoE_track
                self.info['Q_depoList'][3] += depoQ
                self.info['Q_depoList_th_75keV'][3] += depoQ_75keV
                self.info['Q_depoList_th_500keV'][3] += depoQ_500keV
                self.info['Q_depoList_dots_th_75keV'][3] += depoQ_dots
                self.info['Q_depoList_MBox'][3] += depoQ_MBox
                self.info['Q_depoList_MBox_th_75keV'][3] += depoQ_MBox_75keV
                self.info['Q_depoList_MBox_th_500keV'][3] += depoQ_MBox_500keV
                self.info['L_depoList_avg_220PEpMeV'][3] += depoL_220PEpMeV
                self.info['L_depoList_avg_180PEpMeV'][3] += depoL_180PEpMeV
                self.info['L_depoList_avg_140PEpMeV'][3] += depoL_140PEpMeV
                self.info['L_depoList_avg_100PEpMeV'][3] += depoL_100PEpMeV
                self.info['L_depoList_avg_35PEpMeV'][3] += depoL_35PEpMeV
                self.info['L_depoList_MBox_avg_220PEpMeV'][3] += depoL_MBox_220PEpMeV
                self.info['L_depoList_MBox_avg_180PEpMeV'][3] += depoL_MBox_180PEpMeV
                self.info['L_depoList_MBox_avg_140PEpMeV'][3] += depoL_MBox_140PEpMeV
                self.info['L_depoList_MBox_avg_100PEpMeV'][3] += depoL_MBox_100PEpMeV
                self.info['L_depoList_MBox_avg_35PEpMeV'][3] += depoL_MBox_35PEpMeV
                self.info['N_parList'][3] += 1
            elif (pdg == 111):  # pi0
                self.info['E_avail'] += (KE + mass)
                self.info['E_availList'][4] += (KE + mass)
                self.info['E_depoList'][4] += depoE
                self.info['E_depoList_track'][4] += depoE_track
                self.info['Q_depoList'][4] += depoQ
                self.info['Q_depoList_th_75keV'][4] += depoQ_75keV
                self.info['Q_depoList_th_500keV'][4] += depoQ_500keV
                self.info['Q_depoList_dots_th_75keV'][4] += depoQ_dots
                self.info['Q_depoList_MBox'][4] += depoQ_MBox
                self.info['Q_depoList_MBox_th_75keV'][4] += depoQ_MBox_75keV
                self.info['Q_depoList_MBox_th_500keV'][4] += depoQ_MBox_500keV
                self.info['L_depoList_avg_220PEpMeV'][4] += depoL_220PEpMeV
                self.info['L_depoList_avg_180PEpMeV'][4] += depoL_180PEpMeV
                self.info['L_depoList_avg_140PEpMeV'][4] += depoL_140PEpMeV
                self.info['L_depoList_avg_100PEpMeV'][4] += depoL_100PEpMeV
                self.info['L_depoList_avg_35PEpMeV'][4] += depoL_35PEpMeV
                self.info['L_depoList_MBox_avg_220PEpMeV'][4] += depoL_MBox_220PEpMeV
                self.info['L_depoList_MBox_avg_180PEpMeV'][4] += depoL_MBox_180PEpMeV
                self.info['L_depoList_MBox_avg_140PEpMeV'][4] += depoL_MBox_140PEpMeV
                self.info['L_depoList_MBox_avg_100PEpMeV'][4] += depoL_MBox_100PEpMeV
                self.info['L_depoList_MBox_avg_35PEpMeV'][4] += depoL_MBox_35PEpMeV
                self.info['N_parList'][4] += 1
            elif (pdg == 22):  # gamma
                self.info['E_avail'] += (KE + mass)
                self.info['E_availList'][5] += (KE + mass)
                self.info['E_depoList'][5] += depoE
                self.info['E_depoList_track'][5] += depoE_track
                self.info['Q_depoList'][5] += depoQ
                self.info['Q_depoList_th_75keV'][5] += depoQ_75keV
                self.info['Q_depoList_th_500keV'][5] += depoQ_500keV
                self.info['Q_depoList_dots_th_75keV'][5] += depoQ_dots
                self.info['Q_depoList_MBox'][5] += depoQ_MBox
                self.info['Q_depoList_MBox_th_75keV'][5] += depoQ_MBox_75keV
                self.info['Q_depoList_MBox_th_500keV'][5] += depoQ_MBox_500keV
                self.info['L_depoList_avg_220PEpMeV'][5] += depoL_220PEpMeV
                self.info['L_depoList_avg_180PEpMeV'][5] += depoL_180PEpMeV
                self.info['L_depoList_avg_140PEpMeV'][5] += depoL_140PEpMeV
                self.info['L_depoList_avg_100PEpMeV'][5] += depoL_100PEpMeV
                self.info['L_depoList_avg_35PEpMeV'][5] += depoL_35PEpMeV
                self.info['L_depoList_MBox_avg_220PEpMeV'][5] += depoL_MBox_220PEpMeV
                self.info['L_depoList_MBox_avg_180PEpMeV'][5] += depoL_MBox_180PEpMeV
                self.info['L_depoList_MBox_avg_140PEpMeV'][5] += depoL_MBox_140PEpMeV
                self.info['L_depoList_MBox_avg_100PEpMeV'][5] += depoL_MBox_100PEpMeV
                self.info['L_depoList_MBox_avg_35PEpMeV'][5] += depoL_MBox_35PEpMeV
                self.info['N_parList'][5] += 1
            elif (pdg == 1000020040):  # alpha
                self.info['E_avail'] += KE
                self.info['E_availList'][6] += KE
                self.info['E_depoList'][6] += depoE
                self.info['E_depoList_track'][6] += depoE_track
                self.info['Q_depoList'][6] += depoQ
                self.info['Q_depoList_th_75keV'][6] += depoQ_75keV
                self.info['Q_depoList_th_500keV'][6] += depoQ_500keV
                self.info['Q_depoList_dots_th_75keV'][6] += depoQ_dots
                self.info['Q_depoList_MBox'][6] += depoQ_MBox
                self.info['Q_depoList_MBox_th_75keV'][6] += depoQ_MBox_75keV
                self.info['Q_depoList_MBox_th_500keV'][6] += depoQ_MBox_500keV
                self.info['L_depoList_avg_220PEpMeV'][6] += depoL_220PEpMeV
                self.info['L_depoList_avg_180PEpMeV'][6] += depoL_180PEpMeV
                self.info['L_depoList_avg_140PEpMeV'][6] += depoL_140PEpMeV
                self.info['L_depoList_avg_100PEpMeV'][6] += depoL_100PEpMeV
                self.info['L_depoList_avg_35PEpMeV'][6] += depoL_35PEpMeV
                self.info['L_depoList_MBox_avg_220PEpMeV'][6] += depoL_MBox_220PEpMeV
                self.info['L_depoList_MBox_avg_180PEpMeV'][6] += depoL_MBox_180PEpMeV
                self.info['L_depoList_MBox_avg_140PEpMeV'][6] += depoL_MBox_140PEpMeV
                self.info['L_depoList_MBox_avg_100PEpMeV'][6] += depoL_MBox_100PEpMeV
                self.info['L_depoList_MBox_avg_35PEpMeV'][6] += depoL_MBox_35PEpMeV
                self.info['N_parList'][6] += 1
            else:
                self.info['E_avail'] += KE
                self.info['E_availList'][7] += KE
                self.info['E_depoList'][7] += depoE
                self.info['E_depoList_track'][7] += depoE_track
                self.info['Q_depoList'][7] += depoQ
                self.info['Q_depoList_th_75keV'][7] += depoQ_75keV
                self.info['Q_depoList_th_500keV'][7] += depoQ_500keV
                self.info['Q_depoList_dots_th_75keV'][7] += depoQ_dots
                self.info['Q_depoList_MBox'][7] += depoQ_MBox
                self.info['Q_depoList_MBox_th_75keV'][7] += depoQ_MBox_75keV
                self.info['Q_depoList_MBox_th_500keV'][7] += depoQ_MBox_500keV
                self.info['L_depoList_avg_220PEpMeV'][7] += depoL_220PEpMeV
                self.info['L_depoList_avg_180PEpMeV'][7] += depoL_180PEpMeV
                self.info['L_depoList_avg_140PEpMeV'][7] += depoL_140PEpMeV
                self.info['L_depoList_avg_100PEpMeV'][7] += depoL_100PEpMeV
                self.info['L_depoList_avg_35PEpMeV'][7] += depoL_35PEpMeV
                self.info['L_depoList_MBox_avg_220PEpMeV'][7] += depoL_MBox_220PEpMeV
                self.info['L_depoList_MBox_avg_180PEpMeV'][7] += depoL_MBox_180PEpMeV
                self.info['L_depoList_MBox_avg_140PEpMeV'][7] += depoL_MBox_140PEpMeV
                self.info['L_depoList_MBox_avg_100PEpMeV'][7] += depoL_MBox_100PEpMeV
                self.info['L_depoList_MBox_avg_35PEpMeV'][7] += depoL_MBox_35PEpMeV
                self.info['N_parList'][7] += 1

        # for track in self.tracks:
        #     pdg = track.GetPDGCode()
        #     depoE = track.energy['depoTotal']
        #     if (pdg in [13, -13, 11, -11]):
        #         self.info['E_depoList'][0] += depoE
        #     elif (pdg == 2212):
        #         self.info['E_depoList'][1] += depoE
        #     elif (pdg == 2112):
        #         self.info['E_depoList'][2] += depoE
        #     elif (pdg in [211, -211]):
        #         self.info['E_depoList'][3] += depoE
        #     elif (pdg == 111):
        #         self.info['E_depoList'][4] += depoE
        #     else:
        #         self.info['E_depoList'][5] += depoE

    # this is a primary function i design to find the ancester, but i think there are easier way of doing so
    def loopover(self, pdg, track):
        while pdg != 2112:
            ParentId = track.GetParentId()
            track = self.tracks[ParentId]
            pdg = track.GetPDGCode()
            # Break out of the loop if there is no parent (ParentId == -1)
            if ParentId == -1:
                break
        return pdg, track

    def PrintTracksEnergy(self):
        # Loop over all tracks
        switch = 0
        for i in range(self.tracks.size):
            track = self.tracks[i]
            pdg_original = track.GetPDGCode()
            pdg = track.GetPDGCode()  # Get the PDG code of the particle
            track_energy = track.energy.get('depoTotal', 0)

            # Trace back to the parent particle until we find a neutron (PDG code 2112)
            pdg, track = self.loopover(pdg, track)
            ParentId = track.GetParentId()
            while pdg == 2112 and ParentId != -1:
                track = self.tracks[ParentId]
                pdg = track.GetPDGCode()
                pdg, track = self.loopover(pdg, track)
                ParentId = track.GetParentId()
            if pdg == 2112 and ParentId == -1:
                # Instead of printing, store the tuple (or you can use a dict)
                switch = 1
            # results.append(\n)
        return switch

    # for now this function is totally useless because we are still wondering if neutron can deposit a lot of energy

    # select those tracks with neutron as ancester, we can use this function in antinumu samples because we assume hase on
    # energy deposit topology we are always able to find those edep which belongs to neutron
    #
    def select_the_right_track(self):
        results = []  # Initialize an empty list
        # Loop over all tracks
        right_track = []
        for i in range(self.tracks.size):
            track = self.tracks[i]
            for point in track.Points:
                x = point.GetPosition().X()
                y = point.GetPosition().Y()
                z = point.GetPosition().Z()
                t = point.GetPosition().T()
            pdg_original = track.GetPDGCode()
            pdg = track.GetPDGCode()  # Get the PDG code of the particle
            track_energy = track.energy.get('depoTotal', 0)
            # Trace back to the parent particle until we find a neutron (PDG code 2112)
            pdg, track = self.loopover(pdg, track)
            ParentId = track.GetParentId()
            while pdg == 2112 and ParentId != -1:
                track = self.tracks[ParentId]
                pdg = track.GetPDGCode()
                pdg, track = self.loopover(pdg, track)
                ParentId = track.GetParentId()
            if pdg == 2112 and ParentId == -1 and track_energy > self.energy_threshold:
                # Instead of printing, store the tuple (or you can use a dict)
                results.append((pdg_original, x, y, z, t , i ))
                right_track.append(i)
            # results.append(\n)
        #print(results)
        # print(right_track)
        return right_track


    # def guassian_fit(self):
    #     path = "/home/frunshi/frunshiworking/moneyplot/"
    #     cos,events_number,deg_list = self.collect_cos()
    #     cos_list_1 = [arr[0] for arr in cos]
    #     counts, bin_edges = np.histogram(cos_list_1,bins=180,density=True)
    #     bin_centers = (bin_edges[1:] + bin_edges[:-1])/2
    #     # the gaussian function we want:
    #     # amplitude * np.exp(-((x - mean) / stddev)**2 / 2)
    #     def half_gaussian(x, amplitude, mean, stddev):
    #         return np.where(x <=mean, amplitude*np.exp(-((x-mean)**2)/(2*stddev**2)), 0)
    #     mu = 1.0
    #     sigma = 0.1
    #     amplitude_guess = np.max(cos_list_1)
    #     p0_guess = [amplitude_guess, sigma]
    #     #here we take a guess that amplitude equals to 300, peak at 1.0 and sigma =20
    #     try :
    #         popt, pcov = curve_fit(lambda x, amp, std:half_gaussian(x, amp, mu, std)
    #                                , bin_centers,counts,p0=p0_guess)
    #         amplitude_fit, stddev_fit = popt
    #         mean_fit = mu
    #         neutrino_energy = 0.5
    #         x_fit = np.linspace(bin_edges[0],bin_edges[-1],180)
    #         y_fit = half_gaussian(x_fit, amplitude_fit, mean_fit, stddev_fit)
    #         plt.figure(figsize=(10,6))
    #         plt.hist(cos_list_1,bins=180,range=(-1,1),density=True,color = "black")
    #         plt.plot(x_fit,y_fit,'-',color = 'red',label = f"Gaussian fit(A={amplitude_fit},sigma={stddev_fit},mu = {mean_fit})")
    #         plt.title(f"Gaussian fit of {neutrino_energy}GeV AntiNumu events {events_number}")
    #         plt.xlabel("cos of reco and true directions")
    #         plt.ylabel("number of events(density)")
    #         plt.legend()
    #         plt.savefig(path+f"/gaussian_fit_of_{neutrino_energy}GeVAntiNumu.pdf")
    #     except RuntimeError as e:
    #         print(f"Fit failed: {e}")
    #         # As a fallback, plot the data to see why it might have failed
    #         plt.figure(figsize=(10, 6))
    #         plt.hist(cos_list_1, bins=100, range=(-1, 1), density=True, color="black", label="Data")
    #         plt.title("Data Visualization (Fit Failed)")
    #         plt.savefig(path+"data_fit_failed.pdf")
    # def gaussian_fit(self):


    def PrintDepo(self, i):
        # print(f"{self.depos.size} depo points stored in total")
        depo = self.depos[i]
        contrib = depo.GetContributors()
        primaryId = depo.GetPrimaryId()
        edep = depo.GetEnergyDeposit()
        # edepSecond = depo.GetSecondaryDeposit()
        trkLength = depo.GetTrackLength()
        print(contrib, "%6d %.2f %.2f" % (primaryId, edep, trkLength))



    def Use_PrintTrack(self):
        for i in range(self.tracks.size):
            t = self.PrintTrack(i)
        return
    # -----------------------
    def PrintTrack(self, trkId):
        self.PrintTracks(trkId, trkId + 1)
        track = self.tracks[trkId]
        children = track.association['children']
        # print(f'children: {children}')
        for childId in children:
            child = self.tracks[childId]
            name = child.GetName()
            mom = child.GetInitialMomentum()
            KE = mom.E() - mom.M()
            # print(f"{childId} {name}: {KE:.2f} MeV, {child.energy['depoTotal']:.2f} MeV, {child.association['children']}")

        # print(f"{track.Points.size()} points stored in track {trkId}")
        for point in track.Points:
            x = point.GetPosition().X()
            y = point.GetPosition().Y()
            z = point.GetPosition().Z()
            t = point.GetPosition().T()
            print(f"{point.GetProcess()}, {point.GetSubprocess()}, {x}, {y}, {z}, {t}")

        return t  # last point is capture time


    # ------------------------
    def PrintTracks(self, start=0, stop=-1):
        print(f"{self.tracks.size} trajectories stored", )

        print(f"{'pdg':>8}{'name':>8}{'trkId':>6}{'parId':>6}{'acId':>6}{'KE':>10}{'selfDepo':>10}{'allDepo':>10}")
        print(f"{'':>8}{'':>8}{'':>6}{'':>6}{'':>6}{'[MeV]':>10}{'[MeV]':>10}{'[MeV]':>10}")
        print('-'*(8+8+6+6+6+10+10+10))

        neutrontrkId = -1
        neutronKE = -1

        for track in self.tracks[start:stop]:
            pdg = track.GetPDGCode()
            name = track.GetName()
            trkId = track.GetTrackId()
            parId = track.GetParentId()
            mom = track.GetInitialMomentum()
            mass = mom.M()
            KE = mom.E() - mass
            ancestor = track.association['ancestor']
            selfDepo = track.energy['depoTotal']
            allDepo = self.GetEnergyDepoWithDesendents(trkId)
            print(f"{pdg:>8d}{name:>8s}{trkId:>6d}{parId:>6d}{ancestor:>6d}{KE:>10.2f}{selfDepo:>10.2f}{allDepo:>10.2f}")
            # Take the last neutron and check its KE and capture time
            if pdg == 2112:
                neutrontrkId = trkId
                neutronKE = KE

        print('-'*(8+8+6+6+6+10+10+10))
        return [neutrontrkId, neutronKE]



    def neutron_edep_information(self, energy_threshold):

        info = []
        MuonTrack = []
        Entry = self.currentEntry
        if self.SelectNeutronEvent() != 1:
            return info
        mm2cm = 0.1
        neutron_dir,mom = self.read_neutron_momentum()

        for i in range(self.tracks.size):
            track_origin = self.tracks[i]
            # Trace back to find its neutron ancestor
            pdg, track = self.loopover(track_origin.GetPDGCode(), track_origin)
            anc_idx = track.association.get('ancestor')
            if anc_idx is None or anc_idx >= len(self.tracks):
                continue
            if self.tracks[anc_idx].GetPDGCode() != 2112:
                continue

            for di in track_origin.association.get('depoList', []):
                if not (0 <= di < len(self.depos)):
                    continue
                depo = self.depos[di]
                edep = depo.GetEnergyDeposit()
                if edep < energy_threshold:
                    continue
                ancestor = track.association['ancestor']
                pdg_check = self.tracks[ancestor].GetPDGCode()

                start, stop = depo.GetStart(), depo.GetStop()
                midpoint = np.array([
                    (start.X() + stop.X()) / 2 * mm2cm,
                    (start.Y() + stop.Y()) / 2 * mm2cm,
                    (start.Z() + stop.Z()) / 2 * mm2cm
                ])
                midpointlength = np.linalg.norm(midpoint)
                # Compute angle to true neutron direction
                angle = self.compute_angle_between(neutron_dir, midpoint)
                info.append([midpointlength, edep, angle])
            MuonTrack = np.array(MuonTrack)
        return info

    def proton_withoutthreshold(self):

        def nearest_depo_position(track_id):
            """Return depo midpoint closest to the origin for the given track."""

            if track_id is None or track_id < 0:
                return None, None
            if not hasattr(self, 'depos') or track_id >= len(self.tracks):
                return None, None

            mm2cm = 0.1
            closest_pos = None
            closest_time = None
            min_distance = float('inf')
            track = self.tracks[track_id]
            for depo_index in track.association.get('depoList', []):
                if depo_index < 0 or depo_index >= len(self.depos):
                    continue
                depo = self.depos[depo_index]
                start, stop = depo.GetStart(), depo.GetStop()
                midpoint = np.array([
                    (start.X() + stop.X()) / 2 * mm2cm,
                    (start.Y() + stop.Y()) / 2 * mm2cm,
                    (start.Z() + stop.Z()) / 2 * mm2cm,
                ])
                distance = np.linalg.norm(midpoint)
                if distance < min_distance:
                    min_distance = distance
                    closest_pos = midpoint
                    closest_time = (start.T() + stop.T()) / 2
            return closest_pos, closest_time

        energy_track = []
        proton_number = 0
        for particle in self.vertex.Particles:
            pdg = particle.GetPDGCode()
            if pdg == 2212:
                proton_number += 1
                mom = particle.GetMomentum()
                kinetic_energy = mom.E() - mom.M()
                momentum = np.array([mom.X(), mom.Y(), mom.Z()])

                track_id = particle.GetTrackId() if hasattr(particle, 'GetTrackId') else None
                nearest_position = nearest_depo_position(track_id)

                if kinetic_energy > 0:
                    energy_track.append((True, proton_number, momentum, nearest_position))
                else:
                    energy_track.append((False, proton_number, np.zeros(3), nearest_position))
        return energy_track

    def proton_threshold(self):
        """Classify protons by kinetic energy and return truth momentum vectors.

        Returns a list of tuples ``(is_high_energy, index, momentum_vec,
        nearest_position, nearest_time)``. ``is_high_energy`` is True for KE > 70
        MeV, ``index`` is a 1-based counter, ``momentum_vec`` is the 3-vector in
        MeV/c (``np.zeros(3)`` for low energy protons), ``nearest_position`` is
        the midpoint (in cm) of the proton energy deposition closest to the
        vertex (``None`` when unavailable), and ``nearest_time`` is the average
        time of that deposition (``None`` when unavailable)."""

        def nearest_depo_position(track_id):
            """Return depo midpoint closest to the origin for the given track."""

            if track_id is None or track_id < 0:
                return None, None
            if not hasattr(self, 'depos') or track_id >= len(self.tracks):
                return None, None

            mm2cm = 0.1
            closest_pos = None
            closest_time = None
            min_distance = float('inf')
            track = self.tracks[track_id]
            for depo_index in track.association.get('depoList', []):
                if depo_index < 0 or depo_index >= len(self.depos):
                    continue
                depo = self.depos[depo_index]
                start, stop = depo.GetStart(), depo.GetStop()
                midpoint = np.array([
                    (start.X() + stop.X()) / 2 * mm2cm,
                    (start.Y() + stop.Y()) / 2 * mm2cm,
                    (start.Z() + stop.Z()) / 2 * mm2cm,
                ])
                distance = np.linalg.norm(midpoint)
                if distance < min_distance:
                    min_distance = distance
                    closest_pos = midpoint
                    closest_time = (start.T() + stop.T()) / 2
            return closest_pos, closest_time

        energy_track = []
        proton_number = 0
        for particle in self.vertex.Particles:
            pdg = particle.GetPDGCode()
            if pdg == 2212:
                proton_number += 1
                mom = particle.GetMomentum()
                kinetic_energy = mom.E() - mom.M()
                momentum = np.array([mom.X(), mom.Y(), mom.Z()])

                track_id = particle.GetTrackId() if hasattr(particle, 'GetTrackId') else None
                nearest_position, nearest_time = nearest_depo_position(track_id)

                if kinetic_energy > self.proton_trackable_length:
                    energy_track.append((True, proton_number, momentum, nearest_position, nearest_time))
                else:
                    energy_track.append((False, proton_number, np.zeros(3), nearest_position, nearest_time))
        return energy_track


    def pion_threshold(self):
        """Classify charged pions by kinetic energy and return truth momentum vectors.

        Returns a list of tuples ``(is_high_energy, index, momentum_vec)`` where
        ``is_high_energy`` is True for KE > 50 MeV, ``index`` is a 1-based counter,
        and ``momentum_vec`` is the 3-vector in MeV/c (``np.zeros(3)`` for low energy
        pions)."""
        energy_track = []
        pion_number = 0
        for particle in self.vertex.Particles:
            pdg = particle.GetPDGCode()
            if abs(pdg) == 211:
                pion_number += 1
                mom = particle.GetMomentum()
                kinetic_energy = mom.E() - mom.M()
                momentum = np.array([mom.X(), mom.Y(), mom.Z()])
                if kinetic_energy > 0.00001:
                    energy_track.append((True, pion_number, momentum))
                else:
                    energy_track.append((False, pion_number, np.zeros(3)))
        return energy_track



    def lepton_only_reconstruction_antinue_per_event(self):

        neutrino_direction = np.array([0,0,1])
        electron_direction,_ = self.read_electron_direction()
        electron_direction = np.array(electron_direction)
        angle_deviation = self.compute_angle_between(neutrino_direction, electron_direction)
        return np.degrees(angle_deviation)



    def find_nearest_Eedep(self,cone_size):

        """this function is a replica of find_nearest_discontinuous_medep"""


        electron_direction,_ = self.read_electron_direction()
        mm2cm = 0.1
        shortest_distance = float('inf')
        edep_position = [float('inf'), float('inf'), float('inf')]
        threshold = self.energy_threshold
        for track in self.tracks:
            pdg = track.GetPDGCode()
            ancestor = track.association['ancestor']
            ancestor_pdg = self.tracks[ancestor].GetPDGCode()
            if abs(pdg)==11 or abs(ancestor_pdg)==11:
                for di in track.association.get('depoList',[]):
                    depo = self.depos[di]
                    edep = depo.GetEnergyDeposit()
                    if edep < threshold:
                        continue
                    start, stop = depo.GetStart(), depo.GetStop()
                    midpoint = np.array([
                        (start.X() + stop.X()) / 2 * mm2cm,
                        (start.Y() + stop.Y()) / 2 * mm2cm,
                        (start.Z() + stop.Z()) / 2 * mm2cm
                    ])
                    # Compute angle to true neutron direction
                    angle = self.compute_angle_between(electron_direction, midpoint)
                    angle = np.degrees(angle)
                    if angle <= cone_size:
                        continue
                    MidPointLength = np.linalg.norm(midpoint)
                    if shortest_distance > MidPointLength:
                        shortest_distance = MidPointLength
                        edep_position = midpoint
        return edep_position, shortest_distance


    def find_nearest_nedep_EV(self,cone_angle):
        mm2cm = 0.1
        electron_direction, _ = self.read_electron_direction()
        nearest_position = [float('inf'), float('inf'), float('inf')]
        nearest_length = float('inf')
        for i in range(self.tracks.size):
            track = self.tracks[i]
            pdg = track.GetPDGCode()
            ancestor = track.association['ancestor']
            ancestor_pdg = self.tracks[ancestor].GetPDGCode()
            if pdg ==2112 or ancestor_pdg==2112:
                for di in track.association['depoList']:
                    depo = self.depos[di]
                    Edep = depo.GetEnergyDeposit()
                    if Edep >= self.energy_threshold:
                        start, stop = depo.GetStart(), depo.GetStop()
                        t = (start.T() + stop.T()) / 2
                        midpoint = np.array([
                            (start.X() + stop.X()) / 2 * mm2cm,
                            (start.Y() + stop.Y()) / 2 * mm2cm,
                            (start.Z() + stop.Z()) / 2 * mm2cm
                        ])

                        midpoint_direction = midpoint/np.linalg.norm(midpoint)
                        angle_different = self.compute_angle_between(electron_direction, midpoint_direction)
                        angle = np.degrees(angle_different)
                        if angle <= cone_angle:
                            continue
                        else:
                            midpoint_length = np.linalg.norm(midpoint)
                            if midpoint_length < nearest_length:
                                nearest_length = midpoint_length
                                nearest_position = midpoint
            proton_track = self.proton_threshold()
            if proton_track:
                for is_track, _, _, proton_position, _ in proton_track:
                    if is_track or proton_position is None:
                        continue
                    proton_position = np.array(proton_position, dtype=float)
                    if not np.all(np.isfinite(proton_position)):
                        continue
                    distance = np.linalg.norm(proton_position)
                    if distance <= 0:
                        continue
                    angle = self.compute_angle_between(electron_direction, proton_position)
                    if np.degrees(angle) > cone_angle and distance < nearest_length:
                        nearest_length = distance
                        nearest_position = proton_position

        return nearest_position, nearest_length

    def find_nearest_depo_EV(self,electron_energy):

        composite_neutron_ke, composite_neutron_momentum, real_direction = self.composite_neutron()
        electron_direction, _ = self.read_electron_direction()
        if real_direction is not None:
            real_direction = np.array(real_direction, dtype=float)
            if not np.all(np.isfinite(real_direction)) or np.linalg.norm(real_direction) == 0:
                real_direction = None
                print("we cannot find the real direction of neutron from the vertex in the simulation")
        cone_angle_list = Event.cone_applying_for_electron
        n = round(electron_energy * 10)
        n = max(0, min(n, len(cone_angle_list) - 1))
        cone_angle = cone_angle_list[n]

        Eedep_result = self.find_nearest_Eedep(cone_angle)
        pedep_result = self.find_nearest_proton_EV(cone_angle)
        nedep_result = self.find_nearest_nedep_EV(cone_angle)
        noise_result = self.find_nearest_noises(
            self.energy_threshold,
            reference_direction=electron_direction,
            cone_angle=cone_angle)

        # the readout should all be in the form of position + length


        def parse_result(result):
            if result is None:
                return  None, float('inf')
            position, length = result
            position = np.array(position, dtype=float)
            if not np.all(np.isfinite(position)):
                return  None, float('inf')
            distance = np.linalg.norm(position)
            if distance <= 0:
                return  None, float('inf')
            return position, distance

        position_e, distance_e = parse_result(Eedep_result)
        position_p, distance_p = parse_result(pedep_result)
        position_n, distance_n = parse_result(nedep_result)
        position_noise, distance_noise = parse_result(None)
        if noise_result is not None:
            _, noise_position = noise_result
            if noise_position is not None:
                noise_position = np.array(noise_position, dtype=float)
                if np.all(np.isfinite(noise_position)):
                    noise_distance = np.linalg.norm(noise_position)
                    if noise_distance > 0:
                        position_noise = noise_position
                        distance_noise = noise_distance
        # print("position_e", position_e)
        # print("position_p", position_p)
        # print("distance_n", distance_n)
        candidates = []
        if position_e is not None:
            candidates.append({
                'label': 'electron',
                'position': position_e,
                'distance': distance_e,
                'track_pdg': -11,
                'ancestor_pdg': -11,
                'origin': 'Eedep'
            })
        if position_p is not None:
            candidates.append({
                'label': 'proton',
                'position': position_p,
                'distance': distance_p,
                'track_pdg': 2212,
                'ancestor_pdg': 2212,
                'origin': 'proton'
            })
        if position_n is not None:
            candidates.append({
                'label': 'neutron',
                'position': position_n,
                'distance': distance_n,
                'track_pdg': 2112,
                'ancestor_pdg': 2112,
                'origin': 'nedep'
            })
        if position_noise is not None:
            candidates.append({
                'label': 'noise',
                'position': position_noise,
                'distance': distance_noise,
                'track_pdg': 20110715,
                'ancestor_pdg': 20110715,
                'origin': 'noise'
            })

        first_edep = None
        first_edep_source = None
        if candidates:
            winner = min(candidates, key=lambda c: c['distance'])
            first_edep = winner['position']
        else:
            print('no nearest depo outside cone candidates')
            return None
        reconstructed_direction = first_edep/np.linalg.norm(first_edep)
            # print("reconstructed direction is :", reconstructed_direction)
            # print("real direction is :", real_direction)
        if real_direction is None:
            return None, reconstructed_direction
        cos = np.dot(real_direction, reconstructed_direction)/(np.linalg.norm(real_direction)* np.linalg.norm(reconstructed_direction))
        angle_deg = np.degrees(np.arccos(
            np.clip(np.dot(real_direction, reconstructed_direction) / (np.linalg.norm(real_direction) * np.linalg.norm(reconstructed_direction)), -1.0, 1.0)))

        #consider one question, when cos of two direction co1uld be zero?
        #first scenario, real_direction = 0
        # print("cos:", cos)

        return angle_deg,reconstructed_direction






























    def find_discontinuous_muons(self):
        """this is an improvement from find_muon_position function, as it use the function[find_michel_electron]first
        to make sure only the discontinuous part of the muon track remains. In the real detector we need to find an algorithm
        to figure out how to remove the continuous part of the muon track, but here in the reconstruction part it is not necessary
        and time-consuming"""

        MuonTrack = []
        mm2cm = 0.1
        threshold = self.energy_threshold
        midpoint, electron_direction, counter, michel_electron_id, midpoint_end = self.find_Michel_Electron()
        if counter == 0 or michel_electron_id is None:
            return MuonTrack
        track = self.tracks[michel_electron_id]
        initial_children = track.association['children']
        all_descendants = set()
        tracks_to_visit = list(initial_children)  # Use a list as a queue, starting with direct children

        while tracks_to_visit:
            current_track_id = tracks_to_visit.pop(0)  # Get the next track to check
            if current_track_id not in all_descendants:
                all_descendants.add(current_track_id)
                descendant_track = self.tracks[current_track_id]
                # Find the next generation of children and add them to the queue
                grandchildren = descendant_track.association['children']
                if grandchildren:
                    tracks_to_visit.extend(grandchildren)

        # print(f"All descendants of the Michel_Electron are {all_descendants}")
        depoList = track.association['depoList']
        if len(depoList) <=2:
            return MuonTrack
        else:
            depo_end_ME = self.depos[depoList[-1]]
            t_ME = (depo_end_ME.GetStart().T() + depo_end_ME.GetStop().T()) / 2
            for i in all_descendants:
                track = self.tracks[i]
                depoList = track.association['depoList']
                # print(f"this is children {i} \n")
                for di in depoList:
                    depo_start = self.depos[di]
                    edep = depo_start.GetEnergyDeposit()  # MeV
                    # print(f"the energy of energy deposition of children particle is: {edep}")
                    t = (depo_start.GetStart().T() + depo_start.GetStop().T()) / 2
                    # print(f"the time of energy deposition of children particle is: {t}")
                    if edep >= threshold and t > t_ME:
                        x = (depo_start.GetStart().X() + depo_start.GetStop().X()) / 2 * mm2cm
                        y = (depo_start.GetStart().Y() + depo_start.GetStop().Y()) / 2 * mm2cm
                        z = (depo_start.GetStart().Z() + depo_start.GetStop().Z()) / 2 * mm2cm
                        MuonTrack.append([x,y,z])
            return MuonTrack



    def find_nearest_discontinuous_medep(self, energy_threshold: Optional[float] = None) :
        threshold = self._resolve_energy_threshold(energy_threshold)
        MuonTrack = self.find_discontinuous_muons()
        if len(MuonTrack) == 0:
            return None
        nearest_distance = 10000
        i = 0
        nearest_position = np.array([0,0,0])
        for points in MuonTrack:
            x =points[0]
            y = points[1]
            z = points[2]
            distance = np.linalg.norm(points)
            i +=1
            if distance < nearest_distance:
                nearest_distance = distance
                nearest_position = np.array([x,y,z])

        muon_direction,mom = self.read_muon_direction_and_amplitude()
        muon_direction = np.array(muon_direction)
        angle_degrees = np.degrees(
        np.arccos(np.clip(np.dot(muon_direction, nearest_position) / (np.linalg.norm(muon_direction) * np.linalg.norm(nearest_position)), -1.0, 1.0)))
        if nearest_distance < self.minimum_length or nearest_distance > self.inspection_length_direction_side*100:
            angle_degrees = 180
            nearest_position = np.array([100, 10000, 10000])
        if nearest_position is None or angle_degrees is None:
            print("can't find nearest muon position, it is right on the face, not that there is no discontinuous muon")
            angle_degrees = 180
            nearest_position = np.array([100,10000,10000])

        return angle_degrees, nearest_position


    def find_Michel_Electron(self):

        """this function is used to find the starting point of the Michel Electron, and print out its position
        this point, as previously discussed, is what we called "turning point"."""

        electron_direction = np.array([0, 0, 0])
        midpoint = np.array([0, 0, 0])
        midpoint_end = np.array([0, 0, 0])
        counter = 0
        trig = 0
        KE_michelelectron = 0
        michel_electron_id = None
        for i in range(self.tracks.size):
            track = self.tracks[i]
            mm2cm = 0.1
            pdg = track.GetPDGCode()  # Get the PDG code of the particle
            ancestor = track.association['ancestor']
            parId = track.GetParentId()
            track_parent = self.tracks[parId]
            pdg_parent = track_parent.GetPDGCode()
            mom = track.GetInitialMomentum()
            mass = mom.M()
            KE_max = mom.E() - mass

            # here we demand the michel electron shall be a positron
            # in antinumu event since its parent shall be a mu+

            if pdg == 11 and pdg_parent in (13, -13) and KE_max > 20:
                depo_indices = track.association.get('depoList', [])
                if len(depo_indices) < 2:
                    # Need at least two depositions to define the Michel start/end
                    continue
                di_start = depo_indices[1]
                di_end = depo_indices[-1]
                if not (0 <= di_start < len(self.depos)):
                    continue
                if not (0 <= di_end < len(self.depos)):
                    continue
                michel_electron_id = i
                trig += 1
                KE_michelelectron = mom.E() - mass
                depo = self.depos[di_start]
                depo_end = self.depos[di_end]
                # Compute midpoint in cm
                start, stop = depo.GetStart(), depo.GetStop()
                start_end ,start_stop = depo_end.GetStart(), depo_end.GetStop()

                midpoint = np.array([
                    (start.X() + stop.X()) / 2 * mm2cm,
                    (start.Y() + stop.Y()) / 2 * mm2cm,
                    (start.Z() + stop.Z()) / 2 * mm2cm
                ])
                midpoint_end =  np.array([
                    (start_end.X() + start_stop.X()) / 2 * mm2cm,
                    (start_end.Y() + start_stop.Y()) / 2 * mm2cm,
                    (start_end.Z() + start_stop.Z()) / 2 * mm2cm
                ])
                    # print("value given to midpoint",midpoint)
                px = track.GetInitialMomentum().X()
                py = track.GetInitialMomentum().Y()
                pz = track.GetInitialMomentum().Z()
                norm = np.linalg.norm([px,py,pz])
                electron_direction = np.array([px/norm, py/norm, pz/norm])
                # print("value given to electron",electron_direction)
        if trig >=2 :
            print("there is more than one michel electron!")
            print(f"trig = {trig}")
        # elif trig == 0 :
        #     print("No Michel electron!")
        else:
            ###this case there is only one michel electron
            #in principle we only want this to happen
            if np.array_equal(midpoint, [0, 0, 0])  or np.array_equal(electron_direction, [0, 0, 0]):
                # print("unable to find a Michel Electron here.")
                fjaljfla=1
            else:


                if KE_michelelectron > 10:
                    # print(f"there is a Michel Electron in this event, with the position{midpoint}(cm). And the direction {electron_direction}")
                    counter = 1
                # else:
                #     print("Michel electron of too low energy, should ignore")
        return midpoint,electron_direction,counter,michel_electron_id,midpoint_end

    def read_vertex_process(self):
        process = self.vertex.GetReaction()
        print(process)

    def count_michel_electron(self,neutrino_energy):
        total_michel_per_thousand_event = 0
        for i in range(0,1000,1):
            self.Jump(i)
            midpoint,electron_direction,counter,michel_electron_id,midpoint_end = self.find_Michel_Electron()
            total_michel_per_thousand_event += counter
        print(f"in 1000 events there are {total_michel_per_thousand_event} events with michel electron,{total_michel_per_thousand_event/10}% percent")




    def earliest_interaction_time(self):
        trig = self.SelectNeutronEvent()
        mm2m = 0.001
        if trig != 1:
            # print(f"Event {self.currentEntry}: Not a single-neutron event, skipping.")
            return None  # Return a value indicating it was skipped
        earliest_interaction_time = 1000000000000
        earliest_interaction_point = [0,0,0]
        for i in range(self.tracks.size):
            track = self.tracks[i]
            ancestor_id = track.association.get('ancestor')
            ancestor = self.tracks[ancestor_id]
            pdg_original = ancestor.GetPDGCode()
            pdg = track.GetPDGCode()
            depoList = track.association['depoList']
            for di in depoList:
                if di < 0 or di >= len(self.depos):
                    continue
                else:
                    depo = self.depos[di]

                    x = (depo.GetStart().X() + depo.GetStop().X()) / 2 * mm2m
                    y = (depo.GetStart().Y() + depo.GetStop().Y()) / 2 * mm2m
                    z = (depo.GetStart().Z() + depo.GetStop().Z()) / 2 * mm2m
                    t = (depo.GetStart().T() + depo.GetStop().T()) / 2
                    edep = depo.GetEnergyDeposit()
                    if pdg_original != 2112 and pdg!=2112:
                        continue
                    if  t < earliest_interaction_time and edep >self.energy_threshold :
                        earliest_interaction_time = t
                        earliest_interaction_point = [x,y,z]
        if earliest_interaction_time == 1000000000000:
            print("filtered out this event")
        return earliest_interaction_time,earliest_interaction_point




    def nearest_vs_earliest(self, plotter_object):
        """

        Compares the nearest and earliest neutron depositions for the CURRENTLY loaded event.
        The main event loop should be handled by the calling script.

        """
        trig = self.SelectNeutronEvent()
        if trig != 1:
            # print(f"Event {self.currentEntry}: Not a single-neutron event, skipping.")
            return None  # Return a value indicating it was skipped

        result = self.find_nearest_discontinuous_neutron(plotter_object)

        # Check if the function returned valid data
        if result is None or result[0] is None:
            return None  # Propagate the skip

        nearest_position, position_angle, nearest_t = result
        earliest_interaction_time,earliest_interaction_position = self.earliest_interaction_time()
        if nearest_t == earliest_interaction_time:
            earliest_nearest_angle_distribution = 0
        else :
            earliest_nearest_angle_distribution =  np.degrees(np.arccos(
            np.clip(np.dot(earliest_interaction_position, nearest_position) / (np.linalg.norm(earliest_interaction_position) *
                                                                               np.linalg.norm(nearest_position)), -1.0, 1.0)))



        # Return a dictionary or tuple with the results for this single event
        return {
            "is_earliest": earliest_interaction_time >= nearest_t,
            "nearest_time": nearest_t,
            "earliest_time": earliest_interaction_time ,
            "angle_distribution": earliest_nearest_angle_distribution
        }
    def find_nearest_proton_EV(self,cone_angle):

        proton_list = self.proton_threshold()

        if not proton_list:
            return None

        results = self.read_electron_direction()
        if results is None:
            return None

        electron_direction,electron_momentum_amplitude = results
        if electron_momentum_amplitude <= 0 :
            return None

        best_position = None
        best_distance = float('inf')
        for is_track, proton_id, momentum, nearest_position, nearest_time in proton_list:
            if is_track:
                continue
            if nearest_position is None:
                continue

            position = np.array(nearest_position, dtype=float)

            if not np.all(np.isfinite(position)):
                continue
            distance = np.linalg.norm(position)
            if distance <= 0:
                continue
            pos_dir = position / distance

            angle = np.degrees(np.arccos(np.clip(np.dot(electron_direction, pos_dir), -1.0, 1.0)))
            if angle < cone_angle:
                continue
            if distance < best_distance:
                best_distance = distance
                best_position = position

        if best_position is None:
            return None
        if best_distance <= self.minimum_length or best_distance >= self.inspection_length_direction_side*100:
            return None

        return best_position,np.linalg.norm(best_position)
    def find_nearest_discontinuous_proton_EV(self):
        proton_list = self.proton_threshold()
        if not proton_list:
            return None

        electron_direction,_ = self.read_electron_direction()
        if electron_direction is None:
            return None



        best_position = None
        best_distance = float('inf')
        for is_track, proton_id, momentum, nearest_position, nearest_time in proton_list:
            if is_track:
                continue
            if nearest_position is None:
                continue
            position = np.array(nearest_position, dtype=float)
            if not np.all(np.isfinite(position)):
                continue
            distance = np.linalg.norm(position)
            if distance <= 0:
                continue
            if distance < best_distance:
                best_distance = distance
                best_position = position

        if best_position is None:
            return None
        if best_distance <= self.minimum_length or best_distance >= self.inspection_length_direction_side*100:
            return None

        pos_dir = best_position / np.linalg.norm(best_position)
        angle = np.degrees(np.arccos(np.clip(np.dot(electron_direction, pos_dir), -1.0, 1.0)))
        return angle, best_position





    def find_nearest_discontinuous_proton(self,threshold):
        """Return the closest low-energy proton deposition to the vertex.

        The search ignores high-energy protons (``is_track`` is True in
        ``proton_threshold``) and candidates without an associated deposition.
        The result is ``(angle_deg, position_cm)`` where ``angle_deg`` is the
        opening angle relative to the muon direction and ``position_cm`` is the
        midpoint vector in centimeters. Returns ``None`` when no suitable proton
        deposition is available."""

        proton_list = self.proton_threshold()
        if not proton_list:
            return None

        results = self.read_muon_direction_and_amplitude()
        if results is None:
            return None

        muon_direction, _ = results
        muon_direction = np.array(muon_direction, dtype=float)
        mu_norm = np.linalg.norm(muon_direction)
        if mu_norm <= 0 or not np.all(np.isfinite(muon_direction)):
            return None
        mu_dir = muon_direction / mu_norm

        best_position = None
        best_distance = float('inf')
        for is_track, proton_id, momentum, nearest_position, nearest_time in proton_list:
            if is_track:
                continue
            if nearest_position is None:
                continue
            position = np.array(nearest_position, dtype=float)
            if not np.all(np.isfinite(position)):
                continue
            distance = np.linalg.norm(position)
            if distance <= 0:
                continue
            if distance < best_distance:
                best_distance = distance
                best_position = position

        if best_position is None:
            return None
        if best_distance <= self.minimum_length or best_distance >= self.inspection_length_direction_side*100:
            return None
        pos_dir = best_position / np.linalg.norm(best_position)
        angle = np.degrees(np.arccos(np.clip(np.dot(mu_dir, pos_dir), -1.0, 1.0)))
        return angle, best_position
    def find_nearest_noises(self,energy_threshold, reference_direction=None, cone_angle=None):
        nearest_position = None
        nearest_distance = 100000
        noises_list = self.generate_random_energy_from_spectrum(
            number_of_noises_point=self.noise_points,
            radius_of_inspection=self.inspection_length_direction_side,
            energy_threshold=energy_threshold)
        angle_degrees = 0
        if noises_list is None or noises_list == []:
            return None

        if reference_direction is None or cone_angle is None:
            cone_method = self.cone_applying_for_muon
            reference_direction, _ = self.read_muon_direction_and_amplitude()
            muon_energy = self.read_muon_energy()
            index = int(round(muon_energy / 100, 0))
            index = max(0, min(index, len(cone_method) - 1))
            cone_angle = cone_method[index]
        else:
            reference_direction = np.array(reference_direction, dtype=float)
            if not np.all(np.isfinite(reference_direction)):
                reference_direction = np.zeros(3)

        for i in noises_list:

            distance = np.linalg.norm(np.array(i[:3]))
            if distance > nearest_distance:
                continue
            noise_direction = np.array(i[:3]) / distance
            angle = self.compute_angle_between(reference_direction, noise_direction)
            deviation_angle = np.degrees(angle)
            if cone_angle > deviation_angle:
                continue
            nearest_position = np.array(i[:3])
            angle_degrees = deviation_angle

        return angle_degrees, nearest_position




    def find_nearest_depo(self,muon_energy,cone_method,energy_threshold: Optional[float] = None):

        self._last_first_edep_source = None
        threshold = self._resolve_energy_threshold(energy_threshold)
        composite_neutron_ke, composite_neutron_momentum, real_direction = self.composite_neutron()
        if real_direction is not None:
            real_direction = np.array(real_direction, dtype=float)
            if not np.all(np.isfinite(real_direction)) or np.linalg.norm(real_direction) == 0:
                real_direction = None
        cone_angle_list = cone_method
        n = round(muon_energy * 10)
        n = max(0, min(n, len(cone_angle_list) - 1))
        cone_angle = cone_angle_list[n]

        medep_result = self.find_nearest_discontinuous_medep(threshold)
        pedep_result = self.find_nearest_discontinuous_proton(threshold)
        nedep_result = self.find_nearest_discontinuous_nedep(threshold)
        noise_result = self.find_nearest_noises(threshold)

        def parse_result(result):
            if result is None:
                return 180.0, None, float('inf')
            angle, position = result
            position = np.array(position, dtype=float)
            if not np.all(np.isfinite(position)):
                return 180.0, None, float('inf')
            distance = np.linalg.norm(position)
            if distance <= 0:
                return angle, None, float('inf')
            return angle, position, distance

        medep_angle, position_m, distance_m = parse_result(medep_result)
        pedep_angle, position_p, distance_p = parse_result(pedep_result)
        nedep_angle, position_n, distance_n = parse_result(nedep_result)
        noises_angle,position_noises,distance_noise = parse_result(noise_result)

        candidates = []
        if position_m is not None:
            candidates.append({
                'label': 'muon',
                'angle': medep_angle,
                'position': position_m,
                'distance': distance_m,
                'track_pdg': 13,
                'ancestor_pdg': 13,
                'origin': 'medep'
            })
        if position_p is not None:
            candidates.append({
                'label': 'proton',
                'angle': pedep_angle,
                'position': position_p,
                'distance': distance_p,
                'track_pdg': 2212,
                'ancestor_pdg': 2212,
                'origin': 'proton'
            })
        if position_n is not None:
            candidates.append({
                'label': 'neutron',
                'angle': nedep_angle,
                'position': position_n,
                'distance': distance_n,
                'track_pdg': 2112,
                'ancestor_pdg': 2112,
                'origin': 'nedep'
            })
        if position_noises is not None:
            candidates.append({
                'label': 'noise',
                'angle': noises_angle,
                'position': position_noises,
                'distance': distance_noise,
                'track_pdg': 20110715,
                'ancestor_pdg': 20110715,
                'origin': 'noise'
            })

        first_edep = None
        first_edep_source = None

        outside_candidates = [c for c in candidates if c['angle'] > cone_angle]
        if outside_candidates:
            winner = min(outside_candidates, key=lambda c: c['distance'])
            first_edep = winner['position']
            first_edep_source = self._make_source_info(winner['label'], first_edep,
                                                       track_pdg=winner['track_pdg'],
                                                       ancestor_pdg=winner['ancestor_pdg'],
                                                       origin=winner['origin'])
        elif candidates:
            fallback_position = self.find_nearest_depo_outside(muon_energy)
            if fallback_position is not None:
                first_edep = np.array(fallback_position, dtype=float)
                first_edep_source = getattr(self, '_last_depo_outside_source', None)
                if first_edep_source is None:
                    first_edep_source = self._make_source_info('unknown', first_edep, origin='outside_cone')
            else:
                winner = min(candidates, key=lambda c: c['distance'])
                first_edep = winner['position']
                first_edep_source = self._make_source_info(winner['label'], first_edep,
                                                           track_pdg=winner['track_pdg'],
                                                           ancestor_pdg=winner['ancestor_pdg'],
                                                           origin=winner['origin'])
        else:
            fallback_position = self.find_nearest_depo_outside(muon_energy)
            if fallback_position is not None:
                first_edep = np.array(fallback_position, dtype=float)
                first_edep_source = getattr(self, '_last_depo_outside_source', None)
                if first_edep_source is None:
                    first_edep_source = self._make_source_info('unknown', first_edep, origin='outside_cone')

        reconstructed_direction = first_edep
        if reconstructed_direction is None:
            print("reconstructed direction is None")
            if self.last_reco_failure_reason is None:
                self.last_reco_failure_reason = 'no_candidate_neutron_edep'
            return None
        if first_edep_source is None:
            first_edep_source = self._make_source_info('unknown', reconstructed_direction, origin='undefined')
        self._last_first_edep_source = first_edep_source
        reconstructed_direction = np.array(reconstructed_direction, dtype=float)
        reco_norm = np.linalg.norm(reconstructed_direction)
        if reco_norm <= 0 or not np.all(np.isfinite(reconstructed_direction)):
            if self.last_reco_failure_reason is None:
                self.last_reco_failure_reason = 'invalid_reconstructed_direction'
            return None
        reconstructed_direction = reconstructed_direction / reco_norm

        cos = None
        angle_deg = None
        if real_direction is not None:
            real_norm = np.linalg.norm(real_direction)
            if real_norm > 0:
                dot = np.dot(real_direction, reconstructed_direction)
                denom = real_norm  # reconstructed_direction already normalized
                cos = dot / denom
                angle_deg = np.degrees(np.arccos(np.clip(dot / denom, -1.0, 1.0)))

        return cos, angle_deg, reconstructed_direction, first_edep
    def pairing_muon_energy_with_muon_cone(self,muon_energy):
        """input a muon energy in mev and out put its cone"""
        # print(muon_energy)x

        angle=self.cone_applying_for_muon[round(muon_energy/100)]
        return angle

    def collect_cos(self, cone_method, energy_threshold: Optional[float] = None):
        cos_list = []
        angle_deg_list = []
        n_events = 1000
        counter = 0

        for idx in range(n_events):
            self.Jump(idx)
            # if self.selectmuonneutronevent(idx)==0:
            #     continue
            if idx%100==0:
                print(f"read events {idx}/1000 ")
            muon_energy = self.read_muon_energy()
            muon_energy = muon_energy / 1000
            nearest = self.find_nearest_depo(muon_energy, cone_method, energy_threshold)
            if nearest is None:
                continue
            cos, angle_deg, reco_dir, _ = nearest
            cos_list.append(cos)
            angle_deg_list.append(angle_deg)
            counter += 1

        return cos_list,counter,angle_deg_list
    def collect_cos_EV_pureblood(self):
        angle_deg_list = []
        for idx in range(1000):
            self.Jump(idx)
            # if not self.selectelectronneutronevent(idx):
            #     continue
            if idx % 100 == 0:
                print(f"read events {idx}/1000 ")
            electron_energy = self.read_electron_energy()
            electron_energy = electron_energy / 1000
            result = self.find_nearest_depo_EV_pureblood(electron_energy)
            if result is None or result == [None,None]:
                # If reconstruction fails, assume neutron aligns with electron
                angle_deg = 0.0
            else:
                angle_deg, _ = result
                if angle_deg is None:
                    angle_deg = 0.0
            angle_deg_list.append(angle_deg)
        print(angle_deg_list)
        return angle_deg_list

    def find_nearest_depo_EV_pureblood(self,electron_energy):
        composite_neutron_ke, composite_neutron_momentum, real_direction = self.composite_neutron()
        electron_direction, _ = self.read_electron_direction()
        if real_direction is not None:
            real_direction = np.array(real_direction, dtype=float)
            if not np.all(np.isfinite(real_direction)) or np.linalg.norm(real_direction) == 0:
                real_direction = None
                print("we cannot find the real direction of neutron from the vertex in the simulation")
        cone_angle_list = Event.cone_applying_for_electron
        n = round(electron_energy * 10)
        n = max(0, min(n, len(cone_angle_list) - 1))
        cone_angle = cone_angle_list[n]

        Eedep_result = self.find_nearest_Eedep(cone_angle)
        pedep_result = self.find_nearest_proton_EV(cone_angle)
        nedep_result = self.find_nearest_nedep_EV(cone_angle)
        noise_result = self.find_nearest_noises(
            self.energy_threshold,
            reference_direction=electron_direction,
            cone_angle=cone_angle)

        # the readout should all be in the form of position + length

        def parse_result(result):
            if result is None:
                return None, float('inf')
            position, length = result
            position = np.array(position, dtype=float)
            if not np.all(np.isfinite(position)):
                return None, float('inf')
            distance = np.linalg.norm(position)
            if distance <= 0:
                return None, float('inf')
            return position, distance

        position_e, distance_e = parse_result(Eedep_result)
        position_p, distance_p = parse_result(pedep_result)
        position_n, distance_n = parse_result(nedep_result)
        position_noise, distance_noise = parse_result(None)
        if noise_result is not None:
            _, noise_position = noise_result
            if noise_position is not None:
                noise_position = np.array(noise_position, dtype=float)
                if np.all(np.isfinite(noise_position)):
                    noise_distance = np.linalg.norm(noise_position)
                    if noise_distance > 0:
                        position_noise = noise_position
                        distance_noise = noise_distance
        # print("position_e", position_e)
        # print("position_p", position_p)
        # print("distance_n", distance_n)
        candidates = []
        if position_e is not None:
            candidates.append({
                'label': 'electron',
                'position': position_e,
                'distance': distance_e,
                'track_pdg': -11,
                'ancestor_pdg': -11,
                'origin': 'Eedep'
            })
        if position_p is not None:
            candidates.append({
                'label': 'proton',
                'position': position_p,
                'distance': distance_p,
                'track_pdg': 2212,
                'ancestor_pdg': 2212,
                'origin': 'proton'
            })
        if position_n is not None:
            candidates.append({
                'label': 'neutron',
                'position': position_n,
                'distance': distance_n,
                'track_pdg': 2112,
                'ancestor_pdg': 2112,
                'origin': 'nedep'
            })
        if position_noise is not None:
            candidates.append({
                'label': 'noise',
                'position': position_noise,
                'distance': distance_noise,
                'track_pdg': 20110715,
                'ancestor_pdg': 20110715,
                'origin': 'noise'
            })

        first_edep = None
        first_edep_source = None
        if candidates:
            winner = min(candidates, key=lambda c: c['distance'])
            first_edep = winner['position']
        else:
            print('no nearest depo outside cone candidates')
            return None
        reconstructed_direction = first_edep / np.linalg.norm(first_edep)
        # print("reconstructed direction is :", reconstructed_direction)
        # print("real direction is :", real_direction)
        if real_direction is None:
            return None, reconstructed_direction
        cos = np.dot(real_direction, reconstructed_direction) / (
                    np.linalg.norm(real_direction) * np.linalg.norm(reconstructed_direction))
        if winner['ancestor_pdg'] ==2112:
            angle_deg = np.degrees(np.arccos(
                np.clip(np.dot(real_direction, reconstructed_direction) / (
                            np.linalg.norm(real_direction) * np.linalg.norm(reconstructed_direction)), -1.0, 1.0)))

        # consider one question, when cos of two direction co1uld be zero?
        # first scenario, real_direction = 0
        # print("cos:", cos)

            return angle_deg, reconstructed_direction
        else:
            return None,None


    def collect_cos_EV(self):
        """this is the collect cos function for the electron neutrino to see
        how the electron cone improve the neutron direction reconstruction
         in electron neutrino events"""
        angle_deg_list = []
        for idx in range(1000):
            self.Jump(idx)
            # if not self.selectelectronneutronevent(idx):
            #     continue
            if idx % 100 == 0:
                    print(f"read events {idx}/1000 ")
            electron_energy = self.read_electron_energy()
            electron_energy = electron_energy / 1000
            result = self.find_nearest_depo_EV(electron_energy)
            if result is None:
                    # If reconstruction fails, assume neutron aligns with electron
                    angle_deg = 0.0
            else:
                    angle_deg,_ = result
                    if angle_deg is None:
                            angle_deg = 0.0
            angle_deg_list.append(angle_deg)
        print(angle_deg_list)
        return angle_deg_list

    def selectelectronneutronevent(self, idx):
        self.Jump(idx)
        electron_count = 0
        neutron_count = 0
        other_count = 0
        gamma_count = 0
        for particle in self.vertex.Particles:
            pdg = particle.GetPDGCode()
            if pdg == 2112:
                neutron_count += 1
            elif abs(pdg) == 11:
                electron_count += 1
            elif pdg == 22:
                gamma_count += 1
            else:
                other_count += 1

        if electron_count != 1:
            return 0
        if neutron_count < 1:
            return 0
        if other_count != 0:
            return 0
        return 1




    def select_neutron_event_in_antinue(self,idx):
        self.Jump(idx)
        neutron_count = 0
        for particle in self.vertex.Particles:
            pdg = particle.GetPDGCode()
            if pdg == 2112:
                neutron_count += 1
        if neutron_count < 1:
            return 0
        return 1


    def has_michel_electron(self):
        trig = 0
        for i in range(self.tracks.size):
            track = self.tracks[i]
            pdg = track.GetPDGCode()  # Get the PDG code of the particle
            parId = track.GetParentId()
            track_parent = self.tracks[parId]
            pdg_parent = track_parent.GetPDGCode()
            mom = track.GetInitialMomentum()
            mass = mom.M()
            KE_max = mom.E() - mass

            # here we demand the michel electron shall be a positron
            # in antinumu event since its parent shall be a mu+

            if pdg == 11 and pdg_parent in (13, -13) and KE_max > 20:
                depo_indices = track.association.get('depoList', [])
                if len(depo_indices) < 2:
                    # Need at least two depositions to define the Michel start/end
                    continue
                di_start = depo_indices[1]
                di_end = depo_indices[-1]
                if not (0 <= di_start < len(self.depos)):
                    continue
                if not (0 <= di_end < len(self.depos)):
                    continue
                michel_electron_id = i
                trig += 1
        if trig == 1:
            return True
        else:
            return False


    def identify_neutron_scattering_type(self):
        vertex_neutron = []
        for track in self.tracks:
            id = track.GetPDGCode()
            if id ==2112:
                track_parent_id = track.GetParentId()
                if track_parent_id == -1:
                    trackid = track.GetTrackId()
                    vertex_neutron.append(trackid)
            # print("the neutron tracks are:",vertex_neutron)
        sub_process = None
        for track in self.tracks:
            trackid = track.GetTrackId()
            if trackid in vertex_neutron:
                # track.Points[0] is the creation point, so we need the first
                # subsequent entry to learn how the neutron interacts.
                if len(track.Points) < 2:
                    print('neutron track has no recorded scatter beyond creation')
                    continue
                point = track.Points[1]
                process = point.GetProcess()
                sub_process = point.GetSubprocess()
                break
        return sub_process




    def count_first_edep_from_proton(self,cone_method, n_events: int = 1000):
        """Return how often the first candidate edep originates from a proton."""

        max_events = min(n_events, self.nEntry)
        source_counter = Counter()
        failure_counter = Counter()
        total_considered = 0

        for idx in range(max_events):
            self.Jump(idx)
            muon_energy = self.read_muon_energy()
            if muon_energy is None:
                failure_counter['no_muon_energy'] += 1
                continue
            if muon_energy <= 0:
                failure_counter['nonpositive_muon_energy'] += 1
                continue
            muon_energy = muon_energy / 1000.0

            self.last_reco_failure_reason = None
            result = self.find_nearest_depo(muon_energy,cone_method)
            if result is None:
                reason = self.last_reco_failure_reason or 'unknown_failure'
                failure_counter[reason] += 1
                continue

            total_considered += 1
            source_info = getattr(self, '_last_first_edep_source', None)
            label = 'unknown'
            if source_info and source_info.get('label'):
                label = source_info['label']
            source_counter[label] += 1

        return {
            'total_events': max_events,
            'total_considered': total_considered,
            'proton_first_edep': source_counter.get('proton', 0),
            'source_counts': dict(source_counter),
            'failure_counts': dict(failure_counter)
        }


    def find_energy_of_earliest_neutron_energy_deposition(self):
        earliest_timing = 1000000
        mm2cm = 0.1
        earliest_energy = 0
        if self.SelectNeutronEvent()!=1:
            return None
        for track in self.tracks:
            depoList = track.association['depoList']
            for di in depoList:
                depo = self.depos[di]
                edep = depo.GetEnergyDeposit()
                start, stop = depo.GetStart(), depo.GetStop()
                time = (start.T()+stop.T())/2 *mm2cm
                if time <= earliest_timing:
                    earliest_energy = edep
                    earliest_timing = time
        return earliest_energy

    def find_distance_of_earliest_neutron_energy_deposition(
        self,
        energy_threshold: Optional[float] = None
    ):
        """Return the vertex distance of the earliest neutron deposition above threshold."""

        if self.SelectNeutronEvent() != 1:
            return None

        threshold = self._resolve_energy_threshold(energy_threshold)
        mm2m = 0.001
        earliest_time = float('inf')
        earliest_distance = None

        for track in self.tracks:
            ancestor_id = track.association.get('ancestor')
            ancestor_pdg = None
            if ancestor_id is not None and 0 <= ancestor_id < self.tracks.size:
                ancestor_pdg = self.tracks[ancestor_id].GetPDGCode()

            track_pdg = track.GetPDGCode()
            depo_list = track.association.get('depoList', [])

            for di in depo_list:
                if di < 0 or di >= len(self.depos):
                    continue

                if ancestor_pdg != 2112 and track_pdg != 2112:
                    continue

                depo = self.depos[di]
                edep = depo.GetEnergyDeposit()
                if edep < threshold:
                    continue

                start, stop = depo.GetStart(), depo.GetStop()
                time = (start.T() + stop.T()) / 2.0
                if time > earliest_time:
                    continue

                point = np.array([
                    (start.X() + stop.X()) / 2 * mm2m,
                    (start.Y() + stop.Y()) / 2 * mm2m,
                    (start.Z() + stop.Z()) / 2 * mm2m
                ], dtype=float)

                earliest_time = time
                earliest_distance = float(np.linalg.norm(point))

        return earliest_distance

    def collect_low_energy_proton_deposition_distances(
        self,
        energy_threshold: Optional[float] = None
    ):
        """Return all deposition distances from low-energy (<=70 MeV) protons."""

        threshold = self._resolve_energy_threshold(energy_threshold)
        mm2m = 0.001
        distances = []

        for particle in self.vertex.Particles:
            if particle.GetPDGCode() != 2212:
                continue

            mom = particle.GetMomentum()
            kinetic_energy = mom.E() - mom.M()
            if kinetic_energy > self.proton_trackable_length:
                continue

            track_id = particle.GetTrackId() if hasattr(particle, 'GetTrackId') else None
            if track_id is None or track_id < 0 or track_id >= self.tracks.size:
                continue

            track = self.tracks[track_id]
            depo_list = track.association.get('depoList', [])
            for di in depo_list:
                if di < 0 or di >= len(self.depos):
                    continue

                depo = self.depos[di]
                edep = depo.GetEnergyDeposit()
                if edep < threshold:
                    continue

                start, stop = depo.GetStart(), depo.GetStop()
                point = np.array([
                    (start.X() + stop.X()) / 2 * mm2m,
                    (start.Y() + stop.Y()) / 2 * mm2m,
                    (start.Z() + stop.Z()) / 2 * mm2m,
                ], dtype=float)

                distance_m = float(np.linalg.norm(point))
                distances.append(distance_m)

        return distances






    def find_nearest_depo_outside(self,muon_energy):
        """this function is used when both the n-edep and m-edep falls into the cone, in this case we need to find the
        nearest one outside the cone"""
        results= self.read_muon_direction_and_amplitude()
        if results is not None:
            muon_direction,mom =results
        else:
            return
        mm2cm = 0.1
        threshold = self.energy_threshold
        cone_angle_list = self.cone_applying_for_muon #      # for 0.1,0.2,0.3,0.4,0.5,0.6,0.7 GeV muon57 43 41 30 31 32 18
        n= round(muon_energy * 10)
        n = max(0, min(n, len(cone_angle_list) - 1))
        cone_angle = cone_angle_list[n]
        shortest_distance = 100000000
        nearest_position = None
        nearest_source = None
        for track in self.tracks:
            depoList = track.association['depoList']
            for di in depoList:
                depo = self.depos[di]
                edep = depo.GetEnergyDeposit()
                if edep >= threshold:
                    start, stop = depo.GetStart(), depo.GetStop()
                    midpoint = np.array([
                        (start.X() + stop.X()) / 2 * mm2cm,
                        (start.Y() + stop.Y()) / 2 * mm2cm,
                        (start.Z() + stop.Z()) / 2 * mm2cm
                    ], dtype=float)
                    distance = np.linalg.norm(np.array(midpoint))
                    midpoint = midpoint / distance

                    angle = np.degrees(np.arccos(np.clip(np.dot(midpoint, muon_direction), -1.0, 1.0)))
                    if cone_angle < angle and distance < shortest_distance:
                        nearest_position = midpoint
                        shortest_distance = distance
                        nearest_source = self._build_depo_source_info(track, midpoint, origin='outside_cone')
        self._last_depo_outside_source = nearest_source
        if nearest_position is None or abs(np.linalg.norm(nearest_position)) < 0.0001:
            print("here we are not able to find the defined nearest and who ever is the nearest is right on the face")
            if self.last_reco_failure_reason is None:
                self.last_reco_failure_reason = 'outside_cone_search_failed'
        return nearest_position

    def _build_depo_source_info(self, track, position, origin='unknown'):
        track_pdg = track.GetPDGCode() if hasattr(track, 'GetPDGCode') else None
        ancestor_pdg = None
        track_id = track.GetTrackId() if hasattr(track, 'GetTrackId') else None

        ancestor_id = None
        if hasattr(track, 'association'):
            ancestor_id = track.association.get('ancestor')
        if ancestor_id is not None and ancestor_id >= 0 and ancestor_id < self.tracks.size:
            ancestor_track = self.tracks[ancestor_id]
            ancestor_pdg = ancestor_track.GetPDGCode()

        label = self._classify_pdg(track_pdg)
        if label == 'other' and ancestor_pdg is not None:
            ancestor_label = self._classify_pdg(ancestor_pdg)
            if ancestor_label != 'other':
                label = ancestor_label

        return self._make_source_info(label, position, track_pdg=track_pdg,
                                      ancestor_pdg=ancestor_pdg, origin=origin, track_id=track_id)

    def _make_source_info(self, label, position, track_pdg=None, ancestor_pdg=None, origin=None, track_id=None):
        pos_array = None
        if position is not None:
            pos_array = np.array(position, dtype=float)
        return {
            'label': label,
            'track_pdg': track_pdg,
            'ancestor_pdg': ancestor_pdg,
            'origin': origin,
            'position': pos_array,
            'track_id': track_id
        }

    @staticmethod
    def _classify_pdg(pdg_code):
        if pdg_code is None:
            return 'unknown'
        if pdg_code == 2212:
            return 'proton'
        if pdg_code == 2112:
            return 'neutron'
        if abs(pdg_code) == 13:
            return 'muon'
        if abs(pdg_code) == 11:
            return 'electron'
        if pdg_code == 22:
            return 'gamma'
        if abs(pdg_code) in (211, 111):
            return 'pion'
        return 'other'




    def neutron_edep_chromo(self):
        NeutronTrack = []
        mm2cm = 0.1
        threshold = self.energy_threshold

        for track in self.tracks:
            ancestor = track.association['ancestor']
            ancestor_id = self.tracks[ancestor].GetPDGCode()
            if ancestor_id==2112:
                midpoint = [0, 0, 0,0]
                for di in track.association['depoList']:
                    depo = self.depos[di]
                    edep = depo.GetEnergyDeposit()
                    if edep >= threshold:
                        start, stop = depo.GetStart(), depo.GetStop()
                        t = (start.T()+stop.T())/2
                        midpoint = np.array([
                            (start.X() + stop.X()) / 2 * mm2cm,
                            (start.Y() + stop.Y()) / 2 * mm2cm,
                            (start.Z() + stop.Z()) / 2 * mm2cm,
                            t
                        ])
                    if not np.array_equal(midpoint, [0, 0, 0,0]):
                        NeutronTrack.append(midpoint)

        if len(NeutronTrack) == 0:
            return NeutronTrack
        else:
            sorted_neutron_track = sorted(NeutronTrack, key=lambda point: point[3])





    def neutrino_direction_reconstruction(self, n_event, cone_method,energy_threshold: Optional[float] = None,
                                          ):
        """the reconstruction of neutrino shall be decided upon the
        neutrino energy, as it is in the current DUNE collaboration data,
        we input the energy of neutrino, use intrictic function to figure out
        the momentum and amplitude of muon, and then we use the function
        reconstruct_neutron_momentum to read the amplitude and the direction
        of neutron"""

        self.last_reco_failure_reason = None
        # 1.muon read out and neutron reconstruction
        result= self.get_momentum_direction_and_amplitude(n_event)
        if result is  None:
            # print("no neutron momentum found")
            if self.last_reco_failure_reason is None:
                self.last_reco_failure_reason = 'selectmuonneutronevent_failed'
            self.neutrino_direction_reconstruction_result_with_noise = None
            return None
        muon_mom,neutron_mom = result
        muon_mom = np.array(muon_mom)
        muon_direction = muon_mom/np.linalg.norm(muon_mom)
        neutron_direction, neutron_amplitude = self.neutron_direction_and_amplitude_reconstruction(cone_method,energy_threshold)
        # print(f'neutron momentum is {neutron_direction*neutron_amplitude} and muon momentum is {muon_mom}')
        if neutron_direction is None or neutron_amplitude is None:
            print("The neutron direction and/or amplitude were not found. Cannot reconstruct neutrino.")
            if self.last_reco_failure_reason is None:
                self.last_reco_failure_reason = 'neutron_direction_missing'
            self.neutrino_direction_reconstruction_result_with_noise = None
            return
        else:
            neutron_direction,neutron_magnitude = neutron_direction,neutron_amplitude

        # 2.proton read out
        pion_list  = self.pion_threshold()
        proton_list = self.proton_threshold()
        total_proton_momentum = np.zeros(3)
        total_pion_momentum = np.zeros(3)
        if  proton_list:
            for is_track, proton_id, momentum, position, position_time in proton_list:
                if is_track:
                    continue
                total_proton_momentum += momentum

        if pion_list:
            for is_track, pion_id, momentum in pion_list:
                if (is_track):
                    continue
                total_pion_momentum += momentum
        neutrino_direction = np.array(muon_mom + neutron_direction*neutron_amplitude+total_proton_momentum+total_pion_momentum)
        neutrino_direction = neutrino_direction / np.linalg.norm(neutrino_direction)
        angle_deviation = np.degrees(np.arccos(neutrino_direction[2]))
        angle_deviation_lepton = np.degrees(np.arccos(muon_direction[2]))

        self.last_reco_failure_reason = 'success'
        result_tuple = (angle_deviation, angle_deviation_lepton)
        self.neutrino_direction_reconstruction_result_with_noise = result_tuple
        return result_tuple

    def find_nearest_depo_without_noise(self, muon_energy, cone_method, energy_threshold):

            self._last_first_edep_source = None
            threshold = self._resolve_energy_threshold(energy_threshold)
            real_direction, mom = self.read_neutron_momentum()
            if real_direction is not None:
                real_direction = np.array(real_direction, dtype=float)
                if not np.all(np.isfinite(real_direction)) or np.linalg.norm(real_direction) == 0:
                    real_direction = None
            cone_angle_list = cone_method
            n = round(muon_energy * 10)
            n = max(0, min(n, len(cone_angle_list) - 1))
            cone_angle = cone_angle_list[n]

            medep_result = self.find_nearest_discontinuous_medep(threshold)
            pedep_result = self.find_nearest_discontinuous_proton(threshold)
            nedep_result = self.find_nearest_discontinuous_nedep(threshold)

            def parse_result(result):
                if result is None:
                    return 180.0, None, float('inf')
                angle, position = result
                position = np.array(position, dtype=float)
                if not np.all(np.isfinite(position)):
                    return 180.0, None, float('inf')
                distance = np.linalg.norm(position)
                if distance <= 0:
                    return angle, None, float('inf')
                return angle, position, distance

            medep_angle, position_m, distance_m = parse_result(medep_result)
            pedep_angle, position_p, distance_p = parse_result(pedep_result)
            nedep_angle, position_n, distance_n = parse_result(nedep_result)

            candidates = []
            if position_m is not None:
                candidates.append({
                    'label': 'muon',
                    'angle': medep_angle,
                    'position': position_m,
                    'distance': distance_m,
                    'track_pdg': 13,
                    'ancestor_pdg': 13,
                    'origin': 'medep'
                })
            if position_p is not None:
                candidates.append({
                    'label': 'proton',
                    'angle': pedep_angle,
                    'position': position_p,
                    'distance': distance_p,
                    'track_pdg': 2212,
                    'ancestor_pdg': 2212,
                    'origin': 'proton'
                })
            if position_n is not None:
                candidates.append({
                    'label': 'neutron',
                    'angle': nedep_angle,
                    'position': position_n,
                    'distance': distance_n,
                    'track_pdg': 2112,
                    'ancestor_pdg': 2112,
                    'origin': 'nedep'
                })
            first_edep = None
            first_edep_source = None

            outside_candidates = [c for c in candidates if c['angle'] > cone_angle]
            if outside_candidates:
                winner = min(outside_candidates, key=lambda c: c['distance'])
                first_edep = winner['position']
                first_edep_source = self._make_source_info(winner['label'], first_edep,
                                                           track_pdg=winner['track_pdg'],
                                                           ancestor_pdg=winner['ancestor_pdg'],
                                                           origin=winner['origin'])
            elif candidates:
                fallback_position = self.find_nearest_depo_outside(muon_energy)
                if fallback_position is not None:
                    first_edep = np.array(fallback_position, dtype=float)
                    first_edep_source = getattr(self, '_last_depo_outside_source', None)
                    if first_edep_source is None:
                        first_edep_source = self._make_source_info('unknown', first_edep, origin='outside_cone')
                else:
                    winner = min(candidates, key=lambda c: c['distance'])
                    first_edep = winner['position']
                    first_edep_source = self._make_source_info(winner['label'], first_edep,
                                                               track_pdg=winner['track_pdg'],
                                                               ancestor_pdg=winner['ancestor_pdg'],
                                                               origin=winner['origin'])
            else:
                fallback_position = self.find_nearest_depo_outside(muon_energy)
                if fallback_position is not None:
                    first_edep = np.array(fallback_position, dtype=float)
                    first_edep_source = getattr(self, '_last_depo_outside_source', None)
                    if first_edep_source is None:
                        first_edep_source = self._make_source_info('unknown', first_edep, origin='outside_cone')

            reconstructed_direction = first_edep
            if reconstructed_direction is None:
                print("reconstructed direction is None")
                if self.last_reco_failure_reason is None:
                    self.last_reco_failure_reason = 'no_candidate_neutron_edep'
                return None
            if first_edep_source is None:
                first_edep_source = self._make_source_info('unknown', reconstructed_direction, origin='undefined')
            self._last_first_edep_source = first_edep_source
            reconstructed_direction = np.array(reconstructed_direction, dtype=float)
            reco_norm = np.linalg.norm(reconstructed_direction)
            if reco_norm <= 0 or not np.all(np.isfinite(reconstructed_direction)):
                if self.last_reco_failure_reason is None:
                    self.last_reco_failure_reason = 'invalid_reconstructed_direction'
                return None
            reconstructed_direction = reconstructed_direction / reco_norm

            cos = None
            angle_deg = None
            if real_direction is not None:
                real_norm = np.linalg.norm(real_direction)
                if real_norm > 0:
                    dot = np.dot(real_direction, reconstructed_direction)
                    denom = real_norm
                    cos = dot / denom
                    angle_deg = np.degrees(np.arccos(np.clip(dot / denom, -1.0, 1.0)))

            return cos, angle_deg, reconstructed_direction, first_edep
    def neutron_direction_and_amplitude_reconstruction_without_noise(self,cone_method,energy_threshold):


        threshold = self._resolve_energy_threshold(energy_threshold)
        muon_energy = self.read_muon_energy()
        muon_energy = muon_energy / 1000
        nearest = self.find_nearest_depo_without_noise(muon_energy, cone_method, threshold)
        if nearest is None:
            print("invalidate neutron direction information")
            if self.last_reco_failure_reason is None:
                self.last_reco_failure_reason = 'find_nearest_depo_failed'
            self.neutron_direction_reconstruction_result_without_noise = None
            return None, None
        cos, angle_deg, reco_dir, _ = nearest
        reconstructed_neutron_direction = reco_dir
        reconstructed_neutron_magnitude = 0
        energy_sum = 0
        for i in range(self.tracks.size):
            track = self.tracks[i]
            ancestor = track.association['ancestor']
            pdg_ancestor = self.tracks[ancestor].GetPDGCode()
            if pdg_ancestor == 2112:
                depoList = track.association['depoList']
                for di in depoList:
                    depo = self.depos[di]
                    edep = depo.GetEnergyDeposit()
                    if edep >= threshold:
                        energy_sum += edep
        reconstructed_neutron_energy = self.pairing_summed_energy_with_initial_KE(self.neutrino_type, energy_sum, self.inspection_length_energy_side)
        reconstructed_neutron_magnitude = np.sqrt(
            reconstructed_neutron_energy ** 2 + 2 * reconstructed_neutron_energy * 939.58)

        result = (reconstructed_neutron_direction, reconstructed_neutron_magnitude)
        self.neutron_direction_reconstruction_result_without_noise = result
        return result


    def reconstruct_neutron_only(self,idx):
        cone_method =self.cone_applying_for_muon
        energy_threshold = self._resolve_energy_threshold(idx)
        neutron_direction, neutron_amplitude = self.neutron_direction_and_amplitude_reconstruction_without_noise(
            cone_method,
            energy_threshold)





    def pairing_summed_energy_with_initial_KE(self,neutrino_type, summed_neutron_energy, length):

        path = Path('~/frunshiworking/energy_reconstruction_multiplier').expanduser()
        file_path = path/ f'neutron_peak_curve_points_{neutrino_type}_{length}m.txt'
        data = np.loadtxt(file_path, comments = '#')

        initial_KE = data[:,0]
        summed_energy = data[:,1]

        order = np.argsort(summed_energy)
        summed_energy = summed_energy[order]
        initial_KE = initial_KE[order]

        method = 'interp'
        if summed_neutron_energy > summed_energy[-1]:
            return 200
        else:
            return np.interp(summed_neutron_energy, summed_energy, initial_KE)

    def neutrino_direction_reconstruction_without_noises_antinumu(self,n_event,cone_method,energy_threshold: Optional[float] = None):
            self.last_reco_failure_reason = None
            # 1.muon read out and neutron reconstruction
            result = self.get_momentum_direction_and_amplitude(n_event)
            if result is None:
                # print("no neutron momentum found")
                if self.last_reco_failure_reason is None:
                    self.last_reco_failure_reason = 'selectmuonneutronevent_failed'
                self.neutrino_direction_reconstruction_result_without_noise = None
                return None
            muon_mom, neutron_mom = result
            muon_mom = np.array(muon_mom)
            muon_direction = muon_mom / np.linalg.norm(muon_mom)
            neutron_direction, neutron_amplitude = self.neutron_direction_and_amplitude_reconstruction_without_noise(cone_method,
                                                                                                       energy_threshold)
            # print(f'neutron momentum is {neutron_direction*neutron_amplitude} and muon momentum is {muon_mom}')
            if neutron_direction is None or neutron_amplitude is None:
                print("The neutron direction and/or amplitude were not found. Cannot reconstruct neutrino.")
                if self.last_reco_failure_reason is None:
                    self.last_reco_failure_reason = 'neutron_direction_missing'
                self.neutrino_direction_reconstruction_result_without_noise = None
                return
            else:
                neutron_direction, neutron_magnitude = neutron_direction, neutron_amplitude

            # 2.proton read out
            pion_list = self.pion_threshold()
            proton_list = self.proton_threshold()
            total_proton_momentum = np.zeros(3)
            total_pion_momentum = np.zeros(3)
            if proton_list:
                for is_track, proton_id, momentum, position, position_time in proton_list:
                    if is_track:
                        continue
                    total_proton_momentum += momentum

            if pion_list:
                for is_track, pion_id, momentum in pion_list:
                    if (is_track):
                        continue
                    total_pion_momentum += momentum
            neutrino_direction = np.array(
                muon_mom + neutron_direction * neutron_amplitude + total_proton_momentum + total_pion_momentum)
            neutrino_direction = neutrino_direction / np.linalg.norm(neutrino_direction)
            angle_deviation = np.degrees(np.arccos(neutrino_direction[2]))
            angle_deviation_lepton = np.degrees(np.arccos(muon_direction[2]))

            self.last_reco_failure_reason = 'success'
            result_tuple = (angle_deviation, angle_deviation_lepton)
            self.neutrino_direction_reconstruction_result_without_noise = result_tuple
            return result_tuple



    def neutrino_direction_reconstruction_without_neutron(self, n_event):


        self.last_reco_failure_reason = None

            # 1.muon read out and neutron reconstruction
        result = self.get_momentum_direction_and_amplitude(n_event)
        if result is None:
                # print("no neutron momentum found")
            if self.last_reco_failure_reason is None:
                self.last_reco_failure_reason = 'selectmuonneutronevent_failed'
            return None
        muon_mom, neutron_mom = result
        muon_mom = np.array(muon_mom)
        muon_direction = muon_mom / np.linalg.norm(muon_mom)
        neutron_direction, neutron_amplitude = self.neutron_direction_and_amplitude_reconstruction(
            cone_method=self.cone_applying_for_muon ,
            energy_threshold=self.energy_threshold,
        )
            # print(f'neutron momentum is {neutron_direction*neutron_amplitude} and muon momentum is {muon_mom}')
        if neutron_direction is None or neutron_amplitude is None:
            print("The neutron direction and/or amplitude were not found. Cannot reconstruct neutrino.")
            if self.last_reco_failure_reason is None:
                    self.last_reco_failure_reason = 'neutron_direction_missing'
            return
        else:
            neutron_direction, neutron_magnitude = neutron_direction, neutron_amplitude

            # 2.proton read out
        proton_list = self.proton_threshold()
        pion_list = self.pion_threshold()
        total_pion_momentum = np.zeros(3)
        total_proton_momentum = np.zeros(3)
        if proton_list:
            for is_track, proton_id, momentum, position, position_time in proton_list:
                    if is_track:
                        continue
                    total_proton_momentum += momentum
        if pion_list:
            for is_track, pion_id, momentum in pion_list:
                if is_track:
                    continue
                    total_pion_momentum += momentum
        neutrino_direction = np.array(muon_mom + total_proton_momentum+total_pion_momentum)
        neutrino_direction = neutrino_direction/np.linalg.norm(neutrino_direction)
        angle_deviation = np.degrees(np.arccos(neutrino_direction[2]))
        angle_deviation_lepton = np.degrees(np.arccos(muon_direction[2]))

        self.last_reco_failure_reason = 'success'
        return angle_deviation, angle_deviation_lepton

    def simulation_readout(self):
        """Return unit vector of true neutrino momentum from primary daughter sums.

        The primary muon, neutron, proton, and pion momenta are summed directly
        without any energy thresholds before normalizing."""

        muon_total = np.zeros(3)
        neutron_total = np.zeros(3)
        proton_total = np.zeros(3)
        pion_total = np.zeros(3)

        for particle in self.vertex.Particles:
            mom = particle.GetMomentum()
            if not mom:
                continue
            momentum = np.array([mom.X(), mom.Y(), mom.Z()])
            pdg = particle.GetPDGCode()
            if pdg in (-13, 13):
                muon_total += momentum
            elif pdg == 2112:
                neutron_total += momentum
            elif pdg == 2212:
                proton_total += momentum
            elif abs(pdg) in (211, 111):
                pion_total += momentum

        neutrino_mom = muon_total + neutron_total + proton_total + pion_total
        norm = np.linalg.norm(neutrino_mom)
        if norm == 0:
            return neutrino_mom
        return neutrino_mom / norm


    def neutron_direction_and_amplitude_reconstruction_EV(self, energy_threshold: Optional[float] = None,):


        threshold = self._resolve_energy_threshold(energy_threshold)
        electron_energy = self.read_electron_energy()
        depo_result = self.find_nearest_depo_EV(electron_energy)
        if depo_result is None:
            reconstructed_neutron_direction = np.zeros(3)
        else:
            _, reconstructed_neutron_direction = depo_result
        neutron_reconstructed_energy,_ = self.sum_neutron_energy_deposition_EV(threshold)
        neutron_reconstructed_energy = self.pairing_summed_energy_with_initial_KE(self.neutrino_type,neutron_reconstructed_energy,self.inspection_length_energy_side)
        reconstructed_neutron_magnitude = np.sqrt(
            neutron_reconstructed_energy ** 2 + 2 * neutron_reconstructed_energy * 939.58)
        reconstructed_neutron_momentum = reconstructed_neutron_direction*reconstructed_neutron_magnitude
        return reconstructed_neutron_momentum





    def neutron_direction_and_amplitude_reconstruction(self,cone_method,energy_threshold: Optional[float] = None):

        """this function is a sum up of neutron direction and amplitude reconstruction"""
        threshold = self._resolve_energy_threshold(energy_threshold)
        muon_energy = self.read_muon_energy()
        muon_energy = muon_energy / 1000
        nearest = self.find_nearest_depo(muon_energy,cone_method, threshold)
        if nearest is None:
            print("invalidate neutron direction information")
            if self.last_reco_failure_reason is None:
                self.last_reco_failure_reason = 'find_nearest_depo_failed'
            self.neutron_direction_reconstruction_result_with_noise = None
            return None, None
        cos, angle_deg, reco_dir,_ = nearest
        reconstructed_neutron_direction= reco_dir
        reconstructed_neutron_magnitude = 0
        energy_sum = 0
        for i in range(self.tracks.size):
            track = self.tracks[i]
            ancestor = track.association['ancestor']
            pdg_ancestor = self.tracks[ancestor].GetPDGCode()
            if pdg_ancestor == 2112:
                depoList = track.association['depoList']
                for di in depoList:
                    depo = self.depos[di]
                    edep = depo.GetEnergyDeposit()
                    if edep >= threshold:
                        energy_sum += edep
        reconstructed_neutron_energy = self.pairing_summed_energy_with_initial_KE(self.neutrino_type, energy_sum,self.inspection_length_energy_side)
        reconstructed_neutron_magnitude = np.sqrt(reconstructed_neutron_energy**2+2*reconstructed_neutron_energy*939.58)

        result = (reconstructed_neutron_direction, reconstructed_neutron_magnitude)
        self.neutron_direction_reconstruction_result_with_noise = result
        return result

    def find_Energy_and_momentum(self):
        print("--- Starting find_Energy_and_momentum ---")

        momentum_magnitude_list = []
        energy_list = []
        edep_list = []

        n_events_to_process = 9000  # Or self.event.nEntry to run on all events
        counter = 0
        threshold = self.energy_threshold
        for idx in range(n_events_to_process):
            self.Jump(idx)
            if idx%100==0:
                counter += 1
                print(counter*100,"/9000 events processed")

            # Check if this event has only one neutron and one muon
            # Note: Your original check_N_M_only_event() might need similar defensive checks inside it.
            if self.check_N_M_only_event() != 1:
                continue

            # print(f"Event {idx} is a valid candidate. Analyzing particles...")

            # --- Defensive Check 1: Is the primary vertex valid? ---
            if not self.vertex:
                # print(f"Warning: Event {idx} has an invalid primary vertex. Skipping.")
                continue
            energy_sum = 0

            found_neutron = False
            for i in range(self.tracks.size):
                track = self.tracks[i]
                ancestor = track.association['ancestor']
                pdg_ancestor = self.tracks[ancestor].GetPDGCode()
                if pdg_ancestor==2112:
                    depoList = track.association['depoList']
                    for di in depoList:
                        depo = self.depos[di]
                        edep = depo.GetEnergyDeposit()
                        if edep >= threshold:
                            energy_sum += edep
            edep_list.append(energy_sum)

            for i, particle in enumerate(self.vertex.Particles):

                # --- Defensive Check 2: Is the particle object itself valid? ---
                # This is the most likely source of your crash.
                if not particle:
                    # print(f"Warning: Event {idx} has an invalid particle object at index {i}. Skipping.")
                    continue

                # --- Now it is safe to use the particle object ---
                pdg = particle.GetPDGCode()

                if pdg == 2112:  # PDG code for a neutron
                    if found_neutron:
                        # This handles cases where check_N_M_only_event might be imperfect
                        # print(f"Warning: Found a second neutron in event {idx}. Skipping particle to ensure 1n1mu.")
                        continue

                    # print(f"  > Found a neutron (PDG 2112) in event {idx}")
                    found_neutron = True

                    # --- Defensive Check 3: Does this particle have momentum data? ---
                    mom = particle.GetMomentum()
                    if not mom:
                        # print(f"Warning: Neutron in event {idx} has no momentum data. Skipping.")
                        continue

                    # --- All checks passed, now it is safe to calculate and store data ---
                    # .P() gets the magnitude of the 3-momentum vector
                    # .E() gets the total energy
                    # .M() gets the invariant mass
                    momentum_magnitude = mom.P()
                    kinetic_energy = mom.E() - mom.M()

                    momentum_magnitude_list.append(momentum_magnitude)
                    energy_list.append(kinetic_energy)
                    # print(
                        # f"    ... successfully collected E={kinetic_energy:.4f} GeV, |p|={momentum_magnitude:.4f} GeV/c")

        print("--- Finished find_Energy_and_momentum ---")
        return momentum_magnitude_list, energy_list,edep_list

    def find_nearest_discontinuous_nedep_EV(self, energy_threshold: Optional[float] = None):
        NeutronTrack = []
        counter = 0
        threshold = self._resolve_energy_threshold(energy_threshold)
        for track in self.tracks:
            ancestor = track.association['ancestor']
            ancestor_id = self.tracks[ancestor].GetPDGCode()
            mm2cm = 0.1
            # trkid = track.GetParentId()
            # print(f"this is track{trkid}, pdg:{pdg}, parentID:{parentID}")
            if ancestor_id == 2112:
                midpoint = [0, 0, 0]
                # print(f"here is a neutron from first vertex")
                for di in track.association['depoList']:
                    depo = self.depos[di]
                    edep = depo.GetEnergyDeposit()
                    # print(f"edep is {edep}")
                    if edep >= threshold:
                        # print(f"here is a neutron exceed the energy threshold")
                        # midpoint in cm
                        start, stop = depo.GetStart(), depo.GetStop()
                        midpoint = np.array([
                            (start.X() + stop.X()) / 2 * mm2cm,
                            (start.Y() + stop.Y()) / 2 * mm2cm,
                            (start.Z() + stop.Z()) / 2 * mm2cm
                        ])
                        counter += 1


                        # print("there is a discontinuous neutron energy deposition")
                if not np.array_equal(midpoint, [0, 0, 0]):
                    NeutronTrack.append(midpoint)

            proton_track = self.proton_threshold()
            if proton_track:
                for is_track, _, _, proton_position, _ in proton_track:
                    if is_track or proton_position is None:
                        continue
                    proton_position = np.array(proton_position, dtype=float)
                    if not np.all(np.isfinite(proton_position)):
                        continue
                    distance = np.linalg.norm(proton_position)
                    if distance < self.minimum_length:
                        continue
                    if distance <= 0:
                        continue
                    NeutronTrack.append(proton_position)
        # if len(NeutronTrack)== 0:
        #     print("no neutron edep detected")
        # print(f"the neutron track is {NeutronTrack}")
        nearest_distance = 10000
        i = 0
        nearest_position = np.array([0, 0, 0])
        # if counter == 0 and trig == 1:
        #     print("this neutron does not have energy deposition")
        for points in NeutronTrack:
            x = points[0]
            y = points[1]
            z = points[2]
            distance = np.linalg.norm(points)
            i += 1
            if nearest_distance > distance > self.minimum_length and distance <= self.inspection_length_direction_side:
                nearest_distance = distance
                nearest_position = np.array([x, y, z])
        electron_direction,_= self.read_electron_direction()
        # 1. First, check for any NaN values in the arrays
        if np.isnan(electron_direction).any():
            print("we don't have a valid muon direction (contains NaN)")
        elif np.isnan(nearest_position).any():
            # print("we don't have a valid neutron position (contains NaN)")
            frs = 1
        # 2. If no NaNs, THEN check for zero-length vectors
        else:
            norm_muon: floating = np.linalg.norm(electron_direction)
            norm_pos: floating = np.linalg.norm(nearest_position)
            if norm_pos == 0:
                # print("we don't have a valid neutron position (contains NaN)")
                return 0, [1000, 1000, 1000]

            if norm_muon > 0.0001 and norm_pos > 0.0001:
                # This block is now safe from NaNs and zeros
                angle_degrees = np.degrees(
                    np.arccos(np.clip(np.dot(electron_direction, nearest_position) / (norm_muon * norm_pos), -1.0, 1.0))
                )
        # print(f"checking angle_degrees \n")
        # print(f"muon_direction = {muon_direction},nearest position = {nearest_position}")

        angle_degrees = np.degrees(
            np.arccos(
                np.clip(np.dot(electron_direction, nearest_position) / (np.linalg.norm(electron_direction) * norm_pos), -1.0,
                        1.0)))

        return angle_degrees, nearest_position

    def find_nearest_discontinuous_nedep(self, energy_threshold: Optional[float] = None):
        NeutronTrack = []
        counter = 0
        threshold = self._resolve_energy_threshold(energy_threshold)
        for track in self.tracks:
            ancestor = track.association['ancestor']
            ancestor_id = self.tracks[ancestor].GetPDGCode()
            mm2cm = 0.1
            # trkid = track.GetParentId()
            # print(f"this is track{trkid}, pdg:{pdg}, parentID:{parentID}")
            if ancestor_id==2112:
                midpoint = [0, 0, 0]
                # print(f"here is a neutron from first vertex")
                for di in track.association['depoList']:
                    depo = self.depos[di]
                    edep = depo.GetEnergyDeposit()
                    # print(f"edep is {edep}")
                    if edep >= threshold:
                        # print(f"here is a neutron exceed the energy threshold")
                        # midpoint in cm
                        start, stop = depo.GetStart(), depo.GetStop()
                        midpoint = np.array([
                            (start.X() + stop.X()) / 2 * mm2cm,
                            (start.Y() + stop.Y()) / 2 * mm2cm,
                            (start.Z() + stop.Z()) / 2 * mm2cm
                        ])
                        counter +=1
                        # print("there is a discontinuous neutron energy deposition")
                if not np.array_equal(midpoint, [0, 0, 0]):
                    NeutronTrack.append(midpoint)
        # if len(NeutronTrack)== 0:
        #     print("no neutron edep detected")
        # print(f"the neutron track is {NeutronTrack}")
        nearest_distance = 10000
        i = 0
        nearest_position = np.array([0,0,0])
        # if counter == 0 and trig == 1:
        #     print("this neutron does not have energy deposition")
        for points in NeutronTrack:
            x = points[0]
            y = points[1]
            z = points[2]
            distance = np.linalg.norm(points)
            i += 1
            if nearest_distance > distance > self.minimum_length and distance <= self.inspection_length_direction_side*100:
                nearest_distance = distance
                nearest_position = np.array([x,y,z])
        returned_values = self.read_muon_direction_and_amplitude()
        muon_direction, mom = returned_values
        muon_direction= np.array(muon_direction)
        # 1. First, check for any NaN values in the arrays
        if np.isnan(muon_direction).any():
            print("we don't have a valid muon direction (contains NaN)")
        elif np.isnan(nearest_position).any():
            # print("we don't have a valid neutron position (contains NaN)")
            frs=1
        # 2. If no NaNs, THEN check for zero-length vectors
        else:
            norm_muon: floating = np.linalg.norm(muon_direction)
            norm_pos : floating = np.linalg.norm(nearest_position)
            if norm_pos==0:
                # print("we don't have a valid neutron position (contains NaN)")
                return 0, [1000,1000,1000]

            if norm_muon > 0.0001 and norm_pos > 0.0001:
                # This block is now safe from NaNs and zeros
                angle_degrees = np.degrees(
                    np.arccos(np.clip(np.dot(muon_direction, nearest_position) / (norm_muon * norm_pos), -1.0, 1.0))
                )
        # print(f"checking angle_degrees \n")
        # print(f"muon_direction = {muon_direction},nearest position = {nearest_position}")

        angle_degrees = np.degrees(
            np.arccos(np.clip(np.dot(muon_direction, nearest_position) / (np.linalg.norm(muon_direction) * norm_pos), -1.0, 1.0)))

        return angle_degrees, nearest_position


    def check_N_M_only_event(self):
            trig = 0
            neutron_number = 0
            muon_number = 0
            other_particle_number = 0
            for i in range(self.tracks.size):
                track = self.tracks[i]
                pdg = track.GetPDGCode()  # Get the PDG code of the particle
                ParentId = track.GetParentId()
                if ParentId == -1:
                    if pdg==2112:
                        neutron_number += 1
                    elif abs(pdg)== 13:
                        muon_number += 1
                    else:
                        other_particle_number += 1

            if neutron_number ==1 and muon_number ==1 and other_particle_number == 0:
                trig = 1


            return trig


    def compute_angle_between(self, u: np.ndarray, v: np.ndarray) -> float:
        """
        Return arccos(u·v) in radians, or 0 if the inputs are None, zero-length,
        or the dot falls outside [-1,1].
        """
        # 1) catch missing vectors right away
        if u is None or v is None:
            return 0.0

        # 2) coerce to arrays (in case someone passed a list)
        u = np.asarray(u)
        v = np.asarray(v)

        # 3) guard zero-length
        nu = np.linalg.norm(u)
        nv = np.linalg.norm(v)
        if nu == 0 or nv == 0:
            return 0.0

        # 4) safe dot and clamp
        dot = np.dot(u, v) / (nu * nv)
        dot = max(min(dot, 1.0), -1.0)

        # 5) final angle
        return np.arccos(dot)
    #here we write another angle comparison, between the real neutron direction and a direction pointing from neutrino vertex to individual edep







    def edep_angle_difference(self) -> float:
        """
        Find the *first* >0.5 MeV depos by a neutron and return the angle
        between its direction and the neutron's. 0 if none or on error.
        """

        neutron_dir,neutron_mom = self.read_neutron_momentum()
        mm2cm = 0.1
        for track in self.tracks:
            if self.tracks[track.association['ancestor']].GetPDGCode() != 2112:
                continue

            for di in track.association['depoList']:
                if not (0 <= di < len(self.depos)):
                    continue
                depo = self.depos[di]
                edep = depo.GetEnergyDeposit()
                if edep < energy_threshold:
                    continue

                # midpoint in cm
                start, stop = depo.GetStart(), depo.GetStop()
                midpoint = np.array([
                    (start.X() + stop.X()) / 2 * mm2cm,
                    (start.Y() + stop.Y()) / 2 * mm2cm,
                    (start.Z() + stop.Z()) / 2 * mm2cm
                ])

                # we have our two vectors – compute and return
                return self.compute_angle_between(neutron_dir, midpoint)
        # no valid deposit found
        return 0.0

    def SelectNeutronEvent(self):
        trig = 0
        for i in range(self.tracks.size):
            track = self.tracks[i]
            pdg = track.GetPDGCode()  # Get the PDG code of the particle
            ParentId = track.GetParentId()
            if pdg  ==2112 and ParentId == -1:
                trig += 1
        return trig


    def SelectMuonEvent(self):
        trig = 0
        for i in range(self.tracks.size):
            track = self.tracks[i]
            pdg = track.GetPDGCode()  # Get the PDG code of the particle
            ParentId = track.GetParentId()
            if pdg  in (13,-13) and ParentId == -1:
                trig += 1
        return trig
    def SelectProtonEvent(self):
        trig = 0
        for i in range(self.tracks.size):
            track = self.tracks[i]
            pdg = track.GetPDGCode()
            ParentId = track.GetParentId()
            if pdg ==2212 and ParentId == -1:
                trig += 1
        return trig
    def read_neutron_momentum(self):
        """Return the combined direction and momentum magnitude of primary neutrons.

        For single-neutron events this matches the previous behaviour. When
        multiple neutrons leave the vertex, their initial three-momenta are
        summed so downstream reconstruction can treat them as one effective
        neutron, matching the requested workflow."""

        total_momentum = np.array([0.0, 0.0, 0.0])
        for track in self.tracks:
            if track.GetParentId() != -1 or track.GetPDGCode() != 2112:
                continue

            mom = track.GetInitialMomentum()
            total_momentum += np.array([mom.X(), mom.Y(), mom.Z()])

        magnitude = np.linalg.norm(total_momentum)
        if magnitude <= 0.0:
            # print("due to the fact neutron momentum is 0")
            if self.last_reco_failure_reason is None:
                self.last_reco_failure_reason = 'zero_neutron_momentum'
            return [0,0,0], 0.0

        direction_vector = total_momentum / magnitude
        return direction_vector, magnitude

    def find_nearest_discontinuous_neutron(self,plotter_object):
        # Check if muon_energy exists and is not empty
        if not hasattr(plotter_object, 'muon_energy') or len(plotter_object.muon_energy) == 0:
            print("Warning: muon_energy not available in plotter object for this event.")
            return None, None, None  # Return gracefully

        currentenergy = plotter_object.muon_energy[0]  # CORRECT: Accesses the passed-in object
        muon_info = self.read_muon_direction_and_amplitude()

        if muon_info is None:
            print("Warning: Muon direction is not available.")
            return None, None, None

        muon_direction, _ = muon_info
        if muon_direction is None or np.linalg.norm(muon_direction) == 0:
            print("Warning: Muon direction vector is invalid.")
            return None, None, None
        muon_direction = np.array(muon_direction, dtype=float)
        muon_norm = np.linalg.norm(muon_direction)
        mm2m = 0.001
        nearest_distance = float('inf')
        nearest_position = None
        nearest_t = 0

        for i, track in enumerate(self.tracks):
            depoList = track.association['depoList']
            ancestor = track.association['ancestor']
            pdg_original = self.tracks[ancestor].GetPDGCode()
            pdg = self.tracks[i].GetPDGCode()
            if pdg != 2112 and pdg_original != 2112:
                continue
            for di in depoList:
                depo = self.depos[di]
                x = (depo.GetStart().X() + depo.GetStop().X()) / 2 * mm2m
                y = (depo.GetStart().Y() + depo.GetStop().Y()) / 2 * mm2m
                z = (depo.GetStart().Z() + depo.GetStop().Z()) / 2 * mm2m
                t = (depo.GetStart().T() + depo.GetStop().T()) / 2
                edep = depo.GetEnergyDeposit()
                if edep <= self.energy_threshold:
                    continue
                distance = np.linalg.norm([x, y, z])
                if distance < nearest_distance:
                    nearest_distance = distance
                    nearest_position = np.array([x, y, z])
                    nearest_t = t

        # Allow low-energy protons to override neutron deposits when closer to the vertex.
        proton_candidates = self.proton_threshold()
        if proton_candidates:
            for is_track, _, _, proton_position, proton_time in proton_candidates:
                if is_track or proton_position is None:
                    continue
                proton_position = np.array(proton_position, dtype=float)
                if not np.all(np.isfinite(proton_position)):
                    continue
                proton_position_m = proton_position / 100.0  # cm -> m
                proton_distance = np.linalg.norm(proton_position_m)
                if proton_distance <= 0:
                    continue
                if proton_distance < nearest_distance:
                    nearest_distance = proton_distance
                    nearest_position = proton_position_m
                    if proton_time is not None:
                        nearest_t = proton_time

        if nearest_position is None or np.linalg.norm(nearest_position) == 0:
            print("there is a event where we aren't able to find the nearest neutron position")
            nearest_position = np.zeros(3)
            position_angle = 0.0
        else:
            muon_dir = muon_direction / muon_norm
            position_dir = nearest_position / np.linalg.norm(nearest_position)
            position_angle = np.degrees(np.arccos(np.clip(np.dot(muon_dir, position_dir), -1.0, 1.0)))

        return nearest_position, position_angle, nearest_t



    #read the direction of a muon, exactly the same as neutron
    def read_muon_direction_and_amplitude(self):

        """here we just read the muon direction AND AMPLITUDE from
        intrisic function """
        direction_vector = []
        muon_mom_amplitude = 0
        trig = self.SelectMuonEvent()
        if trig != 1:
        #    print("this event does not have one neutron neutrino interaction")
            print("Warning: Muon number is not available in plotter object for this event.")
            return [0,0,0],0
        else:
            #print("this event has muon neutrino interaction")
            i = 0
            track = self.tracks[i]
            parId = track.GetParentId()
            pdg_original = track.GetPDGCode()
            while parId != -1 or pdg_original not in (-13,13):
                i += 1
                track = self.tracks[i]
                parId = track.GetParentId()
                pdg_original = track.GetPDGCode()

            px = track.GetInitialMomentum().X()
            py = track.GetInitialMomentum().Y()
            pz = track.GetInitialMomentum().Z()
            p_square = (px ** 2 + py ** 2 + pz ** 2) ** 0.5
            if p_square == 0.0:
                direction_vector = np.array([0.0, 0.0, 0.0])
            else:
                direction_vector.append(px / p_square)
                direction_vector.append(py / p_square)
                direction_vector.append(pz / p_square)
                direction_vector = np.array(direction_vector)
                muon_mom_amplitude = p_square

        #print("the direction of the muon:",direction_vector)

        return direction_vector,muon_mom_amplitude

    def read_electron_direction(self):
        """Return the unit vector of the primary electron/positron at the vertex."""
        direction = np.zeros(3)
        magnitude = 0.0

        for track in self.tracks:
            if track.GetParentId() != -1:
                continue
            if track.GetPDGCode() not in (11, -11):
                continue

            px = track.GetInitialMomentum().X()
            py = track.GetInitialMomentum().Y()
            pz = track.GetInitialMomentum().Z()
            magnitude = math.sqrt(px ** 2 + py ** 2 + pz ** 2)
            if magnitude == 0.0:
                direction = np.zeros(3)
            else:
                direction = np.array([px, py, pz]) / magnitude
            break

        return direction, magnitude
    def read_electron_energy(self):
        electron_energy = 0
        for i, particle in enumerate(self.vertex.Particles):
            pdg = particle.GetPDGCode()
            KE = 0
            if pdg in (-11, 11):
                mom = particle.GetMomentum()
                KE = mom.E() - mom.M()
                if KE == 0:
                    print("this muon here is without Kinetic energy")
                else:
                    electron_energy = KE
        return electron_energy

    def read_muon_energy(self):
        muon_energy = 0
        for i, particle in enumerate(self.vertex.Particles):
            pdg = particle.GetPDGCode()
            KE = 0
            if pdg in (-13,13):
                mom = particle.GetMomentum()
                KE = mom.E() - mom.M()
                if KE == 0:
                    print("this muon here is without Kinetic energy")
                elif KE<0:
                    print(KE,'is muon energy')
                else:
                    muon_energy = KE
        return muon_energy


    def get_momentum_direction_and_amplitude(self, n_event):

        trig = self.select_muon_neutron_proton_event(n_event)
        # if trig != 1:
        #     # print("trig not 1")
        #     self.last_reco_failure_reason = 'select_muon_neutron_proton_event_failed'
        #     self.last_truth_failure_reason = 'select_muon_neutron_proton_event_failed'
        #     return None

        total_mom_mu = np.array([0.0, 0.0, 0.0])
        total_mom_neu = np.array([0.0, 0.0, 0.0])

        # Loop through all particles in the event's primary vertex.
        for particle in self.vertex.Particles:
            pdg = particle.GetPDGCode()

            # If it's a muon or anti-muon, add its momentum to the total.
            if pdg in (-13, 13):
                mom_mu = particle.GetMomentum()
                total_mom_mu += np.array([mom_mu.X(), mom_mu.Y(), mom_mu.Z()])

            # If it's a neutron, add its momentum to the total.
            if pdg == 2112:
                mom_neu = particle.GetMomentum()
                total_mom_neu += np.array([mom_neu.X(), mom_neu.Y(), mom_neu.Z()])
        
        self.last_truth_failure_reason = 'success'
        return total_mom_mu,total_mom_neu



    def read_neutron_energy(self):
        KE =0
        for i, particle in enumerate(self.vertex.Particles):
            pdg = particle.GetPDGCode()
            if pdg == 2112:
                mom = particle.GetMomentum()
                KE = mom.E() - mom.M()
        return KE

    

    def Brutal_selection_muon(self,angle_of_cone):
        # working for muon neutrino events specifically; assumes Jump already called
        select_right = 1
        angle_degrees, nearest_position = self.find_nearest_discontinuous_nedep()
        # print("angle_degrees:",angle_degrees)
        threshold = self.energy_threshold
        mm2cm = 0.1
        neutron_length = np.linalg.norm(nearest_position)
        if angle_degrees < angle_of_cone:
            select_right = 0
            # print("the nearest neutron energy deposition falls inside of the cone")
            return select_right
        muon_direction,mom = self.read_muon_direction_and_amplitude()
        for i in range(self.tracks.size):
            track = self.tracks[i]
            ancestor_id = track.association["ancestor"]
            ancestor_track = self.tracks[ancestor_id]
            ancestor_pdg = ancestor_track.GetPDGCode()
            pdg = track.GetPDGCode()
            if abs(pdg) == 13 or abs(ancestor_pdg) == 13:
                for di in track.association.get('depoList', []):
                    if not (0 <= di < len(self.depos)):
                        continue
                    depo = self.depos[di]
                    edep = depo.GetEnergyDeposit()
                    if edep < threshold:
                        continue
                    start, stop = depo.GetStart(), depo.GetStop()
                    midpoint = np.array([
                        (start.X() + stop.X()) / 2 * mm2cm,
                        (start.Y() + stop.Y()) / 2 * mm2cm,
                        (start.Z() + stop.Z()) / 2 * mm2cm
                    ])
                    midpointlength = np.linalg.norm(midpoint)
                    if midpointlength >= self.minimum_length:
                        continue
                    # Compute angle to true neutron direction
                    angle = self.compute_angle_between(muon_direction, midpoint)
                    angle = np.degrees(angle)
                    if angle > angle_of_cone and midpointlength < neutron_length:
                        select_right = 0
                        return select_right

        return select_right

    def Brutal_selection_electron(self,angle_of_cone):
        # working for electron neutrino events specifically; assumes Jump already called
        select_right = 1
        angle_degrees, nearest_position = self.find_nearest_discontinuous_nedep_EV()
        neutron_length = np.linalg.norm(nearest_position)
        if angle_degrees < angle_of_cone:
            select_right = 0
            return select_right
        threshold = self.energy_threshold
        mm2cm =0.1
        electron_direction, _ = self.read_electron_direction()
        for i in range(self.tracks.size):
            track = self.tracks[i]
            ancestor_id = track.association["ancestor"]
            ancestor_track = self.tracks[ancestor_id]
            ancestor_pdg = ancestor_track.GetPDGCode()
            pdg = track.GetPDGCode()
            if abs(pdg)==11 or abs(ancestor_pdg)==11:
                for di in track.association.get('depoList', []):
                    if not (0 <= di < len(self.depos)):
                        continue
                    depo = self.depos[di]
                    edep = depo.GetEnergyDeposit()
                    if edep < threshold:
                        continue
                    start, stop = depo.GetStart(), depo.GetStop()
                    midpoint = np.array([
                        (start.X() + stop.X()) / 2 * mm2cm,
                        (start.Y() + stop.Y()) / 2 * mm2cm,
                        (start.Z() + stop.Z()) / 2 * mm2cm
                    ])
                    midpointlength = np.linalg.norm(midpoint)
                    if midpointlength <= self.minimum_length:
                        continue
                    # Compute angle to true neutron direction
                    angle= self.compute_angle_between(electron_direction, midpoint)
                    angle = np.degrees(angle)
                    if angle> angle_of_cone and midpointlength < neutron_length:
                        select_right = 0
                        return select_right

        return select_right

    def particle_energy_per_event(self,particle_pdg_code):
        total_energy=0
        for particle in self.vertex.Particles:
            pdg = particle.GetPDGCode()
            if pdg == particle_pdg_code:
                mom = particle.GetMomentum()
                mass = mom.M()
                KE = mom.E() - mass
                total_energy += KE
        return total_energy
    def particle_momentum_amplitude_per_event(self,particle_pdg_code):
        total_momentum_amplitude=0
        allowed_pdgs = {11, 13, 211,-211,111,-11,-13, 111, 2112, 2212}  # e, μ, π±/π0, n, p
        if abs(particle_pdg_code) not in allowed_pdgs:
            return 0.0
        for particle in self.vertex.Particles:
            pdg = particle.GetPDGCode()
            if pdg == particle_pdg_code:
                mom = particle.GetMomentum()
                mass = mom.M()
                KE = mom.E() - mass
                # print("particle mass",mass,"particle KE",KE)
                momentum = np.sqrt(KE ** 2 + 2 * KE * mass )
                # print("particle momentum",momentum)
                total_momentum_amplitude += momentum
        return total_momentum_amplitude

    def muon_decay_indicator(self):
        muon_decay = 0
        for i in range(self.tracks.size):
                track = self.tracks[i]
                pdg = track.GetPDGCode()
                if abs(pdg)==11:
                    mom = track.GetInitialMomentum()
                    KE = mom.E() - mom.M()
                    if KE >=15:
                        #assume a michel electron is found
                        print("a michel electron is found")
                        muon_decay = 1


        return muon_decay

    def print_neutron_process(self):
        """here we do not constraint to the N-M only events but
        considering all the events with a neutron in it"""
        # print("\n ---Analyzing Secondary Vertices ---")
        primary_vtx_pos = self.vertex.GetPosition()
        process_list = []
        sub_process_list = []
        for i in range(self.tracks.size):
            track = self.tracks[i]
            pdg= self.tracks[i].GetPDGCode()
            if  track.association['children'] and pdg ==2112:
                for point in track.Points:
                    process = point.GetProcess()
                    sub_process = point.GetSubprocess()
                # print(f"  - Interaction:     '{process}' (Sub-process: '{sub_process}')") # Printing the result
                    process_list.append(process)
                    sub_process_list.append(sub_process)
        return process_list, sub_process_list



    def read_particle_direction(self):
        """here we read out the direction of particles from vertex, including muon neutron and what ever other particles we have"""
        particle_direction_list = []
        for i, particle in enumerate(self.vertex.Particles):
            mom = particle.GetMomentum()
            particle_mom = np.array([mom.X(), mom.Y(), mom.Z()])
            particle_mom/=np.linalg.norm(particle_mom)
            pdg= particle.GetPDGCode()
            if pdg !=22:
                particle_direction_list.append([particle_mom[0],particle_mom[1],particle_mom[2], pdg])
        return particle_direction_list




    def selectmuonneutronevent(self,nevent):
        self.Jump(nevent)
        mu_count = 0
        neutron_count = 0
        other_count = 0
        gamma_count = 0
        for particle in self.vertex.Particles:
            pdg = particle.GetPDGCode()
            if pdg == 2112:
                neutron_count += 1
            elif abs(pdg) == 13:
                mu_count += 1
            elif pdg == 22:
                gamma_count += 1
            else:
                other_count += 1

        if mu_count != 1:
            return 0
        if neutron_count < 1:
            return 0
        if other_count != 0:
            return 0
        return 1

    def select_muon_neutron_proton_event(self, nevent):
        """Return 1 when the primary vertex has 1 muon, >=1 neutron, >=1 proton and no other hadrons."""

        self.Jump(nevent)
        mu_count = 0
        neutron_count = 0
        proton_count = 0
        other_count = 0
        gamma_count = 0

        for particle in self.vertex.Particles:
            pdg = particle.GetPDGCode()
            if pdg == 2112:
                neutron_count += 1
            elif pdg == 2212:
                proton_count += 1
            elif abs(pdg) == 13:
                mu_count += 1
            elif pdg == 22:
                gamma_count += 1
            else:
                other_count += 1

        if mu_count != 1:
            return 0
        if neutron_count < 1:
            return 0
        if proton_count < 0:
            return 0
        if other_count != 0:
            return 0

        return 1
    def count_primary_particles(self):
        neutron_event = 0
        multi_neutron_event = 0
        proton_event = 0
        multi_proton_event = 0
        neutron_proton_muon_event = 0
        pion_event = 0
        multi_pion_event = 0
        other_particle_event = 0
        M_N_event = 0
        muon_event = 0
        n_events= 1000
        proton_counts=0
        neutron_counts=0
        for idx in range(n_events):
            self.Jump(idx)
            muon_count = 0
            proton_count = 0
            neutron_count = 0
            pion_count = 0
            other_particle_count = 0
            for track in self.tracks:
                pdg = track.GetPDGCode()
                parentId = track.GetParentId()
                if parentId == -1:
                    if pdg == 2112:
                        neutron_count+=1
                        neutron_counts+=1
                    elif abs(pdg) == 13:
                        muon_count+=1
                        muon_event+=1
                    elif pdg ==2212:
                        proton_count+=1
                        proton_counts+=1
                    elif abs(pdg) == 211 or pdg == 111:
                        pion_event+=1
                        pion_count+=1
                    else:
                        other_particle_count+=1
                        # print(pdg,track.GetTrackId())
            if neutron_count > 0:
                neutron_event+=1
            if proton_count > 0:
                proton_event+=1
            if neutron_count > 1:
                multi_neutron_event+=1
            if proton_count > 1:
                multi_proton_event+=1
            if neutron_count > 0 and proton_count > 0 and muon_count > 0:
                neutron_proton_muon_event+=1
            if pion_count > 1:
                multi_pion_event+=1
            if other_particle_count > 0:
                other_particle_event+=1
            if muon_count==1 and neutron_count==1 and proton_count==0:
                M_N_event+=1

        print(f"in 1000 events,{muon_event}has a muon, {neutron_event}has at least one neutron, {pion_event}has at least one pion,{proton_event}has at least one proton")
        print(f"and {other_particle_event}has a other particle in the first vertex not list above")
        print(f"{neutron_proton_muon_event}events has 1n1p1u")
        print(f"{M_N_event}is the primary Muon Neutron event")
        print(f"{multi_pion_event}events has multiple pion")
        print(f"{multi_proton_event}events has multiple proton")
        print(f"{multi_neutron_event}events has multiple neutron")
        print(f"average proton per event{proton_counts/1000}")
        print(f"average neutron per event{neutron_counts/1000}")
        return (muon_event, neutron_event, pion_event, proton_event, other_particle_event,
                M_N_event, multi_neutron_event, multi_pion_event, multi_proton_event,
                neutron_proton_muon_event,proton_counts,neutron_counts)
    def read_vertex_particles(self):
        print("the particles at the vertex are")
        for i, particle in enumerate(self.vertex.Particles):
            print(particle)

    def Next(self):
        if self.currentEntry != self.nEntry - 1:
            self.currentEntry += 1
        else:
            self.currentEntry = 0
        self.Jump(self.currentEntry)
    def print_neutrino_direction(self):
        px = self.info['Px_nu']
        py = self.info['Py_nu']
        pz = self.info['Pz_nu']
        print(f"neutrino direction is {px,py,pz}")
        return px,py,pz

    # ------------------------
    def Prev(self):
        if self.currentEntry != 0:
            self.currentEntry -= 1
        else:
            self.currentEntry = self.nEntry - 1
        self.Jump(self.currentEntry)

    # -------------------------
    def GetEnergyDepoWithDesendents(self, trkId):
        track = self.tracks[trkId]
        energy = track.energy['depoTotal']
        children = track.association['children']
        for childId in children:
            energy += self.GetEnergyDepoWithDesendents(childId)
        return energy

    # -------------------------
    def GetChargeDepoWithDesendents(self, trkId):
        track = self.tracks[trkId]
        charge = track.energy['depoTotal_charge']
        charge_75keV = track.energy['depoTotal_charge_th_75keV']
        charge_500keV = track.energy['depoTotal_charge_th_500keV']
        children = track.association['children']
        for childId in children:
            charge += self.GetChargeDepoWithDesendents(childId)[0]
            charge_75keV += self.GetChargeDepoWithDesendents(childId)[1]
            charge_500keV += self.GetChargeDepoWithDesendents(childId)[2]
        return [charge, charge_75keV, charge_500keV]

    # -------------------------
    def GetChargeMBoxDepoWithDesendents(self, trkId):
        track = self.tracks[trkId]
        charge = track.energy['depoTotal_charge_MBox']
        charge_75keV = track.energy['depoTotal_charge_MBox_th_75keV']
        charge_500keV = track.energy['depoTotal_charge_MBox_th_500keV']
        children = track.association['children']
        for childId in children:
            charge += self.GetChargeMBoxDepoWithDesendents(childId)[0]
            charge_75keV += self.GetChargeMBoxDepoWithDesendents(childId)[1]
            charge_500keV += self.GetChargeMBoxDepoWithDesendents(childId)[2]
        return [charge, charge_75keV, charge_500keV]

    # -------------------------
    def GetLightDepoWithDesendents(self, trkId):
        track = self.tracks[trkId]
        light_avg_220PEpMeV = track.energy['depoTotal_light_avg_220PEpMeV']
        light_avg_180PEpMeV = track.energy['depoTotal_light_avg_180PEpMeV']
        light_avg_140PEpMeV = track.energy['depoTotal_light_avg_140PEpMeV']
        light_avg_100PEpMeV = track.energy['depoTotal_light_avg_100PEpMeV']
        light_avg_35PEpMeV = track.energy['depoTotal_light_avg_35PEpMeV']
        children = track.association['children']
        for childId in children:
            light_avg_220PEpMeV += self.GetLightDepoWithDesendents(childId)[0]
            light_avg_180PEpMeV += self.GetLightDepoWithDesendents(childId)[1]
            light_avg_140PEpMeV += self.GetLightDepoWithDesendents(childId)[2]
            light_avg_100PEpMeV += self.GetLightDepoWithDesendents(childId)[3]
            light_avg_35PEpMeV += self.GetLightDepoWithDesendents(childId)[4]
        return [light_avg_220PEpMeV, light_avg_180PEpMeV, light_avg_140PEpMeV, light_avg_100PEpMeV, light_avg_35PEpMeV]

    # -------------------------
    def GetLightMBoxDepoWithDesendents(self, trkId):
        track = self.tracks[trkId]
        light_avg_220PEpMeV = track.energy['depoTotal_light_avg_220PEpMeV_MBox']
        light_avg_180PEpMeV = track.energy['depoTotal_light_avg_180PEpMeV_MBox']
        light_avg_140PEpMeV = track.energy['depoTotal_light_avg_140PEpMeV_MBox']
        light_avg_100PEpMeV = track.energy['depoTotal_light_avg_100PEpMeV_MBox']
        light_avg_35PEpMeV = track.energy['depoTotal_light_avg_35PEpMeV_MBox']
        children = track.association['children']
        for childId in children:
            light_avg_220PEpMeV += self.GetLightMBoxDepoWithDesendents(childId)[0]
            light_avg_180PEpMeV += self.GetLightMBoxDepoWithDesendents(childId)[1]
            light_avg_140PEpMeV += self.GetLightMBoxDepoWithDesendents(childId)[2]
            light_avg_100PEpMeV += self.GetLightMBoxDepoWithDesendents(childId)[3]
            light_avg_35PEpMeV += self.GetLightMBoxDepoWithDesendents(childId)[4]
        return [light_avg_220PEpMeV, light_avg_180PEpMeV, light_avg_140PEpMeV, light_avg_100PEpMeV, light_avg_35PEpMeV]

    # -------------------------
    """def GetEnergyDepoWithAncestor(self, acId):
        energy = 0
        for track in self.tracks:
            ancestor = track.association['ancestor']
            if ancestor == acId:
                energy += track.energy['depoTotal']
        return energy"""

    # -----------------------
    def GetEnuFromFileName(self):
        # This is used for Marley low energy files
        # No need for Genie, can read from gRooTracker tree directly
        filename = self.GetFileName().split('/')[-1]
        filename = filename.split('_')[-2]
        filename = filename.replace('MeV', '')
        return float(filename)  # Marley file names already in MeV

    # -----------------------
    def GetnuPDGFromFileName(self):
        # This is used for Marley low energy files
        # No need for Genie, can read from gRooTracker tree directly
        filename = self.GetFileName().split('/')[-1]
        filename = filename.split('_')[-3]
        return filename

    # -----------------------
    def GetFileName(self):
        return self.simTree.GetFile().GetName()

    def GetDepoInformation_all_events(self):
        """
        Returns a list of track lengths (in centimeters) for deposits across all events.
        Assumes mm2cm is defined (e.g., mm2cm = 0.1).
        """
        all_track_lengths = []
        # Loop over all events
        for i in range(self.nEntry):
            self.Jump(i)  # Load event i
            # Now self.depos contains deposits for the current event
            for depo in self.depos:
                trkLength = depo.GetTrackLength() * mm2cm
                all_track_lengths.append(trkLength)
        return all_track_lengths


    def read_vertex_particle_momentum_antinue(self):
        """Return summed three-momenta of first-vertex particles in antinu_e CC events."""

        real_electron_momentum = np.zeros(3)
        real_proton_momentum = np.zeros(3)
        real_proton_momentum_trackable = np.zeros(3)
        real_pion_momentum = np.zeros(3)
        real_other_particle_momentum = np.zeros(3)
        real_neutron_momentum = np.zeros(3)

        for particle in self.vertex.Particles:
            pdg = particle.GetPDGCode()
            mom = particle.GetMomentum()
            momentum = np.array([mom.X(), mom.Y(), mom.Z()])
            if pdg == 2112:
                real_neutron_momentum += momentum
                continue

            if pdg == 2212:
                real_proton_momentum += momentum
                KE = mom.E() - mom.M()
                if KE > self.proton_trackable_length:
                    real_proton_momentum_trackable += momentum
                continue

            if abs(pdg) == 11:
                real_electron_momentum += momentum
                continue

            if abs(pdg) in (211, 111):
                real_pion_momentum += momentum
                continue

            real_other_particle_momentum += momentum

        return (
            real_electron_momentum,
            real_proton_momentum,
            real_pion_momentum,
            real_other_particle_momentum,
            real_neutron_momentum,
            real_proton_momentum_trackable,
        )




    def neutrino_direction_reconstruction_antinue(self, i):
        """Compare antinu_e lepton-only vs. lepton+hadron angle to the true neutrino."""

        self.Jump(i)
        self.last_reco_failure_reason = None

        truth_neutrino = np.array([
            self.info.get('Px_nu', 0.0),
            self.info.get('Py_nu', 0.0),
            self.info.get('Pz_nu', 0.0),
        ])
        # print(truth_neutrino)
        truth_norm = np.linalg.norm(truth_neutrino)
        if truth_norm == 0.0:
            self.last_reco_failure_reason = 'missing_truth_neutrino'
            return None
        real_neutrino_direction = truth_neutrino / truth_norm

        (
            real_electron_momentum,
            real_proton_momentum,
            real_pion_momentum,
            real_other_particle_momentum,
            real_neutron_momentum,
            real_proton_momentum_trackable,
        ) = self.read_vertex_particle_momentum_antinue()
        # print("\n electron momentum", real_electron_momentum)
        # print("\n proton momentum", real_proton_momentum)
        # print("\n neutron momentum", real_neutron_momentum)
        # print("\n lepton plus proton momentum", real_electron_momentum+real_proton_momentum)

        reconstructed_neutron_momentum = self.neutron_direction_and_amplitude_reconstruction_EV(self.energy_threshold)
        if reconstructed_neutron_momentum is None or not np.all(np.isfinite(reconstructed_neutron_momentum)):
            self.last_reco_failure_reason = 'neutron_reconstruction_failed'
            return None

        real_nucleus_momentum = np.array(self.read_final_state_nucleus(), dtype=float)*1000
        # print(real_nucleus_momentum)

        SUM_everything = (
            real_neutron_momentum
            + real_nucleus_momentum
            + real_electron_momentum
            + real_proton_momentum_trackable
            + real_pion_momentum
            + real_other_particle_momentum
        )

        # print(SUM_everything,"is sum everything")
        all_particle_including_neutron_momentum = (
            reconstructed_neutron_momentum
            + real_electron_momentum
            + real_proton_momentum_trackable
            + real_pion_momentum
        )
        all_particle_excluding_neutron_momentum = (
            real_electron_momentum + real_proton_momentum_trackable + real_pion_momentum
        )
        all_particle_excluding_neutron_thresholdless_proton =(
            real_electron_momentum + real_proton_momentum + real_pion_momentum
        )
        TRPM = (
            real_pion_momentum
            + real_proton_momentum
            + real_electron_momentum
            + real_neutron_momentum
            + real_other_particle_momentum
        )

        TRPM_angle = np.degrees(self.compute_angle_between(TRPM, real_neutrino_direction))

        lepton_only_angle = np.degrees(
            self.compute_angle_between(real_electron_momentum, real_neutrino_direction)
        )
        neutron_inclusion_angle = np.degrees(
            self.compute_angle_between(all_particle_including_neutron_momentum, real_neutrino_direction)
        )
        neutron_exclusion_angle = np.degrees(
            self.compute_angle_between(all_particle_excluding_neutron_momentum, real_neutrino_direction)
        )

        all_particle_excluding_neutron_thresholdless_proton_angle = np.degrees(
            self.compute_angle_between(all_particle_excluding_neutron_thresholdless_proton,real_neutrino_direction)
        )

        self.last_reco_failure_reason = 'success'
        return (
            TRPM_angle,
            lepton_only_angle,
            neutron_exclusion_angle,
            neutron_inclusion_angle,
            all_particle_excluding_neutron_thresholdless_proton_angle,
            SUM_everything,
        )

    def read_particle_energy_fraction(self):
        # work for both electron neutrino and muon neutrino
        ###
        ###
        ###
        proton_energy = 0
        neutron_energy = 0
        lepton_energy = 0
        pion_energy = 0
        gamma_energy = 0
        other_energy = 0
        total_energy = 0
        for i, particle in enumerate(self.vertex.Particles):
            pdg = particle.GetPDGCode()
            mom = particle.GetMomentum()
            KE = mom.E() - mom.M()
            total_energy += KE
            if pdg ==2212:
                proton_energy += KE
            elif pdg ==2112:
                neutron_energy += KE
            elif pdg in (11,-11,13,-13):
                lepton_energy += KE
            elif pdg == 22:
                gamma_energy += KE
            elif pdg in (111,211,-211):
                pion_energy += KE
            else:
                other_energy += KE

        return proton_energy, neutron_energy, lepton_energy, pion_energy, gamma_energy,other_energy


##########
## this will be the section where we do function for the radiactive noises
##########
    @staticmethod
    def random_point_in_ball(R):
        # 1) random direction (uniform on sphere)
        u = random.random()          # in [0,1)
        v = random.random()          # in [0,1)
        theta = 2.0 * math.pi * u
        z = 2.0 * v - 1.0            # in [-1, 1]
        xy = math.sqrt(max(0.0, 1.0 - z*z))
        x_dir = xy * math.cos(theta)
        y_dir = xy * math.sin(theta)
        z_dir = z

        # 2) random radius (uniform in volume)
        r = R * (random.random() ** (1.0 / 3.0))
        return r * x_dir, r * y_dir, r * z_dir

    def identify_radiactive_noise_influence(self,number_of_noises_point,radius_of_inspection,energy_threshold: Optional[float] = None):

    #### considering antinumu events for now
        indication=0
        threshold = self._resolve_energy_threshold(energy_threshold)
        # print("the energy threshold is:",threshold)
        muon_energy = self.read_muon_energy()
        muon_direction,mom = self.read_muon_direction_and_amplitude()
        cone_method = self.cone_applying_for_muon
        result = self.find_nearest_depo(muon_energy,cone_method,threshold)
        if result is None:
            return indication
        _, _, _, first_edep = result
        nearest_distance = np.linalg.norm(first_edep)
        noises_list = self.generate_random_energy_from_spectrum(number_of_noises_point, radius_of_inspection,threshold)
        # print("the noise list is:",noises_list)
        if noises_list is None or noises_list == []:
            # print("the indication of radiactive noise is: ", indication)
            return indication
        for i in noises_list:
            # print(i)
            distance = np.linalg.norm(np.array(i[:3]))
            # print("the distance is:",distance)
            # print("the nearest_distance is :",nearest_distance)
            if distance> nearest_distance:
                continue
            noise_direction = np.array(i[:3])/distance
            angle = self.compute_angle_between(muon_direction, noise_direction)
            deviation_angle = np.degrees(angle)
            cone_angle= cone_method[int(round(muon_energy/100,0))]
            if cone_angle > deviation_angle:
                continue
            indication = 1
        # print("the indication of radiactive noise is: ", indication)
        return indication




    def generate_random_energy_from_spectrum(self,number_of_noises_point,radius_of_inspection,energy_threshold: Optional[float] = None):
        ## for now we take the time interval 20us and the radius 5m
        #within 5m we are able to seize about 90 percent of the event
        # and 20us is a time interval we select for the noise, it could be changed later


        preferred_path = Path(__file__).resolve().parent / "Ar39_beta_KEspectrum.root"
        fallback_path = Path("~/frunshiworking/data_antinumu/Ar39_beta_KEspectrum.root").expanduser()
        root_path = preferred_path if preferred_path.exists() else fallback_path

        noises_list = []
        threshold = self._resolve_energy_threshold(energy_threshold)
        hist_name = "c1;1"  # replace after printing f.keys()
        if self._radiative_noise_sampler is None:
            self._radiative_noise_sampler = self.noise_generator.HistogramSpectrumSampler(
                str(root_path), hist_name, seed=123
            )
        energy_draws = self._radiative_noise_sampler.sample(max(1, number_of_noises_point))
        for i in range(number_of_noises_point):
            position = np.array(self.random_point_in_ball(radius_of_inspection), dtype=float) * 100.0  # convert m->cm
            x, y, z = position
            energy = energy_draws[i % len(energy_draws)]
            if energy <= threshold :
                continue
            noises_list.append(np.array([x, y, z, energy]))
        # print(noises_list)
        return noises_list


    def composite_neutron(self):
        """Return the total kinetic energy and net momentum of vertex neutrons."""

        composite_neutron_momentum = np.zeros(3)
        composite_neutron_ke = 0.0

        for particle in self.vertex.Particles:
            pdg = particle.GetPDGCode()
            if pdg == 2112:
                mom = particle.GetMomentum()
                composite_neutron_ke += float(mom.E() - mom.M())
                composite_neutron_momentum += np.array([mom.X(), mom.Y(), mom.Z()], dtype=float)
            else:
                continue

        magnitude = np.linalg.norm(composite_neutron_momentum)
        if magnitude > 0.0:
            composite_neutron_direction = composite_neutron_momentum / magnitude
        else:
            composite_neutron_direction = None

        return composite_neutron_ke, composite_neutron_momentum, composite_neutron_direction

if __name__ == "__main__":
    event = Event(sys.argv[1])
    # event.Jump(0)
    # event.PrintVertex()
    # event.Next()
    # event.PrintVertex()
    # event.PrintTracks(0,8)
    # event.PrintTrack(1)
