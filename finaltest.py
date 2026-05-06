 # package imported
from array import array
from pathlib import Path
import ast
import numpy as np
import ROOT
from event import Event
from plotter import Plotter
import os



n_one_flavor_event_e = 9000
#number of electron neutrino event, 9* 1000 from 0.1 to 1 GeV
n_one_flavor_event_u = 8000
##number of muon neutrino event, 8* 1000 from 0.2 to 1 GeV
n_one_flavor_one_energy = 1000
counter =0

neutrino_energy = 0.5



#this line is used to run one specific energy neutrino
event1 = Event(" ~/frunshiworking/data_antinumu/newevent_antinumu/Geniev3_04_00d_tuneAR23/antinumu/edep_antinumu_all.root ",'Genie') #antinumu
# event2 = Event(" ~/frunshiworking/data_antinumu/newevent_antinumu/Geniev3_04_00d_tuneAR23/edep_antinue_all.root ",'Genie') #antinue
# event1= Event(" ~/frunshiworking/data_antinumu/newevent_antinumu/Geniev3_04_00d_tuneAR23/numu/edep_numu_all.root ",'Genie') #numu
# event2 = Event(" ~/frunshiworking/data_antinumu/newevent_antinumu/Geniev3_04_00d_tuneAR23/edep_nue_all.root ",'Genie') #nue
#
   
# event1 = Event(f" ~/frunshiworking/data_antinumu/newevent_antinumu/Geniev3_04_00d_tuneAR23/antinumu/edep_antinumu_{neutrino_energy}GeV_1kevts.root ",'Genie')
# event2 = Event(f" ~/frunshiworking/data_antinumu/newevent_antinumu/Geniev3_04_00d_tuneAR23/antinue/edep_antinue_{neutrino_energy}GeV_1kevts.root ",'Genie')


# event1 = Event(f" ~/frunshiworking/data_antinumu/newevent_antinumu/Geniev3_04_00d_tuneAR23/numu/edep_numu_{neutrino_energy}GeV_1kevts.root ",'Genie')
# event1 = Event(f"/storage/shared/wshi/atmnu/Geniev3_04_00d_tuneAR23/FlatEspectrum/edep_numu_0.2GeV_1GeV_10kevts.root",'Genie')
# event1 = Event(f" ~/frunshiworking/data_antinumu/newevent_antinumu/Geniev3_04_00d_tuneAR23/antinumu/edep_antinumu_{neutrino_energy}GeV_1kevts.root ",'Genie')
# event1 = Event(f"/home/frunshi/frunshiworking/edep_antinumu_0.2to1.0GeV_1kevts.root",'Genie')


# event2 = Event(f"/storage/shared/wshi/atmnu/Geniev3_04_00d_tuneAR23/antinue/edep_antinue_{neutrino_energy}GeV_1kevts.root ","Genie")
# event2 = Event(f"/storage/shared/wshi/atmnu/Geniev3_04_00d_tuneAR23/FlatEspectrum/edep_antinue_0.1GeV_1GeV_10kevts.root",'Genie')
# event2 = Event(f"/storage/shared/wshi/atmnu/Geniev3_04_00d_tuneAR23/nue/edep_nue_{neutrino_energy}GeV_1kevts.root ","Genie")
# event2 = Event(f"/storage/shared/wshi/atmnu/Geniev3_04_00d_tuneAR23/FlatEspectrum/edep_nue_0.1GeV_1GeV_10kevts.root",'Genie')
# event2 = Event(f"/home/frunshi/frunshiworking/data_antinumu/newevent_antinumu/Geniev3_04_00d_tuneAR23/antinue/edep_antinue_{neutrino_energy}GeV_1kevts.root")
# event2 = Event(f"/home/frunshi/frunshiworking/data_antinumu/newevent_antinumu/Geniev3_04_00d_tuneAR23/nue/edep_nue_{neutrino_energy}GeV_1kevts.root")
#

p_u = Plotter(event1)
# p_e = Plotter(event2)

# print("Look at first", n_one_flavor_one_energy, f"events of {neutrino_energy}GeV antinue ")
energy_threshold = event1.energy_threshold # the energy that a deposition is visible in the detector
proton_trackable_length = 70 # 70 cm for current assumption on whether proton is identifiable or not
inspection_length = 2 # 5 meters sphere centered with the neutrino vertex to be the space of interest
sum_energy_vs_ke_antinue = 0.52
sum_energy_vs_ke_antinumu =0.57
# print("looking at all the antinue event for neutron energy reconstruction")
# print(f'Muon anti neutrino ~~ reconstruction at energy {neutrino_energy}, with setup of \n '
#       f'inner sphere 4 centimeter '
#       f'and outer sphere 1.5 meter on the direction side \n'
#       f'and 1.5 meter on the energy side', '\n peak value')

# print(f"looking at{neutrino_energy}GeV Muon anti neutrino for the earliest vs nearest histogram")
# p.Jump(0,energy_threshold) #entryNo, dQ threshold


#####   commands for money plot and significant other plots
# p.antinumu_event_muon_spectrum()
# p.lepton_only_reco_hist(neutrino_energy)
# p_e.lepton_only_reco_hist_antinue(neutrino_energy)
# p.lepton_proton_pion_reco_hist(neutrino_energy)
# p.truth_reco_hist(neutrino_energy)
# p.lepton_neutron_reco_hist(neutrino_energy)
# p_u.neutron_money_plot(
#     neutrino_energy,
#     true_cone_angle_list_muon_brutal_95,
#     thresholds=[0.5, 0.4, 0.3]
# )
# p_u.neutron_money_plot(neutrino_energy,true_cone_angle_list_muon_brutal_95)
# p_e.neutron_money_plot_e_u_overleap(p_u,p_e,neutrino_energy,true_cone_angle_list)
# p_e.neutron_money_plot_EV(neutrino_energy)
# p_u.neutron_money_plot_u_overlap_cone_angle(neutrino_energy,true_cone_angle_list_p_e,true_cone_angle_list)
# print(f'looking at neutrino energy:{neutrino_energy} nue')
# p_e.neutrino_direction_leptononly_vs_neutroninclude_antinue(neutrino_energy)
# p_u.two_reconstruction_method_compare(neutrino_energy,energy_threshold, require_michel_neutron=False)
# p_u.two_reconstruction_method_compare(neutrino_energy)


#### plots on the understanding of properties of vertex particles

# p.plot_particle_energy_distribution(13,-13)
# p.plot_particle_momentum_distribution(2212)
# p.electron_energy_spectrum()
# p.count_oneelectrononeneutronevent_antinue()
# p_e.neutron_scattering_process()

# p_u.plot_earliest_depo_length_spectra()



# p_u.hadronic_energy_fraction(max_events=10000,)


### plots on neutron momentum magnitude reconstruction
# p_u.neutron_E_E_plot()


# p_u.neutron_blip_ancestry_breakdown(max_events=10000)
# p_u.neutron_E_E_plot_2d_uniform(max_events=10000)
p_u.neutron_E_E_plot_2d()

# p_e.neutron_E_E_plot_EV()
# p_e.neutron_E_E_plot_2d_EV()
# p_u.neutron_e_e_one_plus_one_dimensional()
#### two dimensional individuafl displays
# individual_display_entries = range(200)
# individual_thresholds = (0.3, 0.4, 0.5, 0.6, 0.7)
# for entry in individual_display_entries:
#     print(f"Rendering event display for entry {entry}")
#     p_e.Draw(
#         axis='yz',
#         value='time',
#         energy='MeV',
#         entry=entry,
#         thresholds=individual_thresholds,
#     )

#### other verification commands
# event.count_michel_electron(neutrino_energy)
# event.collect_cos_EV()
# lower_energy_limit =600
# higher_energy_limit =800

# # p_u.antinumu_neutron_distance_spectrum(lower_energy_limit, higher_energy_limit)
# p_u.plot_neutron_distance_spectrum(lower_energy_limit,higher_energy_limit)

# event2.generate_random_energy_from_spectrum(20,5)
# p_u.count_radiative_influence(15,5,energy_threshold)
# p_u.nearest_earliest_angle_distribution(
#     neutrino_energy,
#     p_u,
#     thresholds=[0.0, 0.2, 0.4, 0.6]
# )

# p_u.plot_neutron_earliest_depo_energy_spectrum()


#### noise generator

# event2.generate_random_energy_from_spectrum(1500,5)

####
#loop over all 1000 events and the code of that
####

#
#
# p_u.Jump(0,energy_threshold)
# n= n_one_flavor_one_energy
# while counter<=n:
#     event1.PrintTracks(0,20)
#     p_u.Next(energy_threshold)  # dQ threshold
#      # # # #
#     counter += 1




# p_u.Jump(0,energy_threshold)
# n= n_one_flavor_one_energy
# while counter<=n:
#     event2.PrintTracks(0,20)
#     print("!!!")
    # event1.identify_radiactive_noise_influence(150,5)
    # p_u.Next(energy_threshold)  # dQ threshold
    #
     # # #
    # counter += 1







 # #### plot cone angle using brutal method      #antinumu event
# muon_range= 900
# lower_lim = muon_range - 50
# upper_lim = (muon_range + 50)
# print(f"looking at {muon_range} MeV Muon")
# energy_list = []
# data_path = Path(__file__).resolve().parent / "muon_energy_continuous_spectrum"
# with open(data_path, "r", encoding="utf-8") as energy_file:
#     # file stores a Python list literal, so parse once instead of line-by-line
#     energy = ast.literal_eval(energy_file.read())
#
# for i in range(n_one_flavor_event_u):
#     p_u.Jump(i,energy_threshold)
#     if i%100==0:
#         print(i)
#     muon_energy= energy[i]
#     if lower_lim < muon_energy < upper_lim:
#         energy_list.append(i)
# print(f"number of events: {len(energy_list)}")
# # print(energy_list)
# # print("energy of these 9000 events: \n",energy)
# result_path = Path(__file__).resolve().parent / f"muon_cone_results{muon_range}_MeVMuon_numu.txt"
# with open(result_path, "w", encoding="utf-8") as result_file:
#     for angle in range(1,90,1):
#         correct_number = 0
#         counter = 0
#         for i in energy_list:
#             counter += 1
#             if counter% 1000==0:
#                 print(counter)
#             p_u.Jump(i,energy_threshold)
#             correctness = event1.Brutal_selection_muon(angle)
#             correct_number += correctness
#         correct_ratio = correct_number / len(energy_list)
#         print(f"with cone angle set to {angle}, the correct ratio is {correct_ratio}")
#         print(f"with cone angle set to {angle}, the correct ratio is {correct_ratio}")
#         result_file.write(f"energy {muon_range} MeV | Cone {angle} | Correct {correct_ratio}\n")


# #### plot cone angle using brutal method      #electron event

# energy_window = 50
#
# cone_min = 10
# cone_max = 70
#
# output_dir ="electron_cone_result"
# os.makedirs(output_dir,exist_ok=True)
# electron_range= 900
# lower_lim = electron_range-energy_window
# upper_lim = electron_range + 50
# print(f" looking at {electron_range} MeV Electrons")
# energy_list = []
# for i in range(n_one_flavor_event_e):
#     p_e.Jump(i,energy_threshold)
#     if i%100==0:
#         print("processing ",i,"events")
#     #     print("checking the energy of ",i,"event")
#     electron_energy= event2.read_electron_energy()
#     if lower_lim < electron_energy < upper_lim:
#         energy_list.append(i)
# print(f"number of events: {len(energy_list)}")
# output_file =os.path.join(output_dir,f"electron_{electron_range}.root")
# with open(output_file, "w") as f:
#     f.write("# cone_angle correct_ratio\n")
#     for angle in range(1,90,1):
#         correct_number = 0
#         counter = 0
#         for i in energy_list:
#             counter += 1
#             if counter% 100==0:
#                 print("processing ",counter,"events, of ",angle,"cone")
#             p_e.Jump(i,energy_threshold)
#             correctness = event2.Brutal_selection_electron(angle)
#             correct_number += correctness
#         correct_ratio = correct_number / len(energy_list)
#         print(f"with cone angle set to {angle}, the correct ratio is {correct_ratio}")
#         f.write(f"energy {electron_range} MeV | Cone {angle} | Correct {correct_ratio}\n")
