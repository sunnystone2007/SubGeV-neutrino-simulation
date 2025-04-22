from event import Event
from plotter import Plotter
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# Look at first 10 events
n = 152000
counter = 3
event = Event("/pnfs/dune/persistent/users/weishi/FD3/LArBath/Marley_Edepsim_noSecondaryDeposit/nue/edep_nue_*.root", 'Marley')
#event = Event("/exp/dune/app/users/weishi/VDPDSAna/PNSCali/edep_gammas_6.1MeV_10kevts.root", 'Marley')
p = Plotter(event)

print("Look at first", n, "events")

#dx_all_evts = np.array([])
#trklength_all_evts = np.array([])
#maxdr_all_evts = np.array([])
#maxdt_all_evts = np.array([])
neutron_capture_time = np.array([])

# dQ threshold in MeV
p.Jump(0, 0)

while counter == 3:
    # neutron not captured
    #if counter == 6 or counter == 5 or counter == 2 or counter == 51:
    #if counter == 0:
    p.Draw('xz', value='time', energy= 'MeV')
    #p.Draw('xz', value='charge', energy= 'MeV')
    #p.Draw('yz', value='time', energy= 'MeV')
    #p.Draw('yz', value='charge', energy= 'MeV')
    #p.Draw('xz', value='time', energy= 'MeV')
    #p.hist_LowE_dEdx()
    #p.hist_dx()
    #p.hist_trklength()
    #print("max_edep_dr: ", temp_maxdr, "max_edep_dt: ", temp_maxdt)

    # all events single edep length
    #dx_all_evts = np.concatenate((dx_all_evts, p.hist_dx()))
    # all events tracks length
    #trklength_all_evts = np.concatenate((trklength_all_evts, p.hist_trklength()))
    # all events edep max dist
    #maxdr_all_evts = np.concatenate((maxdr_all_evts, np.array( [p.evt_maxdtdr()[0]] ) ))
    #maxdt_all_evts = np.concatenate((maxdt_all_evts, np.array( [p.evt_maxdtdr()[1]] ) ))
    #print("evt: ", counter)
    #event.PrintVertex()
    #event.PrintTracks(0, 10000)
    #event.PrintTrack(2)
    if event.PrintTracks(0, 10000)[0] != -1:
        # capture
        #print("neutron KE: ", event.PrintTracks(0, 10000)[1], ", neutron trkid: ", event.PrintTracks(0, 10000)[0])
        neutron_capture_time = np.concatenate(( neutron_capture_time, np.array( [event.PrintTrack(event.PrintTracks(0, 10000)[0])] ) ))

    p.Next(0) # dQ threshold

    counter += 1


plt.hist(neutron_capture_time, range=(0, 100000), bins=100)
plt.xlabel('time [ns]')
plt.draw()
plt.savefig('plots/neutron_capture_time.pdf')
plt.clf() # important to clear figure
plt.close()
"""plt.hist(dx_all_evts, range=(0,4), bins=100)
plt.xlabel('dx [cm]')
plt.draw()
plt.savefig('plots/single_edep_dx_all_evts.pdf')
plt.clf() # important to clear figure
plt.close()

plt.hist(trklength_all_evts, range=(0,6), bins=100)
plt.xlabel('track length [cm]')
plt.yscale("log")
plt.draw()
plt.savefig('plots/primary_track_length_all_evts.pdf')
plt.clf() # important to clear figure
plt.close()"""
"""plt.hist(maxdr_all_evts, range=(0,200), bins=100)
plt.xlabel('edeps max dr [cm]')
plt.yscale("log")
plt.draw()
plt.savefig('plots/edeps_max_dr_all_evts.pdf')
plt.clf() # important to clear figure
plt.close()

plt.hist(maxdt_all_evts, range=(0,50), bins=50)
plt.xlabel('edeps max dt [ns]')
plt.yscale("log")
plt.draw()
plt.savefig('plots/edeps_max_dt_all_evts.pdf')
plt.clf() # important to clear figure
plt.close()"""
