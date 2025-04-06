import ROOT
from ROOT import TH1F, TFile
from event_neutron import Event  # or from your module if renamed
import sys

def fill_neutron_direction_cos(event):
    # Dummy logic: you will replace this
    # Let's say you pick the last neutron from PrintTracks()
    trkId, KE = event.PrintTracks()
    if trkId < 0:
        return None
    neutron = event.tracks[trkId]
    dir_vec = neutron.GetInitialMomentum()
    return dir_vec.CosTheta()  # or compute manually

if __name__ == "__main__":
    input_file = sys.argv[1]
    evgen = sys.argv[2]  # 'Genie' or 'Marley'
    outfile = "neutron_dircos.root"

    event = Event(input_file, evgen)

    h = TH1F("h_dircos", "Neutron Direction Cosine;cos(#theta);Events", 100, -1, 1)

    for i in range(event.nEntry):
        event.Jump(i)
        val = fill_neutron_direction_cos(event)
        if val is not None:
            h.Fill(val)

    f_out = TFile(outfile, "RECREATE")
    h.Write()
    f_out.Close()
    print(f"Histogram saved to {outfile}")
