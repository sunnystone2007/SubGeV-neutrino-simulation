import ROOT
from event import Event

# 1) Create & label your 1D histogram
hEventClass = ROOT.TH1F("hEventClass",
                        "Event Classification;Category;Entries(0.5GeV numu)",
                        8, 0.5, 8.5)
labels = [
    "others",         # bin 1
    "mu_pro_pion",    # bin 2
    "mu_pro",         # bin 3
    "mu_pi",          # bin 4
    "mu_only",        # bin 5
    "pro_pion",       # bin 6
    "pro_only",       # bin 7
    "pi_only",        # bin 8
]
for i, lbl in enumerate(labels, start=1):
    hEventClass.GetXaxis().SetBinLabel(i, lbl)

# 2) Loop over your events; for each, call Event.trig()
#    (replace `your_event_files` & constructor args with your own)
for filename in your_event_files:
    evt = Event(filename, 'Genie')    # or however you construct it
    trig_mu, trig_pro, trig_pi, trig_others = evt.trig()

    # 3) Decide which bin to fill
    if trig_others == 1:
        bin_idx = 1
    elif trig_mu >= 1 and trig_pro >= 1 and trig_pi >= 1:
        bin_idx = 2
    elif trig_mu >= 1 and trig_pro >= 1 and trig_pi == 0:
        bin_idx = 3
    elif trig_mu >= 1 and trig_pro == 0 and trig_pi >= 1:
        bin_idx = 4
    elif trig_mu >= 1 and trig_pro == 0 and trig_pi == 0:
        bin_idx = 5
    elif trig_mu == 0 and trig_pro >= 1 and trig_pi >= 1:
        bin_idx = 6
    elif trig_mu == 0 and trig_pro >= 1 and trig_pi == 0:
        bin_idx = 7
    elif trig_mu == 0 and trig_pro == 0 and trig_pi >= 1:
        bin_idx = 8
    else:
        # if you want to count “unknowns” you could add an extra bin/histogram
        continue

    hEventClass.Fill(bin_idx)

# 4) Write out the histogram
out = ROOT.TFile("event_classes.root", "RECREATE")
hEventClass.Write()
out.Close()
