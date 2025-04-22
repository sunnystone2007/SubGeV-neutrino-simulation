#!/usr/bin/env python
import sys
import ROOT
from event_alpha1 import Event  # Ensure that your event module is available

def main(input_filename, output_filename):
    # Initialize the event object.
    evt = Event(input_filename, 'Genie')

    # Create a ROOT histogram: 100 bins from -1 to 1
    hist = ROOT.TH1F("AngleDifferenceHist",
                     "AngleDifferenceHist;angle(in degrees);Events",
                     180, 0, 180)

    # Loop over all events in the file
    for i in range(evt.nEntry):
        evt.Jump(i)
        angle = evt.cos_theta()  # your newly added method
        if angle is not None:
            hist.Fill(angle)
        else:
            print(f"Warning: cos_theta() returned None for event {i}")

    # Fit the histogram with a Gaussian
    # "Q" = quiet, "S" = return fit result
    fit_result = hist.Fit("gaus", "QS")
    gaus = hist.GetFunction("gaus")
    sigma = gaus.GetParameter(2)
    print(f"Fitted Gaussian sigma = {sigma:.5f}")

    # Draw & save the fitted histogram
    canvas = ROOT.TCanvas("c1", "", 800, 600)
    hist.Draw()
    gaus.SetLineColor(ROOT.kRed)
    gaus.Draw("same")
    canvas.SaveAs("AngleDifferenceHist.png")

    # Write the raw histogram (without fit) to the output ROOT file
    out_file = ROOT.TFile(output_filename, "RECREATE")
    hist.Write()
    out_file.Close()
    print(f"Histogram (no fit) saved to {output_filename}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python my_program.py <input_file.root> <output_file.root>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
