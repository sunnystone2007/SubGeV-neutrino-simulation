#!/usr/bin/env python
import sys
import ROOT
from event import Event  # Ensure that your event module is available

def main(input_filename, output_filename):
    # Initialize the event object.
    # Adjust the second argument ('Marley' or 'Genie') based on your simulation.
    evt = Event(input_filename, 'Genie')
    
    # Create a ROOT histogram: 100 bins from -1 to 1
    hist = ROOT.TH1F("cosThetaHist", "Cosine Theta Histogram;cos(theta);Events", 100, -1, 1)
    
    # Loop over all events in the file
    n_entries = evt.nEntry
    for i in range(n_entries):
        evt.Jump(i)  # Load the i-th event.
        
        # Get the cosine of theta for the current event.
        cos_between = evt.cos_theta()  # Your method that returns the cosine value.
        
        # Check if the returned value is valid (optional).
        if cos_between is not None:
            hist.Fill(cos_between)
        else:
            print(f"Warning: cos_theta() returned None for event {i}")
    
    # Write the histogram to a ROOT file.
    out_file = ROOT.TFile(output_filename, "RECREATE")
    hist.Write()
    out_file.Close()
    print(f"Histogram saved to {output_filename}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python my_program.py <input_file.root> <output_file.root>")
        sys.exit(1)
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    main(input_file, output_file)
