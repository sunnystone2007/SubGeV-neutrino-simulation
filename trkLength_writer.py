#!/usr/bin/env python
import sys
import ROOT
from event_firstten import Event  # Adjust the import to match your project structure

def main():
    if len(sys.argv) < 3:
        print("Usage: TrackLength_writer.py <input_root_file> <event_generator>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    evgen = sys.argv[2]

    # Create an instance of your event (this might require additional setup in your code)
    event_inst = Event(input_file, evgen)
    # Jump to an event (here, we use the first event; adjust as needed)
    event_inst.Jump(0)
    
    # Get all the track lengths from the deposits
    track_lengths = event_inst.GetDepoInformation()
    
    if not track_lengths:
        print("No track lengths were retrieved from deposits.")
        sys.exit(1)
        
    # Determine histogram range: you can choose the number of bins and range appropriately.
    # For example, here we use 100 bins from 0 to the maximum track length found.
    max_length = max(track_lengths)
    hist = ROOT.TH1F("TrackLengthHist", "Track Length Distribution;Track Length (cm);Counts", 100, 0, max_length)
    
    # Fill the histogram with each track length value
    for tl in track_lengths:
        hist.Fill(tl)
    
    # Write the histogram to a ROOT file
    output_file = ROOT.TFile("TrackLengthHistogram.root", "RECREATE")
    hist.Write()
    output_file.Close()
    
    print("Histogram saved in TrackLengthHistogram.root")

if __name__ == "__main__":
    main()
