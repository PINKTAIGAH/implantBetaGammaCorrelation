import uproot, sys
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) <= 1:
   sys.exit("No arguments were given")

if len(sys.argv) > 2:
   sys.exit("Argument must be string to root file")

filename = sys.argv[1]


# Open root file
file = uproot.open(filename)

dt = file['nt_aida_implant_beta_dt;1']["dt"].array(library="numpy") 

time_scale = 1e9
time_threshold = 50 * time_scale
bin_widths = [0.1*time_scale, 0.2*time_scale, 0.5*time_scale, 1*time_scale]
fig, ax = plt.subplots(2, 2)

for (bin_width, axis) in zip(bin_widths, ax.flatten()):
   bin_number = (2*time_threshold/time_scale)/bin_width
   # prnt contents of root file
   axis.hist(dt/time_scale, bin_number, histtype="step") 
   axis.set_title("Implant-Decay dt")
   axis.set_xlabel("dt")
   axis.set_ylabel("counts")

plt.show()