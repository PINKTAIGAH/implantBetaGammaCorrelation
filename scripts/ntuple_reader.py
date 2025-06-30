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

dt_first_match = file['nt_aida_implant_beta_dt;1']["dt"].array(library="numpy") 
dt_all_match = file['nt_all_candidate_ionbeta_dt;1']["dt"].array(library="numpy") 

dt = dt_all_match

time_scale = 1e6
time_threshold = 150 * time_scale
bin_widths = [0.1*time_scale, 0.2*time_scale, 0.5*time_scale, 1.5*time_scale]
dt = dt[abs(dt)<time_threshold]
fig, ax = plt.subplots(2, 2)

for (bin_width, axis) in zip(bin_widths, ax.flatten()):
   bin_number = int( (2*time_threshold)/bin_width )
   # prnt contents of root file
   axis.hist(dt/time_scale, bin_number, histtype="step") 
   axis.set_title("Implant-Decay dt")
   axis.set_xlabel("dt")
   axis.set_ylabel("counts")

plt.tight_layout()
plt.show()