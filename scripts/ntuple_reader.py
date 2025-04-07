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
time_threshold = 100 * time_scale
binwidth = 50e6

bin_number = int(2*time_threshold/binwidth)

# prnt contents of root file
plt.hist(dt/time_scale, bin_number, histtype="step") 
plt.title("Implant-Decay dt")
plt.xlabel("dt")
plt.ylabel("counts")
plt.show()