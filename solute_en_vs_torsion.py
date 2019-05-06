"""
Receives a .dat containing a set of dihedral angles and plot U_xs and U1_intra
vs. these dihedrals.

Author: Henrique Musseli Cezar
Date: MAY/2019
"""

import argparse
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Receives a .dat containing a set of dihedral angles and plot U_xs and U1_intra vs these dihedrals")
  parser.add_argument("dihfile", help=".dat containing the dihedral angles, one per line")
  parser.add_argument("ienfile", help="DICE's .ien file")
  parser.add_argument("e12file", help="DICE's .e12 file")

  args = parser.parse_args()

  # read data into arrays
  dihedrals = np.loadtxt(args.dihfile)
  confs = np.loadtxt(args.ienfile, skiprows=1, usecols=(0), dtype=int)
  uintra = np.loadtxt(args.ienfile, skiprows=1, usecols=(1))
  uxs = np.array([])
  cnt = 0
  with open(args.e12file,"r") as f:
    for line in f:
      if "pV" in line:
        break
    for line in f:
      if cnt >= confs.shape[0]:
        break
      if line.split()[0].startswith(str(confs[cnt])):
        cnt += 1
        uxs = np.append(uxs, float(line.split()[2]))

  mpl.rcParams.update({'font.size':18, 'text.usetex':True, 'font.family':'serif', 'ytick.major.pad':4})

  fig, ax1 = plt.subplots()
  ax1.scatter(dihedrals, uintra, label="$U_{intra}$", s=1)
  
  ax2 = ax1.twinx()
  ax2.scatter(dihedrals, uxs, label="$U_{xs}$", color="tab:red", s=1)
  
  ax1.set_xlim([-180,180])
  ax1.set_xticks([-180,-120,-60,0,60,120,180])
  ax1.set_ylabel("$U_{intra}$ (kcal/mol)")
  ax2.set_ylabel("$U_{xs}$ (kcal/mol)")
  h1, l1 = ax1.get_legend_handles_labels()
  h2, l2 = ax2.get_legend_handles_labels()
  ax1.legend(h1+h2, l1+l2, bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=2, markerscale=4)
  ax1.set_xlabel(r"$\phi$ ($^\circ$)")
  fig.tight_layout() 

  plt.savefig("energies_"+os.path.splitext(args.dihfile)[0]+".pdf", bbox_inches='tight')
