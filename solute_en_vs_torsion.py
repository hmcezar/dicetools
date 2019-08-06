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
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic
from distutils.spawn import find_executable

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Receives a .dat containing a set of dihedral angles and plot U_xs and U1_intra vs these dihedrals")
  parser.add_argument("dihfile", help=".dat containing the dihedral angles, one per line")
  parser.add_argument("ienfile", help="DICE's .ien file")
  parser.add_argument("e12file", help="DICE's .e12 file")
  parser.add_argument("nbins", nargs="?", help="number of bins used to bin the data and take the averages over the bins (default = 36)", type=int, default=36)

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

  # plot it
  if find_executable('latex') and find_executable('dvipng'):
    mpl.rcParams.update({'font.size':18, 'text.usetex':True, 'font.family':'serif', 'ytick.major.pad':4})
  else:
    mpl.rcParams.update({'font.size':18, 'font.family':'serif', 'ytick.major.pad':4})

  fig1 = plt.figure()
  ax1 = fig1.add_subplot(1,1,1)
  ax1.scatter(dihedrals, uintra, label="$U_{intra}$", s=1)
  
  ax2 = ax1.twinx()
  ax2.scatter(dihedrals, uxs, label="$U_{xs}$", color="tab:red", s=1)
  
  ax1.set_xlim([-180,180])
  ax1.set_xticks([-180,-120,-60,0,60,120,180])
  ax1.set_ylabel("$U_{intra}$ (kcal/mol)")
  ax1.set_xlabel(r"$\phi$ ($^\circ$)")
  ax2.set_ylabel("$U_{xs}$ (kcal/mol)")
  h1, l1 = ax1.get_legend_handles_labels()
  h2, l2 = ax2.get_legend_handles_labels()
  ax1.legend(h1+h2, l1+l2, bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=2, markerscale=4)

  # set both y axis having the same range
  diffax1 = ax1.get_ylim()[1] - ax1.get_ylim()[0]
  diffax2 = ax2.get_ylim()[1] - ax2.get_ylim()[0]
  if diffax1 > diffax2:
    center = ax2.get_ylim()[0] + diffax2/2.
    ax2.set_ylim([center-diffax1/2.,center+diffax1/2.])
  else:
    center = ax1.get_ylim()[0] + diffax1/2.
    ax1.set_ylim([center-diffax2/2.,center+diffax2/2.])

  fig1.savefig("energies_"+os.path.splitext(args.dihfile)[0]+".pdf", bbox_inches='tight')

  # now bin the data taking averages
  mean_uintra = binned_statistic(dihedrals, uintra, statistic='mean', bins=360/args.nbins, range=(-180.0, 180.0))
  stddev_uintra = binned_statistic(dihedrals, uintra, statistic='std', bins=360/args.nbins, range=(-180.0, 180.0))

  mean_uxs = binned_statistic(dihedrals, uxs, statistic='mean', bins=360/args.nbins, range=(-180.0, 180.0))
  stddev_uxs = binned_statistic(dihedrals, uxs, statistic='std', bins=360/args.nbins, range=(-180.0, 180.0))

  # prepare the data to plot
  midbins = (mean_uintra.bin_edges[1:] + mean_uintra.bin_edges[:-1]) / 2.
  mask = np.ma.masked_invalid(mean_uintra.statistic)
  midbins = midbins[~mask.mask]
  uintra_avg = mean_uintra.statistic[~mask.mask]
  uintra_err = stddev_uintra.statistic[~mask.mask]
  uxs_avg = mean_uxs.statistic[~mask.mask]
  uxs_err = stddev_uxs.statistic[~mask.mask]

  fig2 = plt.figure()
  ax1 = fig2.add_subplot(1,1,1)
  ax1.errorbar(midbins, uintra_avg, marker="o", yerr=uintra_err, label=r"$\langle U_{intra} \rangle $", markersize=6, capsize=2)

  ax2 = ax1.twinx()
  ax2.errorbar(midbins, uxs_avg, marker="s", yerr=uxs_err, label=r"$\langle U_{xs} \rangle$", color="tab:red", markersize=6, capsize=2)

  # write the data to screen
  print("# phi (deg)\tavg(u_intra)\tstddev(u_intra)\tavg(u_uxs)\tstddev(u_uxs)")
  for x, ui, uie, ux, uxe in zip(midbins,uintra_avg,uintra_err,uxs_avg,uxs_err):
    print("%f\t%f\t%f\t%f\t%f" % (x, ui, uie, ux, uxe))

  ax1.set_xlim([-180,180])
  ax1.set_xticks([-180,-120,-60,0,60,120,180])
  ax1.set_xlabel(r"$\phi$ ($^\circ$)")
  ax1.set_ylabel(r"$\langle U_{intra} \rangle$ (kcal/mol)")
  ax2.set_ylabel(r"$\langle U_{xs} \rangle$ (kcal/mol)")
  h1, l1 = ax1.get_legend_handles_labels()
  h2, l2 = ax2.get_legend_handles_labels()
  ax1.legend(h1+h2, l1+l2, bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=2, markerscale=1)

  # set both y axis having the same range
  diffax1 = ax1.get_ylim()[1] - ax1.get_ylim()[0]
  diffax2 = ax2.get_ylim()[1] - ax2.get_ylim()[0]
  if diffax1 > diffax2:
    center = ax2.get_ylim()[0] + diffax2/2.
    ax2.set_ylim([center-diffax1/2.,center+diffax1/2.])
  else:
    center = ax1.get_ylim()[0] + diffax1/2.
    ax1.set_ylim([center-diffax2/2.,center+diffax2/2.])

  fig2.savefig("avg_en_"+os.path.splitext(args.dihfile)[0]+".pdf", bbox_inches='tight')
