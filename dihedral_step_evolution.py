"""
Receives a '.dat' file containing dihedrals to plot them across a step axis.

Author: Henrique Musseli Cezar
Date: FEB/2017
"""

import argparse
import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Receives a '.dat' file containing angles to plot the evolution with the steps.")
  parser.add_argument("fangles", help="path to the file containing the angles, one angle per line.")
  parser.add_argument("stepmult", nargs='?', help="step multiplier, the number of steps between each saved configuration (each angle). Default is 1 and can changed to have an x-axis with the total number of steps", default=1)

  args = parser.parse_args()

  # read data and prepare the lists
  angles = []
  with open(args.fangles,'r') as f:
    for line in f:
      angles.append(float(line.strip()))

  stepmult = int(args.stepmult)
  step = [x*stepmult for x in range(1,len(angles)+1)]

  # plot it
  mpl.rcParams.update({'font.size':18, 'text.usetex':True, 'font.family':'serif', 'ytick.major.pad':4})

  plt.scatter(step,angles,s=2)

  plt.xlabel(r"MC Cycle")
  plt.ylabel(r"$\phi$ ($^\circ$)")
  plt.xlim([0,step[-1]])
  plt.ylim([-180,180])
  plt.yticks([-180,-120,-60,0,60,120,180])
  plt.savefig(os.path.splitext(args.fangles)[0]+".pdf", bbox_inches='tight')
  plt.show()
