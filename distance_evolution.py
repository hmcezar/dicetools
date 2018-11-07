"""
Receives a .xyz and two atom numbers to compute the distance between the two points and plot
a temporal evolution.

Author: Henrique Musseli Cezar
Date: FEB/2017
"""

import argparse
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pybel
import openbabel
import os

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Receives a .xyz or .pdb and two atom numbers to compute the distance between the two points and plot a temporal evolution.")
  parser.add_argument("trjfile", help="path to the file")
  parser.add_argument("a1", help="first atom")
  parser.add_argument("a2", help="second atom")
  parser.add_argument("stepmult", nargs='?', help="step multiplier, the number of steps between each saved configuration (each angle). Default is 1 and can changed to have an x-axis with the total number of steps", default=1)

  args = parser.parse_args()

  a1 = int(args.a1)
  a2 = int(args.a2)

  # read all the molecules from file and get the distances (subtract 1 from a1 since Python indexes start with 0)
  distances = []
  for mol in pybel.readfile(os.path.splitext(args.trjfile)[1][1:], args.trjfile):
    dist = mol.atoms[a1-1].OBAtom.GetDistance(a2)
    print(dist)
    distances.append(dist)

  stepmult = int(args.stepmult)
  step = [x*stepmult for x in range(1,len(distances)+1)]

  # plot it
  mpl.rcParams.update({'font.size':18, 'text.usetex':True, 'font.family':'serif', 'ytick.major.pad':4})

  plt.plot(step,distances)

  plt.xlabel(r"MC Cycle")
  plt.ylabel(r"Distance ($\AA$)")
  plt.xlim([0,step[-1]])

  plt.savefig("distance_evolution.pdf", bbox_inches='tight')
  plt.show()
