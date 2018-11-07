"""
From a file containing the molecules, e.g., from get_solute_xyz.py, and from the .txt
containing the charges for each atom, calculat the dipole moment and print to screen.

Author: Henrique Musseli Cezar
Date: NOV/2016
"""

import argparse
import numpy as np
import pybel

def eA_to_D(val):
  return val/0.20819434

def get_dipole_moments(xyzfile, txtfile):
  # read charges from .txt file
  charges = {}
  with open(txtfile,'r') as f:
    f.readline()
    f.readline()
    natoms = int(f.readline().split()[0])
    for i in range(natoms):
      line = f.readline()
      charges[i] = float(line.split()[5])

  # for every configuration calculate the dipole moment
  for mol in pybel.readfile("xyz", xyzfile):
    mol.OBMol.Center()
    tdip = np.zeros(3)
    for i, atom in enumerate(mol):
      # get dipole
      dipole = [charges[i] * x for x in atom.coords]
      # sum the dipoles
      tdip = np.add(tdip, dipole)

    # print to screen
    print("%f" % np.linalg.norm(eA_to_D(tdip)))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Receives an .xyz file containing molecular configurations and a .txt of the DICE input with charges to calculate the dipole moments.")
  parser.add_argument("xyzfile", help="the .xyz file containing the molecules")
  parser.add_argument("txtfile", help="the .txt file containing the DICE input with the charges")
  args = parser.parse_args()

  xyzfile = args.xyzfile
  txtfile = args.txtfile

  get_dipole_moments(xyzfile, txtfile)

