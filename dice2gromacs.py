"""
Script used to convert the DICE inputs .dfr and .txt to a GROMACS topology .itp and a .gro.

Author: Henrique Musseli Cezar
Date: JAN/2019
"""

import os
import sys
import argparse
import openbabel
import tempfile
import pybel

# from https://stackoverflow.com/a/11541495
def extant_file(x):
    """
    'Type' for argparse - checks that file exists but does not open.
    """
    if not os.path.exists(x):
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: x does not exist
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x

def check_parameters(dfrfile, txtfile):
  # convert the .txt to a .xyz to read in pybel mol to perceive all the bonds, angles, dihedrals..
  fd, temp_path = tempfile.mkstemp(suffix=".xyz")
  fxyz = os.fdopen(fd,'w')
  # table to convert atomic number to symbols
  etab = openbabel.OBElementTable()
  with open(txtfile, 'r') as f:
    line = f.readline()
    combrule = line.strip()
    line = f.readline()
    if int(line.strip()) != 1:
      print("Your .txt should have only one molecule, the one present in the .dfr")
      sys.exit(0)
    line = f.readline()
    natoms = int(line.split()[0])
    # write .xyz header
    fxyz.write("%d\nGenerated from %s\n" % (natoms,txtfile))
    for i in range(natoms):
      line = f.readline()
      atnum = int(line.split()[1])
      x, y, z = [float(x) for x in line.split()[2:5]]
      fxyz.write("%s\t%f\t%f\t%f\n" % (etab.GetSymbol(atnum), x, y, z))
  fxyz.close()
  # read filel to a mol
  mol = pybel.readfile("xyz", temp_path).__next__()
  os.remove(temp_path)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Receives the DICE input files .dfr and .txt to generate the GROMACS input files .itp and .gro.")
  parser.add_argument("dfrfile", type=extant_file, help="the DICE .dfr")
  parser.add_argument("txtfile", type=extant_file, help="the DICE .txt")
  # parser.add_argument("--eq-from-geom", "-g", help="get the equilibrium values for bonds and angles from geometry instead of force field values.", action="store_true")

  args = parser.parse_args()

  check_parameters(args.dfrfile, args.txtfile)