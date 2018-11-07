"""
Given the filename a xyz trajectory and the interval between each save,
print a box configuration.

Author: Henrique Musseli Cezar
Date: JUN/2018
"""

import argparse
import os

def print_configs(fname, svint, nconfs):
  with open(fname,'r') as f:
    natoms = int(f.readline())
    boxsize = natoms + 2
    f.seek(0)

    confnum = 1

    # loop to print nconfs configurations
    for i in range(nconfs):
      # skip svint-1 configurations
      for j in range(svint-1):
        for lnum in range(boxsize):
          f.readline()

      for lnum in range(boxsize):
        line = f.readline()
        print(line, end="")

      confnum += 1

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Separates configurations from a xyz trajectory")
  parser.add_argument("filename", help="the xyz file containing the trajectory")
  parser.add_argument("svint", help="save interval (configuration will be saved each this number of steps)")
  parser.add_argument("nconfs", help="number of configuration to be saved")
  args = parser.parse_args()

  print_configs(args.filename, int(args.svint), int(args.nconfs))