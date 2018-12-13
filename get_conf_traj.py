"""
Given the filename of a DICE xyz trajectory and a number of a configuration
returns the configuration.

Author: Henrique Musseli Cezar
Date: DEC/2018
"""

import argparse

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Given the filename of a DICE xyz trajectory and a number of a configuration returns the configuration.")
  parser.add_argument("trajfile", help="the DICE trajectory in .xyz")
  parser.add_argument("confnum", type=int, help="the number of the configuration to be extracted")

  args = parser.parse_args()

  with open(args.trajfile,'r') as f:
    n_atoms = int(f.readline())
    for line in f:
      if ("Configuration number" in line) and (str(args.confnum) in line):
        # print the config and exit
        print(n_atoms)
        print(line,end='')
        for i in range(n_atoms):
          lin = f.readline()
          print(lin,end='')
        break