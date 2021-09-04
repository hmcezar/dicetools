#!/usr/bin/env python3
"""
Given the filename of batch of several xyz molecules, compute the dihedral
angle between a1-a2-a3-a4 (provided by the user) and print them in stdout.

Author: Henrique Musseli Cezar
Date: JUL/2016
"""

import argparse
import sys
try: 
  import pybel
  import openbabel
except:
  from openbabel import pybel
  from openbabel import openbabel
import os

def get_dihedrals(fname, a1, a2, a3, a4):
	# read all the molecules from file
	for mol in pybel.readfile(os.path.splitext(fname)[1][1:], fname):
		print("%f" % mol.OBMol.GetTorsion(a1,a2,a3,a4))


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Receives a file with several xyz molecules and compute the desired torsional angle for each one of them.")
	parser.add_argument("filename", help="the file containing the molecules")
	parser.add_argument("atom1", help="number of the first atom in the dihedral")
	parser.add_argument("atom2", help="number of the second atom in the dihedral")
	parser.add_argument("atom3", help="number of the third atom in the dihedral")
	parser.add_argument("atom4", help="number of the fourth atom in the dihedral")
	args = parser.parse_args()

	get_dihedrals(args.filename, int(args.atom1), int(args.atom2), int(args.atom3), int(args.atom4))