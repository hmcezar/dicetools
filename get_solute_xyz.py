#!/usr/bin/env python3
"""
Give the filename of a DICE .xyz and number of atoms in the solvent to receive
in stdout the xyz of the solvent for every configuration in the initial file.

Author: Henrique Musseli Cezar
Date: APR/2016
"""

import argparse
import sys

def get_solute(fname, natoms):
	with open(fname, 'r') as f:
		print_atoms = False
		printed = 0
		line = f.readline()
		sys.stdout = open('new_traj_solute.xyz', 'w')
		while line:
			if "Configuration" in line:
				print_atoms = True
				sys.stdout.write("%d\n%s\n" % (natoms, line.rstrip()))
			elif ((printed < natoms) and print_atoms):
				sys.stdout.write(line)
				printed += 1
			else:
				print_atoms = False
				printed = 0

			line = f.readline()
		sys.stdout.close()
        


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Receives simulation boxes and returns the solute conformation for each box.")
	parser.add_argument("filename", help="the xyz containing the simulation boxes")
	parser.add_argument("natoms", help="the number of atoms in the solute")
	args = parser.parse_args()

	get_solute(args.filename, int(args.natoms))
