#!/usr/bin/env python3

"""
Script used to parse Gaussian09/16 partial charges determined via Mulliken, CHelp, CHelpG, 
HLY, MK or MKUFF to GROMACS topology file generated from LigParGen.

Author: Rafael Bicudo Ribeiro
Date: SET/2022
"""

import os
import sys
import argparse

def find_natoms(topfile: str) -> int:
	"""
	Return the number of atoms in the molecule.

	PARAMETERS:
	topfile [type: str] -> topology (.itp) file

	OUTPUT:
	Number of atoms [type: int]
	"""

	with open(topfile, "r") as f:
		words = []

		line = f.readline()
		while "[ atoms ]" not in line:
			line = f.readline()

		while "[ bonds ]" not in line:
			if line.strip().startswith(";") or len(line.strip()) == 0:
				line = f.readline()
			else:
				words = line.split()
				line = f.readline()

	return int(words[0])

def top_order(topfile: str) -> list:
	"""
	Return a vector with atoms types ordered by line.

	PARAMETERS:
	topfile [type: str] -> topology (.itp) file

	OUTPUT:
	Atom types [type: list[str]]
	"""

	with open(topfile, "r") as f:
		atoms_type = []

		line = f.readline()
		while "[ atoms ]" not in line:
			line = f.readline()

		line = f.readline()

		while "[ bonds ]" not in line:
			if line.strip().startswith(";") or len(line.strip()) == 0:
				line = f.readline()
			else:
				words = line.split()
				letters = list(words[4])
				atoms_type.append(letters[0])
				line = f.readline()

	return atoms_type

def read_gaussian_charges(topfile: str, gaussianlogfile: str) -> list:
	"""
	Return a vector with Gaussian partial charges.

	PARAMETERS:
	topfile [type: str] -> topology (.itp) file
	gaussianlogfile [type: str] -> Gaussian09/16 output file with partial charges

	OUTPUT:
	qm_charges [type: list[float]] -> List with partial charge 
	gaussian_atoms [type: list[str]] -> List with the corresponding atoms
	"""

	natoms = find_natoms(topfile)

	qm_charges = []
	gaussian_atoms = []

	with open(gaussianlogfile, "r") as f:
		var_bool = False 
		i = 0

		for line in f:
			if line.find("ESP charges:") != -1:
				var_bool = True

			if var_bool and i in range(natoms+2):
				words = line.split()
				i += 1
				if len(words) > 2:
					qm_charges.append(round(float(words[2]), 4))
					gaussian_atoms.append(words[1])

		var_bool = False

	return qm_charges, gaussian_atoms

def read_mulliken_charges(topfile: str, gaussianlogfile: str) -> list:
	"""
	Return a vector with Gaussian partial charges determined via Mulliken population analysis.

	PARAMETERS:
	topfile [type: str] -> topology (.itp) file
	gaussianlogfile [type: str] -> Gaussian09/16 output file with partial charges

	OUTPUT:
	qm_charges [type: list[float]] -> List with partial charge 
	gaussian_atoms [type: list[str]] -> List with the corresponding atoms
	"""

	natoms = find_natoms(topfile)

	qm_charges = []
	gaussian_atoms = []

	with open(gaussianlogfile, "r") as f:
		var_bool = False
		i = 0

		for line in f:
			if line.find("Mulliken charges:") != -1:
				var_bool = True

			if var_bool and i in range(natoms+2):
				words = line.split()
				i += 1
				if len(words) > 2:
					qm_charges.append(round(float(words[2]), 4))
					gaussian_atoms.append(words[1])

		var_bool = False

	return qm_charges, gaussian_atoms

def correct_total_charge(qm_charges: list) -> list:
	"""
	Change the charge of the first atom to ensure quantum charges sum up to zero.

	PARAMETERS:
	qm_charges [type: list[float]] -> list with partial charges.

	OUTPUT:
	qm_charges [type: list[float]] -> list with corrected partial charge summing up to zero.
	"""

	sum = 0

	for i in qm_charges:
		sum = round(sum + i, 4)

	if sum == 0.0000:
		print("Total charge is zero.")
	else:
		qm_charges[0] = round(qm_charges[0] - sum, 4)
		print("Charge of first atom was increased by %s." % -sum)

	return qm_charges

def reorder_top_charges(qm_charges: list, grofile: str) -> list:
	"""
	Change the ordering of atoms to match with the .gro file.

	PARAMETERS:
	qm_charges [type: list[float]] -> list with partial charges.
	grofile [type: str] -> .gro file where atomic ordering matchs the topology file.

	OUTPUT:
	qm_charges [type: list[float]] -> list with partial charges reordering according to the .gro file.
	"""

	qm_charges_reordered = [0] * len(qm_charges)

	with open(grofile, "r") as f:
		for i in range(3):
			line = f.readline()
		
		for j in range(len(qm_charges)):
			words = line.split()
			qm_charges_reordered[j] = qm_charges[int(words[2])-1]
			line = f.readline()

	print("Charges were reordered using {} file.".format(grofile))

	return qm_charges_reordered

def parse_charges(topfile: str, gaussianlogfile: str, method: str, grofile: str):
	"""
	Write the new parsed_* file with quantum mechanical charges.

	PARAMETERS:
	topfile [type: str] -> topology (.itp) file
	gaussianlogfile [type: str] -> Gaussian09/16 output file with partial charges
	method [type: str] -> method used to determine the partial charges
	grofile [type: str] -> .gro file where atomic ordering matchs the topology file.

	OUTPUT:
	A new "parsed_*.itp" file with charges updated.
	"""

	# Read the file and adjust the header accordingly

	fout = open("parsed_" + os.path.basename(topfile), "w")

	with open(topfile, 'r') as f:
		line = f.readline()
		while ";   nr" not in line:
			if "GENERATED BY LigParGen Server" in line:
				if "reordered by reordered_ligpargen" in line:
					line = "; GENERATED BY LigParGen Server, reordered by reordered_ligpargen and charge parsed by parse_gaussian_charges \n"
				else: 
					line = "; GENERATED BY LigParGen Server and charge parsed by parse_gaussian_charges \n"
			fout.write(line)
			line = f.readline()

		fout.write(line)
		line = f.readline()

		# Use the appropriate function to collect QM charges from each method

		ESP_methods = ['hly', 'chelpg', 'mk', 'chelp', 'mkuff', 'hlygat']
		if method.lower() in ESP_methods:
			qm_charges, gaussian_atoms = read_gaussian_charges(topfile, gaussianlogfile)
		elif method.lower() == 'mulliken':
			qm_charges, gaussian_atoms = read_mulliken_charges(topfile, gaussianlogfile)
		else:
			print('Charge population method not supported (try "python3 parse_gaussian_charges.py -h").')
			sys.exit(0)

		# Check if topology and gaussian output atoms are in the same order

		topology_atoms = top_order(topfile)

		diff_order = False

		for i in range(len(topology_atoms)):
			if topology_atoms[i] != gaussian_atoms[i]:
				diff_order = True
				if args.grofile:
					break
				else:
					print("""Atoms in the topology file are ordered differently from the gaussian output file. 
Atom {}: {} (topology file) and {} (gaussian file).\nPlease order both files in the same way or provide the \
correctly ordered .gro file using the -g flag.""".format(i+1, topology_atoms[i], gaussian_atoms[i]))
					sys.exit(0)

		# Correct the total charge

		qm_charges = correct_total_charge(qm_charges)
		i = 0

		# If required, reorder the atomic charges according to the .gro file matching with .itp file

		if args.grofile:
			qm_charges = reorder_top_charges(qm_charges, grofile)

		# Write the new .itp with QM charges

		while ("opls" in line) or ("ppg" in line) or line.strip().startswith(";"):
			if line.strip().startswith(";"):
				fout.write(line)
				f.readline()
				continue

			words = line.split()
			# line = line.replace(words[6], str(qm_charges[i]))
			line = '{:>6}{:>11}{:>7}{:>7}{:>7}{:>7}{:>11.4f}{:>11.4f}\n'.format(words[0], words[1], words[2], words[3],
					words[4], words[5], qm_charges[i], float(words[7]))
			fout.write(line)
			line = f.readline()
			i = i + 1

		fout.write(line)
			
		for line in f:
			fout.write(line)

	print("The charges were sucessfully parsed.")

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Recieves a GROMACS topology generated from LigParGen and Gaussian output file to parse QM partial charges.")
	parser.add_argument("topfile", type=str, help="the topology file (.top file) containing the molecule data for the OPLS-AA force field.")
	parser.add_argument("gaussianlogfile", type=str, help="the gaussian log file with partial charges.")
	parser.add_argument("--method", "-m", type=str, help="the method used to determine the charge populations (chelp, chelpg, mk, mkuff, hly or hlygat).")
	parser.add_argument("--grofile", "-g", type=str, help="the .gro file with the same order as the topology file.")

	args = parser.parse_args()

	parse_charges(args.topfile, args.gaussianlogfile, args.method, args.grofile)
