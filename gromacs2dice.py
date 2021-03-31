#!/usr/bin/env python3
"""
Script used to convert a GROMACS topology file and geometry to the DICE .dfr and .txt files.

Author: Henrique Musseli Cezar
Date: JAN/2016
"""

import os
import sys
import argparse
import tempfile
import shutil
try: 
  import pybel
  import openbabel
  ob3 = False
except:
  from openbabel import pybel
  from openbabel import openbabel
  ob3 = True
from collections import OrderedDict
from fragGen import generate_fragfile
from clean_dof_dfr import clean_dofs

def nm2a(num):
	return round(num*10.0,4)

def j2cal(num):
	return round(num/4.184, 3)

def check_gromacs_path(path):
	if (not os.path.isdir(path)):
		print("GROMACS top directory is invalid (%s)" % (path))
		sys.exit(0)

	flist = os.listdir(FFPATH)
	if ("ffbonded.itp" not in flist) or ("ffnonbonded.itp" not in flist):
		print("Either ffbonded.itp or ffnonbonded.itp are missing in the force field directory")
		sys.exit(0)

def lookup_ljparam(atype, path, dict_types=[]):
	check_gromacs_path(path)

	with open(os.path.join(path,"ffnonbonded.itp")) as f:
		for line in f:
			if line.strip().startswith(";") or line.strip().startswith("#") or len(line.strip()) == 0:
				continue
			elif line.split()[0] == atype:
				return line
	print("Error: the atom type (%s) was not found in ffnonbonded.itp, aborting..." % atype)
	sys.exit(0)

def lookup_ljparam_ifile(atype, alist):
	for line in alist:
		if atype == line.split()[0]:
			return line
	return

def lookup_ffbond(t1, t2, path):
	check_gromacs_path(path)

	with open(os.path.join(path,"ffbonded.itp")) as f:
		while True:
			line = f.readline()
			if not line:
				break			
			if line.strip().startswith(";") or line.strip().startswith("#") or len(line.strip()) == 0:
				continue
			if line.strip() == "[ bondtypes ]":
				break

		while True:
			line = f.readline()
			if not line or line.startswith("[ "):
				break
			if line.strip().startswith(";") or line.strip().startswith("#") or len(line.strip()) == 0:
				continue
			if (line.split()[0] == t1 and line.split()[1] == t2) or (line.split()[0] == t2 and line.split()[1] == t1):
				return line
	print("Warning: the atom types (%s,%s) were not found in ffbonded.itp... Add manually later" % (t1, t2))
	return "not found"

def lookup_ffangle(t1, t2, t3, path):
	check_gromacs_path(path)

	with open(os.path.join(path,"ffbonded.itp")) as f:
		while True:
			line = f.readline()
			if not line:
				break			
			if line.strip().startswith(";") or line.strip().startswith("#") or len(line.strip()) == 0:
				continue
			if line.strip() == "[ angletypes ]":
				break

		while True:
			line = f.readline()
			if not line or line.startswith("[ "):
				break
			if line.strip().startswith(";") or line.strip().startswith("#") or len(line.strip()) == 0:
				continue
			if (line.split()[0] == t1 and line.split()[1] == t2 and line.split()[2] == t3) or (line.split()[0] == t3 and line.split()[1] == t2 and line.split()[2] == t1):
				return line
	print("Warning: the atom types (%s,%s,%s) were not found in ffbonded.itp... Add manually later" % (t1, t2, t3))
	return "not found"

def lookup_ffimproper(itype, path):
	check_gromacs_path(path)

	with open(os.path.join(path,"ffbonded.itp")) as f:
		for line in f:
			if itype in line:
				return line
	print("Error: the improper dihedral type (%s) was not found in ffbonded.itp, aborting..." % itype)
	sys.exit(0)

def lookup_ffdihedral(t1, t2, t3, t4, dtype, ffname, path):
	check_gromacs_path(path)

	fpos = []
	fnd_lines = []

	with open(os.path.join(path,"ffbonded.itp")) as f:
		while True:
			line = f.readline()
			if not line:
				break			
			if line.strip().startswith(";") or line.strip().startswith("#") or len(line.strip()) == 0:
				continue
			if line.strip() == "[ dihedraltypes ]":
				fpos.append(f.tell())

		# for every dihedtraltypes section, look for the dihedrals
		for i in range(len(fpos)):
			f.seek(fpos[i])

			# try looking for the four atoms in the dihedral
			while True:
				line = f.readline()
				if not line or line.strip().startswith("[ "):
					break
				if line.strip().startswith(";") or line.strip().startswith("#") or len(line.strip()) == 0:
					continue

				# check if the line has the right dihedral type
				if int(line.split()[4]) != dtype:
					continue

				if (line.split()[0] == t1 and line.split()[1] == t2 and line.split()[2] == t3 and line.split()[3] == t4) or (line.split()[0] == t4 and line.split()[1] == t3 and line.split()[2] == t2 and line.split()[3] == t1):
					fnd_lines.append(line)

			# if the complete dihedral was not found, try it with one missing atom
			if (len(fnd_lines) == 0):
				f.seek(fpos[i])

				while True:
					line = f.readline()
					if not line or line.strip().startswith("[ "):
						break
					if line.strip().startswith(";") or line.strip().startswith("#") or len(line.strip()) == 0:
						continue	
					# check if it's not an improper dihedral
					if int(line.split()[4]) != dtype:
						continue

					if (line.split()[0] == "X" and line.split()[1] == t2 and line.split()[2] == t3 and line.split()[3] == t4) or (line.split()[0] == t4 and line.split()[1] == t3 and line.split()[2] == t2 and line.split()[3] == "X"):
						fnd_lines.append(line)
					elif (line.split()[0] == t1 and line.split()[1] == t2 and line.split()[2] == t3 and line.split()[3] == "X") or (line.split()[0] == "X" and line.split()[1] == t3 and line.split()[2] == t2 and line.split()[3] == t1):
						fnd_lines.append(line)
		
			# if still not found, try with two missing atoms
			if (len(fnd_lines) == 0):
				f.seek(fpos[i])

				while True:
					line = f.readline()
					if not line or line.strip().startswith("[ "):
						break
					if line.strip().startswith(";") or line.strip().startswith("#") or len(line.strip()) == 0:
						continue

					# check if the line has the right dihedral type
					if int(line.split()[4]) != dtype:
						continue

					if (line.split()[0] == "X" and line.split()[1] == t2 and line.split()[2] == t3 and line.split()[3] == "X") or (line.split()[0] == "X" and line.split()[1] == t3 and line.split()[2] == t2 and line.split()[3] == "X"):
						fnd_lines.append(line)
						continue

					if (line.split()[0] == "X" and line.split()[1] == "X" and line.split()[2] == t3 and line.split()[3] == t4) or (line.split()[0] == t4 and line.split()[1] == t3 and line.split()[2] == "X" and line.split()[3] == "X"):
						fnd_lines.append(line)
						continue

					if (line.split()[0] == "X" and line.split()[1] == "X" and line.split()[2] == t2 and line.split()[3] == t1) or (line.split()[0] == t1 and line.split()[1] == t2 and line.split()[2] == "X" and line.split()[3] == "X"):
						fnd_lines.append(line)
						continue

	# in opls, we should find just one line and return it
	if "opls" in ffname:
		if len(fnd_lines) == 1:
			return fnd_lines[0]
		elif len(fnd_lines) > 1:
			for i in range(len(fnd_lines)-1):
				if fnd_lines[i].split()[5:11] != fnd_lines[i+1].split()[5:11]:
					print("Error: You have more than one dihedral that is suited to describe your torsional in ffbonded.itp, the lines are shown below.")
					print(fnd_lines[i])
					print(fnd_lines[i+1])
					sys.exit()

			return fnd_lines[0]
	# in amber, we need to sum all the contributions to a single dihedral
	else:
		if len(fnd_lines) >= 1:
			params = [0.0] * 8
			for line in fnd_lines:
				n = int(line.split()[7])
				# have to multiply by 2.0 since I use the 0.5*(...) version of the AMBER definition
				params[n-1] = 2.0*float(line.split()[6])
				params[n+3] = float(line.split()[5])

			retline = "%s\t%s\t%s\t%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f" % (t1,t2,t3,t4,dtype,params[0],params[1],params[2],params[3],params[4],params[5],params[6],params[7])				

			return retline

	print("Warning: the atom types (%s,%s,%s,%s) were not found in ffbonded.itp... Add manually later" % (t1, t2, t3, t4))
	return "not found"

def strip_comment(string, token=';'):
	return string.split(token)[0].strip()

def top2dfr(topfile, geomfile, flexfrag, eqgeom, savefrags, topcharges, ffname, path):
	if "amber" in ffname:
		potname = "AMBER"
	else:
		potname = "OPLS"

	# get the atomic positions from the geometry file
	base, ext = os.path.splitext(geomfile)
	if not ob3:
		etab = openbabel.OBElementTable()
	mol = pybel.readfile(ext[1:],geomfile).__next__()
	molxyzinfo = {}
	for i, atom in enumerate(mol,1):
		if ob3:
			molxyzinfo[i] = [openbabel.GetSymbol(atom.atomicnum),atom.coords]
		else:
			molxyzinfo[i] = [etab.GetSymbol(atom.atomicnum),atom.coords]

	# first pass through the topology file to get the data into a dictionary
	tdata = {}
	tdata["[ improper ]"] = []
	with open(topfile) as f:
		while True:
			line = f.readline()

			if not line:
				break		

			if line.strip().startswith(";") or line.strip().startswith("#") or len(line.strip()) == 0:
				continue

			if strip_comment(line).startswith("[ "):
				if strip_comment(line) not in tdata:
					key = strip_comment(line)
					tdata[key] = []
			else:
				# remove comments from the line
				if "improper" not in line:
					line = strip_comment(line)

				if "improper" in line:
					tdata["[ improper ]"].append(line)
				else:
					try:
						tdata[key].append(line)
					except:
						print("You have a line (%s) of data before assigning a type of entry (such as [ atoms ], [ bonds ] ...)" % (line))
						sys.exit(0)

	# check if atoms were found
	if "[ atoms ]" not in tdata:
		print("The [ atoms ] section was not found in your topology, make sure you're using a single file without #include")
		sys.exit(0)

	# get the atoms and its positions and parameters in a dictionary
	atoms = OrderedDict()
	rdfs = {}
	atom_num = 1
	rdf_label = 1
	for line in tdata["[ atoms ]"]:
		atomlbl = line.split()[0]
		atoms[atomlbl] = []
		fromatomtype = True
		# if atomtypes are in the topology, use them. Otherwise get from ffnonbonded.itp
		ffline = ""
		if "[ atomtypes ]" in tdata:
			ffline = lookup_ljparam_ifile(line.split()[1], tdata["[ atomtypes ]"])
		# get the data from the ffnonbonded.itp
		if not ffline:
			fromatomtype = False
			ffline = lookup_ljparam(line.split()[1], path)
		# if the atomtype was not found, stop
		if not ffline:
			print("Atom type %s was not found neither in the .itp or the force field directory" % (line.split()[1]))
			sys.exit(0)

		# append the data to the list in the same order it will be written in the .dfr
		atoms[atomlbl].append(str(atom_num))
		atom_num += 1
		if "amber" in ffname:
			atomsp = ffline.split()[1]
		else:
			atomsp = ffline.split()[2]
		if atomsp not in rdfs:
			rdfs[atomsp] = str(rdf_label)
			rdf_label += 1
		atoms[atomlbl].append(rdfs[atomsp])
		if "[ atomtypes ]" in tdata:
			atoms[atomlbl].append(mol.atoms[int(atomlbl)-1].atomicnum)
		else:
			atoms[atomlbl].append(atomsp)
		x, y, z =  molxyzinfo[int(atomlbl)][1]
		atoms[atomlbl].append(str(x))
		atoms[atomlbl].append(str(y))
		atoms[atomlbl].append(str(z))
		if ("amber" in ffname) or fromatomtype:
			if topcharges:
				atoms[atomlbl].append(line.split()[6])
			else:
				atoms[atomlbl].append(ffline.split()[-4])
			atoms[atomlbl].append(str(j2cal(float(ffline.split()[-1]))))
			atoms[atomlbl].append(str(nm2a(float(ffline.split()[-2]))))
			# the last one will not be printed but is needed to retrieve the force constants
			if ("[ atomtypes ]" in tdata) and ("opls" in ffname):
				atoms[atomlbl].append(ffline.split()[1])
			else:
				atoms[atomlbl].append(ffline.split()[0])
		else:
			if topcharges:
				atoms[atomlbl].append(line.split()[6])
			else:
				atoms[atomlbl].append(ffline.split()[4])
			atoms[atomlbl].append(str(j2cal(float(ffline.split()[7]))))
			atoms[atomlbl].append(str(nm2a(float(ffline.split()[6]))))
			# the last one will not be printed but is needed to retrieve the opls-aa force constants
			atoms[atomlbl].append(ffline.split()[1])
		# print(atoms[atomlbl])

	# now create the fragment connection list from this file and store it in fraginfo

	# creates a temporary xyz using mkstemp (https://www.logilab.org/blogentry/17873)
	# this file will be used to get the fragment data
	fd, temp_path = tempfile.mkstemp(suffix=".xyz")
	fxyz = os.fdopen(fd,'w')
	fxyz.write(mol.write("xyz"))
	fxyz.close()
	generate_fragfile(temp_path, "header")
	base, ext = os.path.splitext(temp_path)
	fraginfo = []
	with open(base+".dfr") as f:
		while True:
			line = f.readline()

			if line.strip() != "$atoms fragments":
				continue
			else:
				while line.strip() != "$end fragment connection":
					if flexfrag:
						fraginfo.append(line.replace("R","M"))
					else:
						fraginfo.append(line)
					line = f.readline()
				fraginfo.append(line)
				break
	# remove the files
	os.remove(temp_path)
	os.remove(base+".dfr")
	os.remove(base+".txt")
	if savefrags:
		shutil.move(base+"_fragments", os.path.join(os.path.dirname(os.path.abspath(geomfile)), os.path.splitext(os.path.basename(topfile))[0]+"_fragments"))
	else:
		shutil.rmtree(base+"_fragments")

	# !!! units should be converted as in: http://chembytes.wikidot.com/oplsaagro2tnk !!!

	# get the bond info
	bonds = []
	for line in tdata["[ bonds ]"]:
		# get parameters from user's .itp
		if (len(line.split()) == 5):
			if eqgeom:
				bonds.append(line.split()[0]+" "+line.split()[1]+"     \t"+str(round(j2cal(float(line.split()[4]))/(200.0),4))+"\t"+str(round(mol.OBMol.GetBond(int(line.split()[0]),int(line.split()[1])).GetLength(),4))+"\n")
			else:
				bonds.append(line.split()[0]+" "+line.split()[1]+"     \t"+str(round(j2cal(float(line.split()[4]))/(200.0),4))+"\t"+str(nm2a(float(line.split()[3])))+"\n")
		# get parameters from ffbonded.itp
		else:
			ffline = lookup_ffbond(atoms[line.split()[0]][9], atoms[line.split()[1]][9], path)
			if ffline == "not found":
				bonds.append(line.split()[0]+" "+line.split()[1]+"     \tXXX\t"+str(round(mol.OBMol.GetBond(int(line.split()[0]),int(line.split()[1])).GetLength(),4))+"\n")
			elif eqgeom:
				bonds.append(line.split()[0]+" "+line.split()[1]+"     \t"+str(round(j2cal(float(ffline.split()[4]))/(200.0),4))+"\t"+str(round(mol.OBMol.GetBond(int(line.split()[0]),int(line.split()[1])).GetLength(),4))+"\n")
			else:
				bonds.append(line.split()[0]+" "+line.split()[1]+"     \t"+str(round(j2cal(float(ffline.split()[4]))/(200.0),4))+"\t"+str(nm2a(float(ffline.split()[3])))+"\n")

	# get the angle info
	angles = []
	for line in tdata["[ angles ]"]:
		# get parameters from user's .itp
		if (len(line.split()) == 6):
			if eqgeom:
				angles.append(line.split()[0]+" "+line.split()[1]+" "+line.split()[2]+"   \tharmonic\t"+str(j2cal(float(line.split()[5]))/(2.0))+"\t"+str(round(mol.OBMol.GetAngle(mol.OBMol.GetAtom(int(line.split()[0])),mol.OBMol.GetAtom(int(line.split()[1])),mol.OBMol.GetAtom(int(line.split()[2]))),4))+"\n")
			else:
				angles.append(line.split()[0]+" "+line.split()[1]+" "+line.split()[2]+"   \tharmonic\t"+str(j2cal(float(line.split()[5]))/(2.0))+"\t"+str(float(line.split()[4]))+"\n")
		# get parameters from ffbonded.itp
		else:
			ffline = lookup_ffangle(atoms[line.split()[0]][9], atoms[line.split()[1]][9], atoms[line.split()[2]][9], path)
			if ffline == "not found":
				angles.append(line.split()[0]+" "+line.split()[1]+" "+line.split()[2]+"   \tharmonic\tXXX\t"+str(round(mol.OBMol.GetAngle(mol.OBMol.GetAtom(int(line.split()[0])),mol.OBMol.GetAtom(int(line.split()[1])),mol.OBMol.GetAtom(int(line.split()[2]))),4))+"\n")
			elif eqgeom:
				angles.append(line.split()[0]+" "+line.split()[1]+" "+line.split()[2]+"   \tharmonic\t"+str(j2cal(float(ffline.split()[5]))/(2.0))+"\t"+str(round(mol.OBMol.GetAngle(mol.OBMol.GetAtom(int(line.split()[0])),mol.OBMol.GetAtom(int(line.split()[1])),mol.OBMol.GetAtom(int(line.split()[2]))),4))+"\n")
			else:
				angles.append(line.split()[0]+" "+line.split()[1]+" "+line.split()[2]+"   \tharmonic\t"+str(j2cal(float(ffline.split()[5]))/(2.0))+"\t"+ffline.split()[4]+"\n")

	# get the dihedrals info
	dihedrals = []
	pline = {}
	ipline = {}
	dih9 = False
	dih4 = False
	for rline in tdata["[ dihedrals ]"]:
		line = strip_comment(rline)
		ffline = ""

		# get parameters from user's .itp
		if (len(line.split()) == 11 and line.split()[4] == '3'):
			V3 = round(-j2cal(float(line.split()[8])/2.0), 3)
			V2 = round(-j2cal(float(line.split()[7])), 3)
			V1 = round(-2.0*j2cal(float(line.split()[6])) + 3.0*V3, 3)
			dihedrals.append(line.split()[0]+" "+line.split()[1]+" "+line.split()[2]+" "+line.split()[3]+"   \t"+potname+"\t"+str(V1)+"\t"+str(V2)+"\t"+str(V3)+"\t0.0\t0.0\t0.0\n")
		elif (len(line.split()) == 8 and line.split()[4] == '9'):
			dih9 = True
			dihline = "%s %s %s %s" % (line.split()[0],line.split()[1],line.split()[2],line.split()[3])
			if dihline in pline:
				pline[dihline].append(line)
			else:
				pline[dihline] = [line]
		elif (len(line.split()) == 8 and (line.split()[4] == '4' or line.split()[4] == '1')):
			dih4 = True
			dihline = "%s %s %s %s" % (line.split()[0],line.split()[1],line.split()[2],line.split()[3])
			if dihline in ipline:
				ipline[dihline].append(line)
			else:
				ipline[dihline] = [line]
		# get parameters from ffbonded.itp
		elif len(line.split()) == 5:
			ffline = lookup_ffdihedral(atoms[line.split()[0]][9], atoms[line.split()[1]][9], atoms[line.split()[2]][9], atoms[line.split()[3]][9], int(line.split()[4]), ffname, path)
		else:
			print("Error: something is wrong in dihedral line (%s) maybe the number of parameters?" % (line))
			sys.exit(0)

		# parameters from ffbonded.itp need to be converted and stored
		if ffline:
			if ffline == "not found":
				dihedrals.append(line.split()[0]+" "+line.split()[1]+" "+line.split()[2]+" "+line.split()[3]+"   \t"+potname+"\tXXX\tXXX\tXXX"+"\t0.0\t0.0\t0.0\n")
				continue

			if "amber" in ffname:
				# parameters are already of Fourier type, just need to convert to cal
				if float(ffline.split()[8]) == 0.:
					dihedrals.append(line.split()[0]+" "+line.split()[1]+" "+line.split()[2]+" "+line.split()[3]+"   \t"+potname+"\t"+str(j2cal(float(ffline.split()[5])))+"\t"+str(j2cal(float(ffline.split()[6])))+"\t"+str(j2cal(float(ffline.split()[7])))+"\t"+str(round(float(ffline.split()[9]),1))+"\t"+str(round(float(ffline.split()[10]),1))+"\t"+str(round(float(ffline.split()[11]),1))+"\n")
				else:
					dihedrals.append(line.split()[0]+" "+line.split()[1]+" "+line.split()[2]+" "+line.split()[3]+"   \t"+potname+"\t"+str(j2cal(float(ffline.split()[5])))+"\t"+str(j2cal(float(ffline.split()[6])))+"\t"+str(j2cal(float(ffline.split()[7])))+"\t"+str(j2cal(float(ffline.split()[8])))+"\t"+str(round(float(ffline.split()[9]),1))+"\t"+str(round(float(ffline.split()[10]),1))+"\t"+str(round(float(ffline.split()[11]),1))+"\t"+str(round(float(ffline.split()[12]),1))+"\n")
			else:
				if float(ffline.split()[9]) != 0.0 or float(ffline.split()[10]) != 0.0:
					print("Parameters for %s - %s - %s - %s dihedrals are undefined, please treat by hand!" % (atoms[line.split()[0]][9], atoms[line.split()[1]][9], atoms[line.split()[2]][9], atoms[line.split()[3]][9]))
				V3 = round(-j2cal(float(ffline.split()[8])/2.0), 3)
				V2 = round(-j2cal(float(ffline.split()[7])), 3)
				V1 = round(-2.0*j2cal(float(ffline.split()[6])) + 3.0*V3, 3)
				dihedrals.append(line.split()[0]+" "+line.split()[1]+" "+line.split()[2]+" "+line.split()[3]+"   \t"+potname+"\t"+str(V1)+"\t"+str(V2)+"\t"+str(V3)+"\t0.0\t0.0\t0.0\n")

	# after reading all the dihedrals, if type 9 was used, we need to convert
	if dih9:
		for kdih in pline:
			params = [0.0] * 6
			for line in pline[kdih]:
				n = int(line.split()[7])
				# have to multiply by 2.0 since I use the 0.5*(...) version of the AMBER definition
				params[n-1] = 2.0*float(line.split()[6])
				params[n+2] = float(line.split()[5])
			retline = "%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f" % (kdih,9,params[0],params[1],params[2],params[3],params[4],params[5])
			dihedrals.append(retline.split()[0]+" "+retline.split()[1]+" "+retline.split()[2]+" "+retline.split()[3]+"   \t"+potname+"\t"+str(j2cal(float(retline.split()[5])))+"\t"+str(j2cal(float(retline.split()[6])))+"\t"+str(j2cal(float(retline.split()[7])))+"\t"+retline.split()[8]+"\t"+retline.split()[9]+"\t"+retline.split()[10]+"\n")

	# finally, add the improper dihedrals described as proper dihedrals
	if dih4:
		for kdih in ipline:
			params = [0.0] * 6
			for line in ipline[kdih]:
				n = int(line.split()[7])
				# have to multiply by 2.0 since I use the 0.5*(...) version of the AMBER definition
				params[n-1] = 2.0*float(line.split()[6])
				params[n+2] = float(line.split()[5])
			retline = "%s\t%d\t%.3f\t%.3f\t%.3f\t%.2f\t%.2f\t%.2f" % (kdih,9,params[0],params[1],params[2],params[3],params[4],params[5])
			dihedrals.append(retline.split()[0]+" "+retline.split()[1]+" "+retline.split()[2]+" "+retline.split()[3]+"   \t"+"AMBER"+"\t"+str(j2cal(float(retline.split()[5])))+"\t"+str(j2cal(float(retline.split()[6])))+"\t"+str(j2cal(float(retline.split()[7])))+"\t"+retline.split()[8]+"\t"+retline.split()[9]+"\t"+retline.split()[10]+"\n")

	# get the improper dihedrals info
	improper = []
	for line in tdata["[ improper ]"]:
		if len(strip_comment(line).split()) >= 7:
			improper.append(line.split()[0]+" "+line.split()[1]+" "+line.split()[2]+" "+line.split()[3]+"   \t"+str(round(j2cal(float(line.split()[6])),3))+"\t"+line.split()[5]+"\n")			
		else:
			if "opls" in ffname:
				ffline = lookup_ffimproper(line.split()[5], path)
				improper.append(line.split()[0]+" "+line.split()[1]+" "+line.split()[2]+" "+line.split()[3]+"   \t"+str(round(j2cal(float(ffline.split()[3])),3))+"\t"+ffline.split()[2]+"\n")
			else:
				ffline = lookup_ffdihedral(atoms[line.split()[0]][9], atoms[line.split()[1]][9], atoms[line.split()[2]][9], atoms[line.split()[3]][9], 4, path)
				improper.append(line.split()[0]+" "+line.split()[1]+" "+line.split()[2]+" "+line.split()[3]+"   \t"+"OPLS"+"\t"+str(j2cal(float(ffline.split()[5])))+"\t"+str(j2cal(float(ffline.split()[6])))+"\t"+str(j2cal(float(ffline.split()[7])))+"\t"+ffline.split()[8]+"\t"+ffline.split()[9]+"\t"+ffline.split()[10]+"\n")

	# print everything to the output file
	base, ext = os.path.splitext(geomfile)
	with open(base+".txt","w") as f:
		f.write("*\n1\n")
		f.write(str(len(atoms))+" \t %s (generated with gromacs2dice)\n" % (os.path.basename(base)))
		for atom, data in atoms.items():
			f.write("%2d %2d \t %7.4f \t %7.4f \t %7.4f \t %7.4f \t %7.4f \t %7.4f\n" % (int(data[1]), int(data[2]), float(data[3]), float(data[4]), float(data[5]), float(data[6]), float(data[7]), float(data[8])))
		f.write("$end\n")
	with open(base+".dfr","w") as f:
		for line in fraginfo:
			f.write(line)
		f.write("\n$bond\n")
		for line in bonds:
			f.write(line)
		f.write("$end bond\n\n$angle\n")
		for line in angles:
			f.write(line)
		f.write("$end angle\n\n$dihedral\n")
		for line in dihedrals:
			f.write(line)
		if dih4:
			for line in improper:
				f.write(line)
		f.write("$end dihedral\n")
		if improper:
			for line in improper:
				print("\n$improper dihedral\n")
				f.write(line)
				f.write("$end improper dihedral\n")

	if not flexfrag:
		withoutba =	clean_dofs(base+".dfr")
		with open(base+".dfr", 'w') as f:
			f.write(withoutba)

	print("The files %s and %s were successfully generated." % (base+".txt",base+".dfr"))
	# print "Don't forget to check the order of the atoms in the improper dihedrals (central atom first)."

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Receives a GROMACS topology and geometry and creates a '.dfr' (Dice FRagment) file.")
	parser.add_argument("topfile", help="the topology file (.top - if created from gmx pdb2gmx this file is 'topol.top') containing the molecule data for the OPLS-AA force field.")
	parser.add_argument("geomfile", help="the geometry file used to generate the topology (the order of the atoms must be the same of the topology!)")
	parser.add_argument("--gromacs-ff-path", "-p", help="specifies the GROMACS top directory (default: /usr/local/gromacs/share/gromacs/top/)", default="/usr/local/gromacs/share/gromacs/top/")
	parser.add_argument("--force-field", "-f", help="specifies the force field from the list: oplsaa, amber94, amber96, amber99, amber99sb, amber99sb-ildn, ambergs (default OPLS-AA).", default="oplsaa")
	parser.add_argument("--save-fragments", "-s", help="save the fragment configurations in .xyz.", action="store_true")
	parser.add_argument("--flexible-fragments", help="if you will perform a simulation with flexible fragments use this option to have a complete .dfr with all the parameters.", action="store_true")
	parser.add_argument("--eq-from-geom", "-g", help="get the equilibrium values for bonds and angles from geometry instead of force field values.", action="store_true")
	parser.add_argument("--charges-from-topology", "-c", help="Use the charges from the topology file instead of the force field file.", action="store_true")

	args = parser.parse_args()

	# check input consistency
	obabel_sup = ["gro", "acr", "adf", "adfout", "alc", "arc", "bgf", "box", "bs", "c3d1", "c3d2", "cac", "caccrt", "cache", "cacint", "can", "car", "ccc", "cdx", "cdxml", "cht", "cif", "ck", "cml", "cmlr", "com", "copy", "crk2d", "crk3d", "csr", "cssr", "ct", "cub", "cube", "dmol", "dx", "ent", "fa", "fasta", "fch", "fchk", "fck", "feat", "fh", "fix", "fpt", "fract", "fs", "fsa", "g03", "g92", "g94", "g98", "gal", "gam", "gamin", "gamout", "gau", "gjc", "gjf", "gpr", "gr96", "gukin", "gukout", "gzmat", "hin", "inchi", "inp", "ins", "jin", "jout", "mcdl", "mcif", "mdl", "ml2", "mmcif", "mmd", "mmod", "mol", "mol2", "molden", "molreport", "moo", "mop", "mopcrt", "mopin", "mopout", "mpc", "mpd", "mpqc", "mpqcin", "msi", "msms", "nw", "nwo", "outmol", "pc", "pcm", "pdb", "png", "pov", "pqr", "pqs", "prep", "qcin", "qcout", "report", "res", "rsmi", "rxn", "sd", "sdf", "smi", "smiles", "sy2", "t41", "tdd", "test", "therm", "tmol", "txt", "txyz", "unixyz", "vmol", "xed", "xml", "xyz", "yob", "zin"]
	geomfile = os.path.realpath(args.geomfile)
	base, ext = os.path.splitext(geomfile)
	if ext[1:] not in obabel_sup:
		print("The extension of the geometry file (%s) is not supported by OpenBabel" % (ext))
		sys.exit(0)

	if (args.force_field.lower() not in ["oplsaa","amber94","amber96","amber99","amber99sb","amber99sb-ildn","ambergs"]):
		print("You have specified an unsupported force field (%s)" % (args.force_field.lower()))
		sys.exit(0)
	else:
		ffname = args.force_field.lower()
		FFPATH = os.path.join(args.gromacs_ff_path,ffname+".ff")

	# raise warning concerning the charges
	if not args.charges_from_topology:
		print("\n!!!! ATTENTION!: using charges from the force field (GROMACS path)                                 !!!!")
		print('!!!! If you wish to use the charges from your topology, use the "--charges-from-topology" option.  !!!!\n')

	# convert the file
	top2dfr(args.topfile, args.geomfile, args.flexible_fragments, args.eq_from_geom, args.save_fragments, args.charges_from_topology, ffname, FFPATH)
