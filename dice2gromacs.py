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


def read_bonds(fd):
  line = fd.readline()

  dfrBonds = {}
  while "$end" not in line:
    bnd = line.split()[0]+" "+line.split()[1]
    if bnd not in dfrBonds:
      dfrBonds[bnd] = line.split()[2:]
    else:
      dfrBonds[bnd].append(line.split()[2:])

    line = fd.readline()

  return dfrBonds


def read_angles(fd):
  line = fd.readline()

  dfrAngles = {}
  while "$end" not in line:
    ang = line.split()[0]+" "+line.split()[1]+" "+line.split()[2]
    if ang not in dfrAngles:
      dfrAngles[ang] = line.split()[3:]
    else:
      dfrAngles[ang].append(line.split()[3:])

    line = fd.readline()

  return dfrAngles


def read_dihedrals(fd):
  line = fd.readline()

  dfrDihedrals = {}
  while "$end" not in line:
    dih = line.split()[0]+" "+line.split()[1]+" "+line.split()[2]+" "+line.split()[3]
    if dih not in dfrDihedrals:
      dfrDihedrals[dih] = line.split()[4:]
    else:
      dfrDihedrals[dih].append(line.split()[4:])

    line = fd.readline()

  return dfrDihedrals


def read_improper(fd):
  line = fd.readline()

  dfrImpDih = {}
  while "$end" not in line:
    idih = line.split()[0]+" "+line.split()[1]+" "+line.split()[2]+" "+line.split()[3]
    if idih not in dfrImpDih:
      dfrImpDih[idih] = line.split()[4:]
    else:
      dfrImpDih[idih].append(line.split()[4:])

    line = fd.readline()

  return dfrImpDih


def read_dfr_dof(dfrfile):
  with open(dfrfile, 'r') as f:
    while True:
      line = f.readline()

      if not line:
        break

      if "$bond" in line:
        dfrBonds = read_bonds(f)
      elif "$angle" in line:
        dfrAngles = read_angles(f)
      elif "$dihedral" in line:
        dfrDihedrals = read_dihedrals(f)
      elif "$improper" in line:
        dfrImpDih = read_improper(f)

  return dfrBonds, dfrAngles, dfrDihedrals, dfrImpDih

def read_txt_to_mol(txtfile):
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

  # read filel to a mol, remove temp and return
  mol = pybel.readfile("xyz", temp_path).__next__()
  os.remove(temp_path)

  return mol

def check_parameters(dfrfile, txtfile):
  # read the degrees of freedom from dfr into dictionaries
  dfrBonds, dfrAngles, dfrDihedrals, dfrImpDih = read_dfr_dof(dfrfile)
  
  # read the txt to a pybel mol object
  mol = read_txt_to_mol(txtfile)

  # check if all bonds are present in the .dfr
  bondIterator = openbabel.OBMolBondIter(mol.OBMol)
  for bond in bondIterator:
    lbl1 = str(bond.GetBeginAtom().GetId()+1)+" "+str(bond.GetEndAtom().GetId()+1)
    lbl2 = str(bond.GetEndAtom().GetId()+1)+" "+str(bond.GetBeginAtom().GetId()+1)
    if (lbl1 not in dfrBonds) and (lbl2 not in dfrBonds):
      print("The bond (%s) is not specified in your .dfr. Aborting..." % lbl1)
      sys.exit(0)

  # check if all angles are present in the .dfr
  angleIterator = openbabel.OBMolAngleIter(mol.OBMol)
  for angle in angleIterator:
    angidx = [str(x+1) for x in angle]
    lbl1 = angidx[1]+" "+angidx[0]+" "+angidx[2]
    lbl2 = angidx[2]+" "+angidx[0]+" "+angidx[1]
    if (lbl1 not in dfrAngles) and (lbl2 not in dfrAngles):
      print("The angle (%s) is not specified in your .dfr. Aborting..." % lbl1)
      sys.exit(0)

  # check if all dihedrals are present in the .dfr
  torsionIterator = openbabel.OBMolTorsionIter(mol.OBMol)
  for torsional in torsionIterator:
    torsidx = [str(x+1) for x in torsional]
    lbl1 = torsidx[0]+" "+torsidx[1]+" "+torsidx[2]+" "+torsidx[3]
    lbl2 = torsidx[3]+" "+torsidx[2]+" "+torsidx[1]+" "+torsidx[0]
    if (lbl1 not in dfrDihedrals) and (lbl2 not in dfrDihedrals):
      print("The dihedral (%s) is not specified in your .dfr. Aborting..." % lbl1)
      sys.exit(0)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Receives the DICE input files .dfr and .txt to generate the GROMACS input files .itp and .gro.")
  parser.add_argument("dfrfile", type=extant_file, help="the DICE .dfr")
  parser.add_argument("txtfile", type=extant_file, help="the DICE .txt")
  # parser.add_argument("--eq-from-geom", "-g", help="get the equilibrium values for bonds and angles from geometry instead of force field values.", action="store_true")

  args = parser.parse_args()

  check_parameters(args.dfrfile, args.txtfile)