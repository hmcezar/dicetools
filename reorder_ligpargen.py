"""
Receives original .pdb submitted to LigParGen and the output of LigParGen (.gro and .itp)
to reorder the output to be in the same order of the original input

Author: Henrique Musseli Cezar
Date: MAY/2019
"""

import os
import sys
import argparse
import pybel
import openbabel
import rmsd
import numpy as np
import re

def get_mol_info(mol):
  # table to convert atomic number to symbols
  etab = openbabel.OBElementTable()

  q_atoms = []
  q_all = []
  for atom in mol:
    q_atoms.append(etab.GetSymbol(atom.atomicnum))
    q_all.append(atom.coords)

  return np.asarray(q_atoms), np.asarray(q_all)  

def get_atom_correspondence(pdb, gro):
  # read pdb
  pdbmol = pybel.readfile("pdb",pdb).__next__()
  p_atoms, p_all = get_mol_info(pdbmol)

  # read gro
  gromol = pybel.readfile("gro",gro).__next__()
  q_atoms, q_all = get_mol_info(gromol)

  return rmsd.reorder_hungarian(p_atoms, q_atoms, p_all, q_all)

def reorder_gro(gro, mapp):
  fout = open("reordered_"+os.path.basename(gro), "w")
  with open(gro, "r") as f:
    fout.write(f.readline())
    fout.write(f.readline())
    
    # read the lines
    line = f.readline()
    lines = []
    while "UNL" in line:
      lines.append(line)
      line = f.readline()

    # write reordering
    for i in range(len(lines)):
      fout.write(lines[mapp[i]])

    fout.write(line)
  fout.close()

def reorder_itp(itp, mapp):
  fout = open("reordered_"+os.path.basename(itp), "w")
  with open(itp, "r") as f:
    line = f.readline()
    while ";   nr" not in line:
      fout.write(line)
      line = f.readline()

    fout.write(line)
    line = f.readline()

    # read atoms
    lines = []
    while "opls" in line:
      lines.append(line)
      line = f.readline()

    # write atoms in the right order
    for i in range(len(lines)):
      fout.write(re.sub(r'\d+',str(i+1),lines[mapp[i]],1))

    while "[ bonds ]" not in line:
      fout.write(line)
      line = f.readline()

    # write the bonds using the right labels
    while "[ angles ]" not in line:
      a1 = int(line.split()[0])-1
      a2 = int(line.split()[1])-1

      constants = "%s\t%s\t%s" % tuple(line.split()[1:])

      fout.write("%5d %5d     %s" % (mapp[a1],mapp[a2],constants))

  fout.close()


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Receives the original .pdb sent to LigParGen to reorder the output of LigParGen to have the atoms in the same order.")
  parser.add_argument("originalpdb", help="the original pdb uploaded to LigParGen")
  parser.add_argument("outgro", help="the .gro generated by LigParGen")
  parser.add_argument("outitp", help="the .itp generated by LigParGen")
  args = parser.parse_args()

  # map the indexes from one geometry to the others
  mapping = get_atom_correspondence(args.originalpdb, args.outgro)

  # based on the map return the reordered .gro
  reorder_gro(args.outgro, mapping)

  # now reorder the .itp
  reorder_itp(args.outitp, mapping)
