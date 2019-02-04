"""
Script used to convert the DICE inputs .dfr and .txt to a GROMACS topology .itp and a .gro.

Author: Henrique Musseli Cezar
Date: JAN/2019
"""

import os
import sys
import argparse
import tempfile
import openbabel
import pybel
from math import sqrt
try:
  from Queue import Queue
except:
  from queue import Queue

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


def a2nm(num):
  return round(num/10.0, 4)


def cal2j(num):
  return round(num*4.184, 4)


def read_bonds(fd):
  line = fd.readline()

  dfrBonds = {}
  while "$end" not in line:
    bnd = line.split()[0]+" "+line.split()[1]
    if bnd not in dfrBonds:
      dfrBonds[bnd] = [float(x) for x in line.split()[2:]]
    else:
      dfrBonds[bnd].append([float(x) for x in line.split()[2:]])

    line = fd.readline()

  return dfrBonds


def read_angles(fd):
  line = fd.readline()

  dfrAngles = {}
  while "$end" not in line:
    ang = line.split()[0]+" "+line.split()[1]+" "+line.split()[2]
    if ang not in dfrAngles:
      dfrAngles[ang] = [float(x) for x in line.split()[4:]]
    else:
      dfrAngles[ang].append([float(x) for x in line.split()[4:]])

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

  # check if everything was found and raise error if not
  try:
    dfrBonds
  except:
    print("Error: The bonds section was not found in the .dfr, please check your topology.")
    sys.exit(0)
  try:
    dfrAngles
  except:
    print("Error: The angles section was not found in the .dfr, please check your topology.")
    sys.exit(0)
  try:
    dfrDihedrals
  except:
    print("Error: The dihedrals section was not found in the .dfr, please check your topology.")
    sys.exit(0)

  return dfrBonds, dfrAngles, dfrDihedrals, dfrImpDih

def read_txt_to_mol(txtfile):
  # lists to store the charges and LJ parameters
  q = []
  eps = []
  sig = []

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
      # write the atomic symbol and coordinates
      line = f.readline()
      atnum = int(line.split()[1])
      x, y, z = [float(x) for x in line.split()[2:5]]
      fxyz.write("%s\t%f\t%f\t%f\n" % (etab.GetSymbol(atnum), x, y, z))

      # store the parameters
      qv, epsv, sigv = [float(x) for x in line.split()[5:]]
      q.append(qv)
      eps.append(epsv)
      sig.append(sigv)
  
  fxyz.close()

  # read filel to a mol, remove temp and return
  mol = pybel.readfile("xyz", temp_path).__next__()
  os.remove(temp_path)

  return mol, q, eps, sig


def get_pairs(mol):
  infty = -1

  distPairs = []
  # run BFS exploration for each atom to find the distance to others
  for atom in mol.atoms:
    q = Queue()

    vis = [False] * len(mol.atoms)
    dist = [infty] * len(mol.atoms)

    vis[atom.idx-1] = True
    dist[atom.idx-1] = 0
    q.put(atom)

    while not q.empty():
      k = q.get()
      # mark children as visited and queue them
      children = [pybel.Atom(x) for x in openbabel.OBAtomAtomIter(k.OBAtom)]
      for child in children:
        if not vis[child.idx-1]:
          vis[child.idx-1] = True
          dist[child.idx-1] = dist[k.idx-1] + 1
          q.put(child)

    for i in range(len(mol.atoms)):
      if (dist[i] == 3):
        pair1 = "%d %d" % (atom.idx, i+1)
        pair2 = "%d %d" % (i+1, atom.idx)
        if (pair1 not in distPairs) and (pair2 not in distPairs):
          distPairs.append(pair1)

  return distPairs


def read_parameters(dfrfile, txtfile):
  # read the degrees of freedom from dfr into dictionaries
  dfrBonds, dfrAngles, dfrDihedrals, dfrImpDih = read_dfr_dof(dfrfile)
  
  # read the txt to a pybel mol object
  mol, q, eps, sig = read_txt_to_mol(txtfile)

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

  return mol, q, eps, sig, dfrBonds, dfrAngles, dfrDihedrals, dfrImpDih


def itp_from_params(mol, q, eps, sig, dfrBonds, dfrAngles, dfrDihedrals, dfrImpDih):
  # table to convert atomic number to symbols
  etab = openbabel.OBElementTable()

  # !!! units are converted as the reverse of: http://chembytes.wikidot.com/oplsaagro2tnk and based on GROMACS manual

  # write header
  fcontent = """
;
; Generated by dice2gromacs
; https://github.com/hmcezar/dicetools
;
[ atomtypes ]
;name  bond_type  mass    charge   ptype   sigma         epsilon
"""
  # write the atomtypes
  for i, atom in enumerate(mol.atoms):
    fcontent += "att_%03d   %s%03d   %7.4f   0.000    A    %.5e    %.5e\n" % (i+1, etab.GetSymbol(atom.atomicnum), i+1, atom.atomicmass, a2nm(sig[i]), cal2j(eps[i]))

  fcontent += """
[ moleculetype ]
;name            nrexcl
UNL              3

[ atoms ] 
;   nr    type   resi  res    atom  cgnr   charge    mass       
"""

  # write the atoms
  for i, atom in enumerate(mol.atoms):
    fcontent += "%6d   att_%03d   1   UNL    %s%03d   1    %.4f   %7.4f\n" % (i+1, i+1, etab.GetSymbol(atom.atomicnum), i+1, q[i], atom.atomicmass)

  # write the bonds
  fcontent += """
[ bonds ]
;   ai     aj funct     r          k
"""
  for bnd in dfrBonds:
    ai, aj = [int(x) for x in bnd.split()]
    fcontent += "%6d %6d   1    %.4f   %.4f\n" % (ai, aj, a2nm(dfrBonds[bnd][1]), cal2j(dfrBonds[bnd][0])*200.0)

  # write the angles
  fcontent += """
[ angles ]
;   ai     aj     ak funct    theta       cth
"""
  for ang in dfrAngles:
    ai, aj, ak  = [int(x) for x in ang.split()]
    fcontent += "%6d %6d %6d   1    %.4f   %.4f\n" % (ai, aj, ak, dfrAngles[ang][1], cal2j(dfrAngles[ang][0])*2.0)

  # write the proper dihedrals
  fcontent += """
[ dihedrals ]
; proper dihedrals - converted to the RB form from Fourier type if OPLS
;   ai     aj     ak     al   func    params
"""

  fimp = ""
  for dih in dfrDihedrals:
    ai, aj, ak, al = [int(x) for x in dih.split()]
    if dfrDihedrals[dih][0].lower() == "amber":
      # check if it's a proper or improper dihedral
      bondIterator = openbabel.OBMolBondIter(mol.OBMol)
      cnt = 0
      for bond in bondIterator:
        if ((ai == bond.GetBeginAtom().GetId()+1) and (aj == bond.GetEndAtom().GetId()+1)) or ((aj == bond.GetBeginAtom().GetId()+1) and (ai == bond.GetEndAtom().GetId()+1)):
          cnt += 1
        elif ((aj == bond.GetBeginAtom().GetId()+1) and (ak == bond.GetEndAtom().GetId()+1)) or ((ak == bond.GetBeginAtom().GetId()+1) and (aj == bond.GetEndAtom().GetId()+1)):
          cnt += 1
        elif ((ak == bond.GetBeginAtom().GetId()+1) and (al == bond.GetEndAtom().GetId()+1)) or ((al == bond.GetBeginAtom().GetId()+1) and (ak == bond.GetEndAtom().GetId()+1)):
          cnt += 1

      if cnt == 3:
        func = 9
        fparam = [cal2j(float(x))/2.0 for x in dfrDihedrals[dih][1:4]]
        for i, term in enumerate(fparam,1):
          if term != 0.0:
            fcontent += "%6d %6d %6d %6d      %1d  %6.2f   %9.5f   %d\n" % (ai, aj, ak, al, func, float(dfrDihedrals[dih][i+3]), term, i)
      else:
        func = 4
        fparam = [cal2j(float(x))/2.0 for x in dfrDihedrals[dih][1:4]]
        for i, term in enumerate(fparam,1):
          if term != 0.0:
            fimp += "%6d %6d %6d %6d      %1d  %6.2f   %9.5f   %d\n" % (ai, aj, ak, al, func, float(dfrDihedrals[dih][i+3]), term, i)

    elif dfrDihedrals[dih][0].lower() == "opls":
      fparam = [cal2j(float(x)) for x in dfrDihedrals[dih][1:4]]
      c0 = fparam[1] + 0.5*(fparam[0]+fparam[2])
      c1 = 0.5*(-fparam[0]+3.0*fparam[2])
      c2 = -fparam[1]
      c3 = -2.0*fparam[2]
      c4 = 0.0
      c5 = 0.0
      fcontent += "%6d %6d %6d %6d      3  %9.5f   %9.5f   %9.5f   %9.5f   %9.5f   %9.5f\n" % (ai, aj, ak, al, c0, c1, c2, c3, c4, c5)
    else:
      print("Error: Dihedral type (%s) found for dihedral %s in .dfr is not valid." % (dfrDihedrals[dih][0], dih))
      sys.exit(0)

  # write the improper dihedrals
  if fimp or dfrImpDih:
    fcontent += """
[ dihedrals ]
; improper dihedrals
;   ai     aj     ak     al   func    params
"""
  if fimp:
    fcontent += fimp

  for idih in dfrImpDih:
    ai, aj, ak, al = [int(x) for x in idih.split()]
    fcontent += "%6d %6d %6d %6d      4  %6.2f   %9.5f   %d\n" % (ai, aj, ak, al, float(dfrImpDih[idih][1]), cal2j(float(dfrImpDih[idih][0])), 2)

  # write the pairs
  fcontent += """
[ pairs ]
"""
  pairs = get_pairs(mol)
  for pair in pairs:
    ai, aj = [int(x) for x in pair.split()]
    fcontent += "%6d %6d    1\n" % (ai, aj)

  return fcontent


def gen_top(ffname, itpname):
  fcontent = """
; topology generated with dice2gromacs
; https://github.com/hmcezar/dicetools
; using the defaults for %s
""" % ffname

  if ffname == "amber":
    fcontent += """
[ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
1               2               yes             0.5     0.8333

; include the molecule .itp
#include "%s"

[ system ]
 %s simulation

[ molecules ]
; compound   nmol
 UNL         1
""" % (itpname, os.path.splitext(os.path.basename(itpname))[0])
  elif ffname == "opls":
    fcontent += """
[ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
1               3               yes             0.5     0.5

; include the molecule .itp
#include "%s"

[ system ]
 %s simulation

[ molecules ]
; compound   nmol
 UNL         1
""" % (itpname, os.path.splitext(os.path.basename(itpname))[0])
  else:
    print("Error: Invalid force field (%s). Select opls or amber." % ffname)
    sys.exit(0)

  return fcontent


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Receives the DICE input files .dfr and .txt to generate the GROMACS input files .itp and .gro.")
  parser.add_argument("dfrfile", type=extant_file, help="the DICE .dfr")
  parser.add_argument("txtfile", type=extant_file, help="the DICE .txt")
  parser.add_argument("force_field", help='select either "opls" or "amber" to generate the inputs accordingly')

  args = parser.parse_args()

  if args.force_field.lower() not in ["opls", "amber"]:
    print("Error: Invalid force field (%s). Select opls or amber." % args.force_field)
    sys.exit(0)

  # read the geometry into a pybel mol, the LJ+Coulomb into lists and the dfr paramaters into dictionaries
  mol, q, eps, sig, dfrBonds, dfrAngles, dfrDihedrals, dfrImpDih = read_parameters(args.dfrfile, args.txtfile)

  # generate .itp and write
  fcontent = itp_from_params(mol, q, eps, sig, dfrBonds, dfrAngles, dfrDihedrals, dfrImpDih)
  with open(os.path.splitext(args.dfrfile)[0]+".itp", "w") as f:
    f.write(fcontent)

  # generate .top and write
  fcontent = gen_top(args.force_field, os.path.splitext(args.dfrfile)[0]+".itp")
  with open("topology.top", "w") as f:
    f.write(fcontent)

  # get the molecule diameter
  mol.OBMol.Center()
  maxd = 0.0
  for atom in mol:
    x, y, z = atom.coords
    r = sqrt(x*x + y*y + z*z)
    if r > maxd:
      maxd = r

  # modify the unit cell and translate the molecule to the center of box
  boxsize = maxd + 100.0
  cell = openbabel.OBUnitCell()
  cell.SetData(boxsize, boxsize, boxsize, 90.0, 90.0, 90.0)
  mol.OBMol.CloneData(cell)

  disp = boxsize / 2.0
  arr = openbabel.vector3(openbabel.double_array([disp, disp, disp]))
  mol.OBMol.Translate(arr)

  # write .gro
  mol.write("gro", os.path.splitext(args.txtfile)[0]+".gro", overwrite=True)
