"""
Receives a molecular configuration file and create a ".dfr" (Dice FRagment) and ".txt" DICE input files.
The ".txt" contains the atoms (species, xyz coordinates and place for lj params) and the ".dfr" contains the fragment connections and force field parameters.

Author: Henrique Musseli Cezar
Date: JUN/2015
"""

import openbabel
import pybel
import os
import argparse
import sys
import itertools

def split_mol_fragments_daylight(imol):
  # create an OBMol object identical to the original to modificate
  # how-to clone from: http://blueobelisk.shapado.com/questions/how-can-i-clone-copy-an-obmol-or-pybel-molecule
  clone = pybel.ob.OBMol(imol)
  mol = pybel.Molecule(clone)

  # set the Daylight SMARTS pattern to find the rotatable bonds
  rotBondsPattern = pybel.Smarts("[!$(*#*)&!D1]-!@[!$(*#*)&!D1]")

  # find the bonds
  rbonds = rotBondsPattern.findall(mol)

  # store the OBAtom objects of each broken bond
  connectFragsAtomPairs = []
  
  # dictionary with the true atom of each dummy (dummies are the keys)
  dummyAtomCorrespondence = {}

  # for each rotatable bond, add dummy atoms in the ends bind them 
  # to the respective atoms and them break the rotatable bond
  for atom1, atom2 in rbonds:
    currBond = mol.OBMol.GetBond(atom1,atom2)

    dummy1 = mol.OBMol.NewAtom()
    dummy1.SetAtomicNum(currBond.GetBeginAtom().GetAtomicNum()) # dummy atom
    dummy1.SetVector(currBond.GetBeginAtom().GetX(), currBond.GetBeginAtom().GetY(), currBond.GetBeginAtom().GetZ()) # set dummy position

    dummy2 = mol.OBMol.NewAtom()
    dummy2.SetAtomicNum(currBond.GetEndAtom().GetAtomicNum()) # dummy atom
    dummy2.SetVector(currBond.GetEndAtom().GetX(), currBond.GetEndAtom().GetY(), currBond.GetEndAtom().GetZ()) # set dummy position

    # add bonds to the dummy atoms
    mol.OBMol.AddBond(currBond.GetEndAtomIdx(),dummy1.GetIdx(),1)
    mol.OBMol.AddBond(currBond.GetBeginAtomIdx(),dummy2.GetIdx(),1)

    # store the OBAtoms in this bond, makes the dummy correspondence and finally break it
    connectFragsAtomPairs.append([currBond.GetBeginAtom(),currBond.GetEndAtom()])
    dummyAtomCorrespondence[dummy1.GetId()] = currBond.GetBeginAtom().GetId()
    dummyAtomCorrespondence[dummy2.GetId()] = currBond.GetEndAtom().GetId()
    mol.OBMol.DeleteBond(currBond)

  # draw the fragments
  # mol.draw()

  # split the fragments
  fragments = mol.OBMol.Separate()

  # set the fragment indexes
  for i, frag in enumerate(fragments, start=1):
    frag.SetTitle(str(i))

  # store the index of connected pairs
  fragmentConnections = []
  for atom1, atom2 in connectFragsAtomPairs:
    connection = []
    for frag in fragments:
      fragAtomIterator = openbabel.OBMolAtomIter(frag)
      for atom in fragAtomIterator:
        if atom.GetId() == atom1.GetId() or atom.GetId() == atom2.GetId():
          connection.append(frag.GetTitle())
          break
    fragmentConnections.append(connection)

  # return disconnected fragments and the connections list (index of connected fragments) and the dummy to atom correspondence (by index)
  return fragments, fragmentConnections, dummyAtomCorrespondence


def generate_fragfile(filename, outtype, ffparams=None, eqgeom=False):
  # check outtype
  if outtype not in ["flex", "header", "min"]:
    sys.exit('Invalid argument indicating verbosity of .dfr (%s). Use "flex", "header" or "min".' % outtype)

  # get basename and file extension
  base, ext = os.path.splitext(filename)

  # set openbabel file format
  obConversion = openbabel.OBConversion()
  obConversion.SetInAndOutFormats(ext[1:],"xyz")

  # read molecule to OBMol object
  mol = openbabel.OBMol()
  obConversion.ReadFile(mol, filename)

  if ffparams:
    # get atomic labels from pdb
    idToAtomicLabel = {}
    for res in openbabel.OBResidueIter(mol):
      for atom in openbabel.OBResidueAtomIter(res):
        idToAtomicLabel[atom.GetId()] = res.GetAtomID(atom).strip()

    # read force field parameters and store into dictionaries
    labelToSLabel = {}
    charges = {}
    epsilons = {}
    sigmas = {}
    bonds = {}
    angles = {}
    dihedrals = {}
    impropers = {}
    with open(ffparams, 'r') as f:
      line = f.readline()
      # read nb params
      while "$bond" not in line:
        if line.strip().startswith("#") or not line.strip():
          line = f.readline()
          continue

        lbl = line.split()[0]
        charges[lbl] = line.split()[1]
        epsilons[lbl] = line.split()[2]
        sigmas[lbl] = line.split()[3]
        labelToSLabel[lbl] = line.split()[4]

        line = f.readline()

      # read bond params
      line = f.readline()
      while "$angle" not in line:
        if line.strip().startswith("#") or "$end" in line or not line.strip():
          line = f.readline()
          continue

        line = line.replace("–", "-")

        # store the constants for the order of the input and the inverse order
        consts = "\t".join(line.split()[1:])
        bonds[line.split()[0]] = consts
        bonds["-".join(line.split()[0].split("-")[::-1])] = consts

        line = f.readline()

      # read angle params
      line = f.readline()
      while "$dihedral" not in line:
        if line.strip().startswith("#") or "$end" in line or not line.strip():
          line = f.readline()
          continue

        line = line.replace("–", "-")

        # store the constants for the order of the input and the inverse order
        consts = "\t".join(line.split()[1:])
        angles[line.split()[0]] = consts
        angles["-".join(line.split()[0].split("-")[::-1])] = consts

        line = f.readline()

      # read dihedrals
      line = f.readline()
      while "$improper" not in line:
        if line.strip().startswith("#") or "$end" in line or not line.strip():
          line = f.readline()
          continue

        line = line.replace("–", "-")

        # store the constants for the order of the input and the inverse order
        consts = "\t".join(line.split()[1:])
        dihedrals[line.split()[0]] = consts
        dihedrals["-".join(line.split()[0].split("-")[::-1])] = consts

        line = f.readline()

      # read impropers
      line = f.readline()
      while line:
        if line.strip().startswith("#") or "$end" in line or not line.strip():
          line = f.readline()
          continue

        line = line.replace("–", "-")

        # store the constants for the order of the input and the inverse order
        consts = "\t".join(line.split()[1:])
        impropers[line.split()[0]] = consts
        impropers["-".join(line.split()[0].split("-")[::-1])] = consts

        line = f.readline()

    # check if there are unused labels
    for lbl in charges.keys():
      fnd = False
      for i in idToAtomicLabel:
        if lbl == idToAtomicLabel[i]:
          fnd = True
          break
      if not fnd:
        print("!!! WARNING: There are unused atoms in your parameter file (%s) !!!" % lbl)


  # split the molecule
  fragments, fragConnection, dummyToAtom = split_mol_fragments_daylight(mol)

  # dummy atoms ids
  dummyAtoms = dummyToAtom.keys()

  # write molecule to .txt file (passed as ljname to DICE)
  with open(base+".txt","w") as f:
    f.write("*\n1\n")   
    atomToPrint = []
    for frag in fragments:
      fragAtomIterator = openbabel.OBMolAtomIter(frag)
      for atom in fragAtomIterator:
        if atom.GetId() not in dummyAtoms:
          atomToPrint.append(atom)
    # print number of atoms
    f.write(str(len(atomToPrint))+" \t %s (generated with fragGen)\n"%os.path.basename(base))
    # dictionary associating Atomic number with rdf label
    rdfs = {}
    rdf_label = 1
    if ffparams:
      # sort atoms by index and print (this prints the atoms, e.g., in the same order of the xyz input)
      for atom in sorted(atomToPrint, key=lambda atom: atom.GetId()):
        if atom.GetAtomicNum() not in rdfs.keys():
          rdfs[atom.GetAtomicNum()] = str(rdf_label)
          rdf_label += 1
        f.write(rdfs[atom.GetAtomicNum()]+" "+str(atom.GetAtomicNum())+"  \t"+str(atom.GetX())+"      \t"+str(atom.GetY())+"      \t"+str(atom.GetZ())+"      \t"+charges[idToAtomicLabel[atom.GetId()]]+"\t"+epsilons[idToAtomicLabel[atom.GetId()]]+"\t"+sigmas[idToAtomicLabel[atom.GetId()]]+"\n") 
      f.write("$end\n")
    else:
      # sort atoms by index and print (this prints the atoms, e.g., in the same order of the xyz input)
      for atom in sorted(atomToPrint, key=lambda atom: atom.GetId()):
        if atom.GetAtomicNum() not in rdfs.keys():
          rdfs[atom.GetAtomicNum()] = str(rdf_label)
          rdf_label += 1
        f.write(rdfs[atom.GetAtomicNum()]+" "+str(atom.GetAtomicNum())+"  \t"+str(atom.GetX())+"      \t"+str(atom.GetY())+"      \t"+str(atom.GetZ())+"      \t"+"q"+"\t"+"epsilon"+"\t"+"sigma\n") 
      f.write("$end\n")

  # write info to dfr file
  with open(base+".dfr","w") as f:
    # fragments and fragments connections are printed to every outtype
    f.write("$atoms fragments\n")
    fragslst = []
    for frag in fragments:
      f.write(frag.GetTitle()+"\t[ ")
      fragAtomIterator = openbabel.OBMolAtomIter(frag)
      atomlst = [str(dummyToAtom[x.GetId()]+1) if x.GetId() in dummyAtoms else str(x.GetId()+1) for x in fragAtomIterator]
      fragslst.append(atomlst)
      for atom in atomlst:
        f.write(atom+"\t")                    
      if outtype == "min":
        f.write("] R\n")
      else:
        f.write("] F\n")
    f.write("$end atoms fragments\n") 

    f.write("\n$fragment connection\n")
    for frag1, frag2 in fragConnection:
      f.write(frag1+"\t"+frag2+"\n")
    f.write("$end fragment connection\n")

    # bonds are printed to every outtype that is not header, since we need a connection matrix
    if outtype != "header":
      f.write("\n$bond\n")
      bondIterator = openbabel.OBMolBondIter(mol)
      if ffparams:
        for bond in bondIterator:
          try:
            if eqgeom:
              f.write(str(bond.GetBeginAtom().GetId()+1)+" "+str(bond.GetEndAtom().GetId()+1)+"  \t"+bonds[labelToSLabel[idToAtomicLabel[bond.GetBeginAtom().GetId()]]+"-"+labelToSLabel[idToAtomicLabel[bond.GetEndAtom().GetId()]]].split()[0]+"\t"+str("%.6f" % bond.GetLength())+"\n")
            else:
              f.write(str(bond.GetBeginAtom().GetId()+1)+" "+str(bond.GetEndAtom().GetId()+1)+"  \t"+bonds[labelToSLabel[idToAtomicLabel[bond.GetBeginAtom().GetId()]]+"-"+labelToSLabel[idToAtomicLabel[bond.GetEndAtom().GetId()]]]+"\n")
          except KeyError as e:
            print("The parameters for atoms %d %d (%s) was not found in the bonds list\n" % (bond.GetBeginAtom().GetId()+1,bond.GetEndAtom().GetId()+1,e))
            raise                        
      else:
        for bond in bondIterator:
          f.write(str(bond.GetBeginAtom().GetId()+1)+" "+str(bond.GetEndAtom().GetId()+1)+"  \t0.0\t"+str("%.6f" % bond.GetLength())+"\n")
      f.write("$end bond\n")

    # angles are only printed for outtype flex
    if outtype == "flex":
      f.write("\n$angle\n")
      angleIterator = openbabel.OBMolAngleIter(mol)
      if ffparams:
        for angle in angleIterator:
          try:
            if eqgeom:
              atom2 = mol.GetAtomById(angle[0])
              atom1 = mol.GetAtomById(angle[1])
              atom3 = mol.GetAtomById(angle[2])
              aparams = angles[labelToSLabel[idToAtomicLabel[angle[1]]]+"-"+labelToSLabel[idToAtomicLabel[angle[0]]]+"-"+labelToSLabel[idToAtomicLabel[angle[2]]]].split()
              f.write(str(angle[1]+1)+" "+str(angle[0]+1)+" "+str(angle[2]+1)+"   \t"+aparams[0]+"\t"+aparams[1]+"\t"+str("%.6f" % mol.GetAngle(atom1,atom2,atom3))+"\n")
            else:
              f.write(str(angle[1]+1)+" "+str(angle[0]+1)+" "+str(angle[2]+1)+"   \t"+angles[labelToSLabel[idToAtomicLabel[angle[1]]]+"-"+labelToSLabel[idToAtomicLabel[angle[0]]]+"-"+labelToSLabel[idToAtomicLabel[angle[2]]]]+"\n")
          except KeyError as e:
            print("The parameters for atoms %d %d %d (%s) was not found in the angles list\n" % (angle[1]+1,angle[0]+1,angle[2]+1,e))
            raise            
      else:
        for angle in angleIterator:
          # carefully select the atoms to find the angle
          atom2 = mol.GetAtomById(angle[0])
          atom1 = mol.GetAtomById(angle[1])
          atom3 = mol.GetAtomById(angle[2])
          f.write(str(angle[1]+1)+" "+str(angle[0]+1)+" "+str(angle[2]+1)+"   \tharmonic\tK\t"+str("%.6f" % mol.GetAngle(atom1,atom2,atom3))+"\n")
      f.write("$end angle\n")

    # all the dihedrals are printed to outtype flex, but only connection between fragments are printed if outtype is min
    if outtype == "flex":
      f.write("\n$dihedral\n")
      torsionIterator = openbabel.OBMolTorsionIter(mol)
      if ffparams:
        for torsional in torsionIterator:
          # Need to sum 1: http://forums.openbabel.org/Rotable-bonds-tp957795p957798.html
          torsidx = [str(x+1) for x in torsional]
          try:
            f.write(torsidx[0]+" "+torsidx[1]+" "+torsidx[2]+" "+torsidx[3]+"   \t"+dihedrals["-".join([labelToSLabel[idToAtomicLabel[x]] for x in torsional])]+"\n")
          except KeyError as e:
            print("The parameters for atoms %s %s %s %s (%s) was not found in the dihedrals list\n" % (torsidx[0],torsidx[1],torsidx[2],torsidx[3],e))
            raise
      else:
        for torsional in torsionIterator:
          # Need to sum 1: http://forums.openbabel.org/Rotable-bonds-tp957795p957798.html
          torsional = [str(x+1) for x in torsional]
          f.write(torsional[0]+" "+torsional[1]+" "+torsional[2]+" "+torsional[3]+"   \tTYPE\tV1\tV2\tV3\tf1\tf2\tf3\n")
      f.write("$end dihedral\n")

      # improper dihedral = carbon with only 3 atoms connected to it (SP2 hybridization)
      # angle found following this definition --> http://cbio.bmt.tue.nl/pumma/index.php/Theory/Potentials
      f.write("\n$improper dihedral\n")
      atomIterator = openbabel.OBMolAtomIter(mol)
      for atom in atomIterator:
        # print(atom.GetHyb(), atom.GetAtomicNum(), atom.GetValence())
        # if atom.GetAtomicNum() == 6 and atom.GetValence() == 3:
        if atom.GetHyb() == 2 and atom.GetValence() == 3:
          bondIterator = atom.BeginBonds()
          nbrAtom = atom.BeginNbrAtom(bondIterator)
          connectedAtoms = []
          connectedAtoms.append(nbrAtom)
          for i in range(2):
            nbrAtom = atom.NextNbrAtom(bondIterator)
            connectedAtoms.append(nbrAtom)
          if ffparams:
            torsional = [atom.GetId(), connectedAtoms[0].GetId(), connectedAtoms[1].GetId(), connectedAtoms[2].GetId()]
            # create all the permutations to check if one is found
            perms = list(itertools.permutations(torsional[1:]))
            nfound = 0
            for perm in perms:
              try:
                joined = "-".join([labelToSLabel[idToAtomicLabel[torsional[0]]]]+[labelToSLabel[idToAtomicLabel[x]] for x in perm])
                f.write(str(torsional[0]+1)+" "+str(torsional[1]+1)+" "+str(torsional[2]+1)+" "+str(torsional[3]+1)+"    \t"+impropers[joined]+"\n")
              except:
                nfound += 1

            if nfound == len(perms):
              joined = "-".join([labelToSLabel[idToAtomicLabel[torsional[0]]]]+[labelToSLabel[idToAtomicLabel[x]] for x in perms[0]])
              raise KeyError("The key %s (or its permutations) were not found in the improper dihedrals list\n" % (joined))
          else:
            torsional = [atom.GetId()+1, connectedAtoms[0].GetId()+1, connectedAtoms[1].GetId()+1, connectedAtoms[2].GetId()+1]
            torsionAngle = mol.GetTorsion(torsional[0],torsional[1],torsional[2],torsional[3])
            f.write(str(torsional[0])+" "+str(torsional[1])+" "+str(torsional[2])+" "+str(torsional[3])+"    \tV2\t"+str("%.6f" % torsionAngle)+"\n")
      f.write("$end improper dihedral\n")
    
    elif outtype == "min":
      torsionIterator = openbabel.OBMolTorsionIter(mol)
      # tjf = torsionals that join fragments
      tjf = []
      # find the tjfs by checking if all the atoms of a torsional belong to the same fragment
      for tors in torsionIterator:
        tors = [str(x+1) for x in tors]
        istjf = True
        for atomlst in fragslst:
          if (tors[0] in atomlst) and (tors[1] in atomlst) and (tors[2] in atomlst) and (tors[3] in atomlst):
            istjf = False
            break
        if istjf:
          tjf.append(tors)

      f.write("\n$dihedral\n")
      if ffparams:
        for torsidx in tjf:
          torsional = [int(x)-1 for x in torsidx]
          try:
            f.write(torsidx[0]+" "+torsidx[1]+" "+torsidx[2]+" "+torsidx[3]+"   \t"+dihedrals["-".join([labelToSLabel[idToAtomicLabel[x]] for x in torsional])]+"\n")
          except KeyError as e:
            print("The parameters for atoms %s %s %s %s (%s) was not found in the dihedrals list\n" % (torsidx[0],torsidx[1],torsidx[2],torsidx[3],e))
            raise
      else:
        for torsional in tjf:
          f.write(torsional[0]+" "+torsional[1]+" "+torsional[2]+" "+torsional[3]+"   \tTYPE\tV1\tV2\tV3\tf1\tf2\tf3\n")
      f.write("$end dihedral\n")

  # create directory to store the fragments
  if not os.path.exists(base+"_fragments"):
    os.makedirs(base+"_fragments")

  # write framents to the cml files
  for frag in fragments:
    obConversion.WriteFile(frag, os.path.join(base+"_fragments",os.path.basename(filename).split(".")[0]+"_fragment"+frag.GetTitle()+".xyz"))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Receives a molecular structure in an OpenBabel supported file format and creates both the '.dfr' (Dice FRagment) and the '.txt' DICE input files.")
  parser.add_argument("filename", help="file containing the molecular structure")
  parser.add_argument("-p","--force-field-parameters", help="file containing the force field parameters with labels from pdb")
  parser.add_argument("-e","--equilibrium-from-geom", help="when the -p option is used, use the equilibrium values from the geometry and not from the file", action="store_true")

  io_group = parser.add_mutually_exclusive_group()
  io_group.add_argument("--rigid-frags", help="consider the fragments rigid and defines only the fragment connections as flexible (default option)", action="store_true")
  io_group.add_argument("--flexible", help="consider the whole molecule flexible", action="store_true")
  io_group.add_argument("--header", help="print just the header concerning the fragments and fragment connections", action="store_true")  

  args = parser.parse_args()

  if args.equilibrium_from_geom and not args.force_field_parameters:
    sys.exit("The -e option should only be used with -p.")

  obabel_sup = ["acr", "adf", "adfout", "alc", "arc", "bgf", "box", "bs", "c3d1", "c3d2", "cac", "caccrt", "cache", "cacint", "can", "car", "ccc", "cdx", "cdxml", "cht", "cif", "ck", "cml", "cmlr", "com", "copy", "crk2d", "crk3d", "csr", "cssr", "ct", "cub", "cube", "dmol", "dx", "ent", "fa", "fasta", "fch", "fchk", "fck", "feat", "fh", "fix", "fpt", "fract", "fs", "fsa", "g03", "g92", "g94", "g98", "gal", "gam", "gamin", "gamout", "gau", "gjc", "gjf", "gpr", "gr96", "gukin", "gukout", "gzmat", "hin", "inchi", "inp", "ins", "jin", "jout", "mcdl", "mcif", "mdl", "ml2", "mmcif", "mmd", "mmod", "mol", "mol2", "molden", "molreport", "moo", "mop", "mopcrt", "mopin", "mopout", "mpc", "mpd", "mpqc", "mpqcin", "msi", "msms", "nw", "nwo", "outmol", "pc", "pcm", "pdb", "png", "pov", "pqr", "pqs", "prep", "qcin", "qcout", "report", "res", "rsmi", "rxn", "sd", "sdf", "smi", "smiles", "sy2", "t41", "tdd", "test", "therm", "tmol", "txt", "txyz", "unixyz", "vmol", "xed", "xml", "xyz", "yob", "zin"]

  filename = os.path.realpath(args.filename)
  base, ext = os.path.splitext(filename)

  # based on input option define the 'verbosity' of the .dfr
  if args.flexible:
    opt = "flex"
  elif args.header:
    opt = "header"
  else:
    opt = "min"

  if ext[1:].lower() not in obabel_sup:
    sys.exit("The extension of the provided file is not supported by OpenBabel.")

  if (ext[1:].lower() != "pdb") and args.force_field_parameters:
    sys.exit("When the force field parameters are given with '-p' you should have your structure in .pdb with the correct labels.")


  generate_fragfile(filename, opt, args.force_field_parameters, args.equilibrium_from_geom)