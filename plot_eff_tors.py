"""
Receives a DICE ".txt" and a DICE ".dfr" to calculate the energy with a dihedral specified by the reference atoms a1 a2 a3 a4.

Author: Henrique Musseli Cezar
Date: 08/JUN/2018
"""

import argparse
import sys
import os
import numpy as np
from numpy import cos
from numpy import sin
from numpy import sqrt
import matplotlib as mpl
# Force matplotlib to not use any Xwindows backend.
mpl.use('Agg')
import matplotlib.pyplot as plt
try:
  from Queue import Queue
except:
  from queue import Queue

# Coulomb constant with charge in AKMA units
CT_e = 18.2257

# Dictionary to convert between atomic number and symbols
atomsymbols = {
    1:' H', 2:'He', 3:'Li', 4:'Be', 5:' B', 6:' C', 7:' N', 8:' O', 9:' F', 10:'Ne', 11:'Na', 12:'Mg', 
   13:'Al', 14:'Si', 15:' P', 16:' S', 17:'Cl', 18:'Ar', 19:' K', 20:'Ca', 21:'Sc', 22:'Ti', 23:' V', 
   24:'Cr', 25:'Mn', 26:'Fe', 27:'Co', 28:'Ni', 29:'Cu', 30:'Zn', 31:'Ga', 32:'Ge', 33:'As', 34:'Se', 
   35:'Br', 36:'Kr', 37:'Rb', 38:'Sr', 39:' Y', 40:'Zr', 41:'Nb', 42:'Mo', 43:'Tc', 44:'Ru', 45:'Rh', 
   46:'Pd', 47:'Ag', 48:'Cd', 49:'In', 50:'Sn', 51:'Sb', 52:'Te', 53:' I', 54:'Xe', 55:'Cs', 56:'Ba', 
   57:'La', 58:'Ce', 59:'Pr', 60:'Nd', 61:'Pm', 62:'Sm', 63:'Eu', 64:'Gd', 65:'Tb', 66:'Dy', 67:'Ho', 
   68:'Er', 69:'Tm', 70:'Yb', 71:'Lu', 72:'Hf', 73:'Ta', 74:' W', 75:'Re', 76:'Os', 77:'Ir', 78:'Pt', 
   79:'Au', 80:'Hg', 81:'Tl', 82:'Pb', 83:'Bi', 84:'Po', 85:'At', 86:'Rn', 87:'Fr', 88:'Ra', 89:'Ac', 
   90:'Th', 91:'Pa', 92:' U', 93:'Np', 94:'Pu', 95:'Am', 96:'Cm', 97:'Bk', 98:'Cf', 99:'Es', 100:'Fm', 
  101:'Md', 102:'No', 103:'Lw', 104:'XX'
     }

# Dictionary containing the atomic mass of each species used to put mol in the CM for dipole moment calculation
atomicmass= {
  1:1.008, 2:4.003, 3:6.939, 4:9.012, 5:10.811, 6:12.011, 7:14.007, 8:15.999, 9:18.998, 10:20.183, 11:22.989, 12:24.312,
  13:26.982, 14:28.086, 15:30.974, 16:32.064, 17:35.453, 18:39.948, 19:39.102, 20:40.080, 21:44.956, 22:47.900, 23:50.942, 
  24:51.996, 25:54.938, 26:55.847, 27:58.933, 28:58.710, 29:63.540, 30:65.370, 31:69.720, 32:72.590, 33:74.922, 34:78.960,
  35:79.909, 36:83.800, 37:85.470, 38:87.620, 39:88.905, 40:91.220, 41:92.906, 42:95.940, 43:98.000, 44:101.070, 45:102.905,
  46:106.400, 47:107.870, 48:112.400, 49:114.820, 50:118.690, 51:121.750, 52:127.600, 53:126.904, 54:131.300, 55:132.905, 
  56:137.340, 57:138.910, 58:140.120, 59:140.907, 60:144.240, 61:147.000, 62:150.350, 63:151.960, 64:157.250, 65:158.924, 
  66:162.500, 67:164.930, 68:167.260, 69:168.934, 70:173.040, 71:174.970, 72:178.490, 73:180.948, 74:183.850, 75:186.200,
  76:190.200, 77:192.200, 78:195.090, 79:196.967, 80:200.590, 81:204.370, 82:207.190, 83:208.980, 84:210.000, 85:210.000, 
  86:222.000, 87:223.000, 88:226.000, 89:227.000, 90:232.038, 91:231.000, 92:238.030, 93:237.000, 94:242.000, 95:243.000,
  96:247.000, 97:247.000, 98:249.000, 99:254.000, 100:253.000, 101:256.000, 102:254.000, 103:257.000, 104:0.000
}

def rotate_point(atom, pt1, pt2, dphi):
  p1 = np.asarray(pt1)
  v = np.subtract(pt2,pt1)
  a = np.asarray(atom)

  norm = np.linalg.norm(v)
  sqrnorm = norm*norm
  dotprod = np.dot(v,a)

  x = ((p1[0]*(v[1]*v[1]+v[2]*v[2])-v[0]*(p1[1]*v[1]+p1[2]*v[2]-dotprod))*(1.0-cos(dphi))+sqrnorm*a[0]*cos(dphi)+norm*(-p1[2]*v[1]+p1[1]*v[2]-v[2]*a[1]+v[1]*a[2])*sin(dphi))/sqrnorm
  y = ((p1[1]*(v[0]*v[0]+v[2]*v[2])-v[1]*(p1[0]*v[0]+p1[2]*v[2]-dotprod))*(1.0-cos(dphi))+sqrnorm*a[1]*cos(dphi)+norm*(p1[2]*v[0]-p1[0]*v[2]+v[2]*a[0]-v[0]*a[2])*sin(dphi))/sqrnorm
  z = ((p1[2]*(v[0]*v[0]+v[1]*v[1])-v[2]*(p1[0]*v[0]+p1[1]*v[1]-dotprod))*(1.0-cos(dphi))+sqrnorm*a[2]*cos(dphi)+norm*(-p1[1]*v[0]+p1[0]*v[1]-v[1]*a[0]+v[0]*a[1])*sin(dphi))/sqrnorm

  return [x,y,z]

def get_phi(a1, a2, a3, a4):
  rij = np.subtract(a2,a1)
  rjk = np.subtract(a3,a2)
  rlk = np.subtract(a3,a4)

  m = np.cross(rij, rjk)
  n = np.cross(rlk, rjk)

  normm = np.linalg.norm(m)
  normn = np.linalg.norm(n)
  normrjk = np.linalg.norm(rjk)

  return np.arctan2((np.dot(n,rij)*normrjk)/(normm*normn),np.dot(m,n)/(normm*normn))

def energy_tors_amber(params, phi):
  return 0.5*(params[0]*(1.0 + cos(phi)) + params[1]*(1.0 + cos(2.0*phi)) + params[2]*(1.0 + cos(3.0*phi)))

def energy_tors(params, phi):
  return 0.5*(params[0]*(1.0 + cos(phi)) + params[1]*(1.0 - cos(2.0*phi)) + params[2]*(1.0 + cos(3.0*phi)))

def energy_nonbonded(sigma,epsilon,q1,q2,fclb,flj,r):
  sigr = sigma/r
  sigrthree = sigr*sigr*sigr
  sigrsix = sigrthree*sigrthree
  return fclb*((q1*q2)/r) + flj*epsilon*(sigrsix*sigrsix-sigrsix)

def calculate_dipole(atomSp, atomsCoord, atomsNB):
  # put in the center of mass
  xcm = 0.0
  ycm = 0.0
  zcm = 0.0
  tmas = 0.0
  for i, atomr in atomsCoord.items():
    mass = atomicmass[int(atomSp[i])]
    tmas += mass
    xcm += mass * atomr[0]
    ycm += mass * atomr[1]
    zcm += mass * atomr[2]
  xcm /= tmas
  ycm /= tmas
  zcm /= tmas
  centeredCoord = []
  for i, atomr in atomsCoord.items():
    centeredCoord.append([atomr[0]-xcm, atomr[1]-ycm, atomr[2]-zcm])

  # get the dipole moment
  tdip = np.zeros(3)
  for i, atomr in enumerate(centeredCoord,1):
    dip = [atomsNB[i][0]/CT_e * x for x in atomr]
    tdip = np.add(tdip, dip)

  # return the value already in Debyes
  return np.linalg.norm(tdip/0.20819434)

def get_fnb(connInfo, natoms, useamber):
  infty = -1

  fclb = np.zeros((natoms,natoms))
  flj = np.zeros((natoms,natoms))

  # run BFS for each atom, finding the distance with the others
  for atom in connInfo.keys():
    q = Queue()

    vis = [False] * natoms
    dist = [infty] * natoms

    vis[atom-1] = True
    dist[atom-1] = 0
    q.put(atom)

    while not q.empty():
      k = q.get()
      # mark children as visited and queue them
      for child in connInfo[k]:
        if not vis[child-1]:
          vis[child-1] = True
          dist[child-1] = dist[k-1] + 1
          q.put(child)

    for i in range(natoms):
      if (dist[i] == 3):
        if useamber:
          fclb[atom-1][i] = 1.0/1.2
        else:
          fclb[atom-1][i] = 0.5
        flj[atom-1][i] = 0.5
      elif (dist[i] > 3):
        fclb[atom-1][i] = 1.0
        flj[atom-1][i] = 1.0
      elif (dist[i] == infty):
        print("Warning: disconnected atom in .dfr input")
        fclb[atom-1][i] = 1.0
        flj[atom-1][i] = 1.0

  return fclb, flj

def str_to_bool(s):
  if s.lower() == 'true':
    return True
  elif s.lower() == 'false':
    return False
  elif s.lower() == 't':
    return True
  elif s.lower() == 'f':
    return False
  else:
    raise ValueError("Cannot covert {} to a bool".format(s))

def get_dist(p1, p2):
  a = np.array(p1)
  b = np.array(p2)
  return np.linalg.norm(a-b)

def find_rigid_parts(frags, fconn, a1, a2, a3, a4):
  # check wich fragments have a2 and a3
  fragsParts = []
  for idx, frag in frags.items():
    if (a2 in frag) and (a3 in frag):
      fragsParts.append(idx)

  if len(fragsParts) > 2:
    print("You should have atoms %d and %d belonging just to two fragments" % (a2, a3))
    sys.exit(0)

  # based on the fragment connections, build the two parts running a BFS exploration
  pt1 = frags[fragsParts[0]]
  pt2 = frags[fragsParts[1]]

  visited1 = [fragsParts[0], fragsParts[1]]
  visited2 = [fragsParts[0], fragsParts[1]]

  q = Queue()

  # connected to frag1
  q.put(fragsParts[0])

  while not q.empty():
    fref = q.get()
    for frag in fconn[fref]:
      if frag in visited1:
        continue
      else:
        q.put(frag)
        visited1.append(frag)


      for atom in frags[frag]:
        if atom not in pt1:
          pt1.append(atom)


  # connected to frag2
  q.put(fragsParts[1])

  while not q.empty():
    fref = q.get()
    for frag in fconn[fref]:
      if frag in visited2:
        continue
      else:
        q.put(frag)
        visited2.append(frag)

      for atom in frags[frag]:
        if atom not in pt2:
          pt2.append(atom)

  if a1 in pt2:
    aux = pt2
    pt2 = pt1
    pt1 = aux

  return pt1, pt2

def get_potential_curve(txtfile, dfrfile, ab1, ab2, ab3, ab4, npoints, base, printxyz, useamber, gausstop, gaussbot):

  # put gausstop file contents into a string
  if gausstop:
    # check the number of lines to append to the comment the dihedral angle
    topfile = ""
    with open(gausstop, 'r') as f:
      for i, line in enumerate(f,1):
        pass

      f.seek(0)

      for j, line in enumerate(f,1):
        if (j == i-2):
          topfile+=line.rstrip()+" dihedral = ANGLEPLACEHOLDER\n"
        else:
          topfile+=line
      
  # put gaussbot file into a string
  if gaussbot:
    with open(gaussbot, 'r') as f:
      botfile = f.read()

  # read the molecule and nonbonded parameters
  atomSp = {}
  atomsCoord = {}
  nbParams = {}
  with open(txtfile,'r') as f:
    line = f.readline()
    if (line.strip() == '*'):
      mult = True
    else:
      mult = False
    line = f.readline()
    if (int(line.strip()) > 1):
      print("Your .txt file is supposed to contain only one molecule. Aborting..")
      sys.exit(0)
    natoms = int(f.readline().split()[0])
    anum = 1
    for line in f:
      if (line.strip().lower() == "$end"): break
      atomSp[anum] = line.split()[1]
      atomsCoord[anum] = [float(x) for x in line.split()[2:5]]
      nbParams[anum] = [float(x) for x in line.split()[5:]]
      anum += 1

  # read the intramolecular parameters, connection info (bonds) and frag info
  potentialDict = {}
  connInfo = {}
  fragInfo = {}
  fconnInfo = {}
  with open(dfrfile,'r') as f:
    # read the fragments
    line = f.readline()
    while (line.strip().lower() != "$atoms fragments"):
      line = f.readline()
    fnum = 1
    line = f.readline()
    while (line.strip().lower() != "$end atoms fragments"):
      fragInfo[fnum] = [int(x) for x in line.split()[1:-1]]
      fnum += 1
      line = f.readline()

    # read the fragment connections
    while (line.strip().lower() != "$fragment connection"):
      line = f.readline()
    fnum = 1
    line = f.readline()
    while (line.strip().lower() != "$end fragment connection"):
      frag1 = int(line.split()[0])
      frag2 = int(line.split()[1])

      if (frag1 not in fconnInfo.keys()):
        fconnInfo[frag1] = [frag2]
      else:
        fconnInfo[frag1].append(frag2)

      if (frag2 not in fconnInfo.keys()):
        fconnInfo[frag2] = [frag1]
      else:
        fconnInfo[frag2].append(frag1)

      line = f.readline()

    # read the bonds to get the connections later (needed for nonbonded interaction)
    line = f.readline()
    while (line.strip().lower() != "$bond"):
      line = f.readline()
    line = f.readline()
    while (line.strip().lower() != "$end bond"):
      atom1 = int(line.split()[0])
      atom2 = int(line.split()[1])
      if (atom1 not in connInfo.keys()):
        connInfo[atom1] = [atom2]
      else:
        connInfo[atom1].append(atom2)

      if (atom2 not in connInfo.keys()):
        connInfo[atom2] = [atom1]
      else:
        connInfo[atom2].append(atom1)

      line = f.readline()

    # based on the fragments, atomic connections and reference dihedral find the two rigid parts to be moved
    fpt1, fpt2 = find_rigid_parts(fragInfo, fconnInfo, ab1, ab2, ab3, ab4)

    # read the dihedral constants
    line = f.readline()
    while (line.strip().lower() != "$dihedral"):
      line = f.readline()
    dnum = 1
    line = f.readline()
    while (line.strip().lower() != "$end dihedral"):
      if (int(line.split()[1]) == ab2 and int(line.split()[2]) == ab3) or (int(line.split()[1]) == ab3 and int(line.split()[2]) == ab2):
        potentialDict[dnum] = [int(line.split()[0]),int(line.split()[1]),int(line.split()[2]),int(line.split()[3]),float(line.split()[5]),float(line.split()[6]),float(line.split()[7]),float(line.split()[8]),float(line.split()[9]),float(line.split()[10])]
      dnum += 1
      line = f.readline()

  # get the nonbonded factors (0.0, 0.5 or 1.0), which interactions will be evaluated and adjust the constants
  fclb, flj = get_fnb(connInfo, natoms, useamber)

  # the nonzero terms are the ones which we will use to calculate an interaction
  # using a diagonal matrix to not account twice during energy calculations
  ints = np.zeros((natoms,natoms))
  for i in range(natoms-1):
    for j in range(i+1,natoms):
      ints[i][j] = flj[i][j]
  nzl, nzc = ints.nonzero()

  for atom in nbParams.keys():
    nbParams[atom][0] = nbParams[atom][0] * CT_e
    nbParams[atom][1] = sqrt(nbParams[atom][1]) * 2.0
    if (mult):
      nbParams[atom][2] = sqrt(nbParams[atom][2])
    else:
      nbParams[atom][2] = nbParams[atom][2]/2.0

  # get the common bond atom ids and coordinates
  abcoord1 = atomsCoord[ab2]
  abcoord2 = atomsCoord[ab3]

  # angle of the first dihedral, to use as reference
  cphi = get_phi(atomsCoord[ab1], abcoord1, abcoord2, atomsCoord[ab4])
  print("# Initial angle (degrees): %f" % (180.0*cphi/np.pi))

  # calculate the angle increment based on the total number of points to be calculated
  dphi = 2.0*np.pi/npoints
  
  # open xyz output for the trajectory if needed
  if (printxyz):
    fxyz = open(base+'_rotations.xyz','w')

  # open gjf output if needed
  if (gausstop):
    fgjf = open(base+'_scan.gjf','w')

  # now loop changing the angles, calculating energies and storing them properly
  angles = []
  died_energies = []
  nb_energies = []
  dipoles = []
  for i in range(npoints):
    # rotate the atoms of the second fragment
    for atom in fpt2:
      if (atom == ab2) or (atom == ab3):
        continue
      atomsCoord[atom] = rotate_point(atomsCoord[atom], abcoord1, abcoord2, dphi)

    en_tors = 0.0
    # calculate the torsional
    for idx, died in potentialDict.items():
      phi = get_phi(atomsCoord[died[0]],atomsCoord[died[1]],atomsCoord[died[2]],atomsCoord[died[3]])
      if useamber:
        en_tors += energy_tors_amber(died[4:],phi)        
      else:
        en_tors += energy_tors(died[4:],phi)

    # calculate the nonbonded contributions
    en_nb = 0.0
    for i in range(len(nzl)):
      atom1 = nzl[i] + 1
      atom2 = nzc[i] + 1
      # get the distance
      r = get_dist(atomsCoord[atom1], atomsCoord[atom2])
      # determine the interaction constants
      epsi = nbParams[atom1][1] * nbParams[atom2][1]
      if mult:
        sigm = nbParams[atom1][2] * nbParams[atom2][2] 
      else:
        sigm = nbParams[atom1][2] + nbParams[atom2][2] 

      # sum the energy
      en_nb += energy_nonbonded(sigm, epsi, nbParams[atom1][0], nbParams[atom2][0], fclb[atom1-1][atom2-1], flj[atom1-1][atom2-1], r)

    # calculate dipole moment
    dipm = calculate_dipole(atomSp, atomsCoord, nbParams)

    # add angle change to initial angle
    cphi += dphi

    # store the data
    angles.append(cphi)
    died_energies.append(en_tors)
    nb_energies.append(en_nb)
    dipoles.append(dipm)

    # print rotations if needed
    if (printxyz):
      fxyz.write("%d\nDihedral = %f\n"%(natoms,shift_angle(180.0*cphi/np.pi)))
      for i in range(1,natoms+1):
        fxyz.write("%s\t%f\t%f\t%f\n"%(atomsymbols[int(atomSp[i])],atomsCoord[i][0],atomsCoord[i][1],atomsCoord[i][2]))

    # print to .gjf
    if (gausstop):
      fgjf.write(topfile.replace("ANGLEPLACEHOLDER",str(shift_angle(180.0*cphi/np.pi))))
      for i in range(1,natoms+1):
        fgjf.write(" %s\t%f\t%f\t%f\n"%(atomsymbols[int(atomSp[i])],atomsCoord[i][0],atomsCoord[i][1],atomsCoord[i][2]))
      if gaussbot:
        fgjf.write(botfile)
      fgjf.write("\n--link1--\n")

  if (printxyz):
    fxyz.close()

  if (gausstop):
    fgjf.close()

  return angles, died_energies, nb_energies, dipoles

def shift_angle_pos(tetha):
  if tetha < 0.0:
    return tetha+360.0
  elif tetha >= 360.0:
    return tetha-360.0
  else:
    return tetha

def shift_angle(tetha):
  if tetha < 0.0:
    return tetha
  elif tetha >= 180.0:
    return tetha-360.0
  else:
    return tetha

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Receives a DICE ".txt" and a DICE ".dfr" to calculate the energy with a dihedral specified by the reference atoms a1 a2 a3 a4.')
  parser.add_argument("txtfile", help="DICE's .txt file containing the molecule")
  parser.add_argument("dfrfile", help="DICE's .dfr file containing the force field constants fragment information")
  parser.add_argument("a1", help="first atom defining the reference dihedral")
  parser.add_argument("a2", help="second atom defining the reference dihedral")
  parser.add_argument("a3", help="third atom defining the reference dihedral")
  parser.add_argument("a4", help="fourth atom defining the reference dihedral")
  parser.add_argument("npoints", nargs='?', help="number of points used to build the potential curve - default is 50", default=50)
  parser.add_argument("-o", "--output", help="base name for output files")
  parser.add_argument("--printxyz", help="print the trajectory of the rotations", action="store_true")
  parser.add_argument("--amber", help="use AMBER rule to 1-4 interactions and torsional energy", action="store_true")
  parser.add_argument("--gausstop", help="generate a .gjf input with each configuration using topfile passed as argument to this option")
  parser.add_argument("--gaussbot", help="uses the file passed as argument to this option in the end of .gjf before linking the next input")
  parser.add_argument("--shiftangles", help="shift angles to [0,360)", action="store_true")
  parser.add_argument("--shift-min", help="find the minimum of the total energy and shift it to zero. The nonbonded and torsional are shifted based on the angle of the total energy", action="store_true")

  args = parser.parse_args()

  if args.gaussbot and not args.gausstop:
    print("A Gaussian bottom file should always be used with a Gaussian top file.")
    sys.exit(0)

  # define names based on the basename
  if args.output:
    base = args.output
  else:
    base = "tors_"+args.a1+"-"+args.a2+"-"+args.a3+"-"+args.a4

  fout = open(base+".out", "w")
  fdihed = base+"_torsion.pdf"
  fnb = base+"_nonbonded.pdf"
  ftotal = base+"_total.pdf"

  command = ""
  for i in range(len(sys.argv)):
    command += sys.argv[i]+" "
  print("# Running command: %s\n#" % command)
  fout.write("# Running command: %s\n#\n" % command)

  phi, tors_v, nb_v, dip = get_potential_curve(args.txtfile, args.dfrfile, int(args.a1), int(args.a2), int(args.a3), int(args.a4), int(args.npoints), base, args.printxyz, args.amber, args.gausstop, args.gaussbot)

  # convert to degrees and put it in [0,360) or in [-180,180)
  degphi = [180.0*x/np.pi for x in phi]
  if (args.shiftangles):
    shiftphi = [shift_angle_pos(x) for x in degphi]
  else:
    shiftphi = [shift_angle(x) for x in degphi]

  # sort to avoid bugs during plot
  osphi, ostorsen = (list(t) for t in zip(*sorted(zip(shiftphi, tors_v))))
  osphi, osnben = (list(t) for t in zip(*sorted(zip(shiftphi, nb_v))))
  osphi, osdip = (list(t) for t in zip(*sorted(zip(shiftphi, dip))))

  # total energy
  toten = []
  for ten, nben in zip(ostorsen,osnben):
    toten.append(ten+nben)

  # shift energies if needed
  if args.shift_min:
    # get the angle from total energy
    min_idx = np.argmin(toten)
    # shift energies based on this angle
    minval = toten[min_idx]
    toten = [x-minval for x in toten]
    minval = ostorsen[min_idx]
    ostorsen = [x-minval for x in ostorsen]
    minval = osnben[min_idx]
    osnben = [x-minval for x in osnben]

  # print output to screen
  print("# Angle in (degrees), energies in (kcal/mol) and dipole moment in (Debye)")
  fout.write("# Angle in (degrees), energies in (kcal/mol) and dipole moment in (Debye)")
  print("# Angle\t\tTotal en\tTors en\t\tNB en\t\tDip mom")
  fout.write("# Angle\t\tTotal en\tTors en\t\tNB en\t\tDip mom\n")
  for ang, ten, nben, dipm in zip(osphi,ostorsen,osnben,osdip):
    print("%f\t%f\t%f\t%f\t%f"%(ang,ten+nben,ten,nben,dipm))
    fout.write("%f\t%f\t%f\t%f\t%f\n"%(ang,ten+nben,ten,nben,dipm))

  fout.close()

  # plotting options
  mpl.rcParams.update({'font.size':18, 'text.usetex':True, 'font.family':'serif', 'ytick.major.pad':4})

  # plot torsional energy
  plt.plot(osphi,ostorsen)
  plt.xlabel(r"$\phi$ ($^\circ$)")
  if (args.shiftangles):
    plt.xlim([0.0,360.0])
    plt.xticks([0,60,120,180,240,300,360])
  else:
    plt.xlim([-180,180])
    plt.xticks([-180,-120,-60,0,60,120,180])
  plt.ylabel(r"Total dihedral energy (kcal/mol)")
  plt.savefig(fdihed, bbox_inches='tight')
  plt.gcf().clear()

  # plot nonbonded energy
  plt.plot(osphi,osnben)
  plt.xlabel(r"$\phi$ ($^\circ$)")
  if (args.shiftangles):
    plt.xlim([0.0,360.0])
    plt.xticks([0,60,120,180,240,300,360])
  else:
    plt.xlim([-180,180])
    plt.xticks([-180,-120,-60,0,60,120,180])
  plt.ylabel(r"Nonbonded energy (kcal/mol)")    
  plt.savefig(fnb, bbox_inches='tight')
  plt.gcf().clear()

  # plot total energy 
  plt.plot(osphi,toten)
  plt.xlabel(r"$\phi$ ($^\circ$)")
  if (args.shiftangles):
    plt.xlim([0.0,360.0])
    plt.xticks([0,60,120,180,240,300,360])
  else:
    plt.xlim([-180,180])
    plt.xticks([-180,-120,-60,0,60,120,180])
  plt.ylabel(r"Nonbonded and dihedral energy (kcal/mol)")    
  plt.savefig(ftotal, bbox_inches='tight')
