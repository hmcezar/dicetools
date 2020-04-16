#!/usr/bin/env python3
"""
Give the filename of data file to plot a curve of the angular distribution of the data.
This is generated from a histogram, from which the number of bins can be passed as an optional argument.

Author: Henrique Musseli Cezar
Date: OCT/2016
"""

import argparse
import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline
from distutils.spawn import find_executable

def get_pdf(data, nbins):
  p, x = np.histogram(data, density = True, bins = nbins)
  x = x[:-1] + (x[1] - x[0])/2   # convert bin edges to centers
  f = UnivariateSpline(x, p, s=0)
  return x, f

def shift_angle(tetha,shift):
  tetha = float(tetha)
  if shift:
    if tetha < 0.0:
      return tetha+360.0
    elif tetha >= 360.0:
      return tetha-360.0
    else:
      return tetha
  else:
    return tetha

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Receives raw data, bin it and plot the angular distribution.")
  parser.add_argument("filename", help="the filename containing the data with each entry in a line")
  parser.add_argument("nbins", nargs='?', help="the number of bins used by the histogram - if no number is given, the default (100) is used", default=100)
  parser.add_argument("--shiftangles", help="shift the angles to [0,360)", action="store_true")
  args = parser.parse_args()

  # put the data in a list
  # data = np.loadtxt(args.filename,converters={0: shift_angle})
  data = []
  with open(args.filename, 'r') as f:
    for line in f:
      data.append(shift_angle(float(line),args.shiftangles))

  # make a name to the pdf file
  basename = os.path.splitext(args.filename)[0]
  if not os.path.exists(basename+".pdf"):
    pdfname = basename+".pdf"
  else:
    n = 2
    while os.path.exists(basename+"_%02d.pdf"%n):
      n += 1
    pdfname = basename+"_%02d.pdf"%n

  # calculate the pdf (from a histogram and interpolating)
  x, pdf = get_pdf(data, int(args.nbins))

  # write data to file for further use
  with open('pdf.dat','w') as f:
    for i, v in enumerate(x):
      f.write("%f\t%f\n" % (v, pdf(v)))

  # plot
  if find_executable('latex') and find_executable('dvipng'):
    mpl.rcParams.update({'font.size':18, 'text.usetex':True, 'font.family':'serif', 'ytick.major.pad':4})
  else:
    mpl.rcParams.update({'font.size':18, 'font.family':'serif', 'ytick.major.pad':4})    

  plt.plot(x,pdf(x))
  if args.shiftangles:
    plt.xlim([0.0,360.0])
    plt.xticks([0,60,120,180,240,300,360])
  else:
    plt.xlim([-180.0,180.0])
    plt.xticks([-180,-120,-60,0,60,120,180])    

  plt.ylim([0.0,0.06])
  plt.xlabel(r"$\phi$ ($^\circ$)")
  plt.ylabel(r"Probability density function")
  plt.savefig(pdfname, bbox_inches='tight')
