"""
Receives a .dat file contaning a temporal series and calculate the probability
of getting values in a interval for each window of the whole set.
Based on this probabilities, a variance is estimated.

Author: Henrique Musseli Cezar
Date: MAR/2018
"""

import argparse
import numpy as np
import sys

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Receives a file with a temporal series, the number of window and interval to consider a value of the category.")
  parser.add_argument("filename", help="the .dat with the temporal series")
  parser.add_argument("nwindows", help="the number windows")
  parser.add_argument("min", help="minimal value to classify as in the category")
  parser.add_argument("max", help="maximal value to classify as in the category")
  parser.add_argument("min2", nargs='?', help="minimal value to classify as in the category")
  parser.add_argument("max2", nargs='?', help="maximal value to classify as in the category")
  args = parser.parse_args()

  if args.min2 and not args.max2:
    print("You should provide either none or both minimum and maximum of the second interval")
    sys.exit(0)

  twoints = False
  if args.min2:
    twoints = True
    min2 = float(args.min2)
    max2 = float(args.max2)

  min1 = float(args.min)
  max1 = float(args.max)

  series = np.loadtxt(args.filename)

  nwin = int(args.nwindows)

  size_win = len(series)/nwin

  nconf_window = []
  for i in range(nwin):
    beg = int(i*size_win)
    end = int(beg+size_win-1)
    array = np.array(series[beg:end])
    total_classified = 0
    total_classified = np.sum(np.logical_and(array>=min1, array<=max1))
    if twoints:
      total_classified += np.sum(np.logical_and(array>=min2, array<=max2))

    nconf_window.append(total_classified)

  probs = [x/size_win for x in nconf_window]
  print("Average of averages: %f and variance: %f" % (np.average(probs), np.var(probs)))
