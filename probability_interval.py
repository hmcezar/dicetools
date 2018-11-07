"""
Given a file with raw data (a number per line), search for
values in the interval [minval, maxval] and give the probability
of getting a number in it.

Author: Henrique Musseli Cezar
Date: JUL/2016
"""

import argparse

def get_prob_interval(data, minval, maxval):
	# get the how many values are in the interval
	count = 0
	for value in data:
		if ((value >= minval) and (value <= maxval)):
			count += 1
	return float(count)/len(data)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Receive file with raw data and get the probability of having one number in the interval [minval, maxval].")
	parser.add_argument("filename", help="the data file with one number in each line")
	parser.add_argument("minval", help="minimum of the interval")
	parser.add_argument("maxval", help="maximum of the interval")
	parser.add_argument("modify", nargs='?', help="transform data to make torsionals between [0,360)? Default is False.", default=False)
	args = parser.parse_args()

	data = []
	transform = bool(args.modify)
	with open(args.filename, 'r') as f:
		for line in f.readlines():
			value = float(line)
			if (transform and value < 0.0):
				data.append(value+360.0)
			else:
				data.append(value)

	prob = get_prob_interval(data, float(args.minval), float(args.maxval))

	print("The probability of getting one value in the interval [%f, %f] is %.2f %%" % (float(args.minval), float(args.maxval), prob*100.0))