#!/usr/bin/env python
# coding=utf8

"""
%prog <INT> [options]

Creates a histogram of <INT> fragments count distribution

@author:  Thomas Eder
"""

from __future__ import print_function
from __future__ import division

import sys
import numpy as np
import matplotlib.pyplot as plt

from builtins import range
from future.utils import lrange
from optparse import OptionParser
from random import randint
from numpy.random import multivariate_normal, choice

MAX_TEST = 10000 #number of maximal test iterations, can influence runtime

#classes
class HelpfulOptionParser(OptionParser):
	"""An OptionParser that prints full help on errors."""
	def error(self, msg):
		self.print_help(sys.stderr)
		self.exit(2, "\n%s: error: %s\n" % (self.get_prog_name(), msg))


#helper functions
def _callback_list(option, opt, value, parser): #get a list of options and save it
	try:
		tmplist = list(map(lambda x: int(x), value.split(',')))
		setattr(parser.values, option.dest, tmplist)  
	except Exception as e:
		print("got exception: %r, please check input parameters" % (e,))
		exit()

def _callback_matrix(option, opt, value, parser): #get a matrix of options and save it
	try:
		final_option=[]
		for entry in value.split(';'):
			final_option.append(list(map(lambda x: int(x), entry.split(','))))
		setattr(parser.values, option.dest, final_option)
	except Exception as e:
		print("got exception: %r, please check input parameters" % (e,))
		exit()

def _callback_list_float(option, opt, value, parser): #get a list of float options and save it
	try:
		tmplist = list(map(lambda x: float(x), value.split(',')))
		setattr(parser.values, option.dest, tmplist)  
	except Exception as e:
		print("got exception: %r, please check input parameters" % (e,))
		exit()

#main
if __name__ == '__main__':
	#handle input arguments
	parser = HelpfulOptionParser(usage=__doc__)
	
	parser.add_option("--frag-count-sh", default=2.2, dest="frag_count_sh", type="float", help="Shape of gamma distribution for fragment counts [default: %default]")
	parser.add_option("--frag-count-sc", default=20.1, dest="frag_count_sc", type="float", help="Scale of gamma distribution for fragment counts [default: %default]")
	
	parser.add_option("--frag-count-op", default=0.01, dest="frag_count_op", type="float", help="Probability for fragment counts being outliers [default: %default]")
	parser.add_option("--frag-count-om", default=6., dest="frag_count_om", type="float", help="Mean of lognormal distribution for fragment counts of outliers [default: %default]")
	parser.add_option("--frag-count-os", default=0.5, dest="frag_count_os", type="float", help="Sigma of lognormal distribution for fragment counts of outliers[default: %default]")
	
	(options, args) = parser.parse_args()
	num_args = 1
	
	#read in required input arguments
	nfragments = int(args[0])
	
	binwidth=2
	
	s2 = np.random.gamma(shape=options.frag_count_sh,scale=options.frag_count_sc,size=nfragments)
	is_outlier = np.random.choice([True, False], size=nfragments,p=[options.frag_count_op,1-options.frag_count_op])
	s3 = np.random.lognormal(mean=options.frag_count_om,sigma=options.frag_count_os,size=is_outlier.sum())
	sc = []
	sc.extend(s2)
	sc.extend(s3)
	plt.hist(sc, bins=range(min(sc), max(sc) + binwidth, binwidth))
	plt.show()

