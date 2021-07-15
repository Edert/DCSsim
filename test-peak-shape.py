#!/usr/bin/env python
# coding=utf8

"""
%prog <INT> [options]

Creates histograms of <INT> fragment distances and peak shapes

@author:  Thomas Eder
"""

from __future__ import print_function
from __future__ import division

import sys
import numpy as np
import matplotlib.pyplot as plt

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


class Fragment(object):
	"""A fragment in a protein, it has the following attributes:
	
	Attributes:
		domain_name: name or id of the domain
		protein_name:  name or id of the protein
		frag_length: int of the length of the fragment
		rand_shift: int of the shift inside the protein
		position: int of the final position of the fragment
		sequence: string of the sequence of the fragment
		frag_type: string of the fragment type
	"""
	
	def __init__(self, frag_len_mean, frag_len_dev, frag_len_max, frag_dist_on, frag_dist_mn_mean, frag_dist_mn_cov, frag_dist_prob, protein_size):
		"""Return a fragment object"""
		#self.domain_name = domain_name
		#self.protein_name = protein_name
		
		self.frag_length = self._sample_frag_length(frag_len_mean, frag_len_dev, frag_len_max) #get fragment length
		
		if not (frag_dist_on):
			self.rand_shift = randint(-abs(frag_len_mean - protein_size)//2, abs(frag_len_mean - protein_size)//2) #add random shift inside the protein size
		else :
			shift = self._sample_multivariate_normal(protein_size, frag_dist_mn_mean, frag_dist_mn_cov, frag_dist_prob)
			self.rand_shift = int(-abs(frag_len_mean - protein_size)//2) + shift #from most left position on add shift
		
		#self.position = prot_position + self.rand_shift #final position of the fragment
		
		#get sequence
		#self.sequence = genome_seq[self.position - self.frag_length//2 : self.position + self.frag_length//2]
		#self.frag_type = frag_type

	def _sample_frag_length(self, mean, dev, frag_len_max):
		"""Return length of fragment for peaks"""
		return min(abs(int(np.random.normal(loc = mean, scale = dev))), frag_len_max)

	def _sample_multivariate_normal(self, protein_size, frag_dist_mn_mean, frag_dist_mn_cov, frag_dist_prob):
		"""Return a sample from a multivariate normal distribution """
		res = min(max(int(choice(multivariate_normal(frag_dist_mn_mean, frag_dist_mn_cov), p = frag_dist_prob)), 1), protein_size)
		#choose a random element from the n distributions, min = 1 and max = protein_size
		return int(res)

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
	
	parser.add_option("--prot-size", default=150, dest="protein_size", type="int", help="Protein size (this parameter limits the fragment shifts) [default: %default]")
	
	parser.add_option("--frag-len-mean", default=200, dest="frag_len_mean", type="float", help="Set mean of fragments' length [default: %default]")
	parser.add_option("--frag-len-dev", default=20, dest="frag_len_dev", type="float", help="Set deviation of fragments' length [default: %default]")
	parser.add_option("--frag-len-max", default=1000, dest="frag_len_max", type="int", help="Set maximum of fragments' length [default: %default]")
	
	parser.add_option("--frag-dist-on", default=False, action="store_true", dest="frag_dist_on", help="Use multivariate normal distribution for fragment shifts to create peak chapes, shifts are limited by prot-size. The final shift is: postion of peak - prot_size + sampling from distribution  [default: %default]")
	
	parser.add_option("--frag-dist-prob", default=[0.5, 0.5], dest="frag_dist_prob", type="string", action='callback', callback=_callback_list_float,\
	help="Probability fo each of the multivariate normal distributions to be chosen [default: %default]")
	parser.add_option("--frag-dist-muno-mean", default=[20, 100], dest="frag_dist_mn_mean", type="string", action='callback', callback=_callback_list,\
	help="Means of multivariate normal distribution for the shifts of fragments, separator: ',' eg. \"300, 1800\" [default: %default]")
	parser.add_option("--frag-dist-muno-cov", default=[[100,0],[0,500]], dest="frag_dist_mn_cov", type="string", action='callback', callback=_callback_matrix,\
	help="Covariance of multivariate normal distribution for the shifts of fragments, separator: ',' and ';'  eg. \"1000,0;0,5000\" [default: %default]")
	
	(options, args) = parser.parse_args()
	num_args = 1
	
	#read in required input arguments
	n_fragments = int(args[0])
	r_fragshift = []
	r_fraglen = []
	
	for _ in lrange(n_fragments):
		fragment = Fragment(options.frag_len_mean, options.frag_len_dev, options.frag_len_max, options.frag_dist_on, options.frag_dist_mn_mean, options.frag_dist_mn_cov, options.frag_dist_prob, options.protein_size)
		r_fragshift.append(fragment.rand_shift)
		r_fraglen.append(fragment.frag_length)
		
	plt.figure(figsize=(20,10))
	
	plt.subplot(1, 2, 1)
	plt.hist(r_fragshift, bins=50, color = 'g')
	plt.title("Histogram of fragment position shifts")
	plt.ylabel("Count of fragment position shifts")
	plt.xlabel("Range of fragment position shifts")
	
	plt.subplot(1, 2, 2)
	plt.hist(r_fraglen, bins=50, color = 'g')
	plt.title("Histogram of fragment length")
	plt.ylabel("Count of fragment length")
	plt.xlabel("Range of fragment length")
	plt.show()
