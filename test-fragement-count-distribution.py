#!/usr/bin/env python
# coding=utf8

"""
%prog <INT> <INT> [options]

Creates a histogram of <INT> fragments count distribution of <INT> replicates for two samples, for all samples, for sample1 and sample2 and creates an MA-plot based on the supplied beta-values

@author:  Thomas Eder
"""

from __future__ import print_function
from __future__ import division

import sys
import math
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from builtins import range
from future.utils import lrange
from optparse import OptionParser
from random import random, randint
from scipy.stats import beta, laplace
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

def _scale_laplace( nfrags, read_frac, scale):
	"""Scale beta result based on Laplace distribution"""
	loc = 0.5 #fixed
	
	top_nfrgs_scale = laplace.pdf(0.5, loc, scale)#make sure at postion 0.5 is 100%
	scale_factor = 1 / top_nfrgs_scale
	nfrags_scale = laplace.pdf(read_frac, loc, scale) *scale_factor #the bigger the difference from beta the smaller the sample...
	
	if(nfrags_scale == 0):
		nfrags = nfrags 
	else:
		nfrags = int(nfrags *nfrags_scale) #new fragment scaling based on beta result
		
	return nfrags

def _scale_lognorm( nfrags, read_frac, sigma, loc, scale):
	"""Scale number_frags result based on Lognorm distribution"""
	#loc=scale*-1
	#mu=1
	#pos=100
	
	top_beta_scale = sp.stats.lognorm.pdf(1,s=sigma,loc=loc, scale=scale)#to get to 100% at 1
	scale_factor=1/top_beta_scale #this is the factor to set beta_frac to 100% at 1
	beta_scale = sp.stats.lognorm.pdf(nfrags,s=sigma,loc=loc, scale=scale) *scale_factor
	
	#sigma and scale
	#beta_scale = sp.stats.lognorm.pdf(nfrags,s=sigma,loc=loc, scale=scale)
	#mu and sigma
	#beta_scale = sp.stats.lognorm.pdf(nfrags,s=sigma,loc=loc, scale=math.exp(mu))
	
	read_frac_old = read_frac
	
	if(beta_scale == 0):
		read_frac = read_frac 
	else:
		read_frac = ((read_frac-0.5)*beta_scale)+0.5 #center and scale it based on beta_scale from lognorm distribution, we do this to reduce the read_frac towards 0.5 for high values
	
	#print("nfrags: %s read_frac_in: %s scaling_factor: %s read_frac_out: %s" % (nfrags ,read_frac_old, beta_scale, read_frac))
	return read_frac

def _scale_exp(nfrags, read_frac, loc, scale):
	"""Scale number_frags result based on Exponential distribution"""
	#mlambda = scale #jus tot test lambda, rename input variable if we use it...
	
	top_beta_scale = sp.stats.expon.pdf(loc,loc=loc, scale=scale)#to get to 100% at size loc
	scale_factor=1/top_beta_scale #this is the factor to set beta_frac to 100% at 1
	beta_scale = sp.stats.expon.pdf(nfrags,loc=loc, scale=scale) *scale_factor
	
	#beta_scale = sp.stats.expon.pdf(nfrags,loc=loc, scale=scale)
	#beta_scale = sp.stats.expon.pdf(nfrags,loc=loc, scale=1/mlambda)
	
	read_frac_old = read_frac
	
	if(beta_scale == 0):
		read_frac = read_frac 
	else:
		read_frac = ((read_frac-0.5)*beta_scale)+0.5 #center and scale it based on beta_scale from exponential distribution, we do this to reduce the read_frac towards 0.5 for high values
	
	#print("nfrags: %s read_frac_in: %s scaling_factor: %s read_frac_out: %s" % (nfrags ,read_frac_old, beta_scale, read_frac))
	return read_frac

#main
if __name__ == '__main__':
	#handle input arguments
	parser = HelpfulOptionParser(usage=__doc__)
	
	parser.add_option("--frag-count-sh", default=2.2, dest="frag_count_sh", type="float", help="Shape of gamma distribution for fragment counts [default: %default]")
	parser.add_option("--frag-count-sc", default=20.1, dest="frag_count_sc", type="float", help="Scale of gamma distribution for fragment counts [default: %default]")
	
	parser.add_option("--frag-count-op", default=0.01, dest="frag_count_op", type="float", help="Probability for fragment counts being outliers [default: %default]")
	parser.add_option("--frag-count-om", default=6., dest="frag_count_om", type="float", help="Mean of lognormal distribution for fragment counts of outliers [default: %default]")
	parser.add_option("--frag-count-os", default=0.5, dest="frag_count_os", type="float", help="Sigma of lognormal distribution for fragment counts of outliers [default: %default]")
	
	parser.add_option("--frag-count-scaling", default="none", dest="frag_count_scaling", type="string", help="Scaling of fragment distribution, no scaling, scaling of beta result based on fragment counts (with exp) or scaling of fragment counts based on beta result (with laplace) : none , frag , beta [default: %default]")
	parser.add_option("--frag-count-lp-scale", default=0.1, dest="frag_count_lp_sc", type="float", help="Scale for Laplace distribution if frag-count-scaling is frag [default: %default]")
	parser.add_option("--frag-count-ex-loc", default=10, dest="frag_count_ex_lo", type="float", help="Loc for exponential distribution if frag-count-scaling is beta [default: %default]")
	parser.add_option("--frag-count-ex-scale", default=100, dest="frag_count_ex_sc", type="float", help="Scale for exponential distribution if frag-count-scaling is beta [default: %default]")
	
	parser.add_option("--beta", default=[0.5, 0.5], dest="beta_values", type="string", action='callback', callback=_callback_list_float, help="Alpha and Beta of Beta-distribution [default: %default]")
	
	(options, args) = parser.parse_args()
	num_args = 2
	
	#read in required input arguments
	if len(args) != num_args:
		parser.error("At minimum %s parameters are requred to run this tool" % num_args)
	nproteins = int(args[0])
	count_rep = int(args[1])
	
	binwidth = 2
	n_replicates_sample1=count_rep
	n_replicates_sample2=count_rep
	
	fragment_counts=[]
	sample_counts=[]
	sumcounts=[]
	
	for _ in lrange(nproteins):
		
		#get number of fragments per protein
		if(np.random.random() <= options.frag_count_op):
			number_frags = np.random.lognormal(mean=options.frag_count_om,sigma=options.frag_count_os,size=1) #if yes, it is taken from the gamma
		else:
			number_frags = np.random.gamma(shape=options.frag_count_sh,scale=options.frag_count_sc,size=1) #if no, from the lognormal distribution
		
		number_frags = max(1, int(number_frags)) #set to 1 or more
		
		#get read fraction for sample1 (sample2 = 1-x)
		read_frac = np.random.beta(options.beta_values[0], options.beta_values[1]) 
		
		#now choose which way we want to select the number of fragments
		if( options.frag_count_scaling == "none"):
			#no scaling, no changes in nfrags and beta
			number_frags = number_frags
			read_frac = read_frac
		elif (options.frag_count_scaling == "beta") :#nfrags scaling via Laplace, beta not changed
			number_frags = _scale_laplace(number_frags, read_frac, options.frag_count_lp_sc)
		elif (options.frag_count_scaling == "frag") :#nfrags scaling via lognorm distribution, number_frags not changed
			#read_frac = _scale_lognorm(number_frags, read_frac, options.frag_count_ln_si,options.frag_count_ln_lo, options.frag_count_ln_sc)
			read_frac = _scale_exp(number_frags, read_frac, options.frag_count_ex_lo, options.frag_count_ex_sc)
		else:
			print("Unknown scaling method, %s, please choose 'none','frag' or 'beta', exiting now" % (frag_count_scaling))
			exit(1)
		
		#we will distribute the fragments to two samples with x replicates
		number_frags = number_frags * (n_replicates_sample1 + n_replicates_sample2) 
		
		sample1_count = 0
		sample2_count = 0
		fragment_count = 0
		
		for _ in lrange(number_frags):
			
			if random() <= read_frac: #add to sample 1 or sample 2 and to list
				sample1_count += 1
			else:
				sample2_count += 1
			
			fragment_count += 1
			
		#save result for this protein
		sample_counts.append([sample1_count,sample2_count])
		sumcounts.append(number_frags)

	
	#MA-plot
	xma=[]
	yma=[]
	counts1=[]
	counts2=[]
	for entry in sample_counts:
		counts1.append(entry[0])
		counts2.append(entry[1])
		if(entry[1] > 0 and entry[0] > 0):
			x=math.log(entry[0]+entry[1]/2,2)#log2 mean
			xma.append(x) 
			y=math.log(entry[0]/entry[1],2)#log2 fold-change
			yma.append(y) 
	
	print("max: sum %s s1 %s s2 %s" % (max(sumcounts),max(counts1), max(counts2)))
	
	#hist and ma-plot
	fixed_range=[0,200]
	fixed_bins=500
	fig, axs = plt.subplots(4)
	
	axs[0].set_title('fragments all samples')
	axs[0].hist(sumcounts,density=False, bins=fixed_bins, range=fixed_range)
	axs[1].set_title('fragments sample1')
	axs[1].hist(counts1,density=False, bins=fixed_bins, range=fixed_range)
	axs[2].set_title('fragments sample2')
	axs[2].hist(counts2,density=False, bins=fixed_bins, range=fixed_range)
	axs[3].set_title('MA-plot')
	axs[3].scatter(xma, yma)
	
	plt.subplots_adjust(hspace = 0.4)
	plt.show()

	#plt.scatter(xma, yma)
	#plt.show()

