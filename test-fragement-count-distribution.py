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

#main
if __name__ == '__main__':
	#handle input arguments
	parser = HelpfulOptionParser(usage=__doc__)
	
	parser.add_option("--frag-count-sh", default=2.2, dest="frag_count_sh", type="float", help="Shape of gamma distribution for fragment counts [default: %default]")
	parser.add_option("--frag-count-sc", default=20.1, dest="frag_count_sc", type="float", help="Scale of gamma distribution for fragment counts [default: %default]")
	
	parser.add_option("--frag-count-op", default=0.01, dest="frag_count_op", type="float", help="Probability for fragment counts being outliers [default: %default]")
	parser.add_option("--frag-count-om", default=6., dest="frag_count_om", type="float", help="Mean of lognormal distribution for fragment counts of outliers [default: %default]")
	parser.add_option("--frag-count-os", default=0.5, dest="frag_count_os", type="float", help="Sigma of lognormal distribution for fragment counts of outliers [default: %default]")
	
	parser.add_option("--frag-lp-scale", default=0.1, dest="frag_lp_sc", type="float", help="Scale for Laplace distribution [default: %default]")
	parser.add_option("--frag-ln-sigma", default=0.9, dest="frag_ln_si", type="float", help="Sigma for lognorm distribution [default: %default]")
	parser.add_option("--frag-ln-scale", default=100, dest="frag_ln_sc", type="float", help="Scale for lognorm distribution [default: %default]")
	
	parser.add_option("--beta", default=[0.5, 0.5], dest="beta_values", type="string", action='callback', callback=_callback_list_float, help="Alpha and Beta of Beta-distribution [default: %default]")
	
	(options, args) = parser.parse_args()
	num_args = 1
	
	#read in required input arguments
	nproteins = int(args[0])
	count_rep = int(args[1])
	
	binwidth = 2
	n_replicates_sample1=count_rep
	n_replicates_sample2=count_rep
	
	
	#scaling of number of fragments based on beta result
	fragment_counts=[]
	sample_counts=[]
	sumcounts=[]
	
	for _ in lrange(nproteins):
		if(np.random.random() <= options.frag_count_op): #select lognorm or gamma distribution
			nfrags = np.random.lognormal(mean=options.frag_count_om,sigma=options.frag_count_os,size=1) #if yes, it is taken from the gamma
		else:
			nfrags = np.random.gamma(shape=options.frag_count_sh,scale=options.frag_count_sc,size=1) #if no, from the lognormal distribution
		nfrags = max(1, int(nfrags)) 
		
		read_frac = np.random.beta(options.beta_values[0], options.beta_values[1]) #get read fraction for sample1 (sample2 = 1-x)
		
		#laplace distri
		loc = 0.5 #fixed
		scale = options.frag_lp_sc
		
		#make sure at postion 0.5 is 100%
		top_nfrgs_scale = laplace.pdf(0.5, loc, scale)
		scale_factor = 1 / top_nfrgs_scale
		
		nfrags_scale = laplace.pdf(read_frac, loc, scale) *scale_factor #the bigger the differnce from beta the smaller the sample...
		if(nfrags_scale == 0):
			nfrags = nfrags 
		else:
			nfrags = int(nfrags *nfrags_scale) #new fragment scaling based on beta result
		
		nfrags = nfrags *(n_replicates_sample1 + n_replicates_sample2) #we will distribute the fragments to two samples with x replicates
		
		sample1_count=0
		sample2_count=0
		
		for _ in lrange(nfrags):
			
			if random() <= read_frac: #add to sample 1 or sample 2 and to list
				sample1_count += 1
			else:
				sample2_count += 1
		
		#save result for this protein
		sample_counts.append([sample1_count,sample2_count])
		sumcounts.append(nfrags)

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


	#scaling of beta result based on number of fragments
	new_fragment_counts=[]
	new_sample_counts=[]
	new_sumcounts=[]
	
	for _ in lrange(nproteins):
		if(np.random.random() <= options.frag_count_op): #select lognorm or gamma distribution
			nfrags = np.random.lognormal(mean=options.frag_count_om,sigma=options.frag_count_os,size=1) #if yes, it is taken from the gamma
		else:
			nfrags = np.random.gamma(shape=options.frag_count_sh,scale=options.frag_count_sc,size=1) #if no, from the lognormal distribution
		nfrags = max(1, int(nfrags))
		
		read_frac = np.random.beta(options.beta_values[0], options.beta_values[1]) #get read fraction for sample1 (sample2 = 1-x)
		
		#lognorm distri
		sigma=options.frag_ln_si
		scale=options.frag_ln_sc
		loc=scale*-1
		
		
		#to get to 100% at size 10
		top_beta_scale = sp.stats.lognorm.pdf(1,s=sigma,loc=loc, scale=scale)
		scale_factor=1/top_beta_scale #this is the factor to set beta_frac to 100% at 1
		
		beta_scale = sp.stats.lognorm.pdf(nfrags,s=sigma,loc=loc, scale=scale) *scale_factor
		
		#center and scale it based on beta_scale from exponential distri, so reduce to 0.5 for high values
		if(beta_scale == 0):
			read_frac = read_frac 
		else:
			read_frac = ((read_frac-0.5)*beta_scale)+0.5
		
		nfrags = nfrags * (n_replicates_sample1 + n_replicates_sample2) #we will distribute the fragments to two samples with x replicates
		
		sample1_count=0
		sample2_count=0
		
		for _ in lrange(nfrags):
			
			if random() <= read_frac: #add to sample 1 or sample 2 and to list
				sample1_count += 1
			else:
				sample2_count += 1
		
		#save result for this protein
		new_sample_counts.append([sample1_count,sample2_count])
		new_sumcounts.append(nfrags)

	#MA-plot
	new_xma=[]
	new_yma=[]
	new_counts1=[]
	new_counts2=[]
	for entry in new_sample_counts:
		new_counts1.append(entry[0])
		new_counts2.append(entry[1])
		if(entry[1] > 0 and entry[0] > 0):
			x=math.log(entry[0]+entry[1]/2,2)#log2 mean
			new_xma.append(x) 
			y=math.log(entry[0]/entry[1],2)#log2 fold-change
			new_yma.append(y) 


	#no scaling
	old_fragment_counts=[]
	old_sample_counts=[]
	old_sumcounts=[]
	
	for _ in lrange(nproteins):
		if(np.random.random() <= options.frag_count_op): #select lognorm or gamma distribution
			nfrags = np.random.lognormal(mean=options.frag_count_om,sigma=options.frag_count_os,size=1) #if yes, it is taken from the gamma
		else:
			nfrags = np.random.gamma(shape=options.frag_count_sh,scale=options.frag_count_sc,size=1) #if no, from the lognormal distribution
		nfrags = max(1, int(nfrags))
		
		read_frac = np.random.beta(options.beta_values[0], options.beta_values[1]) #get read fraction for sample1 (sample2 = 1-x)
		
		nfrags = nfrags * (n_replicates_sample1 + n_replicates_sample2) #we will distribute the fragments to two samples with x replicates
		
		sample1_count=0
		sample2_count=0
		
		for _ in lrange(nfrags):
			
			if random() <= read_frac: #add to sample 1 or sample 2 and to list
				sample1_count += 1
			else:
				sample2_count += 1
		
		#save result for this protein
		old_sample_counts.append([sample1_count,sample2_count])
		old_sumcounts.append(nfrags)

	#MA-plot
	old_xma=[]
	old_yma=[]
	old_counts1=[]
	old_counts2=[]
	for entry in old_sample_counts:
		old_counts1.append(entry[0])
		old_counts2.append(entry[1])
		if(entry[1] > 0 and entry[0] > 0):
			x=math.log(entry[0]+entry[1]/2,2)#log2 mean
			old_xma.append(x) 
			y=math.log(entry[0]/entry[1],2)#log2 fold-change
			old_yma.append(y) 


	#hist and ma-plot
	fixed_range=[0,200]
	fixed_bins=500
	fig, axs = plt.subplots(3, 4)

	axs[0,0].set_title('No scaling: fragments all samples')
	axs[0,0].hist(old_sumcounts,density=False, bins=fixed_bins, range=fixed_range)
	axs[0,1].set_title('No scaling: fragments sample1')
	axs[0,1].hist(old_counts1,density=False, bins=fixed_bins, range=fixed_range)
	axs[0,2].set_title('No scaling: fragments sample2')
	axs[0,2].hist(old_counts2,density=False, bins=fixed_bins, range=fixed_range)
	axs[0,3].set_title('No scaling: MA-plot')
	axs[0,3].scatter(old_xma, old_yma)
	
	axs[1,0].set_title('Scaling nfragments: fragments all samples')
	axs[1,0].hist(sumcounts,density=False, bins=fixed_bins, range=fixed_range)
	axs[1,1].set_title('Scaling nfragments: fragments sample1')
	#axs[1,1].hist(fragment_counts, bins=range(min(fragment_counts), max(fragment_counts) + binwidth, binwidth))
	axs[1,1].hist(counts1,density=False, bins=fixed_bins, range=fixed_range)
	axs[1,2].set_title('Scaling nfragments: fragments sample2')
	axs[1,2].hist(counts2,density=False, bins=fixed_bins, range=fixed_range)
	axs[1,3].set_title('Scaling nfragments: MA-plot')
	axs[1,3].scatter(xma, yma)
	
	axs[2,0].set_title('Scaling beta: fragments all samples')
	axs[2,0].hist(new_sumcounts,density=False, bins=fixed_bins, range=fixed_range)
	axs[2,1].set_title('Scaling beta: fragments sample1')
	axs[2,1].hist(new_counts1,density=False, bins=fixed_bins, range=fixed_range)
	axs[2,2].set_title('Scaling beta: fragments sample2')
	axs[2,2].hist(new_counts2,density=False, bins=fixed_bins, range=fixed_range)
	axs[2,3].set_title('Scaling beta: MA-plot')
	axs[2,3].scatter(new_xma, new_yma)
	
	plt.subplots_adjust(hspace = 0.4)
	plt.show()


def _test_lognorm():
	sigma=0.9
	loc=-100
	scale=100
	
	x = np.linspace(1,1000, 10000)
	plt.plot(x, sp.stats.lognorm.pdf(x, s=sigma,loc=loc, scale=scale),'r-', lw=5, alpha=0.6, label='lognorm pdf')
	plt.show()
	

def _test_laplace():
	loc= 0.5
	scale = 0.1
	x = np.linspace(0,1, 1000)
	
	plt.plot(x, laplace.pdf(x, loc, scale),'r-', lw=5, alpha=0.6, label='lognorm pdf')
	plt.show()
	
