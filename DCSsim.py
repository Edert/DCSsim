#!/usr/bin/env python
# coding=utf8

"""
%prog <FASTA> <BED> <INT> [options]

Based on the provided <FASTA> sequence differential peaks (DPs) are simulated, restricted to if <--is_white_list> or constrained by the specified <BED> file.
<INT> defines the number of simulated replicates of the two samples.

@author: Thomas Eder
"""

from __future__ import print_function
from __future__ import division

import sys
import os
import time
import datetime
import HTSeq
import numpy as np
import scipy as sp
import multiprocessing as mp
import matplotlib.pyplot as plt

from builtins import range
from future.utils import lrange
from matplotlib.backends.backend_pdf import PdfPages
from optparse import OptionParser
from pybedtools import BedTool, Interval
from random import random, randint
from scipy.stats import beta, laplace
from numpy.random import random_sample, rand, multivariate_normal, choice, gamma, dirichlet
from multiprocessing.managers import BaseManager

MAX_TEST = 10000 #number of maximal test iterations, can influence runtime

#classes
class HelpfulOptionParser(OptionParser):
	"""An OptionParser that prints full help on errors."""
	def error(self, msg):
		self.print_help(sys.stderr)
		self.exit(2, "\n%s: error: %s\n" % (self.get_prog_name(), msg))


class Timer:
	"""Timer to track time to run."""
	def __init__(self):
		self.start = time.time()

	def restart(self):
		self.start = time.time()

	def get_time_hhmmss(self):
		end = time.time()
		m, s = divmod(end - self.start, 60)
		h, m = divmod(m, 60)
		time_str = "%02d:%02d:%02d" % (h, m, s)
		return time_str


class Parameters:
	"""Class to store the main init parameters."""
	def __init__(self, sequence, valid_regions, bins, n_reps_sample1, n_reps_sample2, prot_count_n, prot_count_p, protein_size,\
	frag_len_max, frag_dist_on, frag_dist_mn_mean, frag_dist_mn_cov, frag_dist_prob, prot_dist_mn_mean, prot_dist_mn_cov,\
	frag_count_sh, frag_count_sc, frag_count_op, frag_count_om, frag_count_os, frag_count_scaling, frag_count_lp_sc, frag_count_ln_si, frag_count_ln_sc, \
	beta_values, frag_len_mean, frag_len_dev, skewness, valid_regions_array):
		
		self.sequence = sequence
		self.valid_regions = valid_regions
		self.bins = bins
		self.n_reps_sample1 = n_reps_sample1
		self.n_reps_sample2 = n_reps_sample2
		self.prot_count_n = prot_count_n
		self.prot_count_p = prot_count_p
		self.protein_size = protein_size
		self.frag_len_max = frag_len_max
		self.prot_dist_mn_mean = prot_dist_mn_mean
		self.prot_dist_mn_cov = prot_dist_mn_cov
		self.frag_count_sh = frag_count_sh
		self.frag_count_sc = frag_count_sc
		self.frag_count_op = frag_count_op
		self.frag_count_om = frag_count_om
		self.frag_count_os = frag_count_os
		self.frag_count_scaling = frag_count_scaling
		self.frag_count_lp_sc = frag_count_lp_sc
		self.frag_count_ln_si = frag_count_ln_si
		self.frag_count_ln_sc = frag_count_ln_sc
		self.beta_values = beta_values
		self.frag_len_mean = frag_len_mean
		self.frag_len_dev = frag_len_dev
		self.skewness = skewness
		self.valid_regions_array = valid_regions_array
		self.frag_dist_on = frag_dist_on
		self.frag_dist_mn_mean = frag_dist_mn_mean
		self.frag_dist_mn_cov = frag_dist_mn_cov
		self.frag_dist_prob = frag_dist_prob

	def get_all(self):
		""""Returns all parameters stored"""
		return self.sequence, self.valid_regions, self.bins ,self.n_reps_sample1, self.n_reps_sample2, self.prot_count_n, self.prot_count_p, \
		self.protein_size, self.frag_len_max, self.frag_dist_on, self.frag_dist_mn_mean, self.frag_dist_mn_cov, self.frag_dist_prob, self.prot_dist_mn_mean, \
		self.prot_dist_mn_cov, self.frag_count_sh, self.frag_count_sc, self.frag_count_op, self.frag_count_om, self.frag_count_os, \
		self.frag_count_scaling, self.frag_count_lp_sc, self.frag_count_ln_si, self.frag_count_ln_sc, self.beta_values, self.frag_len_mean, \
		self.frag_len_dev, self.skewness, self.valid_regions_array

class Report_data:
	"""Class to store report data"""
	def __init__(self, n_reps_sample1, n_reps_sample2):
		
		self.r_prot_pos = []
		self.r_dom_nprot = []
		self.r_prot_offset_positions = [] 
		self.r_nfrag = []
		self.r_n_frag1 = []
		self.r_n_frag2 = []
		self.r_fragshift = []
		self.r_fraglen = []
		self.r_n_frag1_reps = [[] for i in lrange(n_reps_sample1)] 
		self.r_n_frag2_reps = [[] for i in lrange(n_reps_sample2)]
		self.db_peaks_sample1 = 0
		self.db_peaks_sample2 = 0
		self.n_domains = 0

	def _extend(self, n_domains, r_prot_pos, r_dom_nprot, r_prot_offset_positions, r_nfrag, r_n_frag1, r_n_frag2, r_fragshift, r_fraglen, r_n_frag1_reps, r_n_frag2_reps, db_peaks_sample1, db_peaks_sample2):
		"""Extend report data lists"""
		self.n_domains += n_domains
		self.r_prot_pos.extend(r_prot_pos)
		self.r_dom_nprot.extend(r_dom_nprot)
		self.r_prot_offset_positions.extend(r_prot_offset_positions)
		self.r_nfrag.extend(r_nfrag)
		self.r_n_frag1.extend(r_n_frag1)
		self.r_n_frag2.extend(r_n_frag2)
		self.r_fragshift.extend(r_fragshift)
		self.r_fraglen.extend(r_fraglen)
		
		rep_index = 0
		for rep in self.r_n_frag1_reps:
			rep.extend(r_n_frag1_reps[rep_index])
			rep_index += 1
		rep_index = 0
		for rep in self.r_n_frag2_reps:
			rep.extend(r_n_frag2_reps[rep_index])
			rep_index += 1
		
		self.db_peaks_sample1 += db_peaks_sample1
		self.db_peaks_sample2 += db_peaks_sample2


class Counter(object):
	"""Counter for the number of domains to create, locks for multiprocessing"""
	def __init__(self, initval=0, maxval=50):
		self.val = mp.Value('i', initval)
		self.lock = mp.Lock()
		self.maxval = maxval

	def increment(self):
		with self.lock:
			if(self.val.value >= self.maxval):
				return -1
			else:
				self.val.value += 1
				return self.val.value

	def get_value(self):
		with self.lock:
			return self.val.value

	def is_max(self):
		with self.lock:
			if(self.val.value >= self.maxval):
				return True
			else:
				return False


class Chromosome(object):
	"""The chromosome were we want to simulate differential binding peaks on, it has the following attributes:
	
	Attributes:
		sequence: str of the genomic sequence
		name: str of the chromosome name
		length: integer of seqeunce length
		domain_list: list of domains
		valid_regions: regions which are not blacklisted
		bins: probabilites of the valid_regions based on their length
		n_reps_sample1: number of replicates in sample1
		n_reps_sample2: number of replicates in sample2
		is_white_list: if the supplied bed file is treated as white-list
		threads: number of processes used to create domains, proteins and fragments
	"""

	def __init__(self, sequence, chrom_name, bedfile_path, n_reps_sample1, n_reps_sample2, is_white_list, threads):
		"""Returns a chromosome object"""
		self.sequence = sequence
		self.name = chrom_name
		self.length = len(sequence)
		self.domain_list = []
		self.valid_regions, self.bins = self._get_valid_regions(bedfile_path, self.length, is_white_list) #define valid regions based on provided bed file
		self.n_reps_sample1 = n_reps_sample1
		self.n_reps_sample2 = n_reps_sample2
		self.is_white_list=is_white_list
		self.threads = threads

	def init_domains(self, n_domains_start, n_domains, prot_count_n, prot_count_p, protein_size, frag_len_max, frag_dist_on, frag_dist_mn_mean, \
	frag_dist_mn_cov, frag_dist_prob, prot_dist_mn_mean, prot_dist_mn_cov, frag_count_sh, frag_count_sc, frag_count_op, frag_count_om, frag_count_os, \
	frag_count_scaling, frag_count_lp_sc, frag_count_ln_si, frag_count_ln_sc, beta_values, frag_len_mean, frag_len_dev, skewness):
		
		"""Creates n_domains on the sequence"""
		j = 0
		work = []
		valid_regions_array = []

		if(len(self.domain_list) == 0):
			
			#build valid regions array for quick checks...
			valid_regions_array = np.zeros(self.length+1,bool) #array of false values
			for region in self.valid_regions:
				for i in lrange(region.start,region.stop+1):
					valid_regions_array[i] = True #set region valid

			#create obj to store all parameters
			parameterobj = Parameters(self.sequence, self.valid_regions, self.bins, self.n_reps_sample1, self.n_reps_sample2, prot_count_n, prot_count_p,\
			protein_size, frag_len_max, frag_dist_on, frag_dist_mn_mean, frag_dist_mn_cov, frag_dist_prob, prot_dist_mn_mean, prot_dist_mn_cov, \
			frag_count_sh, frag_count_sc, frag_count_op, frag_count_om, frag_count_os, frag_count_scaling, frag_count_lp_sc, frag_count_ln_si, frag_count_ln_sc, \
			beta_values, frag_len_mean, frag_len_dev, skewness, valid_regions_array)
			
			if (self.threads > 1): #use  multiple processes
				
				print("Initializing domains %s to %s in %s processes" %(n_domains_start + 1, n_domains_start + n_domains , self.threads))
				
				workers = self.threads #number of workers
				counter = Counter(n_domains_start, n_domains_start + n_domains) #process safe counter from n_domains_start to the number of domains, it is used to name the domains too
				process_list = []
				pipe_list = []
				
				for _ in range(workers): #each worker creates new domains
					recv_end, send_end = mp.Pipe(False) #create unidirectional pipe
					p = mp.Process(target=_init_add_domains_mp, args=(counter, parameterobj, send_end))
					pipe_list.append(recv_end)
					
					p.daemon = True
					p.start() #start process
					process_list.append(p) #append to process list
				
				while True: #fill progress bar
					complete_count = counter.get_value() - n_domains_start 
					_progress_bar(complete_count, n_domains )
					if(complete_count == n_domains): #wait till all workers are done...
						break
					time.sleep(0.1)
				
				print("\n\nCollecting results from %s threads" % self.threads)
				collect_count = 0
				
				for pipe in pipe_list: #go through pipe list and get results from the end points
					domain_list = pipe.recv()
					
					self.domain_list.extend(domain_list)
					collect_count += len(domain_list)
					
					_progress_bar(collect_count, n_domains)
					pipe.close()
				
				for p in process_list:
					p.join()
				
				print("\n")
			
			else:	#use just a single process
				self.domain_list = _init_add_domains(parameterobj, n_domains_start, n_domains)
			
		else:
			print("Domains already initialized")

	def _get_valid_regions(self, bedfile_path, chromosome_length, is_white_list):
		"""Return valid regions and accumulated probabilities to peak one region"""
		if(is_white_list):
			regions = BedTool(bedfile_path)
			return self.__get_valid_region(regions)
		else:
			black_regions = BedTool(bedfile_path)
			valid_regions = BedTool('%s 0 %s' %(self.name, chromosome_length), from_string=True) #create whole chromosome region
			
			return self.__subtract_from_valid_region(valid_regions, black_regions) #subtracted blacklisted regions from valid regions in this case the whole chromosome

	def __subtract_from_valid_region(self, valid_regions, black_regions):
		"""Subtract <black_regions> from <valid_regions> and return bins and <valid_regions>"""
		remaining_regions = valid_regions.subtract(black_regions) #get only valid regions without overlapping blacklisted regions
		
		if(len(remaining_regions) == 0):
			parser.error("No valid sequence region left, please check provided bed-file (blacklist/whitelist) and add valid/remove blacklisted regions or reduce \
			the number of domains '-d'.") 
		
		total_length = 0
		relative_lengths = []
		
		for region in remaining_regions:
			total_length += len(region) #sum up regions
			relative_lengths.append(len(region))
		prob = list(map(lambda x: x/total_length, relative_lengths)) #get probability based on region length
		bins = np.add.accumulate(prob)
		
		return remaining_regions, bins

	def __get_valid_region(self, valid_regions):
		"""Return bins and <valid_regions>"""
		total_length = 0
		relative_lengths = []
		
		for region in valid_regions:
			total_length += len(region)
			relative_lengths.append(len(region))
		prob = list(map(lambda x: x/total_length, relative_lengths))
		bins = np.add.accumulate(prob)
		
		return valid_regions, bins
		
	def wipe_data(self):
		"""Removes all Domains, Proteins and Fragments"""
		for domain in self.domain_list :
			domain.wipe_data()
			del(domain)
		
		self.domain_list = []


class SpikeInChromosome(Chromosome):
	"""The chromosome of the spike-in species, it has the following attributes:
	
	Attributes:
		sequence: str of the genomic sequence
		name: str of the chromosome name
		length: integer of sequence length
	"""

	def __init__(self, sequence, chrom_name):
		"""Return a chromosome object"""
		self.sequence = sequence
		self.name = chrom_name
		self.length = len(sequence)


class Domain(object):
	"""A domain in the genome sequence, it has the following attributes:
	
	Attributes:
		name: name or number of the domain
		protein_list: list of proteins
		n_proteins: integer of the number of proteins
		protein_positions: list of protein positions
		self.protein_offset_positions: list of the protein offset positions from domain start
		n_reps_sample1: number of replicates in sample1
		n_reps_sample2: number of replicates in sample2
	"""

	def __init__(self, n_reps_sample1, n_reps_sample2, name):
		"""Return a domain object"""
		self.name = name
		self.protein_list = []
		self.n_proteins = 0
		self.protein_positions = []
		self.protein_offset_positions = []
		self.n_reps_sample1 = n_reps_sample1
		self.n_reps_sample2 = n_reps_sample2

	def init_proteins(self, genome_sequence, valid_regions, bins, prot_count_n, prot_count_p, protein_size, frag_len_max, frag_dist_on, frag_dist_mn_mean,\
	frag_dist_mn_cov, frag_dist_prob, prot_dist_mn_mean, prot_dist_mn_cov, frag_count_sh, frag_count_sc, frag_count_op, frag_count_om, frag_count_os,\
	frag_count_scaling, frag_count_lp_sc, frag_count_ln_si, frag_count_ln_sc, beta_values, frag_len_mean, frag_len_dev, skewness, valid_regions_array):
		"""Creates n_proteins of proteins in the domain range"""
		if(len(self.protein_list) == 0):
			
			#get number of proteins in this domain
			tmp_n_proteins = self._sample_neg_bin(prot_count_n, prot_count_p) #number will be change if we hit a blacklisted region
			
			#get positions in valid regions...
			self.protein_positions, self.n_proteins, self.protein_offset_positions = self._get_protein_positions(tmp_n_proteins, protein_size, prot_dist_mn_mean, \
			prot_dist_mn_cov, genome_sequence, valid_regions, bins,frag_len_mean, frag_len_max, valid_regions_array)
			
			prot_count = 1
			for pos in self.protein_positions:
				min_protein_pos = 0 #save min and max of fragments
				max_protein_pos = 0
				counts_frag_sample1, counts_frag_sample2 = 0, 0
				tmp_fragments_sample1, tmp_fragments_sample2 = [], []
				
				#init one protein
				proteinobj = Protein(pos, self.n_reps_sample1, self.n_reps_sample2, self.name, prot_count)
				
				#init all fragments in that protein
				proteinobj.init_fragments(genome_sequence, frag_count_sh, frag_count_sc, frag_count_op, frag_count_om, frag_count_os, \
				frag_count_scaling, frag_count_lp_sc, frag_count_ln_si, frag_count_ln_sc, beta_values, \
				frag_len_mean, frag_len_dev, frag_len_max, frag_dist_on, frag_dist_mn_mean, frag_dist_mn_cov, frag_dist_prob, protein_size, skewness)
				
				#add to list
				self.protein_list.append(proteinobj)
				
				prot_count += 1
		else:
			print("Proteins already initialized")

	def _get_protein_positions(self, tmp_n_proteins, protein_size, prot_dist_mn_mean, prot_dist_mn_cov, genome_sequence, valid_regions, bins, frag_len_mean, frag_len_max, \
	valid_regions_array):
		"""Get protein positions in valid regions"""
		count = 0
		while(count < MAX_TEST):
			count += 1

			#get protein distances in this domain
			first_protein_offset_positions = self._sample_multivariate_normal(tmp_n_proteins, protein_size, prot_dist_mn_mean, prot_dist_mn_cov)

			#get most left protein position
			protein_left_position = self._get_left_protein_position(genome_sequence, tmp_n_proteins, int(max(first_protein_offset_positions)), valid_regions, bins, frag_len_mean, frag_len_max, protein_size)

			#check positions 1 to n for valid regions, if not just omit them and add offsets to left_position
			protein_positions, n_proteins, protein_offset_positions = self._get_check_offset_positions(protein_left_position, first_protein_offset_positions, valid_regions_array, protein_size, frag_len_mean)
			
			if(n_proteins > 0): #we found at least one valid region
				return protein_positions, n_proteins, protein_offset_positions
			
		parser.error("Could not find at least one valid protein position please check the bed-file (blacklist/whitelist) and the parameters regarding protein size and number.") 

	def _sample_multivariate_normal(self, k, protein_size, prot_dist_mn_mean, prot_dist_mn_cov):
		"""Return k samples from a multivariate normal distribution """
		res = list(map(lambda x: max(1,int(choice(list(x)))), multivariate_normal(prot_dist_mn_mean, prot_dist_mn_cov, k))) 
		res.insert(0,int(protein_size/2)) #make sure first element is more than half of protein size
		res.pop() #remove last element to be sure that we did not change number of proteins
		#res = [int(x) for x in res] #all elements to int
		return res

	def _get_left_protein_position(self, genome_seq, number_proteins, max_offset, valid_regions, bins, frag_len_mean, frag_len_max, protein_size):
		"""Return position of leftmost protein"""
		offset = (number_proteins - 1) * max_offset #get theoretical maximal offset
		
		count = 0
		while (count < MAX_TEST): # test only MAX_TEST regions after that exit..
			count += 1
			valid = True

			pos = self._get_random_position(valid_regions, bins) #get one position in an randomly selected valid_region

			#check the potential sequence for illegal characters... also do not accept N´s
			tmp_seq = genome_seq[pos - frag_len_max : pos + frag_len_max + offset]

			if tmp_seq.count('A') + tmp_seq.count('C') + tmp_seq.count('G') + tmp_seq.count('T') == len(tmp_seq):
				return pos
				
		parser.error("Could not find a sequence with valid nucleotides, please check the provided fasta file and or bed-file (blacklist/whitelist).") 

	def _sample_neg_bin(self, n, p):
		"""Sample negative binomial distribution"""
		return int(max(np.random.negative_binomial(n, p),1))

	def _get_random_position(self, valid_regions, bins):
		"""Return random start and end positions of region (+/- fragment length mean), based on accumulated probability (based on region size) in bins"""
		
		chosen_region = np.digitize(random_sample(1), bins)[0].item() #get one random bin and convert numpyint to int
		return randint((valid_regions[chosen_region].start ), (valid_regions[chosen_region].stop ))

	def _get_check_offset_positions(self, protein_left_position, protein_offset_positions, valid_regions_array, protein_size, frag_len_mean):
		"""Checks if the offset positions realtive to the left position are in valid regions, if not removes them"""
		new_protein_offset_positions = []
		approx_dif = abs(frag_len_mean - protein_size)//2 #estimate shifts of fragments
		n_proteins = len(protein_offset_positions)

		protein_absolute_positions = protein_left_position + np.cumsum(protein_offset_positions) #add offset positions
		
		rm_list = []
		rm_index = 0
		for offset_pos in protein_absolute_positions: 
			isgood = True
			start = int(offset_pos - approx_dif)
			stop = int(offset_pos + approx_dif)
			
			for i in lrange(start,stop+1):
				if(valid_regions_array[i] == False):#found overlap with blacklisted region
					isgood = False
					n_proteins -= 1
					rm_list.append(rm_index)
					break
			
			if(isgood):#all positions are valid
				new_protein_offset_positions.append(offset_pos)
			rm_index += 1

		protein_offset_positions = np.delete(protein_offset_positions, rm_list) #remove also from offset list
		
		return new_protein_offset_positions, n_proteins, protein_offset_positions

	def wipe_data(self):
		"""Removes all Proteins and Fragments"""
		for protein in self.protein_list :
			protein.wipe_data()
			del(protein)

class Protein(object):
	"""A protein in a domain, it has the following attributes:
	
	Attributes:
		start = start position (most left (5') fragment)
		stop = stop position (most right (3') fragment)
		number_frags = number of fragments
		sample1_mean_fragment_count = mean fragment count in sample1
		sample2_mean_fragment_count = mean fragment count in sample2
		fragment_count = fragment count in this object
		name = name or id of the protein
		domain_name = name or if of its domain
		position = position of the protein (summit of the fragments)
		n_reps_sample1: number of replicates in sample1
		n_reps_sample2: number of replicates in sample2	
		sample1_replicate_counts = counts of fragments per replicate in sample1
		sample2_replicate_counts = counts of fragments per replicate in sample2
		sample1_replicate_list = list of the fragment-objects in the respecitve replicates of sample1
		sample2_replicate_list = list of the fragment-objects in the respecitve replicates of sample2
	"""
	
	def __init__(self, position, n_reps_sample1, n_reps_sample2, domain_name, name):
		"""Return a protein object"""
		self.start = 0
		self.stop = 0
		self.number_frags = 0
		self.sample1_mean_fragment_count = 0
		self.sample2_mean_fragment_count = 0
		self.fragment_count = 0
		self.name = name
		self.domain_name = domain_name
		
		self.position = position
		self.n_replicates_sample1 = n_reps_sample1
		self.n_replicates_sample2 = n_reps_sample2
		
		self.sample1_replicate_counts = [0] * self.n_replicates_sample1
		self.sample2_replicate_counts = [0] * self.n_replicates_sample2
		self.sample1_replicate_list = [[] for i in lrange(self.n_replicates_sample1)]
		self.sample2_replicate_list = [[] for i in lrange(self.n_replicates_sample2)]

	def init_fragments(self, genome_sequence, frag_count_sh, frag_count_sc, frag_count_op, frag_count_om, frag_count_os, \
	frag_count_scaling, frag_count_lp_sc, frag_count_ln_si, frag_count_ln_sc, beta_values, \
	frag_len_mean, frag_len_dev, frag_len_max, frag_dist_on, frag_dist_mn_mean, frag_dist_mn_cov, frag_dist_prob, protein_size, skewness):
		"""Creates n_fragments in this protein"""
		if(self.fragment_count == 0 ):
			min_protein_pos = 0
			max_protein_pos = 0
			sample1_list = []
			sample2_list = []
			sample1_count = 0
			sample2_count = 0
			
			#get number of fragments per protein
			if(np.random.random() <= frag_count_op):
				self.number_frags = np.random.lognormal(mean=frag_count_om,sigma=frag_count_os,size=1) #if yes, it is taken from the gamma
			else:
				self.number_frags = np.random.gamma(shape=frag_count_sh,scale=frag_count_sc,size=1) #if no, from the lognormal distribution
			
			self.number_frags = max(1, int(self.number_frags)) #set to 1 or more
			
			#get read fraction for sample1 (sample2 = 1-x)
			read_frac = np.random.beta(options.beta_values[0], options.beta_values[1]) 
			
			
			#now choose which way we want to select the number of fragments
			if( frag_count_scaling == "none"):
				#no scaling, no changes in nfrags and beta
				self.number_frags = self.number_frags
				read_frac = read_frac
			elif (frag_count_scaling == "beta") :#nfrags scaling via Laplace, beta not changed
				self.number_frags = self._scale_laplace(self.number_frags, read_frac, frag_count_lp_sc)
			elif (frag_count_scaling == "frag") :#nfrags scaling via lognorm distribution, number_frags not changed
				read_frac = self._scale_lognorm(self.number_frags, read_frac, frag_count_ln_si, frag_count_ln_sc)
			else:
				print("Unknown scaling method, %s, please choose 'none','frag' or 'beta', exiting now" % (frag_count_scaling))
				exit(1)
			
			#we will distribute the fragments to two samples with x replicates
			self.number_frags = self.number_frags * (self.n_replicates_sample1 + self.n_replicates_sample2) 
			
			for _ in lrange(self.number_frags):
				
				#init one fragment with position, sequence and type
				fragmentobj=Fragment(self.position, genome_sequence, "p", frag_len_mean, frag_len_dev, frag_len_max, frag_dist_on, frag_dist_mn_mean, \
				frag_dist_mn_cov, frag_dist_prob, protein_size, self.domain_name, self.name )
				
				#to get peak start and stop we need to save min and max of fragment positions
				if(min_protein_pos == 0 or (fragmentobj.position - fragmentobj.frag_length//2) < min_protein_pos):
					min_protein_pos = fragmentobj.position - fragmentobj.frag_length//2
				if(max_protein_pos == 0 or (fragmentobj.position + fragmentobj.frag_length//2) > max_protein_pos):
					max_protein_pos = fragmentobj.position + fragmentobj.frag_length//2
				
				if random() <= read_frac: #add to sample 1 or sample 2 and to list
					sample1_count += 1
					sample1_list.append(fragmentobj)
				else:
					sample2_count += 1
					sample2_list.append(fragmentobj)
				
				self.fragment_count += 1
			
			#save positions
			self.start = min_protein_pos
			self.stop = max_protein_pos
			
			#save get distribution into replicates
			self.sample1_replicate_list, self.sample1_replicate_counts = self._fragments_to_replicates(sample1_list, skewness, self.n_replicates_sample1)
			self.sample2_replicate_list, self.sample2_replicate_counts = self._fragments_to_replicates(sample2_list, skewness, self.n_replicates_sample2)
			
			#save mean counts per sample
			self.sample1_mean_fragment_count = sample1_count / self.n_replicates_sample1
			self.sample2_mean_fragment_count = sample2_count / self.n_replicates_sample2
		
		else:
			print("Fragments already initialized")

	def _fragments_to_replicates(self, fragment_list, skewness, n_replicates):
		"""puts fragments in a randomly selected replicate"""
		#get fragments into replicates
		bins = np.cumsum([0] + list(dirichlet(np.array([skewness] * n_replicates), 1)[0])) #sample dirichlet dist, get array, make list, add 0 to beginning and sum up all values to 1, so we end up with an array of bins with probabilities
		replicate_counts = [0] * n_replicates
		replicate_list = [[] for i in lrange(n_replicates)]
		
		for fragmentobj in fragment_list: #go through fragments
			index = np.digitize([random()], bins)[0].item() - 1 #get index-1 of random number (between 0 or 1) fitting to the right bin
			replicate_list[index].append(fragmentobj) #add element to replicate via choosen index
			replicate_counts[index] += 1
		
		return replicate_list, replicate_counts #return number of fragments per replicate

	def _scale_laplace(self, nfrags, read_frac, scale):
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
		
	def _scale_lognorm(self, nfrags, read_frac, sigma, scale):
		"""Scale number_frags result based on Lognorm distribution"""
		loc=scale*-1
		
		top_beta_scale = sp.stats.lognorm.pdf(1,s=sigma,loc=loc, scale=scale)#to get to 100% at size 10
		scale_factor=1/top_beta_scale #this is the factor to set beta_frac to 100% at 1
		beta_scale = sp.stats.lognorm.pdf(nfrags,s=sigma,loc=loc, scale=scale) *scale_factor
		
		if(beta_scale == 0):
			read_frac = read_frac 
		else:
			read_frac = ((read_frac-0.5)*beta_scale)+0.5 #center and scale it based on beta_scale from exponential distri, so reduce to 0.5 for high values
		return read_frac
		
	def _get_factor(self, x, loc=0.5, scale=0.2):
		"""laplace distribution for distribution into replicates"""
		return laplace.pdf(x, loc, scale)

	def __lt__(self, other): #for sorting
		return self.start < other.start

	def wipe_data(self):
		"""Removes all Fragments"""
		self.sample1_replicate_counts = []
		self.sample2_replicate_counts = []
		self.sample1_replicate_list = []
		self.sample2_replicate_list = []

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
	
	def __init__(self, prot_position, genome_seq, frag_type, frag_len_mean, frag_len_dev, frag_len_max, frag_dist_on, frag_dist_mn_mean, frag_dist_mn_cov,\
	frag_dist_prob, protein_size, domain_name, protein_name):
		"""Return a fragment object"""
		self.domain_name = domain_name
		self.protein_name = protein_name
		
		self.frag_length = self._sample_frag_length(frag_len_mean, frag_len_dev, frag_len_max) #get fragment length
		
		if not (frag_dist_on): #simply random choosen shift
			self.rand_shift = randint(-abs(frag_len_mean - protein_size)//2, abs(frag_len_mean - protein_size)//2) #add random shift inside the protein size
		else : #shift from multivariate distribution defined by user
			shift = self._sample_multivariate_normal(protein_size, frag_dist_mn_mean, frag_dist_mn_cov, frag_dist_prob)
			self.rand_shift = int(-abs(frag_len_mean - protein_size)//2) + shift #from most left position on add shift

		self.position = prot_position + self.rand_shift #final position of the fragment
		
		#get sequence
		self.sequence = genome_seq[self.position - self.frag_length//2 : self.position + self.frag_length//2]
		self.frag_type = frag_type

	def _sample_frag_length(self, mean, dev, frag_len_max):
		"""Return length of fragment for peaks"""
		return min(abs(int(np.random.normal(loc = mean, scale = dev))), frag_len_max)

	def _sample_multivariate_normal(self, protein_size, frag_dist_mn_mean, frag_dist_mn_cov, frag_dist_prob):
		"""Return a sample from a multivariate normal distribution """
		res = min(max(int(choice(multivariate_normal(frag_dist_mn_mean, frag_dist_mn_cov), p = frag_dist_prob)), 1), protein_size)
		#choose a random element from the n distributions, min = 1 and max = protein_size
		return int(res)


class NoiseFragment(Fragment):
	"""A noise fragment it has the following attributes:
	
	Attributes:
		position: integer of the final position of the fragment
		frag_length: integer of th length of the fragment
		sequence: string of the sequence of the fragment
		frag_type: string of the fragment type
	"""
	
	def __init__(self, position, genome_seq, frag_type, frag_len_mean, frag_len_dev, frag_len_max):
		"""Return a fragment object"""
		self.position = position
		self.frag_length = self._sample_frag_length(frag_len_mean, frag_len_dev, frag_len_max) #get fragment length
		self.sequence = genome_seq[self.position - self.frag_length//2 : self.position + self.frag_length//2]
		self.frag_type = frag_type


class Sample(object):
	"""A sample which holds all replikates with fragments per replikate, it has the following attributes:
	
	Attributes:
		replicate_list: list of replicate-objects
		n_replicates: number of replicates
		sample_id: string of the id/name
		chromosome: str of the chromosome name
		weighted_bins: list of weights of the bins of the chromosome
		noise_regions: list of the regions the weighted bins are there for
		myinput: input-object for this sample
	"""

	def __init__(self, chromosome, n_replicates, sample_id, weighted_bins, noise_regions):
		"""Return a sample object"""
		self.replicate_list = []
		self.n_replicates = n_replicates
		self.sample_id = sample_id
		self.chromosome = chromosome
		self.weighted_bins = []
		self.noise_regions = []
		self.weighted_bins = weighted_bins
		self.noise_regions = noise_regions
		self.myinput=Input()
		
		rep_id = 1 #start with 1 and increase per replicate
		for _ in lrange(n_replicates): #create replicates
			replicateobj = Replicate(rep_id)
			self.replicate_list.append(replicateobj)
			rep_id += 1

	def get_fragment_number(self):
		"""Returns the number of fragments in this sample"""
		n_fragments = 0
		for replicate in self.replicate_list:
			n_fragments += replicate.get_fragment_number()
		return n_fragments

	def add_db_fragments(self, fragment_list, index):
		"""Adds a fragment-list from a db peak to the chosen index of the replicates"""
		self.replicate_list[index].add_fragments(fragment_list) #add fragment

	def add_noise(self, back_avg, read_length, frag_len_mean, frag_len_dev, frag_len_max):
		"""Adds noise to all replicates of this sample"""
		for rep in self.replicate_list:
			print("Adding noise fragments to sample %s replicate %s" % (self.sample_id, rep.replicate_id))
			rep.add_noise(self.chromosome, self.weighted_bins, self.noise_regions, back_avg, read_length, frag_len_mean, frag_len_dev, frag_len_max)

	def write_fasta(self, prefix, read_length, read_count, append = False):
		"""Write reads of fragments in the replicates to a fasta file"""
		for rep in self.replicate_list:#write fasta per replicate
			read_count = rep.write_fasta(self.sample_id, prefix, read_length, read_count, append)
		
		self.myinput.write_fasta(self.sample_id, prefix, read_length)#write input fasta
		
		return read_count

	def add_input(self, frag_len_mean, frag_len_dev, frag_len_max, back_avg, read_length):
		"""Adds input to this sample"""
		print("Adding Input fragments to sample %s " % (self.sample_id))
		self.myinput.create_fragments(self.chromosome, self.weighted_bins, self.noise_regions, frag_len_mean, frag_len_dev, frag_len_max, back_avg, read_length)

	def add_spike_in(self, spike_in_chrom, si_weighted_bins, si_regions, spike_in_cov, read_length, frag_len_mean, frag_len_dev, frag_len_max):
		"""Adds spike-in fragments to this sample"""	
		for rep in self.replicate_list:
			print("Adding spike-in fragments to sample %s replicate %s" % (self.sample_id, rep.replicate_id))
			rep.add_spike_in(spike_in_chrom, si_weighted_bins, si_regions, spike_in_cov, read_length, frag_len_mean, frag_len_dev, frag_len_max)
			
		print("Adding spike-in fragments to Input of sample %s " % (self.sample_id))
		self.myinput.add_spike_in(spike_in_chrom, si_weighted_bins, si_regions, spike_in_cov, read_length, frag_len_mean, frag_len_dev, frag_len_max)


class Replicate(object):
	"""A replikate with fragments per replikate, it has the following attributes:
	
	Attributes:
		fragment_list: list of fragment-objects
		n_peak_fragments: integer of number fragments from peaks
		n_noise_fragments: integer of number fragments from noise
		n_spike_in_fragments: integer of number of spike-in fragments
		replicate_id: str of the id/name
	"""

	def __init__(self, replicate_id):
		"""Return a sample object"""
		self.fragment_list= []
		self.n_peak_fragments = 0
		self.n_noise_fragments = 0
		self.n_spike_in_fragments = 0
		self.replicate_id = replicate_id

	def get_fragment_number(self):
		"""Return number of fragments in this replicate"""
		return len(self.fragment_list)

	def add_fragments(self, fragment_list):
		"""Add a fragment to this replicate"""
		self.fragment_list.extend(fragment_list)
		self.n_peak_fragments+=len(fragment_list)

	def add_noise(self, chrom, bins, noise_regions, back_avg, read_length, frag_len_mean, frag_len_dev, frag_len_max):
		"""Add noise to this replicate"""
		number_noise_reads = int(back_avg * chrom.length // read_length) #only the provided genome coverage is taken into account

		s = (int(max(bins))-1) * np.random.random_sample(number_noise_reads) #get probability per noise-fragment
		chosen_regions = np.digitize(s, bins) #get indices of bins where the reads belong
		
		for cri in chosen_regions: #get indexes 
			pos = randint(noise_regions[cri].start, noise_regions[cri].stop) #get random position inside the region
			noisefragmentobj=NoiseFragment(pos, chrom.sequence, "n", frag_len_mean, frag_len_dev, frag_len_max) #create new noise fragment with position, sequence and type
			self.fragment_list.append(noisefragmentobj) #add to fragment list
			self.n_noise_fragments+=1

	def add_spike_in(self, spike_in_chrom, si_weighted_bins, si_regions, spike_in_cov, read_length, frag_len_mean, frag_len_dev, frag_len_max):
		"""Add spike-in fragments to this replicate"""
		number_spike_in_reads = int(spike_in_cov * spike_in_chrom.length // read_length) #only the provided genome coverage is taken into account

		s = (int(max(si_weighted_bins))-1) * np.random.random_sample(number_spike_in_reads) #get probability per noise-fragment
		chosen_regions = np.digitize(s, si_weighted_bins) #get indices of bins where the reads belong
		
		for cri in chosen_regions: #get indexes 
			pos = randint(si_regions[cri].start, si_regions[cri].stop) #get random position inside the region
			sifragmentobj=NoiseFragment(pos,spike_in_chrom.sequence,"s", frag_len_mean, frag_len_dev, frag_len_max) #create spike-in fragment, with position, sequence and type
			self.fragment_list.append(sifragmentobj) #add to fragment list
			self.n_spike_in_fragments+=1

	def write_fasta(self, sample_id, prefix, read_length, read_count, append = False):
		"""Write reads of fragment to a fasta file"""
		if(len(self.fragment_list) != 0):

			name = prefix + "_sample%s-rep%s.fasta" %(sample_id, self.replicate_id)
			print("Writing reads to fasta file: %s" % (name))
			
			if(append):
				fastafile = open(name, "a") #open file stream to write
			else:
				fastafile = open(name, "w")
			
			for frag in self.fragment_list:#randomly choose if forward or reverse
				if random() < 0.5: #forward
					seq = frag.sequence[:read_length]
					if not len(seq) == read_length:
						seq = seq.rjust(read_length,'N') #if it is shorter than the read length fill it up with N´s
					
				else: #reverse
					seq = frag.sequence[len(frag.sequence) - read_length:]
					if not len(seq) == read_length:
						seq = seq.rjust(read_length,'N')
					seq = self._reverse_complement(seq)
					
				#write
				if(frag.frag_type == "p"): #add more domain and protein infos
					fastafile.write(">read%s_%s_d%sp%s\n%s\n" %(read_count, frag.frag_type, frag.domain_name, frag.protein_name, seq))
				else:
					fastafile.write(">read%s_%s\n%s\n" %(read_count, frag.frag_type, seq))
				read_count += 1
				
			return read_count
		else:
			print("No fragments in replicate to write")

	def _reverse_complement(self, sequence, rev=True):
		"""Return the reverse complement of a DNA sequence sequence"""
		_complement = dict(A="T", T="A", C="G", G="C", N="N")
		reverse = reversed(sequence) if rev else sequence
		try:
			reversecomplement = (_complement[x] for x in reverse)
		except KeyError:
			print("Seqeunce %s could not be transfered to referce complement" %sequence)
			return sequence
			
		return "".join(reversecomplement)  # join the elements into a string


class Input(Replicate):
	"""The input with fragments, it has the following attributes:
	
	Attributes:
		fragment_list: list of fragment-objects
		n_fragments: integer of number fragments
		n_spike_in_fragments: integer of spike-in fragments
	"""

	def __init__(self):
		"""Return a sample object"""
		self.fragment_list= []
		self.n_fragments = 0
		self.n_spike_in_fragments = 0

	def create_fragments(self, chrom, bins, noise_regions, frag_len_mean, frag_len_dev, frag_len_max, back_avg, read_length):
		"""Add fragments to this input"""
		number_noise_reads = int(back_avg * chrom.length // read_length) #only the provided genome coverage is taken into account

		s = (int(max(bins))-1) * np.random.random_sample(number_noise_reads) #get probability per noise-fragment
		chosen_regions = np.digitize(s, bins) #get indices of bins where the reads belong
		
		for cri in chosen_regions: #get indexes 
			pos = randint(noise_regions[cri].start, noise_regions[cri].stop) #get random position inside the region
			noisefragmentobj=NoiseFragment(pos,chrom.sequence,"i", frag_len_mean, frag_len_dev, frag_len_max) #create new input fragment with position, sequence and type
			self.fragment_list.append(noisefragmentobj) #add to fragment list
			self.n_fragments+=1

	def write_fasta(self, sample_id, prefix, read_length):
		"""Write reads of fragment to a fasta file"""
		if(len(self.fragment_list) != 0): #if no input is required there are no fragments so do not write a fasta
			name = prefix + "_sample%s-INPUT.fasta" %(sample_id)
			print("Writing reads to fasta file: %s" % (name))
			
			fastafile = open(name, "w") #open file stream to write
			entry_count = 1
			
			for frag in self.fragment_list:
				if random() < 0.5: #randomly choose if forward or reverse
					seq_forward = frag.sequence[:read_length]
					if not len(seq_forward) == read_length:
						seq_forward = seq_forward.rjust(read_length,'N') #if it is shorter than the read length fill it up with N´s
					fastafile.write(">read%s_%s\n%s\n" %(entry_count, frag.frag_type, seq_forward))
				else:
					seq_reverse = frag.sequence[len(frag.sequence) - read_length:]
					if not len(seq_reverse) == read_length:
						seq_reverse = seq_reverse.rjust(read_length,'N')
					seq_reverse = self._reverse_complement(seq_reverse)
					fastafile.write(">read%s_%s\n%s\n" %(entry_count, frag.frag_type, seq_reverse))
				
				entry_count += 1


#helper functions
def _progress_bar(count, total, suffix=''):
	"""Simple progress bar, takes tha actual count of completed tasls and the total tasks that need to be completed"""
	bar_len = 60
	filled_len = int(round(bar_len * count / float(total)))
	
	percents = round(100.0 * count / float(total), 1)
	bar = '=' * filled_len + '-' * (bar_len - filled_len)
	
	if(suffix == ''):
		sys.stdout.write('[%s] %s%s \r' % (bar, percents, '%'))
	else:
		sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', suffix))
	sys.stdout.flush()  

def _init_add_domains_mp(counter, parametersobj, results_list):
	"""Initialize domains for multiprocessing"""
	tmp_results_list = [] #store created object in this list
	while True:
		domain_number = counter.increment() #add 1 to the counter, one more domain done...
		if(domain_number == -1): #reached maximum end here...
			break
		sequence, valid_regions, bins, n_reps_sample1, n_reps_sample2, prot_count_n, prot_count_p, protein_size, frag_len_max, frag_dist_on, frag_dist_mn_mean, \
		frag_dist_mn_cov, frag_dist_prob, prot_dist_mn_mean, prot_dist_mn_cov, frag_count_sh, frag_count_sc, frag_count_op, frag_count_om, frag_count_os, \
		frag_count_scaling, frag_count_lp_sc, frag_count_ln_si, frag_count_ln_sc, beta_values, frag_len_mean, frag_len_dev, skewness, \
		valid_regions_array = parametersobj.get_all()
		
		domainobj = Domain(n_reps_sample1, n_reps_sample2, domain_number) #init one domain
		domainobj.init_proteins(sequence, valid_regions, bins, prot_count_n, prot_count_p, protein_size, frag_len_max, frag_dist_on, frag_dist_mn_mean, \
		frag_dist_mn_cov, frag_dist_prob, prot_dist_mn_mean, prot_dist_mn_cov, frag_count_sh, frag_count_sc, frag_count_op, frag_count_om, frag_count_os, \
		frag_count_scaling, frag_count_lp_sc, frag_count_ln_si, frag_count_ln_sc, beta_values, frag_len_mean, frag_len_dev, skewness, valid_regions_array) #init all proteins in that domain
		
		tmp_results_list.append(domainobj) #save object in tmp list
	
	results_list.send(tmp_results_list) #done so save all results...
	return

def _init_add_domains(parametersobj, n_domains_start, n_domains):
	"""Initialize domains in one process"""
	results_list = []
	domain_count = 1
	print("Initializing domains %s to %s " % (n_domains_start + 1, n_domains_start + n_domains))

	for domain_number in lrange(n_domains_start + 1, n_domains_start + n_domains + 1):
		sequence, valid_regions, bins, n_reps_sample1, n_reps_sample2, prot_count_n, prot_count_p, protein_size, frag_len_max, frag_dist_on, frag_dist_mn_mean, \
		frag_dist_mn_cov, frag_dist_prob, prot_dist_mn_mean, prot_dist_mn_cov, frag_count_sh, frag_count_sc, frag_count_op, frag_count_om, frag_count_os, \
		frag_count_scaling, frag_count_lp_sc, frag_count_ln_si, frag_count_ln_sc, beta_values, frag_len_mean, frag_len_dev, skewness, \
		valid_regions_array = parametersobj.get_all()
		
		domainobj = Domain(n_reps_sample1, n_reps_sample2, domain_number) #init one domain
		domainobj.init_proteins(sequence, valid_regions, bins, prot_count_n, prot_count_p, protein_size, frag_len_max, frag_dist_on, frag_dist_mn_mean, \
		frag_dist_mn_cov, frag_dist_prob, prot_dist_mn_mean, prot_dist_mn_cov, frag_count_sh, frag_count_sc, frag_count_op, frag_count_om, frag_count_os, \
		frag_count_scaling, frag_count_lp_sc, frag_count_ln_si, frag_count_ln_sc, beta_values, frag_len_mean, frag_len_dev, skewness, valid_regions_array) #init all proteins in that domain
		
		results_list.append(domainobj) #save object in tmp list
		_progress_bar(domain_count, n_domains)
		domain_count += 1
	
	print("\n")
	return results_list

def _get_genome(path_to_fasta, chrom):
	"""Return genome-sequence of fasta file """
	seq = None

	#check fasta for desired chromosome
	for fastaseq in HTSeq.FastaReader(path_to_fasta):
		if fastaseq.name == chrom:
			seq = str(fastaseq).upper() #set all bases to upper case
			break
		else:
			continue
	
	if seq is None:
		parser.error("The requested chromosome (%s) is not in the supplied fasta file (%s), exiting" %(chrom, path_to_fasta))
		sys.exit(1)
	else: #found chromosome
		return seq

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

def _get_weight_bins(chrom_length, back_res, back_s, back_c, chrom):
	"""Calculate weights for bins over a chromosome"""
	
	#do not care for white or blacklisted regions this will be used for noise and input only 
	weights = [] #list of weights per interval/bin
	noise_regions = [] #list of bed-tools intervals
	cur_pos = 0 #start with postion 0
	
	while ((cur_pos+back_res+1) < chrom_length):
		noise_regions.append(Interval(chrom, cur_pos, cur_pos+back_res-1, strand='*'))
		cur_pos += back_res
		cur_weight = _sample_gamma(back_s, back_c)
		weights.append(cur_weight)

	noise_regions.append(Interval(chrom, cur_pos, chrom_length, strand='*')) #add last bin
	
	weighted_bins = np.add.accumulate(weights) #weights per bin
	return weighted_bins, noise_regions, weights
	
def _sample_gamma(s, c):
	"""Sample from Gamma for region's weights for noise"""
	return gamma(s, c)

def _write_output(chromosome, min_counts, dp_thres, prefix, append = False):
	"""Compare the two samples and write out only DB peaks"""
	db_peaks_sample1 = []
	db_peaks_sample2 = []
	all_peaks = []
	
	for domain in chromosome.domain_list : #go through all domains
		for protein in domain.protein_list : #go through all proteins
			
			#get mean
			counts1 = protein.sample1_mean_fragment_count
			counts2 = protein.sample2_mean_fragment_count
			
			#check if me meet threshold
			if counts1 > counts2 and counts1 > min_counts and (counts1 / (counts1 + counts2)) >= dp_thres: #up in sample1
				db_peaks_sample1.append(protein)
				all_peaks.append((protein, 1, ("d" + str(protein.domain_name) + "p" + str(protein.name))))
			elif counts2 > counts1 and counts2 > min_counts and (counts2 / (counts1 + counts2)) >= dp_thres: #up in sample2
				db_peaks_sample2.append(protein)
				all_peaks.append((protein, 2, ("d" + str(protein.domain_name) + "p" + str(protein.name))))
			else: #not meeting user defined DP parameters or is not up in either sample
				all_peaks.append((protein, 0, ("d" + str(protein.domain_name) + "p" + str(protein.name))))
	
	#print
	print("Writing peaks to bed file: %s" % (prefix + "_sample1-peaks.bed"))
	_write_bed(chromosome, db_peaks_sample1, prefix + "_sample1-peaks.bed", append)
	
	print("Writing peaks to bed file: %s" % (prefix + "_sample2-peaks.bed"))
	_write_bed(chromosome, db_peaks_sample2, prefix + "_sample2-peaks.bed", append)
	
	#write complete peaks file
	print("Writing peaks to bed file: %s" % (prefix + "_all_peaks.bed"))
	_write_bed_all(chromosome, all_peaks, prefix + "_all_peaks.bed", append)
	
	return len(db_peaks_sample1), len(db_peaks_sample2)

def _write_bed(chromosome, protein_list, filename, append):
	"""Write bedfile with only DB peaks."""
	if(append):
		bedfile = open(filename, "a")
	else:
		bedfile = open(filename, "w")
	
	for prot in protein_list:
		bedfile.write("%s\t%s\t%s\n" %(chromosome.name, prot.start, prot.stop))

def _write_bed_all(chromosome, protein_list_plus, filename, append):
	"""Write bedfile with allpeaks."""
	if(append):
		bedfile = open(filename, "a")
	else:
		bedfile = open(filename, "w")
		
		#header write only for first batch
		bedfile.write("#Chromsome\tStart\tStop\tSample1_mean_counts\tSample2_mean_counts%s%s\tStatus\tPeak_name\n" %( \
		"".join(list(map(lambda x: '\tSample1_replicate%i_counts' % x, lrange(1,chromosome.n_reps_sample1+1)))),\
		"".join(list(map(lambda x: '\tSample2_replicate%i_counts' % x, lrange(1,chromosome.n_reps_sample2+1))))))
	
	for prot, stat, name in protein_list_plus:
		bedfile.write("%s\t%s\t%s\t%.2f\t%.2f\t%s\t%s\t%s\t%s\n" %(chromosome.name, prot.start, prot.stop, prot.sample1_mean_fragment_count, prot.sample2_mean_fragment_count,\
		 "\t".join(list(map(lambda x: str(x), prot.sample1_replicate_counts))),  "\t".join(list(map(lambda x: str(x), prot.sample2_replicate_counts))), stat, name))

def _collect_report_data(chromosome, sample1, sample2, db_peaks_sample1, db_peaks_sample2, report_obj):
	"""Collect report pdf and extend report object."""
	
	#Get numbers:
	r_prot_pos, r_dom_nprot, r_prot_offset_positions, r_nfrag = [], [], [], []
	r_n_frag1, r_n_frag2, r_fragshift, r_fraglen = [], [], [], []
	r_n_frag1_reps, r_n_frag2_reps = [[] for i in lrange(chromosome.n_reps_sample1)], [[] for i in lrange(chromosome.n_reps_sample2)]
	
	print("Collecting data for report")
	domain_count = 0
	n_domains = len(chromosome.domain_list)
	for domain in chromosome.domain_list : #go through all domains
		r_dom_nprot.append(domain.n_proteins)
		r_prot_offset_positions.extend(domain.protein_offset_positions)
		
		for protein in domain.protein_list : #go through all proteins
			r_prot_pos.append(protein.position)
			r_nfrag.append(protein.number_frags)
			r_n_frag1.append(protein.sample1_mean_fragment_count)
			r_n_frag2.append(protein.sample2_mean_fragment_count)
			
			for rep in lrange(chromosome.n_reps_sample1):#go through all replicates
				r_n_frag1_reps[rep].append(len(protein.sample1_replicate_list[rep]))#add fragments to corresonding index of replicate 
				
				for fragment in protein.sample1_replicate_list[rep]: #all fragments
					r_fragshift.append(fragment.rand_shift)
					r_fraglen.append(fragment.frag_length)
					
			for rep in lrange(chromosome.n_reps_sample2):#go through all replicates
				r_n_frag2_reps[rep].append(len(protein.sample2_replicate_list[rep]))#add fragments to corresonding index of replicate
				
				for fragment in protein.sample2_replicate_list[rep]: #all fragments
					r_fragshift.append(fragment.rand_shift)
					r_fraglen.append(fragment.frag_length)
		domain_count += 1
		_progress_bar(domain_count, n_domains)
	

	#store data in report object
	report_obj._extend(n_domains, r_prot_pos, r_dom_nprot, r_prot_offset_positions, r_nfrag, r_n_frag1, r_n_frag2, r_fragshift, r_fraglen, r_n_frag1_reps, r_n_frag2_reps, db_peaks_sample1, db_peaks_sample2)

def _write_report(chromosome, prefix, noise_regions, bin_weights, spike_in_chrom, spike_in_regions, spike_in_bin_weights, sample1, sample2, report_obj):
	"""Write report pdf."""

	#write report:
	print("Writing report to: " +  prefix + "_report.pdf")
	
	with PdfPages(prefix + "_report.pdf") as pdf:
		#start with table
		plt.style.use("ggplot")

		fig, ax = plt.subplots()
		ax.xaxis.set_visible(False) #hide axis
		ax.yaxis.set_visible(False)
		table_data = [["Domains", report_obj.n_domains]]
		table_data.append(["Proteins", len(report_obj.r_prot_pos)])
		table_data.append(["Fragments in DB-peaks", len(report_obj.r_fraglen)])
		table_data.append(["Replicates sample1", chromosome.n_reps_sample1])
		table_data.append(["Replicates sample2", chromosome.n_reps_sample2])
		table_data.append(["Upregulated peaks in sample1", report_obj.db_peaks_sample1])
		table_data.append(["Upregulated peaks in sample2", report_obj.db_peaks_sample2])
		table_data.append(["Input noise-fragments in sample1", sample1.myinput.n_fragments])
		table_data.append(["Input spike-in-fragments in sample1", sample1.myinput.n_spike_in_fragments])
		table_data.append(["Input noise-fragments in sample2", sample2.myinput.n_fragments])
		table_data.append(["Input spike-in-fragments in sample1", sample2.myinput.n_spike_in_fragments])
		
		#check samples for input und spike-in
		rep_count = 0 #need this to get index of replicate in report object
		for rep in sample1.replicate_list:
			
			#table_data.append(["Peak fragments in sample1 replicate" + str(rep.replicate_id), rep.n_peak_fragments])
			table_data.append(["Peak fragments in sample1 replicate" + str(rep.replicate_id), sum(report_obj.r_n_frag1_reps[rep_count])])
			table_data.append(["Noise fragments in sample1 replicate" + str(rep.replicate_id), rep.n_noise_fragments])
			table_data.append(["Spike-in fragments in sample1 replicate" + str(rep.replicate_id), rep.n_spike_in_fragments])
			rep_count += 1
			
		rep_count = 0
		for rep in sample2.replicate_list:
			#table_data.append(["Peak fragments in sample2 replicate" + str(rep.replicate_id), rep.n_peak_fragments])
			table_data.append(["Peak fragments in sample2 replicate" + str(rep.replicate_id), sum(report_obj.r_n_frag2_reps[rep_count])])
			table_data.append(["Noise fragments in sample2 replicate" + str(rep.replicate_id), rep.n_noise_fragments])
			table_data.append(["Spike-in fragments in sample2 replicate" + str(rep.replicate_id), rep.n_spike_in_fragments])
			rep_count += 1
			
		collabel=("Description", "#")
		ax.table(cellText=table_data, colLabels=collabel, loc='center', rowLoc='right')
		pdf.savefig()
		plt.close()
		
		#and now the plots...
		nbins = [] #create bins out of interval stops
		for intervall in noise_regions:
			nbins.append(intervall.stop)
		prot_counts = np.histogram(report_obj.r_prot_pos, nbins)[0] #sum up positions to bins
		nbins = nbins[:-1] #trim last element
		
		plt.plot(nbins, prot_counts, color = 'k')
		plt.title("Histogram of domain positions over %s" % (chromosome.name))
		plt.ylabel("Number of domains")
		plt.xlabel("Positions in %s" % (chromosome.name))
		plt.axis([0,max(nbins),0,max(prot_counts)])
		pdf.savefig()
		plt.close()

		plt.hist(report_obj.r_dom_nprot, bins=50, color = 'b')
		plt.title("Histogram of proteins per domain")
		plt.ylabel("Count of proteins with same range")
		plt.xlabel("Range of proteins")
		pdf.savefig()
		plt.close()

		plt.hist(report_obj.r_prot_offset_positions, bins=50, color = 'r')
		plt.title("Histogram of protein distances inside the domains")
		plt.ylabel("Count of protein distances")
		plt.xlabel("Range of protein distances")
		pdf.savefig()
		plt.close()

		plt.hist(report_obj.r_nfrag, bins=50, color = 'g')
		plt.title("Histogram of fragments per protein")
		plt.ylabel("Count of fragment numbers")
		plt.xlabel("Range of fragment numbers") #neg bin values / 2 samples * number of replicates
		pdf.savefig()
		plt.close()

		plt.hist(report_obj.r_n_frag1, bins=50, histtype='stepfilled', color='y', label='Sample1')
		plt.hist(report_obj.r_n_frag2, bins=50, histtype='stepfilled', color='m', alpha=0.5, label='Sample2')
		legend = plt.legend(loc='upper right')
		plt.title("Histogram of fragments per protein in sample1 and sample2")
		plt.ylabel("Count of fragment numbers")
		plt.xlabel("Range of fragment numbers")
		pdf.savefig()
		plt.close()
		
		plt.hist(report_obj.r_fragshift, bins=50, color = 'g')
		plt.title("Histogram of fragment position shifts")
		plt.ylabel("Count of fragment position shifts")
		plt.xlabel("Range of fragment position shifts")
		pdf.savefig()
		plt.close()
		
		plt.hist(report_obj.r_fraglen, bins=50, color = 'g')
		plt.title("Histogram of fragment length")
		plt.ylabel("Count of fragment length")
		plt.xlabel("Range of fragment length")
		pdf.savefig()
		plt.close()

		plt.plot(nbins, list(bin_weights), color = 'k')
		plt.title("Noise bin positions over %s" % (chromosome.name))
		plt.ylabel("Weights")
		plt.xlabel("Positions in %s" % (chromosome.name))
		plt.axis([0,max(nbins),0,max(list(bin_weights))])
		pdf.savefig()
		plt.close()
		
		if(spike_in_chrom.name != ""):
			sibins = [] #create bins out of interval stops
			for intervall in spike_in_regions:
				sibins.append(intervall.stop)
			sibins = sibins[:-1] #trim last element
		
			plt.plot(sibins, list(spike_in_bin_weights), color = 'k')
			plt.title("Spike-in bin positions over %s" % (spike_in_chrom.name))
			plt.ylabel("Weights")
			plt.xlabel("Positions in %s" % (spike_in_chrom.name))
			plt.axis([0,max(sibins),0,max(list(spike_in_bin_weights))])
			pdf.savefig()
			plt.close()
		
		for rep in lrange(chromosome.n_reps_sample1):
			plt.hist(report_obj.r_n_frag1_reps[rep], bins=50, histtype='stepfilled', alpha=0.5, label='Sample1_rep'+str(rep+1))
		legend = plt.legend(loc='upper right')
		plt.title("Histogram of fragments per protein in replicates of sample1")
		plt.ylabel("Count of fragment numbers")
		plt.xlabel("Range of fragment numbers")
		pdf.savefig()
		plt.close()

		for rep in lrange(chromosome.n_reps_sample2):
			plt.hist(report_obj.r_n_frag2_reps[rep], bins=50, histtype='stepfilled', alpha=0.5, label='Sample2_rep'+str(rep+1))
		legend = plt.legend(loc='upper right')
		plt.title("Histogram of fragments per protein in replicates of sample2")
		plt.ylabel("Count of fragment numbers")
		plt.xlabel("Range of fragment numbers")
		pdf.savefig()
		plt.close()
		
		#pdf metadata via the PdfPages object:
		d = pdf.infodict()
		d['Title'] = 'DB-ChIP-Sim PDF Report'
		d['Author'] = 'Thomas Eder'
		d['Subject'] = 'Report of simulated ChIP-seq data'
		d['Keywords'] = 'in silico ChIP-seq data'
		d['CreationDate'] = datetime.datetime.today()
		d['ModDate'] = datetime.datetime.today()

#main
if __name__ == '__main__':
	#handle input arguments
	parser = HelpfulOptionParser(usage=__doc__)
	
	parser.add_option("-c", "--chrom", default='chr1', dest="chrom", type="string", help="Chromosome used for simulation [default: %default]")
	parser.add_option("-d", "--domain-counts", default=1000, dest="n_domains", type="int", help="Number of domains with DPs [default: %default]")
	parser.add_option("-l", "--length", default=50, dest="read_length", type="int", help="Read length [default: %default]")
	parser.add_option("-p", "--prefix", default='sim', dest="prefix", type="string", help="Prefix for output files [default: %default]")
	
	parser.add_option("--prot-size", default=150, dest="protein_size", type="int", help="Protein size (this parameter limits the fragment shifts) [default: %default]")
	
	parser.add_option("--prot-count-n",default=1, dest="prot_count_n", type="float",help="n of negative binomial distribution for protein counts [default: %default]")
	parser.add_option("--prot-count-p",default=0.9, dest="prot_count_p", type="float",help="p of negative binomial distribution for protein counts [default: %default]")
	
	parser.add_option("--prot-dist-muno-mean", default=[300, 1800, 900], dest="prot_dist_mn_mean", type="string", action='callback', callback=_callback_list,\
	help="Means of multivariate normal distribution for distances between proteins, separator: ',' eg. \"300, 1800\" [default: %default]")
	parser.add_option("--prot-dist-muno-cov", default=[[1000,0,0],[0,5000,0],[0,0,5000]], dest="prot_dist_mn_cov", type="string", action='callback', callback=_callback_matrix,\
	help="Covariance of multivariate normal distribution for distances between proteins, separator: ',' and ';'  eg. \"1000,0;0,5000\" [default: %default]")
	
	parser.add_option("--frag-len-mean", default=200, dest="frag_len_mean", type="float", help="Set mean of fragments' length [default: %default]")
	parser.add_option("--frag-len-dev", default=20, dest="frag_len_dev", type="float", help="Set deviation of fragments' length [default: %default]")
	parser.add_option("--frag-len-max", default=1000, dest="frag_len_max", type="int", help="Set maximum of fragments' length [default: %default]")
	
	parser.add_option("--frag-count-sh", default=2.2, dest="frag_count_sh", type="float", help="Shape of gamma distribution for fragment counts [default: %default]")
	parser.add_option("--frag-count-sc", default=20.1, dest="frag_count_sc", type="float", help="Scale of gamma distribution for fragment counts [default: %default]")
	parser.add_option("--frag-count-op", default=0.01, dest="frag_count_op", type="float", help="Probability for fragment counts being outliers [default: %default]")
	parser.add_option("--frag-count-om", default=6., dest="frag_count_om", type="float", help="Mean of lognormal distribution for fragment counts of outliers [default: %default]")
	parser.add_option("--frag-count-os", default=0.5, dest="frag_count_os", type="float", help="Sigma of lognormal distribution for fragment counts of outliers [default: %default]")
	
	parser.add_option("--frag-count-scaling", default="none", dest="frag_count_scaling", type="string", help="Scaling of fragment distribution, no scaling, scaling of beta result based on fragment counts (with lognorm) or scaling of fragment counts based on beta result (with laplace) : none , frag , beta [default: %default]")
	parser.add_option("--frag-count-lp-scale", default=0.1, dest="frag_count_lp_sc", type="float", help="Scale for Laplace distribution if frag-count-scaling is beta [default: %default]")
	parser.add_option("--frag-count-ln-sigma", default=0.9, dest="frag_count_ln_si", type="float", help="Sigma for lognorm distribution if frag-count-scaling is frag [default: %default]")
	parser.add_option("--frag-count-ln-scale", default=100, dest="frag_count_ln_sc", type="float", help="Scale for lognorm distribution if frag-count-scaling is frag [default: %default]")
	
	parser.add_option("--frag-dist-on", default=False, action="store_true", dest="frag_dist_on", help="Use multivariate normal distribution for fragment shifts to create peak chapes, shifts are limited by prot-size. The final shift is: postion of peak - prot_size + sampling from distribution [default: %default]")
	parser.add_option("--frag-dist-prob", default=[0.5, 0.5], dest="frag_dist_prob", type="string", action='callback', callback=_callback_list_float,\
	help="Probability fo each of the multivariate normal distributions to be chosen [default: %default]")
	parser.add_option("--frag-dist-muno-mean", default=[20, 100], dest="frag_dist_mn_mean", type="string", action='callback', callback=_callback_list,\
	help="Means of multivariate normal distribution for the shifts of fragments, separator: ',' eg. \"300, 1800\" [default: %default]")
	parser.add_option("--frag-dist-muno-cov", default=[[100,0],[0,500]], dest="frag_dist_mn_cov", type="string", action='callback', callback=_callback_matrix,\
	help="Covariance of multivariate normal distribution for the shifts of fragments, separator: ',' and ';'  eg. \"1000,0;0,5000\" [default: %default]")
	
	parser.add_option("--beta", default=[0.5, 0.5], dest="beta_values", type="string", action='callback', callback=_callback_list_float, help="Alpha and Beta of Beta-distribution [default: %default]")
	
	parser.add_option("--dp-thres", default=0.6, dest="dp_thres", type="float", help="Threshold of reads/fragments to define a DB peak [default: %default]")
	parser.add_option("-m", "--min-counts", default=25, dest="min_counts", type="int", help="Minimum number of reads/fragments for a DB peak [default: %default]")
	
	parser.add_option("-s", "--skewness", default=10, dest="skewness", type="float", help="Variance between replicates (the higher, the less variance) [default: %default]")
	
	parser.add_option("--back-avg", default=0.25, dest="back_avg", type="float", help="Average background coverage for noise [default: %default]")
	parser.add_option("--back-res", default=1000, dest="back_res", type="int", help="Resolution for ChIP-seq noise estimates (also used for spike-in and input) [default: %default]")
	parser.add_option("--back-c", default=20, dest="back_c", type="int", help="Gamma distribution scale (theta) for noise model (also used for spike-in and input) [default: %default]")
	parser.add_option("--back-s", default=1, dest="back_s", type="int", help="Gamma distribution shape (k) for noise model (also used for spike-in and input) [default: %default]")
	
	parser.add_option("--is-white-list", default=False, action="store_true", dest="is_white_list", help="Set provided bed-file to white-list and alow only DPs in these regions [default: %default]")
	
	parser.add_option("--no-input", default=False, action="store_true", dest="no_input", help="Create no input-fasta per sample [default: %default]")
	parser.add_option("--no-fasta", default=False, action="store_true", dest="write_no_fasta", help="Do not create fasta files, only the report will be created [default: %default]")
	parser.add_option("--no-report", default=False, action="store_true", dest="write_no_report", help="Do not create pdf report [default: %default]")
	parser.add_option("--no-noise", default=False, action="store_true", dest="no_noise", help="Do not add noise to the simulated reads/peaks [default: %default]")
	
	parser.add_option("--spike-in", default=False, action="store_true", dest="spike_in", help="Add spike-in reads from a defined reference fasta [default: %default]")
	parser.add_option("--si-fasta", default='spike_in.fasta', dest="spike_in_fastafile_path", type="string", help="Path to reference fasta file for spike in [default: %default]")
	parser.add_option("--si-chrom", default='chr1', dest="spike_in_chrom", type="string", help="Chromosome of spike-in to simulate. [default: %default]")
	parser.add_option("--si-cov", default=0.25, dest="spike_in_cov", type="float", help="Average background coverage for spike-in [default: %default]")
	
	parser.add_option("-t", "--threads", default=1, dest="threads", type="int", help="Number of threads to use [default: %default]")
	parser.add_option("--batch-size", default=10000, dest="batch_size", type="int", help="Number of domains calculated in one batch to limit memory usage [default: %default]")
	
	(options, args) = parser.parse_args()
	num_args = 3
	
	#read in required input arguments
	if len(args) != num_args:
		parser.error("At minimum %s parameters are requred to run this tool" % num_args)
	fastafile_path = args[0]
	bedfile_path = args[1]
	count_rep = int(args[2])
	
	#some sanity checks on input arguments
	if not (options.frag_count_op >= 0 and options.frag_count_op <= 1):
		parser.error("Probability for fragment counts being outliers %s needs to be >= 0 and <= 1" % (options.frag_count_op))
	if not (options.frag_count_sh > 0):
		parser.error("Shape of gamma distribution for fragment counts %s needs to be > 0" % (options.frag_count_sh))
	if not (options.frag_count_sc > 0):
		parser.error("Scale of gamma distribution for fragment counts %s needs to be > 0" % (options.frag_count_sc))
	if not (options.frag_count_om > 0):
		parser.error("Mean of lognormal distribution for fragment counts of outliers %s should be positive" % (options.frag_count_om))
	if not (options.frag_count_os > 0):
		parser.error("Sigma of lognormal distribution for fragment counts of outliers %s needs to be > 0" % (options.frag_count_os))
	
	if not (options.frag_count_scaling == "none" or options.frag_count_scaling == "frag" or options.frag_count_scaling == "beta"):
		parser.error("frag_count_scaling %s must ether be 'none', 'frag' or 'beta'" % (options.frag_count_scaling))
	
	if not (options.frag_len_max > options.frag_len_mean): 
		parser.error("Fragment length max %s needs to be bigger than the mean %s" % (options.frag_len_max, options.frag_len_mean))
	if (options.min_counts < 1):
		parser.error("The minimum number of reads/fragments for a DB peak %s needs to be higher than 0" % options.min_counts)
	
	if not (len(options.beta_values) == 2 and options.beta_values[0] > 0  and options.beta_values[1] > 0 ):
		parser.error("Please check your Alpha and Beta parameters for the Beta-distribution %s it should be similar to [0.5, 0.5], strictly two entries" % options.beta_values)
	if not (len(options.prot_dist_mn_mean) >= 2):
		parser.error("Please check your means of multivariate normal distribution for distances between proteins %s it should be similar to [300, 1800, 900], two or more entries" % options.prot_dist_mn_mean)
	if not (len(options.prot_dist_mn_cov) == len(options.prot_dist_mn_mean)):
		parser.error("Please check your covariance of multivariate normal distribution for distances between proteins %s it should be similar to [[1000, 0, 0], [0, 5000,0], [0, 0, 5000]], the covariance lists need to be equal to the number of supplied means of multivariate normal distribution %s" % (options.prot_dist_mn_cov,options.prot_dist_mn_mean))
	for cov in options.prot_dist_mn_cov:
		if not (len(cov) == len(options.prot_dist_mn_mean)):
			parser.error("Please check your covariance of multivariate normal distribution for distances between proteins it should be similar to [[1000, 0, 0], [0, 5000,0], [0, 0, 5000]] covariances for [300, 1800, 900] means, the covariance lists need to be equal to the number of supplied means of multivariate normal distribution, your covariances: %s and means: %s" % (options.prot_dist_mn_cov,options.prot_dist_mn_mean))
	
	if not (len(options.frag_dist_mn_cov) == len(options.frag_dist_mn_mean)):
		parser.error("Please check your covariance of multivariate normal distribution for the shifts of fragments %s it should be similar to [[1000, 0, 0], [0, 5000,0], [0, 0, 5000]], the covariance lists need to be equal to the number of supplied means of multivariate normal distribution %s" % (options.frag_dist_mn_cov,options.frag_dist_mn_mean))
	if not (len(options.frag_dist_prob) == len(options.frag_dist_mn_mean)):
		parser.error("Please check your probability parameters for the shifts of fragments %s the number of probailities should be the same as the number of supplied means of multivariate normal distribution %s" % (options.frag_dist_prob,options.frag_dist_mn_mean))
	if not (sum(options.frag_dist_prob) == 1):
		parser.error("Please check your probability parameters for the shifts of fragments %s the sum of all probabilities needs to be 1" % options.frag_dist_prob)
	for cov in options.frag_dist_mn_cov:
		if not (len(cov) == len(options.frag_dist_mn_mean)):
			parser.error("Please check your covariance of multivariate normal distribution for the shifts of fragments it should be similar to [[1000, 0, 0], [0, 5000,0], [0, 0, 5000]] covariances for [300, 1800, 900] means, the covariance lists need to be equal to the number of supplied means of multivariate normal distribution, your covariances: %s and means: %s" % (options.prot_dist_mn_cov,options.prot_dist_mn_mean))
	
	if not (os.path.isfile(fastafile_path)):
		parser.error("File %s does not exist" % fastafile_path)
	if not (os.path.isfile(bedfile_path)):
		parser.error("File %s does not exist" % bedfile_path)
	if (options.spike_in and not os.path.isfile(options.spike_in_fastafile_path)):
		parser.error("File %s does not exist" % options.spike_in_fastafile_path)
	
	if not (options.spike_in_cov > 0):
		parser.error("Spike-in coverage %s needs to be higher than 0 " % options.spike_in_cov)
	if not (options.back_avg > 0 ):
		parser.error("The average background coverage %s needs to be higher than 0 " % options.back_avg)
	if not (options.dp_thres > 0 and options.dp_thres < 1):
		parser.error("The threshold of reads/fragments to define a DB peak %s needs to be between 0 and 1" % options.dp_thres)
	
	if not (options.threads > 0):
		parser.error("The number of threads %s needs to be bigger than 0" % options.threads)
	if (options.threads > mp.cpu_count()):
		print("\nYou are intending to use more threads (%s) than availible CPUs (%s), this might be disadvantageous for the runtime" % (options.threads, mp.cpu_count()))
	
	prog_timer = Timer() #start timer
	
	try: #now get going..

		#load chromosome
		genome_seq = _get_genome(fastafile_path, options.chrom)
		
		chromosomeobj = Chromosome(genome_seq, options.chrom, bedfile_path, count_rep, count_rep, options.is_white_list, options.threads)
		
		#calculate weighted bins on complete chromosome
		(weighted_bins, noise_regions, bin_weights) = _get_weight_bins(chromosomeobj.length,options.back_res, options.back_s, options.back_c, options.chrom) 
		
		report_obj = Report_data(count_rep, count_rep) #object to store report data
		
		domain_batches = []
		cur_pos = 1
		
		while ((cur_pos + options.batch_size) < options.n_domains):
			domain_batches.append(options.batch_size)
			cur_pos += options.batch_size
		
		domain_batches.append(options.n_domains - cur_pos + 1) #last batch

		is_first = True #first batch to create files and headers
		domain_count = 0 #count total created domains
		batch_count = 1
		read_count1 = 1 #make sure all fasta entries have the right number
		read_count2 = 1


		for n_domains in domain_batches : #run batches for db-peaks to limit maximal memory usage
			strlen1 = len(str(batch_count)) #format the output string somewhat appealing
			strlen2 = len(str(len(domain_batches)))
			strlenf1 = max((70 - (23 + strlen1 + strlen2)) // 2, 1)
			strlenf2 = max(69 - (23 + strlen1 + strlen2 + strlenf1), 1)
			print("\n\n%s Working on batch %s of %s %s \n" %(('*'*strlenf1), batch_count,len(domain_batches),('*'*strlenf2)))
			
			#create domains, proteins and fragments
			chromosomeobj.init_domains(domain_count, n_domains, options.prot_count_n, options.prot_count_p, options.protein_size, options.frag_len_max,\
			options.frag_dist_on, options.frag_dist_mn_mean, options.frag_dist_mn_cov,  options.frag_dist_prob, options.prot_dist_mn_mean, options.prot_dist_mn_cov,\
			options.frag_count_sh, options.frag_count_sc, options.frag_count_op, options.frag_count_om, options.frag_count_os, options.frag_count_scaling, \
			options.frag_count_lp_sc, options.frag_count_ln_si, options.frag_count_ln_sc, options.beta_values, options.frag_len_mean, options.frag_len_dev, \
			options.skewness)
			
			domain_count += n_domains

			#save results (objects of chromosome, proteins and fragments) into sample1 and sample2
			print("Distributing fragments to samples and replicates")
			sample1 = Sample(chromosomeobj, count_rep, "1", weighted_bins, noise_regions)
			sample2 = Sample(chromosomeobj, count_rep, "2", weighted_bins, noise_regions)
			
			#go through chromosome and get proteins to get their counts/replicate distribution
			d_count = 0
			
			for domain in chromosomeobj.domain_list : #go through all proteins
				for protein in domain.protein_list : #for each replicate
					for rep in lrange(count_rep) : #add fragments to corresonding index of replicate
						sample1.add_db_fragments(protein.sample1_replicate_list[rep], rep) 
						sample2.add_db_fragments(protein.sample2_replicate_list[rep], rep)
				d_count += 1
				_progress_bar(d_count, n_domains)
			print("\n")
			
			#append to a file or create a new one...
			if(is_first):
				append = False
				is_first = False
			else:
				append = True
			
			#now write peak fragments to fasta
			if not options.write_no_fasta:
				read_count1 = sample1.write_fasta(options.prefix, options.read_length, read_count1, append)
				read_count2 = sample2.write_fasta(options.prefix, options.read_length, read_count2, append)
				print("") 
			
			#write bed files
			db_peaks_sample1, db_peaks_sample2 = _write_output(chromosomeobj, options.min_counts, options.dp_thres, options.prefix, append)
			print("")
			
			#collect report data
			if not options.write_no_report:
				_collect_report_data(chromosomeobj, sample1, sample2, db_peaks_sample1, db_peaks_sample2, report_obj)
				print("")
			
			#wipe db fragments in chromosome and samples and replicates
			chromosomeobj.wipe_data()
			del(sample1)
			del(sample2)
			
			batch_count += 1
		print("\n\n%s Finishing %s \n" % ("*"*29,"*"*29))

		#empty samples for noise and input
		sample1 = Sample(chromosomeobj, count_rep, "1", weighted_bins, noise_regions)
		sample2 = Sample(chromosomeobj, count_rep, "2", weighted_bins, noise_regions)
		
		#add noise
		if not options.no_noise:
			sample1.add_noise(options.back_avg, options.read_length, options.frag_len_mean, options.frag_len_dev, options.frag_len_max)
			sample2.add_noise(options.back_avg, options.read_length, options.frag_len_mean, options.frag_len_dev, options.frag_len_max)
			print("")
		
		#create input per sample
		if not options.no_input:
			sample1.add_input(options.frag_len_mean, options.frag_len_dev, options.frag_len_max, options.back_avg, options.read_length)
			sample2.add_input(options.frag_len_mean, options.frag_len_dev, options.frag_len_max, options.back_avg, options.read_length)
			print("")
		
		#add spike-in
		if options.spike_in:
			#spike in fragments/reads are created like noise fragments just on a different genome/chromosome 
			spike_in_seq = _get_genome(options.spike_in_fastafile_path, options.spike_in_chrom)
			spike_in_chromobj = SpikeInChromosome(spike_in_seq, options.spike_in_chrom)
			(spike_in_weighted_bins, spike_in_regions, spike_in_bin_weights) = _get_weight_bins(spike_in_chromobj.length, options.back_res, options.back_s, \
			options.back_c, options.chrom) 
			
			sample1.add_spike_in(spike_in_chromobj, spike_in_weighted_bins, spike_in_regions, options.spike_in_cov, options.read_length, options.frag_len_mean,\
			options.frag_len_dev, options.frag_len_max)
			sample2.add_spike_in(spike_in_chromobj, spike_in_weighted_bins, spike_in_regions, options.spike_in_cov, options.read_length, options.frag_len_mean,\
			options.frag_len_dev, options.frag_len_max)
			print("")
		else:
			spike_in_chromobj = SpikeInChromosome("", "") #just a dummy spike-in-chromosome object that _write_report will work
			spike_in_regions = []
			spike_in_bin_weights = []
		
		#write spike-in and noise to fasta
		if not options.write_no_fasta:
			sample1.write_fasta(options.prefix, options.read_length, read_count1, True) #true because we want to append the noise reads
			sample2.write_fasta(options.prefix, options.read_length, read_count2, True)
			print("") 
		
		#write report
		if not options.write_no_report:
			_write_report(chromosomeobj, options.prefix, noise_regions, bin_weights, spike_in_chromobj, spike_in_regions, spike_in_bin_weights, sample1, sample2, report_obj)
			print("")
		
	except KeyboardInterrupt:
		print("Interrupted by user")

	time_hhmmss = prog_timer.get_time_hhmmss()
	print ("Time elapsed for simulation: %s" % (time_hhmmss))
	exit(0)
