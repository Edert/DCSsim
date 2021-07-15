# Simulate differential ChIP-seq regions

DCSsim simulates differential ChIP-seq data for two samples (e.g. treatment and control). The results are sequence reads in fasta format for the two samples with n (user defined) replicates plus a report as pdf. DCSsim supports multiple threads and is able to work in batches to account for limited memory. We strongly encourage new users to start with the two parameter simulation scripts test-peak-shape.py and test-fragment-count-distribution.py. These can be applied to test the parameters for sampling without creating reads but output diverse histograms and plots.


## DCSsim ##
**DCSsim.py** \<FASTA\> \<BED\> \<INT\> [options]

Based on the provided \<FASTA\> sequence differential peaks (DPs) are simulated, restricted to if --is_white_list or constrained by the specified \<BED\>-file. \<INT\> defines the number of simulated replicates of the two samples.

Options:

-c, --chrom	Chromosome used for simulation, default='chr1',

-d, --domain-counts	Number of domains/clusters with DPs, default=1000

-l, --length 	Read length, default=50

-p, --prefix 	Prefix for output files, default='sim'

--prot-size 	Protein size (this parameter limits the fragment shifts), default=150

--prot-count-n	n of negative binomial distribution for protein counts, default=1

--prot-count-p	p of negative binomial distribution for protein counts, default=0.9

--prot-dist-muno-mean	Means of multivariate normal distribution for distances between proteins, separator: ',' eg. "300, 1800", default=300, 1800, 900

--prot-dist-muno-cov	Covariance of multivariate normal distribution for distances between proteins, separator: ',' and ';'  eg. "1000,0;0,5000", default=[[1000,0,0],[0,5000,0],[0,0,5000]]

--frag-len-mean	Set mean of fragments' length, default=200

--frag-len-dev	Set deviation of fragments' length, default=20

--frag-len-max	Set maximum of fragments' length, default=1000

--frag-count-sh	Shape of gamma distribution for fragment counts, default=2.2

--frag-count-sc	Scale of gamma distribution for fragment counts, default=20.1

--frag-count-op	Probability for fragment counts being outliers, default=0.01

--frag-count-om	Mean of lognormal distribution for fragment counts of outliers, default=6.0

--frag-count-os	Sigma of lognormal distribution for fragment counts of outliers, default=0.5

--frag-count-scaling	Scaling of fragment distribution, no scaling, scaling of beta result based on fragment counts (with exponential distribution) or scaling of fragment counts based on beta result (with Laplace distribution): none, frag, beta, default="none"

--frag-count-lp-scale	Scale for Laplace distribution if frag-count-scaling is frag, default=0.1

--frag-count-ex-loc	Loc for exponential distribution if frag-count-scaling is beta, default=10

--frag-count-ex-scale	Scale for exponential distribution if frag-count-scaling is beta, default=100
	
--frag-dist-on	Use multivariate normal distribution for fragment shifts to create peak shapes, shifts are limited by prot-size. The final shift is: postion of peak - prot_size + sampling from distribution, default=False 

--frag-dist-prob, dest	Probability for each of the multivariate normal distributions to be chosen, default=[0.5, 0.5]

--frag-dist-muno-mean	Means of multivariate normal distribution for the shifts of fragments, separator: ',' eg. "300, 1800", default=[20, 100]

--frag-dist-muno-cov	Covariance of multivariate normal distribution for the shifts of fragments, separator: ',' and ';'  eg. "1000,0;0,5000", default=[[100,0],[0,500]]

--beta 	Alpha and Beta of Beta-distribution, default=[0.5, 0.5]	

--dp-thres	Threshold of reads/fragments to define a DB peak, default=0.6

-m, --min-counts	Minimum number of reads/fragments for a DB peak, default=25

-s, --skewness	Variance between replicates (the higher, the less variance), default=10

--back-avg	Average background coverage for noise, default=0.25

--back-res	Resolution for ChIP-seq noise estimates (also used for spike-in and input), default=1000

--back-c		Gamma distribution scale (theta) for noise model (also used for spike-in and input), default=20

--back-s		Gamma distribution shape (k) for noise model (also used for spike-in and input), default=1

--is-white-list	Set provided bed-file to white-list and alow only DPs in these regions ", default=False

--no-input	Create no input/control-fasta per sample, default=False

--no-fasta	Do not create fasta files, only the report will be created, default=False

--no-report	Do not create pdf report, default=False

--no-noise	Do not add noise to the simulated reads/peaks, default=False

--spike-in	Add spike-in reads from a defined reference fasta, default=False

--si-fasta	Path to reference fasta file for spike in, default='spike_in.fasta'

--si-chrom	Chromosome of spike-in reference fasta file to simulate, default='chr1'

--si-cov		Background coverage for spike-in, default=0.25

-t, --threads	Number of threads to use, default=1

--batch-size	Number of domains/clusters calculated in one batch to limit memory usage, default=10000

## test-fragement-count-distribution ##
**test-fragement-count-distribution.py** \<INT\> \<INT\> [options]
	
Creates a histogram of \<INT\> protein interactions sites of \<INT\> replicates for two samples, for all samples, for sample1 and sample2 and creates an MA-plot based on the supplied beta-values
	
Options:
	
--frag-count-sh	Shape of gamma distribution for fragment counts, default=2.2
	
--frag-count-sc	Scale of gamma distribution for fragment counts, default=20.1
	
--frag-count-op	Probability for fragment counts being outliers, default=0.01
	
--frag-count-om	Mean of lognormal distribution for fragment counts of outliers, default=6.0
	
--frag-count-os	Sigma of lognormal distribution for fragment counts of outliers, default=0.5
	
--frag-count-scaling	Scaling of fragment distribution, no scaling, scaling of beta result based on fragment counts (with exponential distribution) or scaling of fragment
counts based on beta result (with Laplace distribution): none, frag, beta, default="none"
	
--frag-count-lp-scale	Scale for Laplace distribution if frag-count-scaling is frag, default=0.1
	
--frag-count-ex-loc	Loc for exponential distribution if frag-count-scaling is beta, default=10
	
--frag-count-ex-scale	Scale for exponential distribution if frag-count-scaling is beta, default=100
	
--beta	Alpha and Beta of Beta-distribution, default=[0.5, 0.5]

## test-peak-shape ##

**test-peak-shape.py** \<INT\> [options]
	
Creates histograms of \<INT\> fragment distances and peak shapes
	
Options:
	
--prot-size	Protein size (this parameter limits the fragment shifts), default=150
	
--frag-len-mean	Set mean of fragments' length, default=200
	
--frag-len-dev	Set deviation of fragments' length, default=20
	
--frag-len-max	Set maximum of fragments' length , default=1000	
	
--frag-dist-on	Use multivariate normal distribution for fragment shifts to create peak shapes, shifts are limited by prot-size. The final shift is: position of peak - prot_size + sampling from distribution, default=False
	
--frag-dist-prob	Probability for each of the multivariate normal distributions to be chosen, default=[0.5, 0.5]
	
--frag-dist-muno-mean	Means of multivariate normal distribution for the shifts of fragments, separator: ',' eg. "300, 1800", default=[20, 100]
	
--frag-dist-muno-cov	Covariance of multivariate normal distribution for the shifts of fragments, separator: ',' and ';'  eg. "1000,0;0,5000", default=[[100,0],[0,500]]


