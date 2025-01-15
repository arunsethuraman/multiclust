# MULTICLUST
Multinomial clustering of multi-allelic, multi-locus population genomic data to assign population structure.

# Installation

MULTICLUST is written in the C programming language, and will therefore need to be compiled prior to execution. To install from source, download the GitHub repository to your machine and install via:

git clone https://github.com/arunsethuraman/multiclust/
cd multiclust
make

This should create an executable named "multiclust" inside the folder.

# Synopsis

	multiclust -k < n > | -1 < n > -2 < n > [-a -b -B < n > -bou < d > -c -C < n > -d < s > -e < d > -E < d > -f < d > -g < d > -h
		-i < n > -I -m < n > --missing < n > -M < n > -n < n > -o < s > -p < n > --projection --plus -R -s < n > -t < n > -T < d > -u < s > -v -w < s > -x -f < s > --format < s >
		--simulate <qfile> <pfile> <ofile>

		where < n > stands for integer, < s > for string, < d > for double
  
# Description

MULTICLUST clusters multivariate discrete data (here multi-locus genomic data, multi-allelic, polyploid) observed on a sample of individuals using the EM algorithm (both accelerated and unaccelerated).  It handles data missing at random.  It assumes loci within an individual are independent (i.e. in linkage equilibrium). It allows the admixture model, where each locus is independently drawn from a cluster, or the mixture model, where each individual (i.e. the entire genotype of an individual) is drawn from a cluster.

# Options

	-a	Choose admixture model (default: no).
 
	-b, --bootstrap	Bootstrap test of H0: K=<k>-1 vs. Ha: K=<k>, where <k> is
		given by -k option.  Specify number of bootstraps as
		argument (default: 0).
  
	--bound	Lower bound for allele and mixing/admixing proportions
		(default: 1.000000e-08).
  
	-B	DEBUG ONLY: turns on opt->test_run
 
	-c	Constrain mixing proportions identical across individuals
		(only enforced with -a; default: no).
  
	-C	The maximum number of iterations to fit (default: 0).
 
	-d	Directory where output files are written (default: ./).
 
	-e	Allowable log likelihood relative error for convergence
		(default: 0.0e+00).
  
	-E	Allowable log likelihood absolute error for convergence
		(default: 1.0e-04).
  
	-f	Name of data file (STRUCTURE format).
 
	--format
		Format of data output file (default: stru).
			stru	STRUCTURE format, the default.
			ped	Plink's ped format.
   
	-g	Adjust step size at most this many times (default: 0)
 
	-h	Display help.
 
	-i	Initial iterations prior to acceleration (default: 0)
 
	-I	Alleles are indices (no sorting, etc.) (default: no)
 
	-I1	Alleles are indices plus 1 (no sorting, etc.) (default: no)
 
	-k	The number of clusters to fit (default: 6).
 
	-1	The minimum number of clusters to fit (default: 6).
 
	-2	The maximum number of clusters to fit (default: 6).
 
	-m	The number of Rand EM initializations, 0 to avoid Rand EM
		(default: 50).
  
	--missing
		Integer value that indicates missing (default: -9).
  
	-M	Parallelization option
		(default: no parallelization). 
  
	-n	Number of initializations to run EM to convergence
		(default: 50).
  
	-o	Option to create unique output file name
 
	-p	The ploidy (default: 2).
 
	--projection
		Turn off simplex projection (default: on).
  
	--plus
		Plus one to alleles when writing data (default: on).
  
	-r	Random number (default: 1234567).
 
	-R	Data file in R format (default: no).
 
	-s	The acceleration scheme (default: No acceleration).
 
		Acceleration options:
  
			0 (default) - no acceleration
   
			1 - SQUAREM version 1
   
			2 - SQUAREM version 2
   
			3 - SQUAREM version 3
   
			4 - Quasi Newton version 1 (1 secant condition)
   
			5 - Quasi Newton version 2 (2 secant conditions)
   
			6 - Quasi Newton version 3 (3 secant conditions)
   
	--simulate <qfile> <pfile> [<ofile>]	Simulate data from admixture <qfile>, <pfile>, and write data to <ofile>.
 
	-u	Iterate until beat target
 
		-u n #: repeat until reach same max # times (default: 0)
  
		-u l #: repeat until reach max log likelihood # (default: 0.000000)
  
	-t	The time (in minutes) to maximize likelihood (default: 0).
		Be sure to check convergence if you set the above!
  
	-v	Level of verbosity (default: 3)
 
	-w	Write best fits to file (default: yes)
 
		Also used to time repeated runs, as follows:
  
		-w n < n >: repeat at least < n > times (default: 1)
  
		-w t < n >: repeat at least < n > minutes (default: 0)
  
		-w m < n >: repeat at most < n > minutes (default: 0)
  
	-x	Use block relaxation algorithm (default: no).

# Citation
 
MULTICLUST - Fast Multinomial Clustering of multiallelic genotypes to infer genetic population structure by Arun Sethuraman, Wei-Chen Chen, Margaret Wanjiku, Karin S. Dorman (forthcoming) Please contact Dr. Sethuraman (asethuraman@sdsu.edu) with questions, or post bug reports to this GitHub page.
