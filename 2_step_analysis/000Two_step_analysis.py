#!/usr/bin/python


"""

Suppose we have an approximately equimolar concentration of amplicons, which we are interested
in sequencing to genotype by GT-seq. This simulation attempts to provide an estimate of the
read depth per individual and the total number of reads needed to ensure allelic variation
is captured. 



Execute file with (in Python3):

import os
working_directory = "/Users/joehe/Documents/Python3_projects/Cost_effective_genotyping/PCR_sampling_simulation/2_step_analysis"
os.chdir (working_directory)
exec (open ("000Two_step_analysis.py").read())

"""


import os
import random
import math


##########################################################################################
# Input parameters
##########################################################################################


number_of_individuals = 1000

number_of_chromosomes = 8
read_depth_per_chromosome = 10 #Suggested to be 10 (Catchen et al. 2013)

individual_replications = 1000 # must be >> 1/individual_p_value
individual_p_value = 0.05


##########################################################################################
# Datsets
##########################################################################################


all_individual_data = {}


##########################################################################################
# Read depth per individual
##########################################################################################


"""

Given the protocol of GT-seq, there is no way to distinguish between homeologous chromosomes
in a polploid individual. Thus, we are interested in estimating the minimum read depth per
individual, such that every (potentially) unique chromosome is probably covered at least some
minimum number of times.

Here, we use a monte-carlo style algorithm to estimate read depth per individual:

	1. 

"""


individual_lengths = []

for replicate in range (1, individual_replications + 1):
	all_individual_data [replicate] = {}
	
	success_chromosomes = set ()
	while success_chromosomes < set (range (1, number_of_chromosomes + 1)):
	
		chromosome_sample = random.randrange(1, number_of_chromosomes + 1)
		
		if chromosome_sample not in all_individual_data [replicate]:
			all_individual_data [replicate] [chromosome_sample] = 1
		else:
			all_individual_data [replicate] [chromosome_sample] += 1
		
		if all_individual_data [replicate] [chromosome_sample] == read_depth_per_chromosome:
			success_chromosomes.add (chromosome_sample)
	
	individual_lengths.append (sum (all_individual_data [replicate].values()))


individual_lengths = sorted (individual_lengths)
individual_estimate = individual_lengths [math.floor((1 - individual_p_value) * individual_replications)]

print ("estimate of read depth per individual needed to sample all " + str (number_of_chromosomes) + 
" chromosomes with probability of failure " + str (individual_p_value) + 
" and read depth per chromosome " + str (read_depth_per_chromosome) + " is " 
+ str (individual_estimate))



































