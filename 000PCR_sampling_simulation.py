#!/usr/bin/python

"""

Script to simulate the random sampling of primers in sequencing to identify number of
reads necessary to be confident that all relavant amplicons are sampled.


Execute file with (in Python3):

import os
working_directory = "/Users/joehe/Documents/Python3_projects/Cost_effective_genotyping/PCR_sampling_simulation"
os.chdir (working_directory)
exec (open ("000PCR_sampling_simulation.py").read())

"""

import os
import random
import math

##########################################################################################
# Input parameters
##########################################################################################


number_of_amplicons = 1600
specific_amplicons_needed = 500 # <= number_of_amplicons
probability_of_failure = 0.05 # 1/probability_of_failure << simulation_repetition

simulation_repetition = 1000


##########################################################################################
# Datasets
##########################################################################################


all_sampled_data = {}
ordered_length_list = []


##########################################################################################
# Simulation Sequence
##########################################################################################


for simulation_count, simulation in enumerate (range (1, simulation_repetition + 1)):
	
	all_sampled_data [simulation] = []
	print ("calculating simulation " + str (simulation_count + 1) + " of " + str (simulation_repetition))
	
	amplicon_count = 0
	while not set (range (1, specific_amplicons_needed +1)) < set (all_sampled_data [simulation]):
		amplicon_count += 1
		random_amplicon = random.randrange (1, number_of_amplicons +1)
		if random_amplicon not in all_sampled_data [simulation]:
			all_sampled_data [simulation].append (random_amplicon)
			
	ordered_length_list.append (amplicon_count)


##########################################################################################
# read depth estimation
##########################################################################################


ordered_length_list = sorted (ordered_length_list)
read_depth_estimate = ordered_length_list [simulation_repetition - math.floor (probability_of_failure * simulation_repetition)]

print ("Estimate of read depth to identify at least " + str (specific_amplicons_needed) + 
" specific amplicons from " + str (number_of_amplicons) + " with probability of failure <= "
 + str (probability_of_failure) + " is " + str (read_depth_estimate))

