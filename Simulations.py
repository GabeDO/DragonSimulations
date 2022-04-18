# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 13:59:20 2022

@author: Gabe
"""

import csv
import time
import random
from simuPOP import *
import simuPOP as sim
from simuPOP.utils import importPopulation, export
from simuPOP.sampling import drawRandomSample


LocNumber = 4691
start = time.time()


#with open('SimulationLoci.csv', 'a', newline='') as myfile:
#     wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
#     wr.writerow(["Heterozygosity", "Loci", "Replicate"])


#now lets start the simulation

#making a output file
with open('Output.csv', 'a', newline='') as myfile:
     wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
     header = []
     for locus in range(LocNumber):
         header.append("locus "+str(locus))
     wr.writerow(header)

#function to record passed time    
def timecount(param):    
    end = time.time()   
    time_elapsed = end - param
    #print(time_elapsed)
    return True


# a little function to export the expected het values for each loci of a population. 
def exportFunc(pop):    
    #stat(pop, alleleFreq = ALL_AVAIL) 
    #export(pop, format='FASTA', output='test'+ str(pop.dvars().gen)+'.txt')
    stat(pop, alleleFreq=ALL_AVAIL)
    H_exp = []
    for i in pop.dvars().alleleFreq:    
        H_exp.append(1 - sum([x*x for x in pop.dvars().alleleFreq[i].values()]))
               
    return H_exp




for rep in range(50):
    #importing the population from a genepop file
    pop = importPopulation('GENEPOP', 'Gen1Genepop.txt')
    timecount(start)
    
    #telling the script what variables I'll need
    stat(pop, alleleFreq=ALL_AVAIL, subPops= ALL_AVAIL)
    
    #exporting gen 0 
    HetGen0 = exportFunc(pop)
    
    
    ##This bit is a lil wack. So I'm looking through every loci of every individual, 
    #and setting all missing values (0) to be whatever the most common value at that 
    #locus is. So:
    #for each dragon
    for dragon in range(129):
        #at each locus
        for locus in range(LocNumber):
            #if the variant is equal to 0
            if pop.individual(dragon).allele(locus) == 0:
                
                #set a random value
                RandomNum = random.uniform(0, 1)
                #now i find out what the freq of the most common allele
                MaxProp = max(pop.dvars().alleleFreq[locus].values())

                #organise the variants and their props into a sorted list
                y = sorted(list(pop.dvars().alleleFreq[locus].items()), key=lambda tup: tup[1])
                              
                
                #Do some RNG to make fill in missing data with a variant, with chance of variant equal to its proportion in the pop
                # If the random value I made is lower than the most common variant, I set x to the most common variant
                if RandomNum <= MaxProp:
                    x = y[-1][0]
                    
                else:
                    #if its not, I fist check to make sure the 2nd most common variant isn't a 0 (eg a missing value)
                    if y[-2][0] == 0:
                        #if it is, I just set x to the most common variant again
                        x = y[-1][0]
                    else:
                        #if its not 0, then it becomes the 2nd most common variant
                        x = y[-2][0]
                        
                
                    
                #then I set the missing data to be x
                pop.individual(dragon).setAllele(x,locus)
            else:
                yeet = 1
    
    #start the evolution
    pop.evolve(
            #setting it to random mating
            matingScheme = RandomMating(),   
            #50/50 sex split at the moment - I can change that tho
            initOps = [
                        InitSex()],
    
            #telling the sim to go for 2 gens (to gen 3)
            gen=2
            )
    #telling the script what variables I'll need again now that the pop has evolved
    stat(pop, alleleFreq=ALL_AVAIL, subPops= ALL_AVAIL)
    
    #export final data
    HetGen3 = exportFunc(pop)
    
    DeltaHet = [x - y for x, y in zip(HetGen0, HetGen3)]
    
    
    #write output to the output file. 
    with open('Output.csv', 'a', newline='') as myfile:
         wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
         dataoutput = []
         for locus in range(LocNumber):
             dataoutput.append(DeltaHet[locus])
         wr.writerow(dataoutput)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
#Old code, dont worry about it
    
    
    
    
    #now this is where I calculate He for each of the 4691 loci
    #for loci in range(4691):
    #    #get allele props for each possible variant
    #    #Zero is missing data, which should always be 0 after I removed all the 0s
    #    Zero = pop.dvars().alleleFreq[loci][0]
    #    A = pop.dvars().alleleFreq[loci][1]
    #    T = pop.dvars().alleleFreq[loci][2]
    #    C = pop.dvars().alleleFreq[loci][3]
    #    G = pop.dvars().alleleFreq[loci][4]
    #    
    #    #calculate proportions of each variant - probs not needed now that 0s are gone
    #    Ap = A/sum([A,T,C,G])
    #    Tp = T/sum([A,T,C,G])
    #    Cp = C/sum([A,T,C,G])
    #   Gp = G/sum([A,T,C,G])
        
    #    #the square of each proportion, to help calc heterozygosity
    #    Ah=Ap**2
    #    Th=Tp**2
    #    Ch=Cp**2
    #    Gh=Gp**2
        
    #    #calculate heterozygosity
    #    He = 1-sum([Ah,Th,Ch,Gh])

    #    #write het to my data file, along with what loci its for, and what replicate
    #    with open('SimulationLoci.csv', 'a', newline='') as myfile:
    #        wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    #        wr.writerow([He, loci, rep])
    #        
    #    #loop it again for the next loci. 