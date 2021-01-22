import numpy as np
from scipy.fft import fft, ifft
from scipy import stats

from os import getcwd, listdir

# Loading in the Data
file_path = getcwd() + "/data/Test1/Riboviria"
listdir(file_path)
#Practice with MLDSP methods


#making a function that changes from Si to Ni
#Integer representation = T:0, C:1, A:2, G:3
dict_of_bases = {"T":0, "C":1, "A":2, "G":3}
numeric = []
def numerical(dna_strand):
    for base in dna_strand:
        numeric.append(dict_of_bases[base])
    return np.array(numeric)

dna_strand_1 = np.array(numerical(file_path[0]))
dna_strand_2 = np.array(numerical(dna_strand2))
type(dna_strand_1)
type(dna_strand_2)
#run strand 1 and then rerun for strand 2

#Calculating discrete numerical representation
𝐹𝑖(𝑘)=∑𝑗=0𝑝−1𝑓(𝑆𝑖(𝑗))⋅𝑒(−2𝜋𝑖/𝑝)𝑘𝑗
#Results  = [6, -1 -3i , 0, 1+3i]
#Logic =
#to get 6 do (1) + (2) + (3)
#to get -1 -3i do 1 + 3*e^(-pi(i)/2) + 2*e^(-pi(i)) + 0
#to get 0 do 1 + 3*e^(-pi(i)) + 2*e^(-2pi(i)) + 0
#to get -1+3i 1 + 3*e^(-3/2(pi(i))) + 2*e^(-3(pi)(i)) + 0

dft_strand_1 = fft(dna_strand_1)
dft_strand_2 = fft(dna_strand_2)


#Finding magnitude
mag_list = []
def mag(list_of_dtf):
    mag_list = abs(list_of_dtf)
    return mag_list


mag_strand_1 = mag(dft_strand_1)
len(mag_strand_1)
mag_strand_2 = mag(dft_strand_2)
len(mag_strand_2)



#Finding Pearson's Correlation Coefficiet
stats.pearsonr(mag_strand_1, mag_strand_2)

#first result = pearson's correlation coefficient so close to 1, higher similarity
#second result = p-value (0 to 1), the lower the p-value the more statistically significant the correlation is
#the lower the p-value the more accurate the results
#comparing identical strands gives a p-values of 0 and a PCC of 1
