import math


#Practice with MLDSP methods

#Making the Si dna_strand
dna_strand = "CGAT"




#Integer representation using if statements
list_numeric = []
def numeric(strand):
    for i in strand:
         if (i == "T"):
             list_numeric.append(0)
         elif (i == "C"):
             list_numeric.append(1)
         elif(i == "A"):
             list_numeric.append(2)
         elif (i == "G"):
             list_numeric.append(3)

numeric(dna_strand)
list_numeric


#Calculating discrete numerical representation
𝐹𝑖(𝑘)=∑𝑗=0𝑝−1𝑓(𝑆𝑖(𝑗))⋅𝑒(−2𝜋𝑖/𝑝)𝑘𝑗
#Results  = [6, -1 -3i , 0, 1+3i]
#Logic =
#to get 6 do (1) + (2) + (3)
#to get -1 -3i do 1 + 3*e^(-pi(i)/2) + 2*e^(-pi(i)) + 0
#to get 0 do 1 + 3*e^(-pi(i)) + 2*e^(-2pi(i)) + 0
#to get -1+3i 1 + 3*e^(-3/2(pi(i))) + 2*e^(-3(pi)(i)) + 0

list_dtf = []
length = len(list_numeric)
length

def dtf(list_of_numbers, length):
    k = 0
    j = 0
    sum = 0
    pi = np.pi
    e = np.e
    i = 1j
    exponent = (-2 * pi * i) / length
    for m in range(k,length):
        for n in list_of_numbers:
            sum = sum + (n * np.exp(exponent * m * j))
            j = j+1
        list_dtf.append(sum)
        sum = 0
        j = 0



dtf(list_numeric, length)

list_dtf


mag_list = []
def mag(list_of_dtf):
    for i in list_of_dtf:
        i = abs(i)
        mag_list.append(i)

list_test = [-1 -4j, 2-3j, -4, 7]
mag(list_test)
print(mag_list)
