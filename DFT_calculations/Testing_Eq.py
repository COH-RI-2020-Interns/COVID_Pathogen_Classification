#Practice with MLDSP methods

#Making the Si dna_strand
dna_strand = "CGAT"


#making a function that changes from Si to Ni
#Integer representation = T:0, C:1, A:2, G:3
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
