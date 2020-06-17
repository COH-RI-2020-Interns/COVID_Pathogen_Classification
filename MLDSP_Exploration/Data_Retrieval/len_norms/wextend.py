import numpy as np
import matplotlib.pyplot as plt



DNA_Seq = [ 1 ,-1 ,-1, 1 ]
len(DNA_Seq)
#Goal = antisymmetric padding equal to a certain length
antisymmetric_DNA = []
def antisymmetric(DNA_sequence, desired_length):
    beg_num = []
    end_num = []
    num_to_add = (int)((desired_length - len(DNA_sequence))/2)
    middle = (int)(len(DNA_sequence)/2 - 1)
    for i in DNA_Seq[middle::-1]:
        beg_num.append(i)
    for j in DNA_Seq[len(DNA_Seq):middle:-1]:
        end_num.append(j)
    for i in range(0,num_to_add):
        antisymmetric_DNA.append(beg_num[i])
    for i in DNA_sequence:
        antisymmetric_DNA.append(i)
    for i in range(0,num_to_add):
        antisymmetric_DNA.append(end_num[i])
    return antisymmetric_DNA



padded = antisymmetric(DNA_Seq, 8)
padded
-1 1 1 -1 -1 1 -1

DNA_Seq = [1,2,3,4]
DNA_Seq[3:0]
