import numpy as np
import pandas as pd
from os import getcwd
import matplotlib.pyplot as plt
%matplotlib inline
from IPython.display import Image


loc = getcwd() + "/data"
human = pd.read_table(loc + "/human_data.txt")
chimp = pd.read_table(loc + '/chimp_data.txt')
dog = pd.read_table(loc + '/dog_data.txt')
#Chimp, Dog, Human all have several fasta sequences
# We must seperate each sequence into kmers
# The sequence and class the sequence belongs to are listed
# The class is the gene family the sequence belongs to
# 0 = g - protein coupled receptors / 1 = Tyrosine Kinase / 2 = Tyrosine phosphatase/ 3 = Synthethase / 4 = Synthase / 5 = Ion Channel / 6 = Transcription Factor
chimp.head()
dog.head()
human.head()


#dog is a more divergent species than a chimp to a human

# function to convert sequence strings into k-mer words, default size = 6 (hexamer words)
def getKmers(sequence, size=6):
    return [sequence[x:x+size].lower() for x in range(len(sequence) - size + 1)]

#Steps to making kmers:
# 1.Within the dataframe, creating a new column called words and applying it to the sequence
# 2.removing the sequence column with one long
human['words'] = human.apply(lambda x: getKmers(x['sequence']), axis=1)
human = human.drop('sequence', axis=1)
chimp['words'] = chimp.apply(lambda x: getKmers(x['sequence']), axis=1)
chimp = chimp.drop('sequence', axis=1)
dog['words'] = dog.apply(lambda x: getKmers(x['sequence']), axis=1)
dog = dog.drop('sequence', axis=1)




#Making the kmers into a sentence
human_texts = list(human['words'])
for item in range(len(human_texts)):
    human_texts[item] = ' '.join(human_texts[item])
y_h = human.iloc[:, 0].values

chimp_texts = list(chimp['words'])
for item in range(len(chimp_texts)):
    chimp_texts[item] = ' '.join(chimp_texts[item])
y_c = chimp.iloc[:, 0].values                       # y_c for chimp

dog_texts = list(dog['words'])
for item in range(len(dog_texts)):
    dog_texts[item] = ' '.join(dog_texts[item])
y_d = dog.iloc[:, 0].values                         # y_d for dog


human_texts[0]
len(human_texts)
chimp_texts[0]
len(chimp_texts)
dog_texts[0]
len(dog_texts)

# Creating the Bag of Words model using CountVectorizer()
# This is equivalent to k-mer counting
# The n-gram size of 4 was previously determined by testing
#converting into uniform length feature vectors
from sklearn.feature_extraction.text import CountVectorizer
cv = CountVectorizer(ngram_range=(4,4))
X = cv.fit_transform(human_texts)
X_chimp = cv.transform(chimp_texts)
X_dog = cv.transform(dog_texts)

print(X.shape)
print(X_chimp.shape)
print(X_dog.shape)
