import numpy as np
import pandas as pd
from os import getcwd
import matplotlib.pyplot as plt
%matplotlib inline
from IPython.display import Image
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score
from sklearn.naive_bayes import MultinomialNB
from sklearn.feature_extraction.text import CountVectorizer




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
y_h
chimp_texts = list(chimp['words'])
for item in range(len(chimp_texts)):
    chimp_texts[item] = ' '.join(chimp_texts[item])
y_c = chimp.iloc[:, 0].values                       # y_c for chimp

dog_texts = list(dog['words'])
for item in range(len(dog_texts)):
    dog_texts[item] = ' '.join(dog_texts[item])
y_d = dog.iloc[:, 0].values        # y_d for dog

human_texts
len(human_texts)
chimp_texts
len(chimp_texts)
dog_texts[0]
len(dog_texts)

# Creating the Bag of Words model using CountVectorizer()
# This is equivalent to k-mer counting
# The n-gram size of 4 was previously determined by testing
#converting into uniform length feature vectors

cv = CountVectorizer(ngram_range=(4,4))
X = cv.fit_transform(human_texts)
X_chimp = cv.transform(chimp_texts)
X_dog = cv.transform(dog_texts)




print(X.shape)
print(X_chimp.shape)
print(X_dog.shape)



#Graphing the class distributions in human, chimp and dog
human['class'].value_counts().sort_index().plot.bar()
chimp['class'].value_counts().sort_index().plot.bar()
dog['class'].value_counts().sort_index().plot.bar()


# Splitting the human dataset into the training set and test set
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X,y_h, test_size = 0.20, random_state=42)

print(X_train.shape)
print(X_test.shape)


### Multinomial Naive Bayes Classifier , save as Naive Bayes###
# The alpha parameter was determined by grid search previously
classifier = MultinomialNB(alpha=0.1)
classifier.fit(X_train, y_train)
y_pred = classifier.predict(X_test)



print("Confusion matrix\n")
print(pd.crosstab(pd.Series(y_test, name='Actual'), pd.Series(y_pred, name='Predicted')))
def get_metrics(y_test, y_predicted):
    accuracy = accuracy_score(y_test, y_predicted)
    precision = precision_score(y_test, y_predicted, average='weighted')
    recall = recall_score(y_test, y_predicted, average='weighted')
    f1 = f1_score(y_test, y_predicted, average='weighted')
    return accuracy, precision, recall, f1
accuracy, precision, recall, f1 = get_metrics(y_test, y_pred)
print("accuracy = %.3f \nprecision = %.3f \nrecall = %.3f \nf1 = %.3f" % (accuracy, precision, recall, f1))



# Predicting the chimp, dog and worm sequences
y_pred_chimp = classifier.predict(X_chimp)
y_pred_dog = classifier.predict(X_dog)

# performance on chimp genes
print("Confusion matrix\n")
print(pd.crosstab(pd.Series(y_c, name='Actual'), pd.Series(y_pred_chimp, name='Predicted')))
accuracy, precision, recall, f1 = get_metrics(y_c, y_pred_chimp)
print("accuracy = %.3f \nprecision = %.3f \nrecall = %.3f \nf1 = %.3f" % (accuracy, precision, recall, f1))

# performance on dog genes
print("Confusion matrix\n")
print(pd.crosstab(pd.Series(y_d, name='Actual'), pd.Series(y_pred_dog, name='Predicted')))
accuracy, precision, recall, f1 = get_metrics(y_d, y_pred_dog)
print("accuracy = %.3f \nprecision = %.3f \nrecall = %.3f \nf1 = %.3f" % (accuracy, precision, recall, f1))

#The model seems to perform well on human data. It also does on Chimpanzee. That might not be a surprise since the chimp and human are so similar genetically. The performance on dog is not quite as good. We would expect this since the dog is more divergent from human than the chimpanze.
