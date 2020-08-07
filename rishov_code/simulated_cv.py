import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

X = np.random.randint(0, 100, (1000))

X2 = 2*X

X3 = 3*X

df = pd.DataFrame.from_dict({'col1': X, 'col2': X2, 'col3': X3})

Y = df['col1']

X = df.drop('col1', axis=1)

sample_dict = {}

count = 0
train_df = []
test_df = []
while count < 10:
    sample = df.sample(800)
    test = df.drop(sample.index)
    train_df.append(sample)
    test_df.append(test)
    count += 1

sample_dict = {'train': train_df, 'test': test_df}

# Stratified Sampling
YZ = np.zeros(700)
YO = np.ones(300)

YZ = pd.DataFrame.from_dict({'col': YZ})
YO = pd.DataFrame.from_dict({'col': YO})

Y = pd.concat([YZ, YO], axis=0)

X = np.random.randint(0, 1000, (1000, 3))

f_df = pd.DataFrame(X, columns=['f1', 'f2', 'f3'])

df = pd.concat([Y.reset_index(), f_df], axis=1)

new_count = 0

cv_dict = {}

X_tr = []
X_te = []
Y_tr = []
Y_te = []

X = df.drop('col', axis=1)
Y = df['col']

while new_count < 10:
     X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, stratify=Y)
     X_tr.append(X_train)
     X_te.append(X_test)
     Y_tr.append(Y_train)
     Y_te.append(Y_test)
     new_count += 1

cv_dict['X Train'] = X_tr
cv_dict['X Test'] = X_te
cv_dict['Y Train'] = Y_tr
cv_dict['Y Test'] = Y_te

cv_dict['X Train'][0]
