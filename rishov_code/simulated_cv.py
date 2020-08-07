import numpy as np
import pandas as pd

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
