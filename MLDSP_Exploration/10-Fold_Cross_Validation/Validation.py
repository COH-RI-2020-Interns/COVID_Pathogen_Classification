from sklearn import cross_validation





data = cross_validation.KFold(len(train_set), n_folds=10, indices=False) 
