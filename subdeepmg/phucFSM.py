import numpy
import pandas
def varianceFilter(X_train,X_test,desiredFeature):
	result_X_train=numpy.zeros((X_train.shape[0], 0))
	result_X_test =numpy.zeros((X_test.shape[0], 0))
	# fill index
	vars=numpy.var(X_train,axis=0)
	# append index
	VARS=[]
	for i in range(len(vars)):
		VARS.append((vars[i],i))
	VARS.sort()
	for _, index in VARS[len(VARS):len(VARS)-desiredFeature-1:-1]:
		result_X_train=numpy.append(result_X_train, X_train[:,index].reshape(X_train.shape[0], 1), axis=1)
		result_X_test =numpy.append(result_X_test , X_test [:,index].reshape(X_test.shape[0], 1), axis=1)
	return result_X_train, result_X_test
def anovaFValue(X_train, X_test, desiredFeature, Y_train):
	from sklearn.feature_selection import SelectKBest
	from sklearn.feature_selection import f_classif
	fvalue_selector = SelectKBest(f_classif, k=desiredFeature)
	X_train=fvalue_selector.fit_transform(X_train, Y_train)
	X_test=fvalue_selector.transform(X_test)
	return X_train, X_test
def perceptronCheat(X_train,X_test,sqrSize,Y_train):
	import sklearn
	from sklearn.linear_model import Perceptron
	from sklearn.neural_network import MLPClassifier
	clf = Perceptron(tol=0, max_iter=2000)
	clf.fit (X_train, Y_train)
	coef = clf.coef_[0]
	tmp = []
	i = 0
	for a in coef: 
		tmp.append ((a, i))
		i+=1
	tmp.sort ()
	totalSize = sqrSize
	head = totalSize // 2
	tail = totalSize - head
	selected = [index for (value, index) in tmp[:head]]
	selected.extend ([index for (value, index) in tmp[len(tmp)-1:len(tmp)-tail-1:-1]])

	result_X_train=numpy.zeros((X_train.shape[0], 0))
	result_X_test =numpy.zeros((X_test.shape[0], 0))
	result_X_train=numpy.append(result_X_train, X_train[:,selected].reshape(X_train.shape[0], len(selected)), axis=1)
	result_X_test =numpy.append(result_X_test , X_test [:,selected].reshape(X_test.shape[0],  len(selected)), axis=1)
	return result_X_train, result_X_test	
if __name__ == "__main__":
	pass