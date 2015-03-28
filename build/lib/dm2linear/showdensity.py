import numpy as np

def eval_cgtb(cgtb,X):
	"""
		evaluate the value of contracted gaussian type basis function 
		on point X
	"""
	(alpha_list,coeff_list,(l,m,n),A) = cgtb
	output = 0.
	AX2 = pow(A[0]-X[0],2)+pow(A[1]-X[1],2)+pow(A[2]-X[2],2)
	for i,alpha in enumerate(alpha_list):
		output += coeff_list[i]*pow(A[0]-X[0],l)*pow(A[1]-X[1],m)*pow(A[2]-X[2],n)*np.exp(-alpha*AX2)
	return output

def eval_bf(bf,X):
	"""
		evaluate the value of each cgtb in bf on point X
		and output as a vector
	"""
	output = np.empty(len(bf))
	for i,cgtb in enumerate(bf):
		output[i] = eval_cgtb(cgtb,X)
	return output

def dm2density(dm,bf,X):
	"""
		density represented by density matrix at point X
	"""
	bfX = eval_bf(bf,X)
	return np.dot(bfX,np.dot(dm,bfX))

def linear2density(c,bf,X):
	"""
		density represented by linear combination 
		of gaussian type basis at point X
	"""
	bfX = eval_bf(bf,X)
	return np.dot(c,bfX)





