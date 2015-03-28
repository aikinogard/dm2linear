from pyints import *
import numpy as np
from scipy.linalg import cho_factor,cho_solve

def contracted_overlap2(cgtb1,cgtb2):
	"""
		Compute the integral of the product of two contracted gaussians
		cgtb1 = (alpha1_list,coeff1_list,(l1,m1,n1),A)
		cgtb2 = (alpha2_list,coeff2_list,(l2,m2,n2),B)

	"""
	(alpha1_list,coeff1_list,(l1,m1,n1),A) = cgtb1
	(alpha2_list,coeff2_list,(l2,m2,n2),B) = cgtb2
	result = 0.
	for mu1,alpha1 in enumerate(alpha1_list):
		for mu2,alpha2 in enumerate(alpha2_list):
			#mu1 and mu2 are all index of gtb in cgtb
			result += coeff1_list[mu1]*coeff2_list[mu2]*overlap(alpha1,(l1,m1,n1),A,alpha2,(l2,m2,n2),B)
	return result

def two_gaussian_overlap(bf1,bf2):
	"""
		return the two gaussian overlap matrix
	"""
	S = np.empty((len(bf1),len(bf2)))

	for i,cgtb1 in enumerate(bf1):
		for j,cgtb2 in enumerate(bf2):
			S[i,j] = contracted_overlap2(cgtb1,cgtb2)

	return S

def gaussian_product(gtb1,gtb2):
	"""
		output a gaussian from the product of two gtb
		Eq.(2.5) in THO66
	"""
	(alpha1,(l1,m1,n1),A) = gtb1
	(alpha2,(l2,m2,n2),B) = gtb2
	P = gaussian_product_center(alpha1,A,alpha2,B)
	gamma = alpha1+alpha2
	exp_pre = exp(-alpha1*alpha2*dist2(A,B)/gamma)
	x_factor = [binomial_prefactor(i,l1,l2,P[0]-A[0],P[0]-B[0]) for i in xrange(l1+l2+1)]
	y_factor = [binomial_prefactor(i,m1,m2,P[1]-A[1],P[1]-B[1]) for i in xrange(m1+m2+1)]
	z_factor = [binomial_prefactor(i,n1,n2,P[2]-A[2],P[2]-B[2]) for i in xrange(n1+n2+1)]
	return P,gamma,exp_pre,x_factor,y_factor,z_factor

def overlap3(gtb1,gtb2,gtb3):
	"""
		Compute the integral of the product of three gaussians
	"""
	# make a new gaussian as a product of gaussian1 and gaussian2
	P,gamma,exp_pre,x_factor,y_factor,z_factor = gaussian_product(gtb1,gtb2)
	lmn12_list = [(l12,m12,n12) for l12 in xrange(len(x_factor)) 
								for m12 in xrange(len(y_factor)) 
								for n12 in xrange(len(z_factor))]
	(alpha3,(l3,m3,n3),C) = gtb3
	result = 0.
	for lmn12 in lmn12_list:
		(l12,m12,n12) = lmn12
		result += exp_pre*x_factor[l12]*y_factor[m12]*z_factor[n12]*\
						overlap(gamma,lmn12,P,alpha3,(l3,m3,n3),C)
	return result

def contracted_overlap3(cgtb1,cgtb2,cgtb3):
	"""
		Compute the integral of the product of three contracted gaussians
	"""
	(alpha1_list,coeff1_list,(l1,m1,n1),A) = cgtb1
	(alpha2_list,coeff2_list,(l2,m2,n2),B) = cgtb2
	(alpha3_list,coeff3_list,(l3,m3,n3),C) = cgtb3
	result = 0.
	for mu1,alpha1 in enumerate(alpha1_list):
		for mu2,alpha2 in enumerate(alpha2_list):
			for mu3,alpha3 in enumerate(alpha3_list):
				#mu1,mu2 and mu3 are all index of gtb in cgtb
				result += coeff1_list[mu1]*coeff2_list[mu2]*coeff3_list[mu3]*\
							overlap3((alpha1,(l1,m1,n1),A),
										(alpha2,(l2,m2,n2),B),
										(alpha3,(l3,m3,n3),C))
	return result

def three_gaussian_overlap(bf1,bf2,bf3):
	"""
		return the three gaussian overlap tensor
	"""
	T = np.empty((len(bf1),len(bf2),len(bf3)))

	for i,cgtb1 in enumerate(bf1):
		for j,cgtb2 in enumerate(bf2):
			for k,cgtb3 in enumerate(bf3):
				T[i,j,k] = contracted_overlap3(cgtb1,cgtb2,cgtb3)
	return T

def dm2linear(dm,dm_bf,lin_bf):
	"""
		dm is the reduced density matrix
		dm_bf is the atomic basis used to represent dm
		lin_bf is the auxiliary basis
	"""
	T = three_gaussian_overlap(lin_bf,dm_bf,dm_bf)
	d = np.tensordot(T,dm,axes=([1,2],[0,1]))
	S = two_gaussian_overlap(lin_bf,lin_bf)
	Q = cho_factor(S)
	return cho_solve(Q, d)












