import subprocess
import dm2linear.cgtbints
import numpy as np


def group_overlap3_MAE(alpha_group,l_group,m_group,n_group,A_group,random_max_num=None):
	"""
		this function will test all possible combinations of alpha,l,m,n,A
		Those groups do not need to have the same length
	"""
	combo = [(alpha,(l,m,n),A) for alpha in alpha_group 
								for l in l_group 
								for m in m_group
								for n in n_group
								for A in A_group]

	if random_max_num!=None:
		randidx = np.random.choice(len(combo), random_max_num)
		combo = [combo[i] for i in randidx]

	c = 0
	tot = (len(combo))**3
	err = 0.				
	for gtb1 in combo:
		for gtb2 in combo:
			for gtb3 in combo:
				c += 1
				result_py,result_math = overlap3_diff(gtb1,gtb2,gtb3)
				aberror = np.abs(result_py-result_math)
				print '%d/%d Abs Error %e'%(c,tot,aberror)
				err += aberror
	return err/tot

def overlap3_diff(gtb1,gtb2,gtb3):
	(alpha1,(l1,m1,n1),A) = gtb1
	(alpha2,(l2,m2,n2),B) = gtb2
	(alpha3,(l3,m3,n3),C) = gtb3
	result_py = dm2linear.cgtbints.overlap3(gtb1,gtb2,gtb3)
	result_math = float(subprocess.check_output(['/Applications/Mathematica.app/Contents/MacOS/MathematicaScript',
									'-script','gaussian3.m',
								'{%f,{%d,%d,%d},{%f,%f,%f},%f,{%d,%d,%d},{%f,%f,%f},%f,{%d,%d,%d},{%f,%f,%f}}'
								%(alpha1,l1,m1,n1,A[0],A[1],A[2],
									alpha2,l2,m2,n2,B[0],B[1],B[2],
									alpha3,l3,m3,n3,C[0],C[1],C[2])]))
	return result_py,result_math

if __name__=='__main__':
	print 'Test the accuracy of overlap integral of product of three gaussians\n'
	print 'Test single input:'
	alpha1 = 0.1
	(l1,m1,n1) = (2,0,0)
	A = (1.,0.,0.)
	gtb1 = (alpha1,(l1,m1,n1),A)

	alpha2 = 0.3
	(l2,m2,n2) = (2,0,0)
	B = (2.,0.,0.)
	gtb2 = (alpha2,(l2,m2,n2),B)

	alpha3 = 0.3
	(l3,m3,n3) = (2,0,0)
	C = (2.,0.,0.)
	gtb3 = (alpha3,(l3,m3,n3),C)

	print 'integral A B C'
	result_py,result_math = overlap3_diff(gtb1,gtb2,gtb3)
	print 'result:%f, benchmark:%f'%(result_py,result_math)
	print 'result_py-result_math=%e\n'%(result_py-result_math)

	print 'integral B A C'
	result_py,result_math = overlap3_diff(gtb1,gtb2,gtb3)
	print 'result:%f, benchmark:%f'%(result_py,result_math)
	print 'result_py-result_math=%e\n'%(result_py-result_math)

	print 'integral C B A'
	result_py,result_math = overlap3_diff(gtb1,gtb2,gtb3)
	print 'result:%f, benchmark:%f'%(result_py,result_math)
	print 'result_py-result_math=%e\n'%(result_py-result_math)

	print 'Test a group of input:'
	MAE = group_overlap3_MAE([0.1,1.0],[0,1,2],[0,1,2],[0,1,2],[(0.,0.,1.),(0.,0.,-1.)],8)
	print 'The mean absolute error (MAE) of the overlap3 function on test set is %e'%(MAE)


	"""
		OUTPUT:
	Test the accuracy of overlap integral of product of three gaussians

	Test single input:
	integral A B C
	result:46.510372, benchmark:46.510372
	result_py-result_math=-2.351896e-12

	integral B A C
	result:46.510372, benchmark:46.510372
	result_py-result_math=-2.351896e-12

	integral C B A
	result:46.510372, benchmark:46.510372
	result_py-result_math=-2.351896e-12

	Test a group of input:
	1/512 Abs Error 0.000000e+00
	2/512 Abs Error 1.626303e-18
	3/512 Abs Error 6.329572e-16
	4/512 Abs Error 0.000000e+00
	5/512 Abs Error 0.000000e+00
	6/512 Abs Error 3.736594e-15
	7/512 Abs Error 0.000000e+00
	8/512 Abs Error 4.336809e-19
	9/512 Abs Error 1.517883e-18
	10/512 Abs Error 0.000000e+00
	11/512 Abs Error 0.000000e+00
	12/512 Abs Error 0.000000e+00
	13/512 Abs Error 0.000000e+00
	14/512 Abs Error 0.000000e+00
	15/512 Abs Error 0.000000e+00
	16/512 Abs Error 0.000000e+00
	17/512 Abs Error 6.329572e-16
	18/512 Abs Error 0.000000e+00
	19/512 Abs Error 0.000000e+00
	20/512 Abs Error 0.000000e+00
	21/512 Abs Error 0.000000e+00
	22/512 Abs Error 0.000000e+00
	23/512 Abs Error 0.000000e+00
	24/512 Abs Error 0.000000e+00
	25/512 Abs Error 0.000000e+00
	26/512 Abs Error 0.000000e+00
	27/512 Abs Error 0.000000e+00
	28/512 Abs Error 0.000000e+00
	29/512 Abs Error 1.998401e-15
	30/512 Abs Error 0.000000e+00
	31/512 Abs Error 0.000000e+00
	32/512 Abs Error 0.000000e+00
	33/512 Abs Error 0.000000e+00
	34/512 Abs Error 0.000000e+00
	35/512 Abs Error 0.000000e+00
	36/512 Abs Error 1.776357e-15
	37/512 Abs Error 0.000000e+00
	38/512 Abs Error 0.000000e+00
	39/512 Abs Error 3.330669e-16
	40/512 Abs Error 0.000000e+00
	41/512 Abs Error 3.736594e-15
	42/512 Abs Error 0.000000e+00
	43/512 Abs Error 0.000000e+00
	44/512 Abs Error 0.000000e+00
	45/512 Abs Error 0.000000e+00
	46/512 Abs Error 0.000000e+00
	47/512 Abs Error 0.000000e+00
	48/512 Abs Error 0.000000e+00
	49/512 Abs Error 0.000000e+00
	50/512 Abs Error 0.000000e+00
	51/512 Abs Error 0.000000e+00
	52/512 Abs Error 0.000000e+00
	53/512 Abs Error 2.775558e-16
	54/512 Abs Error 0.000000e+00
	55/512 Abs Error 0.000000e+00
	56/512 Abs Error 0.000000e+00
	57/512 Abs Error 4.336809e-19
	58/512 Abs Error 0.000000e+00
	59/512 Abs Error 0.000000e+00
	60/512 Abs Error 0.000000e+00
	61/512 Abs Error 0.000000e+00
	62/512 Abs Error 0.000000e+00
	63/512 Abs Error 0.000000e+00
	64/512 Abs Error 0.000000e+00
	65/512 Abs Error 1.517883e-18
	66/512 Abs Error 0.000000e+00
	67/512 Abs Error 0.000000e+00
	68/512 Abs Error 0.000000e+00
	69/512 Abs Error 0.000000e+00
	70/512 Abs Error 0.000000e+00
	71/512 Abs Error 0.000000e+00
	72/512 Abs Error 0.000000e+00
	73/512 Abs Error 0.000000e+00
	74/512 Abs Error 2.602085e-18
	75/512 Abs Error 6.505213e-19
	76/512 Abs Error 0.000000e+00
	77/512 Abs Error 0.000000e+00
	78/512 Abs Error 1.387779e-16
	79/512 Abs Error 0.000000e+00
	80/512 Abs Error 1.009221e-17
	81/512 Abs Error 0.000000e+00
	82/512 Abs Error 6.505213e-19
	83/512 Abs Error 5.963112e-19
	84/512 Abs Error 0.000000e+00
	85/512 Abs Error 0.000000e+00
	86/512 Abs Error 1.040834e-17
	87/512 Abs Error 0.000000e+00
	88/512 Abs Error 3.252607e-19
	89/512 Abs Error 0.000000e+00
	90/512 Abs Error 0.000000e+00
	91/512 Abs Error 0.000000e+00
	92/512 Abs Error 1.532108e-14
	93/512 Abs Error 0.000000e+00
	94/512 Abs Error 0.000000e+00
	95/512 Abs Error 5.551115e-16
	96/512 Abs Error 0.000000e+00
	97/512 Abs Error 0.000000e+00
	98/512 Abs Error 0.000000e+00
	99/512 Abs Error 0.000000e+00
	100/512 Abs Error 0.000000e+00
	101/512 Abs Error 1.110223e-16
	102/512 Abs Error 0.000000e+00
	103/512 Abs Error 0.000000e+00
	104/512 Abs Error 0.000000e+00
	105/512 Abs Error 0.000000e+00
	106/512 Abs Error 1.110223e-16
	107/512 Abs Error 1.040834e-17
	108/512 Abs Error 0.000000e+00
	109/512 Abs Error 0.000000e+00
	110/512 Abs Error 8.881784e-15
	111/512 Abs Error 0.000000e+00
	112/512 Abs Error 3.469447e-18
	113/512 Abs Error 0.000000e+00
	114/512 Abs Error 0.000000e+00
	115/512 Abs Error 0.000000e+00
	116/512 Abs Error 4.440892e-16
	117/512 Abs Error 0.000000e+00
	118/512 Abs Error 0.000000e+00
	119/512 Abs Error 2.220446e-16
	120/512 Abs Error 0.000000e+00
	121/512 Abs Error 0.000000e+00
	122/512 Abs Error 1.009221e-17
	123/512 Abs Error 4.336809e-19
	124/512 Abs Error 0.000000e+00
	125/512 Abs Error 0.000000e+00
	126/512 Abs Error 6.938894e-18
	127/512 Abs Error 0.000000e+00
	128/512 Abs Error 3.469447e-18
	129/512 Abs Error 6.329572e-16
	130/512 Abs Error 0.000000e+00
	131/512 Abs Error 0.000000e+00
	132/512 Abs Error 0.000000e+00
	133/512 Abs Error 0.000000e+00
	134/512 Abs Error 0.000000e+00
	135/512 Abs Error 0.000000e+00
	136/512 Abs Error 0.000000e+00
	137/512 Abs Error 0.000000e+00
	138/512 Abs Error 6.505213e-19
	139/512 Abs Error 5.963112e-19
	140/512 Abs Error 0.000000e+00
	141/512 Abs Error 0.000000e+00
	142/512 Abs Error 1.040834e-17
	143/512 Abs Error 0.000000e+00
	144/512 Abs Error 3.252607e-19
	145/512 Abs Error 0.000000e+00
	146/512 Abs Error 7.047314e-19
	147/512 Abs Error 2.196052e-16
	148/512 Abs Error 0.000000e+00
	149/512 Abs Error 0.000000e+00
	150/512 Abs Error 2.277692e-15
	151/512 Abs Error 0.000000e+00
	152/512 Abs Error 2.439455e-19
	153/512 Abs Error 0.000000e+00
	154/512 Abs Error 0.000000e+00
	155/512 Abs Error 0.000000e+00
	156/512 Abs Error 1.598721e-14
	157/512 Abs Error 0.000000e+00
	158/512 Abs Error 0.000000e+00
	159/512 Abs Error 3.108624e-15
	160/512 Abs Error 0.000000e+00
	161/512 Abs Error 0.000000e+00
	162/512 Abs Error 0.000000e+00
	163/512 Abs Error 0.000000e+00
	164/512 Abs Error 0.000000e+00
	165/512 Abs Error 7.771561e-16
	166/512 Abs Error 0.000000e+00
	167/512 Abs Error 0.000000e+00
	168/512 Abs Error 0.000000e+00
	169/512 Abs Error 0.000000e+00
	170/512 Abs Error 3.469447e-18
	171/512 Abs Error 2.277692e-15
	172/512 Abs Error 0.000000e+00
	173/512 Abs Error 0.000000e+00
	174/512 Abs Error 1.438849e-13
	175/512 Abs Error 0.000000e+00
	176/512 Abs Error 2.602085e-18
	177/512 Abs Error 0.000000e+00
	178/512 Abs Error 0.000000e+00
	179/512 Abs Error 0.000000e+00
	180/512 Abs Error 3.108624e-15
	181/512 Abs Error 0.000000e+00
	182/512 Abs Error 0.000000e+00
	183/512 Abs Error 1.110223e-16
	184/512 Abs Error 0.000000e+00
	185/512 Abs Error 0.000000e+00
	186/512 Abs Error 3.252607e-19
	187/512 Abs Error 1.897354e-19
	188/512 Abs Error 0.000000e+00
	189/512 Abs Error 0.000000e+00
	190/512 Abs Error 1.734723e-18
	191/512 Abs Error 0.000000e+00
	192/512 Abs Error 1.355253e-19
	193/512 Abs Error 0.000000e+00
	194/512 Abs Error 0.000000e+00
	195/512 Abs Error 0.000000e+00
	196/512 Abs Error 0.000000e+00
	197/512 Abs Error 1.998401e-15
	198/512 Abs Error 0.000000e+00
	199/512 Abs Error 0.000000e+00
	200/512 Abs Error 0.000000e+00
	201/512 Abs Error 0.000000e+00
	202/512 Abs Error 0.000000e+00
	203/512 Abs Error 0.000000e+00
	204/512 Abs Error 1.532108e-14
	205/512 Abs Error 0.000000e+00
	206/512 Abs Error 0.000000e+00
	207/512 Abs Error 5.551115e-16
	208/512 Abs Error 0.000000e+00
	209/512 Abs Error 0.000000e+00
	210/512 Abs Error 0.000000e+00
	211/512 Abs Error 0.000000e+00
	212/512 Abs Error 1.421085e-14
	213/512 Abs Error 0.000000e+00
	214/512 Abs Error 0.000000e+00
	215/512 Abs Error 3.552714e-15
	216/512 Abs Error 0.000000e+00
	217/512 Abs Error 0.000000e+00
	218/512 Abs Error 1.754152e-14
	219/512 Abs Error 3.552714e-15
	220/512 Abs Error 0.000000e+00
	221/512 Abs Error 0.000000e+00
	222/512 Abs Error 1.309672e-10
	223/512 Abs Error 0.000000e+00
	224/512 Abs Error 5.059247e-14
	225/512 Abs Error 4.440892e-16
	226/512 Abs Error 0.000000e+00
	227/512 Abs Error 0.000000e+00
	228/512 Abs Error 0.000000e+00
	229/512 Abs Error 0.000000e+00
	230/512 Abs Error 0.000000e+00
	231/512 Abs Error 0.000000e+00
	232/512 Abs Error 0.000000e+00
	233/512 Abs Error 0.000000e+00
	234/512 Abs Error 0.000000e+00
	235/512 Abs Error 0.000000e+00
	236/512 Abs Error 1.018634e-10
	237/512 Abs Error 0.000000e+00
	238/512 Abs Error 0.000000e+00
	239/512 Abs Error 2.000888e-11
	240/512 Abs Error 0.000000e+00
	241/512 Abs Error 0.000000e+00
	242/512 Abs Error 1.776357e-15
	243/512 Abs Error 7.549517e-15
	244/512 Abs Error 0.000000e+00
	245/512 Abs Error 0.000000e+00
	246/512 Abs Error 2.182787e-11
	247/512 Abs Error 0.000000e+00
	248/512 Abs Error 6.106227e-16
	249/512 Abs Error 0.000000e+00
	250/512 Abs Error 0.000000e+00
	251/512 Abs Error 0.000000e+00
	252/512 Abs Error 5.059247e-14
	253/512 Abs Error 0.000000e+00
	254/512 Abs Error 0.000000e+00
	255/512 Abs Error 8.881784e-16
	256/512 Abs Error 0.000000e+00
	257/512 Abs Error 0.000000e+00
	258/512 Abs Error 0.000000e+00
	259/512 Abs Error 0.000000e+00
	260/512 Abs Error 1.776357e-15
	261/512 Abs Error 0.000000e+00
	262/512 Abs Error 0.000000e+00
	263/512 Abs Error 3.330669e-16
	264/512 Abs Error 0.000000e+00
	265/512 Abs Error 0.000000e+00
	266/512 Abs Error 0.000000e+00
	267/512 Abs Error 0.000000e+00
	268/512 Abs Error 0.000000e+00
	269/512 Abs Error 1.110223e-16
	270/512 Abs Error 0.000000e+00
	271/512 Abs Error 0.000000e+00
	272/512 Abs Error 0.000000e+00
	273/512 Abs Error 0.000000e+00
	274/512 Abs Error 0.000000e+00
	275/512 Abs Error 0.000000e+00
	276/512 Abs Error 0.000000e+00
	277/512 Abs Error 7.771561e-16
	278/512 Abs Error 0.000000e+00
	279/512 Abs Error 0.000000e+00
	280/512 Abs Error 0.000000e+00
	281/512 Abs Error 4.440892e-16
	282/512 Abs Error 0.000000e+00
	283/512 Abs Error 0.000000e+00
	284/512 Abs Error 0.000000e+00
	285/512 Abs Error 0.000000e+00
	286/512 Abs Error 0.000000e+00
	287/512 Abs Error 0.000000e+00
	288/512 Abs Error 0.000000e+00
	289/512 Abs Error 0.000000e+00
	290/512 Abs Error 4.996004e-16
	291/512 Abs Error 5.551115e-16
	292/512 Abs Error 0.000000e+00
	293/512 Abs Error 0.000000e+00
	294/512 Abs Error 4.092726e-12
	295/512 Abs Error 0.000000e+00
	296/512 Abs Error 1.774022e-16
	297/512 Abs Error 0.000000e+00
	298/512 Abs Error 0.000000e+00
	299/512 Abs Error 0.000000e+00
	300/512 Abs Error 0.000000e+00
	301/512 Abs Error 4.092726e-12
	302/512 Abs Error 0.000000e+00
	303/512 Abs Error 0.000000e+00
	304/512 Abs Error 0.000000e+00
	305/512 Abs Error 7.771561e-16
	306/512 Abs Error 0.000000e+00
	307/512 Abs Error 0.000000e+00
	308/512 Abs Error 0.000000e+00
	309/512 Abs Error 0.000000e+00
	310/512 Abs Error 0.000000e+00
	311/512 Abs Error 0.000000e+00
	312/512 Abs Error 0.000000e+00
	313/512 Abs Error 0.000000e+00
	314/512 Abs Error 0.000000e+00
	315/512 Abs Error 0.000000e+00
	316/512 Abs Error 0.000000e+00
	317/512 Abs Error 1.774022e-16
	318/512 Abs Error 0.000000e+00
	319/512 Abs Error 0.000000e+00
	320/512 Abs Error 0.000000e+00
	321/512 Abs Error 3.736594e-15
	322/512 Abs Error 0.000000e+00
	323/512 Abs Error 0.000000e+00
	324/512 Abs Error 0.000000e+00
	325/512 Abs Error 0.000000e+00
	326/512 Abs Error 0.000000e+00
	327/512 Abs Error 0.000000e+00
	328/512 Abs Error 0.000000e+00
	329/512 Abs Error 0.000000e+00
	330/512 Abs Error 1.110223e-16
	331/512 Abs Error 1.040834e-17
	332/512 Abs Error 0.000000e+00
	333/512 Abs Error 0.000000e+00
	334/512 Abs Error 8.881784e-15
	335/512 Abs Error 0.000000e+00
	336/512 Abs Error 3.469447e-18
	337/512 Abs Error 0.000000e+00
	338/512 Abs Error 3.469447e-18
	339/512 Abs Error 2.277692e-15
	340/512 Abs Error 0.000000e+00
	341/512 Abs Error 0.000000e+00
	342/512 Abs Error 1.438849e-13
	343/512 Abs Error 0.000000e+00
	344/512 Abs Error 2.602085e-18
	345/512 Abs Error 0.000000e+00
	346/512 Abs Error 0.000000e+00
	347/512 Abs Error 0.000000e+00
	348/512 Abs Error 1.018634e-10
	349/512 Abs Error 0.000000e+00
	350/512 Abs Error 0.000000e+00
	351/512 Abs Error 2.000888e-11
	352/512 Abs Error 0.000000e+00
	353/512 Abs Error 0.000000e+00
	354/512 Abs Error 0.000000e+00
	355/512 Abs Error 0.000000e+00
	356/512 Abs Error 0.000000e+00
	357/512 Abs Error 4.092726e-12
	358/512 Abs Error 0.000000e+00
	359/512 Abs Error 0.000000e+00
	360/512 Abs Error 0.000000e+00
	361/512 Abs Error 0.000000e+00
	362/512 Abs Error 1.421085e-14
	363/512 Abs Error 1.469935e-13
	364/512 Abs Error 0.000000e+00
	365/512 Abs Error 0.000000e+00
	366/512 Abs Error 2.619345e-10
	367/512 Abs Error 0.000000e+00
	368/512 Abs Error 7.105427e-15
	369/512 Abs Error 0.000000e+00
	370/512 Abs Error 0.000000e+00
	371/512 Abs Error 0.000000e+00
	372/512 Abs Error 2.364686e-11
	373/512 Abs Error 0.000000e+00
	374/512 Abs Error 0.000000e+00
	375/512 Abs Error 9.094947e-13
	376/512 Abs Error 0.000000e+00
	377/512 Abs Error 0.000000e+00
	378/512 Abs Error 6.938894e-18
	379/512 Abs Error 6.938894e-18
	380/512 Abs Error 0.000000e+00
	381/512 Abs Error 0.000000e+00
	382/512 Abs Error 6.217249e-15
	383/512 Abs Error 0.000000e+00
	384/512 Abs Error 2.081668e-17
	385/512 Abs Error 0.000000e+00
	386/512 Abs Error 0.000000e+00
	387/512 Abs Error 0.000000e+00
	388/512 Abs Error 0.000000e+00
	389/512 Abs Error 2.775558e-16
	390/512 Abs Error 0.000000e+00
	391/512 Abs Error 0.000000e+00
	392/512 Abs Error 0.000000e+00
	393/512 Abs Error 0.000000e+00
	394/512 Abs Error 0.000000e+00
	395/512 Abs Error 0.000000e+00
	396/512 Abs Error 4.440892e-16
	397/512 Abs Error 0.000000e+00
	398/512 Abs Error 0.000000e+00
	399/512 Abs Error 2.220446e-16
	400/512 Abs Error 0.000000e+00
	401/512 Abs Error 0.000000e+00
	402/512 Abs Error 0.000000e+00
	403/512 Abs Error 0.000000e+00
	404/512 Abs Error 3.108624e-15
	405/512 Abs Error 0.000000e+00
	406/512 Abs Error 0.000000e+00
	407/512 Abs Error 1.110223e-16
	408/512 Abs Error 0.000000e+00
	409/512 Abs Error 0.000000e+00
	410/512 Abs Error 1.776357e-15
	411/512 Abs Error 7.549517e-15
	412/512 Abs Error 0.000000e+00
	413/512 Abs Error 0.000000e+00
	414/512 Abs Error 2.182787e-11
	415/512 Abs Error 0.000000e+00
	416/512 Abs Error 6.106227e-16
	417/512 Abs Error 7.771561e-16
	418/512 Abs Error 0.000000e+00
	419/512 Abs Error 0.000000e+00
	420/512 Abs Error 0.000000e+00
	421/512 Abs Error 0.000000e+00
	422/512 Abs Error 0.000000e+00
	423/512 Abs Error 0.000000e+00
	424/512 Abs Error 0.000000e+00
	425/512 Abs Error 0.000000e+00
	426/512 Abs Error 0.000000e+00
	427/512 Abs Error 0.000000e+00
	428/512 Abs Error 2.364686e-11
	429/512 Abs Error 0.000000e+00
	430/512 Abs Error 0.000000e+00
	431/512 Abs Error 9.094947e-13
	432/512 Abs Error 0.000000e+00
	433/512 Abs Error 0.000000e+00
	434/512 Abs Error 1.332268e-15
	435/512 Abs Error 1.221245e-15
	436/512 Abs Error 0.000000e+00
	437/512 Abs Error 0.000000e+00
	438/512 Abs Error 9.094947e-13
	439/512 Abs Error 0.000000e+00
	440/512 Abs Error 3.885781e-16
	441/512 Abs Error 0.000000e+00
	442/512 Abs Error 0.000000e+00
	443/512 Abs Error 0.000000e+00
	444/512 Abs Error 8.881784e-16
	445/512 Abs Error 0.000000e+00
	446/512 Abs Error 0.000000e+00
	447/512 Abs Error 5.551115e-17
	448/512 Abs Error 0.000000e+00
	449/512 Abs Error 4.336809e-19
	450/512 Abs Error 0.000000e+00
	451/512 Abs Error 0.000000e+00
	452/512 Abs Error 0.000000e+00
	453/512 Abs Error 0.000000e+00
	454/512 Abs Error 0.000000e+00
	455/512 Abs Error 0.000000e+00
	456/512 Abs Error 0.000000e+00
	457/512 Abs Error 0.000000e+00
	458/512 Abs Error 1.009221e-17
	459/512 Abs Error 4.336809e-19
	460/512 Abs Error 0.000000e+00
	461/512 Abs Error 0.000000e+00
	462/512 Abs Error 6.938894e-18
	463/512 Abs Error 0.000000e+00
	464/512 Abs Error 3.469447e-18
	465/512 Abs Error 0.000000e+00
	466/512 Abs Error 3.252607e-19
	467/512 Abs Error 1.897354e-19
	468/512 Abs Error 0.000000e+00
	469/512 Abs Error 0.000000e+00
	470/512 Abs Error 1.734723e-18
	471/512 Abs Error 0.000000e+00
	472/512 Abs Error 1.355253e-19
	473/512 Abs Error 0.000000e+00
	474/512 Abs Error 0.000000e+00
	475/512 Abs Error 0.000000e+00
	476/512 Abs Error 5.059247e-14
	477/512 Abs Error 0.000000e+00
	478/512 Abs Error 0.000000e+00
	479/512 Abs Error 8.881784e-16
	480/512 Abs Error 0.000000e+00
	481/512 Abs Error 0.000000e+00
	482/512 Abs Error 0.000000e+00
	483/512 Abs Error 0.000000e+00
	484/512 Abs Error 0.000000e+00
	485/512 Abs Error 1.774022e-16
	486/512 Abs Error 0.000000e+00
	487/512 Abs Error 0.000000e+00
	488/512 Abs Error 0.000000e+00
	489/512 Abs Error 0.000000e+00
	490/512 Abs Error 6.938894e-18
	491/512 Abs Error 6.938894e-18
	492/512 Abs Error 0.000000e+00
	493/512 Abs Error 0.000000e+00
	494/512 Abs Error 6.217249e-15
	495/512 Abs Error 0.000000e+00
	496/512 Abs Error 2.081668e-17
	497/512 Abs Error 0.000000e+00
	498/512 Abs Error 0.000000e+00
	499/512 Abs Error 0.000000e+00
	500/512 Abs Error 8.881784e-16
	501/512 Abs Error 0.000000e+00
	502/512 Abs Error 0.000000e+00
	503/512 Abs Error 5.551115e-17
	504/512 Abs Error 0.000000e+00
	505/512 Abs Error 0.000000e+00
	506/512 Abs Error 3.469447e-18
	507/512 Abs Error 1.897354e-19
	508/512 Abs Error 0.000000e+00
	509/512 Abs Error 0.000000e+00
	510/512 Abs Error 1.387779e-17
	511/512 Abs Error 0.000000e+00
	512/512 Abs Error 8.935813e-19
	The mean absolute error (MAE) of the overlap3 function on test set is 1.451953e-12
	"""
