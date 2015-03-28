import subprocess
import dm2linear.cgtbints
import numpy as np


def group_overlap2_MAE(alpha_group,l_group,m_group,n_group,A_group,random_max_num=None):
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
	tot = (len(combo))**2
	err = 0.				
	for gtb1 in combo:
		for gtb2 in combo:
			c += 1
			result_py,result_math = overlap2_diff(gtb1,gtb2)
			aberror = np.abs(result_py-result_math)
			print '%d/%d Abs Error %e'%(c,tot,aberror)
			err += aberror
	return err/tot

def overlap2_diff(gtb1,gtb2):
	(alpha1,(l1,m1,n1),A) = gtb1
	(alpha2,(l2,m2,n2),B) = gtb2
	result_py = dm2linear.cgtbints.overlap(alpha1,(l1,m1,n1),A,alpha2,(l2,m2,n2),B)
	result_math = float(subprocess.check_output(['/Applications/Mathematica.app/Contents/MacOS/MathematicaScript',
									'-script','gaussian2.m',
								'{%f,{%d,%d,%d},{%f,%f,%f},%f,{%d,%d,%d},{%f,%f,%f}}'
								%(alpha1,l1,m1,n1,A[0],A[1],A[2],
									alpha2,l2,m2,n2,B[0],B[1],B[2])]))
	return result_py,result_math

if __name__ == '__main__':
	print 'Test the accuracy of overlap integral of product of two gaussians\n'
	print 'Test single input:'
	alpha1 = 0.5
	(l1,m1,n1) = (1,0,0)
	A = (1.,2.,2.)
	gtb1 = (alpha1,(l1,m1,n1),A)

	alpha2 = 0.5
	(l2,m2,n2) = (1,0,0)
	B = (1.,2.,2.)
	gtb2 = (alpha2,(l2,m2,n2),B)

	result_py,result_math = overlap2_diff(gtb1,gtb2)
	print 'result:%f, benchmark:%f'%(result_py,result_math)
	print 'result_py-result_math=%e'%(result_py-result_math)
	
	print 'Test a group of input:'
	MAE = group_overlap2_MAE([0.1,1.0],[0,1,2],[0,1,2],[0,1,2],[(0.,0.,1.),(0.,0.,-1.)],20)
	print 'The mean absolute error (MAE) of the overlap function on test set is %e'%(MAE)



	"""
	OUTPUT:

	Test the accuracy of overlap integral of product of two gaussians

	Test single input:
	result:2.784164, benchmark:2.784164
	result_py-result_math=-1.332268e-15
	Test a group of input:
	1/400 Abs Error 5.115908e-13
	2/400 Abs Error 0.000000e+00
	3/400 Abs Error 0.000000e+00
	4/400 Abs Error 0.000000e+00
	5/400 Abs Error 0.000000e+00
	6/400 Abs Error 0.000000e+00
	7/400 Abs Error 0.000000e+00
	8/400 Abs Error 0.000000e+00
	9/400 Abs Error 1.932609e-15
	10/400 Abs Error 0.000000e+00
	11/400 Abs Error 0.000000e+00
	12/400 Abs Error 5.115908e-13
	13/400 Abs Error 0.000000e+00
	14/400 Abs Error 2.775558e-15
	15/400 Abs Error 1.705303e-13
	16/400 Abs Error 0.000000e+00
	17/400 Abs Error 0.000000e+00
	18/400 Abs Error 3.930731e-16
	19/400 Abs Error 0.000000e+00
	20/400 Abs Error 0.000000e+00
	21/400 Abs Error 0.000000e+00
	22/400 Abs Error 1.818989e-11
	23/400 Abs Error 1.932609e-15
	24/400 Abs Error 0.000000e+00
	25/400 Abs Error 0.000000e+00
	26/400 Abs Error 0.000000e+00
	27/400 Abs Error 0.000000e+00
	28/400 Abs Error 0.000000e+00
	29/400 Abs Error 0.000000e+00
	30/400 Abs Error 6.821210e-13
	31/400 Abs Error 1.818989e-12
	32/400 Abs Error 0.000000e+00
	33/400 Abs Error 9.094947e-13
	34/400 Abs Error 0.000000e+00
	35/400 Abs Error 0.000000e+00
	36/400 Abs Error 0.000000e+00
	37/400 Abs Error 0.000000e+00
	38/400 Abs Error 0.000000e+00
	39/400 Abs Error 6.821210e-13
	40/400 Abs Error 0.000000e+00
	41/400 Abs Error 0.000000e+00
	42/400 Abs Error 1.932609e-15
	43/400 Abs Error 2.220446e-16
	44/400 Abs Error 0.000000e+00
	45/400 Abs Error 0.000000e+00
	46/400 Abs Error 0.000000e+00
	47/400 Abs Error 0.000000e+00
	48/400 Abs Error 0.000000e+00
	49/400 Abs Error 0.000000e+00
	50/400 Abs Error 2.220446e-16
	51/400 Abs Error 8.881784e-16
	52/400 Abs Error 0.000000e+00
	53/400 Abs Error 4.440892e-16
	54/400 Abs Error 0.000000e+00
	55/400 Abs Error 0.000000e+00
	56/400 Abs Error 0.000000e+00
	57/400 Abs Error 0.000000e+00
	58/400 Abs Error 0.000000e+00
	59/400 Abs Error 5.551115e-17
	60/400 Abs Error 0.000000e+00
	61/400 Abs Error 0.000000e+00
	62/400 Abs Error 0.000000e+00
	63/400 Abs Error 0.000000e+00
	64/400 Abs Error 2.273737e-13
	65/400 Abs Error 2.220446e-16
	66/400 Abs Error 0.000000e+00
	67/400 Abs Error 0.000000e+00
	68/400 Abs Error 0.000000e+00
	69/400 Abs Error 0.000000e+00
	70/400 Abs Error 0.000000e+00
	71/400 Abs Error 0.000000e+00
	72/400 Abs Error 0.000000e+00
	73/400 Abs Error 0.000000e+00
	74/400 Abs Error 0.000000e+00
	75/400 Abs Error 0.000000e+00
	76/400 Abs Error 0.000000e+00
	77/400 Abs Error 2.220446e-16
	78/400 Abs Error 0.000000e+00
	79/400 Abs Error 0.000000e+00
	80/400 Abs Error 0.000000e+00
	81/400 Abs Error 0.000000e+00
	82/400 Abs Error 0.000000e+00
	83/400 Abs Error 0.000000e+00
	84/400 Abs Error 2.220446e-16
	85/400 Abs Error 5.551115e-17
	86/400 Abs Error 0.000000e+00
	87/400 Abs Error 0.000000e+00
	88/400 Abs Error 0.000000e+00
	89/400 Abs Error 0.000000e+00
	90/400 Abs Error 0.000000e+00
	91/400 Abs Error 0.000000e+00
	92/400 Abs Error 0.000000e+00
	93/400 Abs Error 0.000000e+00
	94/400 Abs Error 0.000000e+00
	95/400 Abs Error 0.000000e+00
	96/400 Abs Error 0.000000e+00
	97/400 Abs Error 5.551115e-17
	98/400 Abs Error 0.000000e+00
	99/400 Abs Error 0.000000e+00
	100/400 Abs Error 0.000000e+00
	101/400 Abs Error 0.000000e+00
	102/400 Abs Error 0.000000e+00
	103/400 Abs Error 0.000000e+00
	104/400 Abs Error 0.000000e+00
	105/400 Abs Error 0.000000e+00
	106/400 Abs Error 1.364242e-12
	107/400 Abs Error 1.705303e-13
	108/400 Abs Error 1.307951e-12
	109/400 Abs Error 0.000000e+00
	110/400 Abs Error 0.000000e+00
	111/400 Abs Error 0.000000e+00
	112/400 Abs Error 0.000000e+00
	113/400 Abs Error 0.000000e+00
	114/400 Abs Error 0.000000e+00
	115/400 Abs Error 0.000000e+00
	116/400 Abs Error 9.094947e-13
	117/400 Abs Error 0.000000e+00
	118/400 Abs Error 0.000000e+00
	119/400 Abs Error 0.000000e+00
	120/400 Abs Error 1.942890e-16
	121/400 Abs Error 0.000000e+00
	122/400 Abs Error 0.000000e+00
	123/400 Abs Error 0.000000e+00
	124/400 Abs Error 0.000000e+00
	125/400 Abs Error 0.000000e+00
	126/400 Abs Error 1.705303e-13
	127/400 Abs Error 5.115908e-13
	128/400 Abs Error 5.684342e-13
	129/400 Abs Error 0.000000e+00
	130/400 Abs Error 0.000000e+00
	131/400 Abs Error 0.000000e+00
	132/400 Abs Error 0.000000e+00
	133/400 Abs Error 0.000000e+00
	134/400 Abs Error 0.000000e+00
	135/400 Abs Error 0.000000e+00
	136/400 Abs Error 4.185442e-13
	137/400 Abs Error 0.000000e+00
	138/400 Abs Error 0.000000e+00
	139/400 Abs Error 0.000000e+00
	140/400 Abs Error 1.332268e-15
	141/400 Abs Error 0.000000e+00
	142/400 Abs Error 0.000000e+00
	143/400 Abs Error 0.000000e+00
	144/400 Abs Error 0.000000e+00
	145/400 Abs Error 0.000000e+00
	146/400 Abs Error 1.307951e-12
	147/400 Abs Error 5.684342e-13
	148/400 Abs Error 1.000444e-11
	149/400 Abs Error 0.000000e+00
	150/400 Abs Error 0.000000e+00
	151/400 Abs Error 0.000000e+00
	152/400 Abs Error 0.000000e+00
	153/400 Abs Error 0.000000e+00
	154/400 Abs Error 0.000000e+00
	155/400 Abs Error 0.000000e+00
	156/400 Abs Error 5.115908e-13
	157/400 Abs Error 0.000000e+00
	158/400 Abs Error 0.000000e+00
	159/400 Abs Error 0.000000e+00
	160/400 Abs Error 6.938894e-17
	161/400 Abs Error 1.932609e-15
	162/400 Abs Error 0.000000e+00
	163/400 Abs Error 0.000000e+00
	164/400 Abs Error 0.000000e+00
	165/400 Abs Error 0.000000e+00
	166/400 Abs Error 0.000000e+00
	167/400 Abs Error 0.000000e+00
	168/400 Abs Error 0.000000e+00
	169/400 Abs Error 2.255141e-16
	170/400 Abs Error 0.000000e+00
	171/400 Abs Error 0.000000e+00
	172/400 Abs Error 1.932609e-15
	173/400 Abs Error 0.000000e+00
	174/400 Abs Error 9.469018e-17
	175/400 Abs Error 5.551115e-17
	176/400 Abs Error 0.000000e+00
	177/400 Abs Error 0.000000e+00
	178/400 Abs Error 5.551115e-17
	179/400 Abs Error 0.000000e+00
	180/400 Abs Error 0.000000e+00
	181/400 Abs Error 0.000000e+00
	182/400 Abs Error 6.821210e-13
	183/400 Abs Error 2.220446e-16
	184/400 Abs Error 0.000000e+00
	185/400 Abs Error 0.000000e+00
	186/400 Abs Error 0.000000e+00
	187/400 Abs Error 0.000000e+00
	188/400 Abs Error 0.000000e+00
	189/400 Abs Error 0.000000e+00
	190/400 Abs Error 1.136868e-12
	191/400 Abs Error 4.092726e-12
	192/400 Abs Error 0.000000e+00
	193/400 Abs Error 5.115908e-13
	194/400 Abs Error 0.000000e+00
	195/400 Abs Error 0.000000e+00
	196/400 Abs Error 0.000000e+00
	197/400 Abs Error 0.000000e+00
	198/400 Abs Error 0.000000e+00
	199/400 Abs Error 4.185442e-13
	200/400 Abs Error 0.000000e+00
	201/400 Abs Error 0.000000e+00
	202/400 Abs Error 1.818989e-12
	203/400 Abs Error 8.881784e-16
	204/400 Abs Error 0.000000e+00
	205/400 Abs Error 0.000000e+00
	206/400 Abs Error 0.000000e+00
	207/400 Abs Error 0.000000e+00
	208/400 Abs Error 0.000000e+00
	209/400 Abs Error 0.000000e+00
	210/400 Abs Error 4.092726e-12
	211/400 Abs Error 1.818989e-11
	212/400 Abs Error 0.000000e+00
	213/400 Abs Error 2.728484e-12
	214/400 Abs Error 0.000000e+00
	215/400 Abs Error 0.000000e+00
	216/400 Abs Error 0.000000e+00
	217/400 Abs Error 0.000000e+00
	218/400 Abs Error 0.000000e+00
	219/400 Abs Error 1.046361e-12
	220/400 Abs Error 0.000000e+00
	221/400 Abs Error 5.115908e-13
	222/400 Abs Error 0.000000e+00
	223/400 Abs Error 0.000000e+00
	224/400 Abs Error 0.000000e+00
	225/400 Abs Error 0.000000e+00
	226/400 Abs Error 0.000000e+00
	227/400 Abs Error 0.000000e+00
	228/400 Abs Error 0.000000e+00
	229/400 Abs Error 1.932609e-15
	230/400 Abs Error 0.000000e+00
	231/400 Abs Error 0.000000e+00
	232/400 Abs Error 5.115908e-13
	233/400 Abs Error 0.000000e+00
	234/400 Abs Error 2.775558e-15
	235/400 Abs Error 1.705303e-13
	236/400 Abs Error 0.000000e+00
	237/400 Abs Error 0.000000e+00
	238/400 Abs Error 3.930731e-16
	239/400 Abs Error 0.000000e+00
	240/400 Abs Error 0.000000e+00
	241/400 Abs Error 0.000000e+00
	242/400 Abs Error 9.094947e-13
	243/400 Abs Error 4.440892e-16
	244/400 Abs Error 0.000000e+00
	245/400 Abs Error 0.000000e+00
	246/400 Abs Error 0.000000e+00
	247/400 Abs Error 0.000000e+00
	248/400 Abs Error 0.000000e+00
	249/400 Abs Error 0.000000e+00
	250/400 Abs Error 5.115908e-13
	251/400 Abs Error 2.728484e-12
	252/400 Abs Error 0.000000e+00
	253/400 Abs Error 0.000000e+00
	254/400 Abs Error 0.000000e+00
	255/400 Abs Error 0.000000e+00
	256/400 Abs Error 0.000000e+00
	257/400 Abs Error 0.000000e+00
	258/400 Abs Error 0.000000e+00
	259/400 Abs Error 1.743934e-13
	260/400 Abs Error 0.000000e+00
	261/400 Abs Error 2.775558e-15
	262/400 Abs Error 0.000000e+00
	263/400 Abs Error 0.000000e+00
	264/400 Abs Error 0.000000e+00
	265/400 Abs Error 0.000000e+00
	266/400 Abs Error 0.000000e+00
	267/400 Abs Error 0.000000e+00
	268/400 Abs Error 0.000000e+00
	269/400 Abs Error 9.469018e-17
	270/400 Abs Error 0.000000e+00
	271/400 Abs Error 0.000000e+00
	272/400 Abs Error 2.775558e-15
	273/400 Abs Error 0.000000e+00
	274/400 Abs Error 2.914335e-16
	275/400 Abs Error 1.387779e-16
	276/400 Abs Error 0.000000e+00
	277/400 Abs Error 0.000000e+00
	278/400 Abs Error 6.539753e-17
	279/400 Abs Error 0.000000e+00
	280/400 Abs Error 0.000000e+00
	281/400 Abs Error 1.705303e-13
	282/400 Abs Error 0.000000e+00
	283/400 Abs Error 0.000000e+00
	284/400 Abs Error 0.000000e+00
	285/400 Abs Error 0.000000e+00
	286/400 Abs Error 0.000000e+00
	287/400 Abs Error 0.000000e+00
	288/400 Abs Error 0.000000e+00
	289/400 Abs Error 5.551115e-17
	290/400 Abs Error 0.000000e+00
	291/400 Abs Error 0.000000e+00
	292/400 Abs Error 1.705303e-13
	293/400 Abs Error 0.000000e+00
	294/400 Abs Error 1.387779e-16
	295/400 Abs Error 1.818989e-12
	296/400 Abs Error 0.000000e+00
	297/400 Abs Error 0.000000e+00
	298/400 Abs Error 1.110223e-16
	299/400 Abs Error 0.000000e+00
	300/400 Abs Error 0.000000e+00
	301/400 Abs Error 0.000000e+00
	302/400 Abs Error 0.000000e+00
	303/400 Abs Error 0.000000e+00
	304/400 Abs Error 0.000000e+00
	305/400 Abs Error 0.000000e+00
	306/400 Abs Error 9.094947e-13
	307/400 Abs Error 4.185442e-13
	308/400 Abs Error 5.115908e-13
	309/400 Abs Error 0.000000e+00
	310/400 Abs Error 0.000000e+00
	311/400 Abs Error 0.000000e+00
	312/400 Abs Error 0.000000e+00
	313/400 Abs Error 0.000000e+00
	314/400 Abs Error 0.000000e+00
	315/400 Abs Error 0.000000e+00
	316/400 Abs Error 4.092726e-12
	317/400 Abs Error 0.000000e+00
	318/400 Abs Error 0.000000e+00
	319/400 Abs Error 0.000000e+00
	320/400 Abs Error 1.932609e-15
	321/400 Abs Error 0.000000e+00
	322/400 Abs Error 0.000000e+00
	323/400 Abs Error 0.000000e+00
	324/400 Abs Error 2.220446e-16
	325/400 Abs Error 5.551115e-17
	326/400 Abs Error 0.000000e+00
	327/400 Abs Error 0.000000e+00
	328/400 Abs Error 0.000000e+00
	329/400 Abs Error 0.000000e+00
	330/400 Abs Error 0.000000e+00
	331/400 Abs Error 0.000000e+00
	332/400 Abs Error 0.000000e+00
	333/400 Abs Error 0.000000e+00
	334/400 Abs Error 0.000000e+00
	335/400 Abs Error 0.000000e+00
	336/400 Abs Error 0.000000e+00
	337/400 Abs Error 5.551115e-17
	338/400 Abs Error 0.000000e+00
	339/400 Abs Error 0.000000e+00
	340/400 Abs Error 0.000000e+00
	341/400 Abs Error 3.930731e-16
	342/400 Abs Error 0.000000e+00
	343/400 Abs Error 0.000000e+00
	344/400 Abs Error 0.000000e+00
	345/400 Abs Error 0.000000e+00
	346/400 Abs Error 0.000000e+00
	347/400 Abs Error 0.000000e+00
	348/400 Abs Error 0.000000e+00
	349/400 Abs Error 5.551115e-17
	350/400 Abs Error 0.000000e+00
	351/400 Abs Error 0.000000e+00
	352/400 Abs Error 3.930731e-16
	353/400 Abs Error 0.000000e+00
	354/400 Abs Error 6.539753e-17
	355/400 Abs Error 1.110223e-16
	356/400 Abs Error 0.000000e+00
	357/400 Abs Error 0.000000e+00
	358/400 Abs Error 4.163336e-17
	359/400 Abs Error 0.000000e+00
	360/400 Abs Error 0.000000e+00
	361/400 Abs Error 0.000000e+00
	362/400 Abs Error 6.821210e-13
	363/400 Abs Error 5.551115e-17
	364/400 Abs Error 0.000000e+00
	365/400 Abs Error 0.000000e+00
	366/400 Abs Error 0.000000e+00
	367/400 Abs Error 0.000000e+00
	368/400 Abs Error 0.000000e+00
	369/400 Abs Error 0.000000e+00
	370/400 Abs Error 4.185442e-13
	371/400 Abs Error 1.046361e-12
	372/400 Abs Error 0.000000e+00
	373/400 Abs Error 1.743934e-13
	374/400 Abs Error 0.000000e+00
	375/400 Abs Error 0.000000e+00
	376/400 Abs Error 0.000000e+00
	377/400 Abs Error 0.000000e+00
	378/400 Abs Error 0.000000e+00
	379/400 Abs Error 2.728484e-12
	380/400 Abs Error 0.000000e+00
	381/400 Abs Error 0.000000e+00
	382/400 Abs Error 0.000000e+00
	383/400 Abs Error 0.000000e+00
	384/400 Abs Error 0.000000e+00
	385/400 Abs Error 0.000000e+00
	386/400 Abs Error 1.942890e-16
	387/400 Abs Error 1.332268e-15
	388/400 Abs Error 6.938894e-17
	389/400 Abs Error 0.000000e+00
	390/400 Abs Error 0.000000e+00
	391/400 Abs Error 0.000000e+00
	392/400 Abs Error 0.000000e+00
	393/400 Abs Error 0.000000e+00
	394/400 Abs Error 0.000000e+00
	395/400 Abs Error 0.000000e+00
	396/400 Abs Error 1.932609e-15
	397/400 Abs Error 0.000000e+00
	398/400 Abs Error 0.000000e+00
	399/400 Abs Error 0.000000e+00
	400/400 Abs Error 5.551115e-17
	The mean absolute error (MAE) of the overlap function on test set is 2.373332e-13
	"""
