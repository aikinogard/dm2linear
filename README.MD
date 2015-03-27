#dm2linear

## What is dm2linear?

This is a python module expanding electronic density on [gaussian type basis functions](http://en.wikipedia.org/wiki/Gaussian_orbital) from reduced density matrix.

## Is it fast?

No. 
Gaussian integrals are evaluated by pyints.py from Richard P. Muller's [PyQuante](http://pyquante.sourceforge.net). One day I will use [libint](http://sourceforge.net/projects/libint/) to speed up the speed, but not now.

This code is for people who want to do some gaussian integrals in their python code but do not want to spend too much time on finding libint API to python(if there is any).

## Variables

### basis
			read from basis file.
			basis = {'h':[(0,[15.675,3.606,1.208],[0.019,0.063,0.120]),
							(0,[0.473],[0.059]),
							(2,[2.216],[0.003])],
					'o':[.....],
					'n':[.....],
					...}

### bf
		    representation of basis functions in this module.
			list
			bf = [cgto1,cgto2,cgto3,...]

### cgtb
		    contracted gaussian type basis.
			tuple
			(alpha_list,coeff_list,(l,m,n),A)
### gtb 
			single gaussian type basis function.
			tuple
			(alpha,coeff,(l,m,n),A)
