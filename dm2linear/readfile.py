import re
import numpy as np
from cgtbints import contracted_overlap2
from pyints import fact2

def get_lmn(order):
	if order==0:
		return [(0,0,0)]
	elif order==1:
		return [(1,0,0),(0,1,0),(0,0,1)]
	elif order==2:
		return [(2,0,0),(0,2,0),(0,0,2),(1,1,0),(1,0,1),(0,1,1)]
	elif order==3:
		return [(3,0,0),(0,3,0),(0,0,3),(2,1,0),(2,0,1),(1,2,0),(0,2,1),(1,0,2),(0,1,2),(1,1,1)]

def read_coord(filename):
	"""
		read coord file and output in the format:
		coord = [('h',(0.,0.,0.5)),('h',(0.2,0.,0.2)),('o',(0.,0.,0.))]
	"""
	print 'read atomic coordinates from coord file:',filename
	with open(filename,'r') as f:
		fc = f.read()
	m = re.search(r'coord.*?\n(.*?)\$',fc,re.DOTALL)
	lines = m.group(1).strip().split('\n')
	coord = []
	for line in lines:
		line_sp = line.split()
		coord.append((line_sp[3],(float(line_sp[0]),float(line_sp[1]),float(line_sp[2]))))
	return coord


def read_basis_file(filename):
	"""
		read basis from single file
		output in the format:
		basis = {'h':[(0,[15.675,3.606,1.208],[0.019,0.063,0.120]),
						(0,[0.473],[0.059]),
						(2,[2.216],[0.003])],
				'o':[.....],
				'n':[.....],
				...}
	"""
	order_dict = {'s':0,'p':1,'d':2,'f':3}
	print 'input basis from basis file:',filename

	output = {}
	with open(filename,'r') as f:
		fc = f.read()
	raw_block = fc.split('\n*\n')[1:-1]

	for i in xrange(len(raw_block)):
	    if raw_block[i][0]!=' ':
			atom_bftype = raw_block[i].split('\n')[0].split(' ')
			atom = atom_bftype[0]
			bftype = atom_bftype[1]

			print '    input basis for atom %s(%s)'%(atom,bftype)	        
			lines = raw_block[i+1].split('\n')
			atom_basis = []
			for idx,line in enumerate(lines):
				line_sp = line.split()
				if line_sp[1] in 'spdf':
					alpha_list = [float(lines[idx+1+j].split()[0]) for j in xrange(int(line_sp[0]))]
					coeff_list = [float(lines[idx+1+j].split()[1]) for j in xrange(int(line_sp[0]))]
					atom_basis.append((order_dict[line_sp[1]],alpha_list,coeff_list))
			output[atom] = atom_basis
	return output

def read_basis_lib(path,atoms=['h'],bftype='universal'):
	"""
		read basis from library folder
		atoms: ['h','li','ca']
		bftype: type of basis function set, default='universal'
		read basis file and output in the format:
		basis = {'h':[(0,[15.675,3.606,1.208],[0.019,0.063,0.120]),
						(0,[0.473],[0.059]),
						(2,[2.216],[0.003])],
				'o':[.....],
				'n':[.....],
				...}
	"""
	order_dict = {'s':0,'p':1,'d':2,'f':3}
	if path[-1]!='/':
		path+='/'
	print 'input basis from basis library:',path

	output = {}
	for atom in atoms:
		print '    input basis for atom %s(%s)'%(atom,bftype)
		with open(path+atom,'r') as f:
			fc = f.read()
		m = re.search(r'\s'+bftype+'\n.*?\n\s(.*?)\*',fc,re.DOTALL)
		try:
			raw = m.group(1).strip()
		except:
			raise RuntimeError(bftype+' is not available for atom '+atom)
		if 'D' in raw:
			lines = raw.replace('D','E').split('\n')
		else:
			lines = raw.split('\n')
		
		atom_basis = []
		for idx,line in enumerate(lines):
			line_sp = line.split()
			if line_sp[1] in 'spdf':
				alpha_list = [float(lines[idx+1+i].split()[0]) for i in xrange(int(line_sp[0]))]
				coeff_list = [float(lines[idx+1+i].split()[1]) for i in xrange(int(line_sp[0]))]
				atom_basis.append((order_dict[line_sp[1]],alpha_list,coeff_list))
		output[atom] = atom_basis
	return output

def read_dscf(head,filename):
	"""
		read the matrix from turbomole output
	"""

	print 'input matrix (head:%s) from: %s'%(head,filename)
	with open(filename,'r') as f:
		fc = f.read()
	m = re.search(head + r'.*?\n(.*?)END', fc, re.DOTALL)
	lines = m.group(1).strip().split('\n')
	dim = len(lines)
	mat = np.zeros((dim,dim))
	lin_entries = np.sum([line.split() for line in lines])
	lin_mat = np.array([float(entry) for entry in lin_entries])
	mat[np.tril_indices(dim)] = lin_mat
	mat = mat+mat.T-np.diag(mat.diagonal())
	return mat

def get_norm_constant_gtb(gtb):
	"""
		return the normalization constant
		THO66 Eq.(2.2)
	"""
	(alpha,(l,m,n),A) = gtb
	up = pow(2,2*(l+m+n)+1.5)*pow(alpha,l+m+n+1.5)
	down = fact2(2*l-1)*fact2(2*m-1)*fact2(2*n-1)*pow(np.pi,1.5)
	return np.sqrt(up/down)

def norm_cgtb(cgtb):
	"""
		normalize each gtb and combine with coeff and then norm cgtb again
	"""
	(alpha_list,coeff_list,lmn,A) = cgtb
	# normalize each gtb
	coeff_list_N1 = []
	for mu1,alpha in enumerate(alpha_list):
		N1 = get_norm_constant_gtb((alpha,lmn,A))
		coeff_list_N1.append(coeff_list[mu1]*N1)
	# combine
	cgtb_N1 = (alpha_list,coeff_list_N1,lmn,A)
	# normalize cgtb
	N2 = 1./np.sqrt(contracted_overlap2(cgtb_N1,cgtb_N1))
	coeff_list_N2 = [c*N2 for c in coeff_list_N1]	
	return (alpha_list,coeff_list_N2,lmn,A)

def create_basis_list(basis,coord):
	"""
		create basis from atomic basis and atomic coordinates
		basis = {'h':[(0,[15.675,3.606,1.208],[0.019,0.063,0.120]),
						(0,[0.473],[0.059]),
						(2,[2.216],[0.003])],
				'o':[.....],
				'n':[.....],
				...}
		coord = [('h',(0.,0.,0.5)),('h',(0.2,0.,0.2)),('o',(0.,0.,0.))]

		bf = [cgtb1,cgtb2,cgtb3,cgtb4...] 
			= [(alpha1_list,coeff1_list,(l1_order1,m1_order1,n_order1),A1),
					(alpha1_list,coeff1_list,(l1_order1,m1_order1,n_order1),A1),
					(alpha2_list,coeff2_list,(l1_order2,m1_order2,n_order2),A2),...]

		the basis functions will be sorted in turbomole order
		for example, if we have an LiH, each of the atom has 1s, 2s, p
		then the order is 1s(Li) 2s(Li) 1s(H) 2s(H) p(Li) p(Li)
		order > atom > basis

	"""
	
	bf = []
	for order in [0,1,2,3]:
		for atom in coord:
			for single_basis in basis[atom[0]]:
				if single_basis[0] == order:
					for lmn in get_lmn(order):
						cgtb = (single_basis[1],single_basis[2],lmn,atom[1])
						bf.append(norm_cgtb(cgtb))
	print 'Basis set includes %d cgtb functions.'%len(bf)
	return bf

if __name__ == '__main__':
	coord = read_coord('data/coord')
	basis = read_basis_lib('data/auxbasis/',atoms=['h','o'],bftype='def-TZVPP')
	for cgtb in create_basis_list(basis,coord):
		print cgtb
