import dm2linear.readfile
import dm2linear.cgtbints
import dm2linear.showdensity
import numpy as np
import matplotlib.pyplot as plt
import os

def test_density(data_path,auxbasis_path,show_td=False):
	if data_path[-1]!='/':
		data_path = data_path + '/'
	if auxbasis_path[-1]!='/':
		auxbasis_path = auxbasis_path + '/'

	print 'read density matrix from dscf.out'
	dm = dm2linear.readfile.read_dscf('alpha RDM in AO basis',data_path+'dscf.out')
	print dm
	print '\n'

	coord = dm2linear.readfile.read_coord(data_path+'coord')
	print '\n'
	lin_basis = dm2linear.readfile.read_basis_lib(auxbasis_path,atoms=['h'],bftype='universal')
	print '\n'
	dm_basis = dm2linear.readfile.read_basis_file(data_path+'basis')
	print '\n'

	print 'creating basis functions for expansion...'
	lin_bf = dm2linear.readfile.create_basis_list(lin_basis,coord)
	print '\n'
	print 'creating basis functions of density matrix...'
	dm_bf = dm2linear.readfile.create_basis_list(dm_basis,coord)
	print 'evaluating overlap matrix:'
	dm_S = dm2linear.cgtbints.two_gaussian_overlap(dm_bf,dm_bf)
	print dm_S
	print 'trace =',np.trace(np.dot(dm_S,dm))
	print 'read overlap matrix from dscf.out'
	dscf_S = dm2linear.readfile.read_dscf('overlap matrix',data_path+'dscf.out')
	print dscf_S
	print dscf_S.shape
	print dm.shape
	print 'trace =',np.trace(np.dot(dscf_S,dm))
	print '\n'

	c = dm2linear.cgtbints.dm2linear(dm,dm_bf,lin_bf)
	print 'expansion coefficients:',c
	print '\n'

	plt.figure()
	for y in [0,1,2,3]:
		n_dm = np.empty(200)
		for i,z in enumerate(np.linspace(-8,8,200)):
			X = (0,y,z)
			n_dm[i] = dm2linear.showdensity.dm2density(dm,dm_bf,X)
		n_lin = np.empty(200)
		for i,z in enumerate(np.linspace(-8,8,200)):
			X = (0,y,z)
			n_lin[i] = dm2linear.showdensity.linear2density(c,lin_bf,X)
		plt.plot(np.linspace(-8,8,200),n_lin,'-',label='linear y='+str(y))
		plt.plot(np.linspace(-8,8,200),n_dm,'--',label='DM y='+str(y))
	ylim = plt.get(plt.gca(), 'ylim')
	for atom_coord in coord:
		atom_z = atom_coord[1][2]
		plt.vlines(atom_z,ylim[0],ylim[1],colors=u'k', linestyles=u'dashed')
	plt.xlabel(r'$z$')
	plt.ylabel(r'$n(0,y,z)$')

	if show_td:
		vec = np.loadtxt(data_path+'td.vec')
		for y in [0,1,2,3]:
			plt.plot(vec[vec[:,1]==y,0],vec[vec[:,1]==y,2],'k.',label='td.vec y='+str(y))

	plt.legend(loc=0)
	plt.show()
	#plt.savefig(os.path.split(data_path[:-1])[1]+'.pdf',bbox_inches='tight')

	return c

if __name__=='__main__':
## test h2+ on an artificial basis with only one uncontracted gtb on each atom
	test_density(data_path='../data/h2+mini/gtb',auxbasis_path='../data/auxbasis/',show_td=True)
## test h2+ on an artificial basis with only one contracted gtb on each atom
	test_density(data_path='../data/h2+mini/cgtb',auxbasis_path='../data/auxbasis/',show_td=True)
## test h2+ on def2-QZVPPD basis, by different atomic distance
	for b in [1.0,3.0,5.0,7.0,11.0]:
		test_density(data_path='../data/h2+/'+str(b),auxbasis_path='../data/auxbasis/',show_td=True)

