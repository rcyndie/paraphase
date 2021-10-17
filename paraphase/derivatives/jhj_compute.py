import numpy as np

def jhj_compute(jac, jh, alpha):
	"""
	Return approximate to Hessian.

	"""

	#
	n_timint, n_freint, n_ant, n_par, n_ccor = alpha.shape

	##Initialise jhj.
	jhj = np.zeros((n_timint*n_freint*n_ant*n_par*n_ccor, n_timint*n_freint*n_ant*n_par*n_ccor), dtype=jac.dtype)

	n_parccor = n_par*n_ccor

	for k in range(n_timint*n_freint*n_ant):
		jhj[k*n_parccor:(k+1)*n_parccor, k*n_parccor:(k+1)*n_parccor] = np.dot(jh[k*n_parccor:(k+1)*n_parccor, :], 
																			jac[:, k*n_parccor:(k+1)*n_parccor])

	return jhj