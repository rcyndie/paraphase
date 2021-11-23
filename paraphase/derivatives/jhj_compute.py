import numpy as np

def jhj_compute(jac, jh, alpha):
	"""
	Return approximate to Hessian.

	"""

	#
	n_timint, n_freint, n_ant, n_ccor, n_par = alpha.shape
	n_parccor = n_par*n_ccor

	##Initialise jhj.
	jhj = np.zeros((n_timint*n_freint*n_ant*n_parccor, n_timint*n_freint*n_ant*n_parccor), dtype=jac.dtype)

	for k in range(n_timint*n_freint*n_ant):
		jhj[k*n_parccor:(k+1)*n_parccor, k*n_parccor:(k+1)*n_parccor] = np.dot(jh[k*n_parccor:(k+1)*n_parccor, :], 
																			jac[:, k*n_parccor:(k+1)*n_parccor])

	return jhj


def jhj_compute_slow(jac, jh, alpha):
	"""
	Evaluates JHJ explicitly.
	Remember this is not a block diagonal approximation. 

	"""

	#
	n_timint, n_freint, n_ant, n_ccor, n_par = alpha.shape
	n_tim = jac.shape[0]
	n_fre = jac.shape[1]

	#Initialise JHJ.
	jhj = np.zeros((alpha.shape+alpha.shape), dtype=jac.dtype)

	for tt in range(n_timint):
		for ff in range(n_freint):
			for ant in range(n_ant):
				for cc in range(n_ccor):
					for par in range(n_par):
						for t in range(n_tim):
							for f in range(n_fre):
								for p in range(n_ant):
									for q in range(n_ant):
										for c in range(n_ccor):
											jhj[tt, ff, ant, cc, par, tt, ff, ant, cc, par] += jh[tt, ff, ant, cc, par, t, f, p, q, c] * jac[t, f, p, q, c, tt, ff, ant, cc, par]

	return jhj

