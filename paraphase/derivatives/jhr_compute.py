import numpy as np

def jhr_compute(jh, residuals):
	"""
	Return JHr.

	"""

	return np.dot(jh, residuals)


def jhr_compute_slow(jh, residuals, alpha):
	"""
	Evaluate explicit JHr.
	
	"""

	#
	n_timint, n_freint, n_ant, n_ccor, n_par = alpha.shape
	n_tim = residuals.shape[0]
	n_fre = residuals.shape[1]

	#Initialise JHr.
	jhr = np.zeros(alpha.shape, dtype=jh.dtype)

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
											jhr[tt, ff, ant, cc, par] += jh[tt, ff, ant, cc, par, t, f, p, q, c] * residuals[t, f, p, q, c]

	return jhr


def jhr_compute0(residuals):
	"""
	Evaluates JHr at a point.

	"""

	def get_jhr(jhi):
		return jhi.dot(residuals)

	return get_jhr