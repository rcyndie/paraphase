import numpy as np

def jh_compute(jac):
	"""
	Returns complex conjugate transpose of Jacobian.

	"""

	return np.conjugate(jac.T)


def jh_compute_slow(data_arr, model_arr, alpha, lcoord, mcoord, chan_freq, freqf, basis):
	"""
	Evaluates an explicit computation of JH.

	"""

	n_timint, n_freint, n_ant, n_ccor, n_par = alpha.shape
	n_dir = model_arr.shape[0]
	n_tim = data_arr.shape[0]
	n_fre = data_arr.shape[1]

	#Initialise complex conjugate transpose of Jacobian.
	jh = np.zeros(alpha.shape+data_arr.shape, dtype=data_arr.dtype)

	for tt in range(n_timint):
		for ff in range(n_freint):
			for ant in range(n_ant):
				for cc in range(n_ccor):
					for par in range(n_par):
						#end of param axes
						#start of data axes
						for t in range(n_tim):
							for f in range(n_fre):
								fprof = freqf(chan_freq[f])
								for p in range(n_ant):
									for q in range(n_ant):
										for c in range(n_ccor):
											# end of data axes
											model = model_arr[:, t, f, p, q, c]
											for d in range(n_dir):
												#Get partial derivative of the phase.
												if ant==p and p!=q and cc==c and tt==t//n_timint and ff==f//n_freint:
													L = basis(lcoord[d], mcoord[d])
													Lalpha = L.dot(alpha[tt, ff, ant, cc])
													prefact = -1j * L[par] * fprof
													phase = fprof * Lalpha
													gp = np.exp(1j * phase)
													gq = np.exp(1.0j*fprof*L.dot(alpha[tt, ff, q, cc]))
													jh[tt, ff, ant, cc, par, t, f, p, q, c] += prefact * gp * model[d] * np.conj(gq)
												elif ant==q and q!=p and cc==c and tt==t//n_timint and ff==f//n_freint:
													L = basis(lcoord[d], mcoord[d])
													Lalpha = L.dot(alpha[tt, ff, ant, cc])
													prefact = 1j * L[par] * fprof
													phase = fprof * Lalpha
													gq = np.exp(1j * phase)
													gp = np.exp(1.0j*fprof*L.dot(alpha[tt, ff, p, cc]))
													jh[tt, ff, ant, cc, par, t, f, p, q, c] += prefact * gp * model[d] * np.conj(gq)
												
	return jh