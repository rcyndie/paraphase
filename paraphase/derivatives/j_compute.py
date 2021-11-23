import numpy as np
	

def j_compute(data_arr, model_arr, lcoord, mcoord, alpha, chan_freq, freqf, basis, datadiag=None):
	"""
	Evaluates Jacobian and treat 'basis' as an operator.

	"""

	n_dir = model_arr.shape[0]
	n_tim = data_arr.shape[0]
	n_fre = data_arr.shape[1]
	n_timint, n_freint, n_ant, n_ccor, n_par = alpha.shape

	#Remember, t_int is size of the interval and
	#n_timint is the number of intervals.
	t_int = n_tim//n_timint
	f_int = n_fre//n_freint

	#Initialise Jacobian.
	jac = np.zeros(data_arr.shape+alpha.shape, dtype=data_arr.dtype)

	for t in range(n_tim):
		tt = t//t_int
		for f in range(n_fre):
			ff = f//f_int
			fprof = freqf(chan_freq[f])
			for p in range(n_ant):
				for q in range(n_ant):
					for c in range(n_ccor):
						# end of data axes
						if datadiag:
							model = model_arr[:, t, f, p, q, c]
						else:
							model = model_arr[:, t, f, p, q, c, c]
						# start of param axes
						for par in range(n_par):
							for d in range(n_dir):
								#Get partial derivative of the phase.
								L = basis(lcoord[d], mcoord[d])
								Lalphap = L.dot(alpha[tt, ff, p, c])
								Lalphaq = L.dot(alpha[tt, ff, q, c])
								prefact = 1j * L[par] * fprof
								gp = np.exp(1.0j * fprof * Lalphap)
								gq = np.exp(1.0j * fprof * Lalphaq)
								if datadiag:
									jac[t, f, p, q, c, tt, ff, p, c, par] += prefact * gp * model[d] * np.conj(gq)
									jac[t, f, p, q, c, tt, ff, q, c, par] += -prefact * gp * model[d] * np.conj(gq)
								else:
									jac[t, f, p, q, c, c, tt, ff, p, c, par] += prefact * gp * model[d] * np.conj(gq)
									jac[t, f, p, q, c, c, tt, ff, q, c, par] += -prefact * gp * model[d] * np.conj(gq)
												
	return jac


def j_compute_slow(data_arr, model_arr, freq, alpha, fprofile, basis, lcoord, mcoord):
	"""
	Returns the Jacobian.

	"""

	ntime, nchan, nant, _, ncorr = data_arr.shape
	nsource = model_arr.shape[-1]
	ntimeint, nchanint, nant, ncorr, nparam = alpha.shape

	#Initialise Jacobian.
	jac = np.zeros(data_arr.shape+alpha.shape, dtype=data_arr.dtype)

	for t in range(ntime):
		for f in range(nchan):
			fprof = fprofile(freq[f])
			for p in range(nant):
				for q in range(nant):
					for c in range(ncorr):
						# end of data axes
						model = model_arr[:, t, f, p, q, c]
						# start of param axes
						for tt in range(ntimeint):
							for ff in range(nchanint):
								for ant in range(nant):
									for cc in range(ncorr):
										for par in range(nparam):
											for s in range(nsource):
												#Get partial derivative of the phase.
												if ant==p and p!=q and cc==c and tt==t//ntimeint and ff==f//nchanint:
													L = basis(lcoord[s], mcoord[s])
													Lalpha = L.dot(alpha[tt, ff, ant, cc])
													prefact = 1j * L[par] * fprof
													phase = fprof * Lalpha
													gp = np.exp(1j * phase)
													gq = np.exp(1.0j*fprof*L.dot(alpha[tt, ff, q, cc]))
													jac[t, f, p, q, c, tt, ff, ant, cc, par] += prefact * gp * model[s] * np.conj(gq)
												elif ant==q and q!=p and cc==c and tt==t//ntimeint and ff==f//nchanint:
													L = basis(lcoord[s], mcoord[s])
													Lalpha = L.dot(alpha[tt, ff, ant, cc])
													prefact = -1j * L[par] * fprof
													phase = fprof * Lalpha
													gq = np.exp(1j * phase)
													gp = np.exp(1.0j*fprof*L.dot(alpha[tt, ff, p, cc]))
													jac[t, f, p, q, c, tt, ff, ant, cc, par] += prefact * gp * model[s] * np.conj(gq)


	return jac
