import numpy as np

def j_compute(data_arr, model_arr, gains, gparams, alpha, basis, datadiag=None):
	"""
	Returns the Jacobian.

	"""

	n_dir, n_timint, n_fre, n_ant, n_ccor = gains.shape
	n_tim = data_arr.shape[0]
	n_par = gparams["n_par"]
	n_freint = gparams["alpha_shape"][1]

	#Initialise Jacobian.
	jac = np.zeros(data_arr.shape+alpha.shape, dtype=data_arr.dtype)

	for t in range(n_tim):
		tt = t//n_timint
		for f in range(n_fre):
			ff = f//n_freint
			for p in range(n_ant):
				for q in range(p):  #note only doing this for q < p
					for d in range(n_dir):
						for k in range(n_ccor):
							for param in range(n_par):
								#Get partial derivative of the phase.
								dphidalpha = 1.0j * gparams["chan_freq"][f] * basis[d, param]
								#if diagonal data, consider,
								if datadiag:
									jac[t, f, p, q, k, tt, ff, p, param, k] += dphidalpha * gains[d, tt, f, p, k]* model_arr[d, t, f, p, q, k] * np.conj(gains[d, tt, f, q, k].T)
									jac[t, f, p, q, k, tt, ff, q, param, k] += -dphidalpha * gains[d, tt, f, p, k]* model_arr[d, t, f, p, q, k] * np.conj(gains[d, tt, f, q, k].T)
								# else:
								# 	jac[t, f, p, q, k, k, tt, ff, p, param, k] += dphidalpha * gains[d, tt, f, p, k]* model_arr[d, t, f, p, q, k, k] * np.conj(gains[d, tt, f, q, k].T)
								# 	jac[t, f, p, q, k, k, tt, ff, q, param, k] += -dphidalpha * gains[d, tt, f, p, k]* model_arr[d, t, f, p, q, k, k] * np.conj(gains[d, tt, f, q, k].T)
					# Set [q,p] element as conjugate of [p,q].
					jac[t, f, q, p] = np.conj(jac[t, f, p, q])

	return np.reshape(jac, (data_arr.size, alpha.size))

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
						model = model_arr[t, f, p, q, c]
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
