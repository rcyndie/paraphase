import numpy as np

def j_compute(data_arr, model_arr, gains, gparams, alpha, basis):
	"""
	Returns the Jacobian.

	"""
	
	n_timint, n_freint, n_ant, n_par, n_ccor = alpha.shape
	n_dir, _, n_fre, _, _ = gains.shape
	n_tim = data_arr.shape[0]

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
								dphidalpha = 1.0j * gparams["chan_freq"][f] * basis[param, d]
								jac[t, f, p, q, k, k, tt, ff, p, param, k] += dphidalpha * gains[d, tt, f, p, k] * model_arr[d, t, f, p, q, k, k] * np.conj(gains[d, tt, f, q, k]).T
								jac[t, f, p, q, k, k, tt, ff, q, param, k] += -dphidalpha * gains[d, tt, f, p, k] * model_arr[d, t, f, p, q, k, k] * np.conj(gains[d, tt, f, q, k]).T

					#Set [q,p] element as conjugate of [p,q].
					jac[t, f, q, p] = np.conj(jac[t, f, p, q])

	jac = np.reshape(jac, (n_tim*n_fre*n_ant*n_ant*n_ccor*n_ccor, n_timint*n_freint*n_ant*n_par*n_ccor))

	return jac
	