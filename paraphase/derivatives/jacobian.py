import numpy as np

def compute_jacobian_residual(self, data_arr, model_arr, gains):
	"""
	Returns the Jacobian.

	"""

	#import pdb; pdb.set_trace()

	#Initialise Jacobian.
	self.jac_shape = [self.n_tim, self.n_fre, self.n_ant, self.n_ant, self.n_cor, self.n_timint, self.n_freint, self.n_ant, self.n_param, self.n_cor] 
	jac = np.zeros(self.jac_shape, dtype=self.dtype)

	for t in range(self.n_tim):
		tt = t//self.t_int
		for f in range(self.n_fre):
			ff = f//self.f_int
			for p in range(self.n_ant):
				for q in range(p):  #note only doing this for q < p
					for s in range(self.n_dir):
						for k in range(self.n_cor):
							#Get Jacobian.
							for param in range(self.n_param):
								#Get partial derivative of the phase.
								dphidalpha = 1.0j * self.chunk_fs[f] * self.basis[param, s]
								jac[t, f, p, q, k, tt, ff, p, param, k] += dphidalpha * gains[s, tt, f, p, k, k] * model_arr[s, 0, t, f, p, q, k, k] * np.conj(gains[s, tt, f, q, k, k]).T #I do not need to transpose gains_q (scalar).
								jac[t, f, p, q, k, tt, ff, q, param, k] += -dphidalpha * gains[s, tt, f, p, k, k] * model_arr[s, 0, t, f, p, q, k, k] * np.conj(gains[s, tt, f, q, k, k]).T

					#Set [q,p] element as conjugate of [p,q].
					jac[t, f, q, p] = np.conj(jac[t, f, p, q])

	##Reshape the Jacobian to a 2D shape.
	jac = np.reshape(jac, (self.n_tim*self.n_fre*self.n_ant*self.n_ant*self.n_cor, self.n_timint*self.n_freint*self.n_ant*self.n_param*self.n_cor))


	return jac