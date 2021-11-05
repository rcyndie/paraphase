import numpy as np

def jh_compute(jac):
	"""
	Return complex conjugate transpose of Jacobian.

	"""

	return np.conjugate(jac.T)