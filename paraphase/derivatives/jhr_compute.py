import numpy as np

def jhr_compute(jh, residuals):
	"""
	Return JHr.

	"""

	return np.dot(jh, residuals)