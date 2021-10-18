import numpy as np
from paraphase.calibration.basis_compute import basis_compute
from paraphase.calibration.gains_compute import gains_compute
from paraphase.derivatives.j_compute import j_compute
from paraphase.derivatives.jh_compute import jh_compute
from paraphase.derivatives.jhj_compute import jhj_compute
from paraphase.calibration.residuals_compute import residuals_compute
from paraphase.derivatives.jhr_compute import jhr_compute


def calibratewith(data_arr, model_arr, msrcs, bparams, gparams, sparams, datadiag=None):
	"""
	The function calibrates given data with specified 'solver'.

	"""
	
	#
	basis = basis_compute(msrcs, bparams)
	#For unity gains.
	alpha = np.zeros(gparams["alpha_shape"], dtype=float)
	gains = gains_compute(basis, alpha, gparams)
	jac = j_compute(data_arr, model_arr, gains, gparams, alpha, basis, datadiag=datadiag)
	jh = jh_compute(jac)
	jhj = jhj_compute(jh, jac, alpha)
	residuals = residuals_compute(data_arr, model_arr, gains, datadiag=datadiag)
	jhr = jhr_compute(jh, residuals)

	"""
	while True:
		#get J
		#get JHJ, JHr
		#compute update

		if ...

	"""

	return data_arr
		


