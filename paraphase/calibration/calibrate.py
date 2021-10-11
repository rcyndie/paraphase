import numpy as np
from paraphase.calibration.basis_compute import basis_compute
from paraphase.calibration.gains_compute import gains_compute
from paraphase.derivatives.j_compute import j_compute


def calibratewith(data, msrcs, bparams, gparams, tol):
	"""
	The function calibrates given data with specified 'solver'.

	"""

	#compute basis <<< done!
	#get gains (also input gains_shape)
	#compute J
	#compute JH
	#compute JHJ
	#get residuals
	#compute JHr

	#
	basis = basis_compute(msrcs, bparams)
	#For unity gains.
	alpha = np.zeros(gparams["alpha_shape"], dtype=float)
	gains = gains_compute(basis, alpha, gparams)
	jac = j_compute(data, model, gains)

	"""
	while True:
		#get J
		#get JHJ, JHr
		#compute update

		if ...

	"""


	return data
		


