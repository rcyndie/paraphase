import numpy as np
from paraphase.calibration.basis_compute import basis_compute
from paraphase.calibration.gains_compute import gains_compute


def calibratewith(data, msrcs, bparams, alpha_shape, gains_shape, tol, solver):
	"""
	The function calibrates given data with specified 'solver'.

	"""

	#compute basis
	#get gains (also input gains_shape)
	#compute J
	#compute JH
	#compute JHJ
	#get residuals
	#compute JHr

	#
	basis = basis_compute(msrcs, bparams, solver)
	# gains_compute()

	"""
	while True:
		#get J
		#get JHJ, JHr
		#compute update

		if ...

	"""


	return data
		


