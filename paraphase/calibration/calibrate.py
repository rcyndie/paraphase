import numpy as np
import time
from tqdm import tqdm

from paraphase.calibration.basis_compute import basis_compute
from paraphase.calibration.gains_compute import gains_compute
from paraphase.derivatives.j_compute import j_compute
from paraphase.derivatives.jh_compute import jh_compute
from paraphase.derivatives.jhj_compute import jhj_compute
from paraphase.calibration.residuals_compute import residuals_compute
from paraphase.derivatives.jhr_compute import jhr_compute
from paraphase.calibration.plot_generate import chi2_plot


def calibratewith(data_arr, model_arr, msrcs, bparams, gparams, sparams, datadiag=None):
	"""
	The function calibrates given data with specified 'solver'
	(Consider a diagonal block-wise JHJ).

	"""
	
	#
	basis = basis_compute(msrcs, bparams)
	#For unity gains.
	alpha = np.zeros(gparams["alpha_shape"], dtype=float)
	gains = gains_compute(basis, alpha, gparams)

	#Compute an initial Chi2 value.
	residuals = residuals_compute(data_arr, model_arr, gains, datadiag=datadiag)
	dof = data_arr.size - alpha.size
	chi20 = (np.linalg.norm(residuals)) / dof
	chi2_arr = np.array([chi20])

	
	for itern in tqdm(range(sparams["itermax"]), desc="DDCalibrating"):
		jac = j_compute(data_arr, model_arr, gains, gparams, alpha, basis, datadiag=datadiag)
		jh = jh_compute(jac)
		jhj = jhj_compute(jh, jac, alpha)
		jhr = jhr_compute(jh, residuals)

		try:
			#The GN updates differ for the gtype options.
			delta_alpha = delta_compute(alpha, jhj, jhr, gparams, sparams)

		except Exception as e:
		#In case, delta_alpha (or jhj) approaches a singular matrix.
			print(e)
			delta_alpha = np.zeros(gparams["alpha_shape"], dtype=alpha.dtype)

		#Compute update step.
		alpha_old = alpha.copy()
		alpha = alpha + delta_alpha
		alpha = (alpha + alpha_old) / 2.

		gains = gains_compute(basis, alpha, gparams)
		residuals = residuals_compute(data_arr, model_arr, gains, datadiag=datadiag)
		chi2i = (np.linalg.norm(residuals)) / dof
		chi2_arr = np.append(chi2_arr, chi2i)
		
		if chi2i < sparams["deltachi"]:
			print("Reached solution stagnancy at iternation number ", itern)
			break

		#Progress bar?!
		time.sleep(0.01)
	
	np.save(sparams["outputdir"]+"/alpha.npy", alpha)
	chi2_plot(chi2_arr, sparams)

def delta_compute(alpha, jhj, jhr, gparams, sparams):
	"""
	The function computes delta_alpha.

	"""

	n_timint, n_freint, n_ant, n_par, n_ccor = gparams["alpha_shape"]

	#Initialise delta_alpha.
	delta_alpha = np.zeros((n_timint*n_freint*n_ant, n_par*n_ccor), dtype=alpha.dtype)
	x = np.reshape(alpha, delta_alpha.shape)

	for k in range(n_timint*n_freint*n_ant):
		block_jhj = jhj[k*n_par*n_ccor:(k+1)*n_par*n_ccor, k*n_par*n_ccor:(k+1)*n_par*n_ccor]
		if gparams["gtype"] == "ppoly":
			delta_alpha[k] = (np.linalg.solve(block_jhj + sparams["lambda1"]*np.diag(block_jhj), jhr[k*n_par*n_ccor:(k+1)*n_par*n_ccor])).real
		elif gparams["gtype"] == "pcov":
			delta_alpha[k] = (np.linalg.solve(block_jhj + np.eye(block_jhj.shape[0]), jhr[k*n_par*n_ccor:(k+1)*n_par*n_ccor] - x[k])).real
	
	return delta_alpha.reshape(gparams["alpha_shape"])



