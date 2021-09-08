from __future__ import print_function
import numpy as np
from abc import ABCMeta

#Maybe treat this like the MasterMachine class from CubiCal?!
class MasterMachine(metaclass=ABCMeta):
	"""
	Provides base class for parametrised phase-only solver.
	And, lays the basic requirements for all machines.

	"""

	def __init__(self, args):
		"""
		Initialises a phase-only parametrised gain solver.

		Args:
			args (dict or namespace?)
			....

		"""

		self.args = args

	def jhj_compute(self):
		return NotImplementedError


	def update_compute(self):
		return NotImplementedError

	def residual_compute(self):
		return NotImplementedError

	def apply_gains(self):
		return NotImplementedError


class ParametrisedPhase(MasterMachine):
	"""
	This class implements phase-only parametrised gain machine with a polynomial
	basis.

	"""

	def __init__(self, args):
		"""

		"""

		MasterMachine.__init__(self)
		self.args = args
		print(self.args)

		self.n_param = options.get("pphase-nparam", sources.shape[0])

		#Set random seed for alpha.
		np.random.seed(3)
		self.alpha = np.zeros((self.n_timint, self.n_freint, self.n_ant, self.n_param, self.n_cor))
		
		self.basis = get_basis(sources)


	def init_gains(self):
		"""

		"""
		
		self.gains = np.empty(self.gain_shape, dtype=self.dtype)
		


