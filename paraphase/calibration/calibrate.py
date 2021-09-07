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

		self.args = args
		print(self.args)


