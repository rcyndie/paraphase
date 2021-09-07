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

	def compute_jhj(self):
		return NotImplementedError


	def compute_update(self):
		return NotImplementedError

	def compute_residual(self):
		return NotImplementedError

	def apply_gains(self):
		return NotImplementedError

