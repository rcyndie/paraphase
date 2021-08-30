from __future__ import print_function
import numpy as np


#Maybe treat this like the MasterMachine class from CubiCal?!
class ParametrisedPhase:
	"""
	Implements parametrised phase-only solver.

	"""

	def __init__(self, args):
		"""
		Initialises a phase-only parametrised gain solver.

		**about arguments

		"""

		self.args = args

