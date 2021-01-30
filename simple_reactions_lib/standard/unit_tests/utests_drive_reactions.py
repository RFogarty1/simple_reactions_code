
import unittest
import unittest.mock as mock

import simple_reactions_lib.standard.drive_reactions as tCode


class TestDriveUntilSteadyState(unittest.TestCase):

	def setUp(self):
		self.nStepsPerCheck = 400
		self.speciesNames = ["Mg","X"]
		self.concChanges = [ {"Mg":50,"X":30},
		                     {"Mg":3, "X":10},
		                     {"Mg":2, "X":3} ]
		self.tolerancePerStep = 4

		self.createTestObjs()

	def createTestObjs(self):
		self.controller = mock.Mock()
		self.controller.idx = 0

		def _advanceNSteps(n):
			self.controller.idx += 1

		def _getConcDiffs():
			return self.concChanges[self.controller.idx]

		self.controller.doNextNSteps.side_effect = _advanceNSteps
		self.controller._getConcDiffs.side_effect = _getConcDiffs

	def _runTestFunct(self):
		args = [self.controller, self.speciesNames, self.nStepsPerCheck, self.tolerancePerStep]
		tCode.runControllerUntilSteadyState(*args)

	def testTerminatesAtExpectedTime(self):
		expIdx = 2
		self._runTestFunct()
		actIdx = self.controller.idx
		self.assertEqual(expIdx, actIdx)



