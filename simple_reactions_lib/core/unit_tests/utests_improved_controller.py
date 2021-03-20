
import itertools as it
import unittest
import unittest.mock as mock

import simple_reactions_lib.core.core_classes as coreHelp
import simple_reactions_lib.core.improved_controller as tCode


class TestReactionControllerImproved(unittest.TestCase):

	def setUp(self):
		self.startReactants = [coreHelp.ChemSpeciesStd( "Mg", 2 )]
		self.temp = 25
		self.pot = 4
		self.createTestObjs()

	def createTestObjs(self):
		self.propagator = mock.Mock()
		self.testObjA = tCode.ReactionControllerImproved(self.propagator, self.startReactants, temperature=self.temp, potential=self.pot)

	def testResetFunction(self):
		expReactants = self.startReactants
		self.testObjA.currentReactants = [coreHelp.ChemSpeciesStd("O", 4)]
		self.assertNotEqual(expReactants, self.testObjA.currentReactants)
		self.testObjA.reset()
		self.assertEqual(expReactants, self.testObjA.currentReactants)

	#Failed to get assert called stuff to work properly so lazily left it out
	def testMoveForward(self):
		#setup
		def _propagate(inpReactants, *args, **kwargs):
			for reactant in inpReactants:
				if reactant.name=="Mg":
					reactant.conc -= 10
		self.propagator.propagate.side_effect = lambda inpReactants, *args, **kwargs: _propagate(inpReactants, *args, **kwargs)

		#test
		timeStep = 10
		expReactants = [coreHelp.ChemSpeciesStd( "Mg", -8 )]
		self.testObjA.moveForwardByT(timeStep)


		actReactants = self.testObjA.currentReactants
		self.assertEqual(expReactants, actReactants)


class TestPropagatorStandard(unittest.TestCase):

	def setUp(self):
		self.targTimeStep = 1.5
		self.nameA, self.concA = "Mg", 2
		self.stepSizes = [0.1,0.4,1]
		self.concChanges = [{self.nameA:1}, {self.nameA:3}, {self.nameA:4}]
		self.temp = 250
		self.pot = 2.5
		self.minConc = 0
		self.maxConc = 1
		self.createTestObjs()

	def createTestObjs(self):
		self.inpReactants = [coreHelp.ChemSpeciesStd( self.nameA, self.concA )]
		self.concChangesFinder = mock.Mock()
		retVals = it.zip_longest(self.stepSizes, self.concChanges)
		self.concChangesFinder.getConcChangesForNextTimeStep.side_effect = retVals
		self.testObjA = tCode.ConcsPropagatorStandard(self.concChangesFinder)

	def _runTestFunct(self):
		return self.testObjA.propagate(self.inpReactants, self.targTimeStep, temperature=self.temp, potential=self.pot)

	def testExpectedUpdateA(self):
		expReactants = [coreHelp.ChemSpeciesStd( self.nameA, 2+1+3+4 )]
		self._runTestFunct()

		expRequestedTimeSteps = [1.5, 1.4, 1.0]
		self.assertEqual( len(expRequestedTimeSteps), len(self.concChangesFinder.getConcChangesForNextTimeStep.call_args_list) )
		for idx,callArgs in enumerate(self.concChangesFinder.getConcChangesForNextTimeStep.call_args_list):
			args,kwargs = callArgs
			expTimeStep, actTimeStep = expRequestedTimeSteps[idx], args[1]
			self.assertAlmostEqual(expTimeStep, actTimeStep)

		self.assertEqual(expReactants, self.inpReactants)

	@unittest.skip("Bad idea; breaks conservation rules and hides numerical problems that shouldnt be hidden")
	def testUpdateConcsStaysWithinLimits(self):
		self.reactants = [coreHelp.ChemSpeciesStd("Mg",0.15), coreHelp.ChemSpeciesStd("O",0.9)]
		self.concChanges = {"Mg":-0.2, "O":0.15}

		self.testObjA._updateConcs(self.reactants, self.concChanges)

		expUpdatedVals = [coreHelp.ChemSpeciesStd("Mg",0), coreHelp.ChemSpeciesStd("O",1)]
		self.assertEqual( expUpdatedVals,self.reactants )


class TestStandardConcChangesFinder(unittest.TestCase):

	def setUp(self):
		self.nameA, self.nameB = "Mg", "O"
		self.rateA, self.rateB = 30, 60
		self.maxTimeStep = 10
		self.variableConcSpecies = [self.nameA, self.nameB]
		self.inpReactants = mock.Mock()
		self.temperature, self.potential = 300, 0
		self.maxConcChange = 0.1
		self.createTestObjs()

	def createTestObjs(self):
		self.rateDictA = {self.nameA:self.rateA, self.nameB:self.rateB}
		self.rateCalculatorA = mock.Mock()
		self.rateCalculatorA.getRates.side_effect = lambda *args,**kwargs: self.rateDictA

		self.testObjA = tCode.ConcChangesFinderStandard(self.rateCalculatorA, self.variableConcSpecies, maxConcChange=self.maxConcChange)

	def _runTestFunct(self):
		currKwargs = {"temperature":self.temperature, "potential":self.potential}
		return self.testObjA.getConcChangesForNextTimeStep(self.inpReactants, self.maxTimeStep, **currKwargs)

	def testExpectedChangesAndTimestepsA(self):
		expChanges = {self.nameA:0.05, self.nameB:0.1}
		expTimeStep = self.maxConcChange / self.rateB
		actTimeStep, actChanges = self._runTestFunct()

		self.assertAlmostEqual(expTimeStep, actTimeStep)
		self._checkExpAndActConcChangesEqual( expChanges, actChanges )

	def testExpectedChanges_bFixedConc(self):
		self.variableConcSpecies = [self.nameA]
		self.createTestObjs()
		expTimeStep = self.maxConcChange / self.rateA
		expChanges = {self.nameA:0.1}

		actTimeStep, actChanges = self._runTestFunct()
		self.assertAlmostEqual(expTimeStep, actTimeStep)
		self._checkExpAndActConcChangesEqual(expChanges, actChanges)

	def testExpectedUsingMaxTimeStep(self):
		self.variableConcSpecies = [self.nameA]
		self.rateA = 0.01
		self.maxTimeStep = 8
		self.createTestObjs()
		expTimeStep = self.maxTimeStep
		expChanges = {self.nameA:0.08}
		actTimeStep, actChanges = self._runTestFunct()
		self.assertAlmostEqual(expTimeStep, actTimeStep)
		self._checkExpAndActConcChangesEqual(expChanges, actChanges)

	def _checkExpAndActConcChangesEqual(self, expChanges, actChanges):
		self.assertEqual( expChanges.keys(), actChanges.keys() )
		for key in expChanges.keys():
			self.assertAlmostEqual( expChanges[key], actChanges[key] )


class TestRateCalculatorStandard(unittest.TestCase):

	def setUp(self):
		self.nameA, self.concA = "Mg", 2
		self.nameB, self.concB = "O", 3
		self.reactionA_returnA, self.reactionA_returnB = 5, -3 
		self.reactionB_returnA, self.reactionB_returnB = 2, 4

		self.temp = 40
		self.pot = 2

		self.createTestObjs()

	def createTestObjs(self):
		self.chemA = coreHelp.ChemSpeciesStd( self.nameA, self.concA )
		self.chemB = coreHelp.ChemSpeciesStd( self.nameB, self.concB )
		self.inpReactants = {self.nameA:self.chemA, self.nameB:self.chemB}
		self.reactionA, self.reactionB = mock.Mock(), mock.Mock()


		reactionAReturnVals = {self.nameA:self.reactionA_returnA, self.nameB:self.reactionA_returnB}
		reactionBReturnVals = {self.nameA:self.reactionB_returnA, self.nameB:self.reactionB_returnB}
		self.reactionA.getChangesInReactants.side_effect = lambda *args,**kwargs: reactionAReturnVals
		self.reactionB.getChangesInReactants.side_effect = lambda *args,**kwargs: reactionBReturnVals


		self.testObjA = tCode.RateCalculatorStandard([self.reactionA, self.reactionB])

	def _runTestFunct(self):
		return self.testObjA.getRates(self.inpReactants, temperature=self.temp, potential=self.pot)

	def testExpectedRatesA(self):
		expRates = {self.nameA: self.reactionA_returnA+self.reactionB_returnA,
		            self.nameB: self.reactionA_returnB+self.reactionB_returnB}
		actRates = self._runTestFunct()

		#Check expected calls made
		self.reactionA.getChangesInReactants.assert_called_with(1, self.inpReactants, self.temp, potential=self.pot)
		self.reactionB.getChangesInReactants.assert_called_with(1, self.inpReactants, self.temp, potential=self.pot)
		self.assertEqual(expRates, actRates)





