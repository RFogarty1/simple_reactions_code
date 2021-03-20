
import copy
import itertools as it
import math
import unittest
import unittest.mock as mock

import simple_reactions_lib.core.core_units as unitHelp
import simple_reactions_lib.core.core_classes as tCode


class TestChemSpeciesBase(unittest.TestCase):

	def setUp(self):
		self.name = "fake_name"
		self.conc = 4.7
		self.concFmt = "{:.5g}"
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.ChemSpeciesStd(self.name, self.conc, strConcFmt=self.concFmt)

	def testExpectedStrGiven(self):
		expStr = "ChemSpeciesStd: fake_name conc=" + self.concFmt.format(self.conc)
		actStr = str(self.testObjA)
		self.assertEqual(expStr, actStr)

	def testEqualObjsCompareEqual(self):
		objA = copy.deepcopy(self.testObjA)
		self.createTestObjs()
		objB = self.testObjA
		self.assertEqual(objA, objB)
		
	def testUnequalObjsCompareUnequal_diffNames(self):
		objA = copy.deepcopy(self.testObjA)
		self.name = self.name + "_extended"
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)

	def testUnequalObjsCompareUnequal_diffConcs(self):
		objA = copy.deepcopy(self.testObjA)
		self.conc += 1.2
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA,objB)

class TestStandardReactionTemplate(unittest.TestCase):

	def setUp(self):
		self.barrier = 1.5
		self.prefactor = 10
		self.timeStep = 4
		self.temperature = 50
		self.reactantA = tCode.ChemSpeciesStd("Mg", 5)
		self.productA = tCode.ChemSpeciesStd("X", 2)
		self.reactants, self.products = ["Mg"], ["X"]
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.ChemReactionTemplate(self.reactants, self.products, self.barrier, self.prefactor)

	@mock.patch("simple_reactions_lib.core.core_classes.ChemReactionTemplate._getReactantConcRateFactor")
	@mock.patch("simple_reactions_lib.core.core_classes.ChemReactionTemplate._getTafelFactor")
	@mock.patch("simple_reactions_lib.core.core_classes.ChemReactionTemplate._getk0")
	def testExpectedRateSimple(self, mockedGetK0, mockedGetTafel, mockedGetReactantConcFactor):
		#1) mock things
		inputReactants, temp, pH, pot = mock.Mock(), mock.Mock(), mock.Mock(), mock.Mock()
		expk0, expkTafel, expReactantConcRateFactor = 2, 4, 6
		mockedGetK0.side_effect = lambda *args,**kwargs: expk0
		mockedGetTafel.side_effect = lambda *args,**kwargs: expkTafel
		mockedGetReactantConcFactor.side_effect = lambda *args,**kwargs: expReactantConcRateFactor

		expRate = expk0*expkTafel*expReactantConcRateFactor
		actRate = self.testObjA.getReactionRate(inputReactants, temp, pH=pH, potential=pot)

		mockedGetK0.assert_called_with(temp)
		mockedGetTafel.assert_called_with(inputReactants, temp, pH, pot)
		mockedGetReactantConcFactor.assert_called_with(inputReactants)
		self.assertEqual(expRate, actRate)

	@mock.patch("simple_reactions_lib.core.core_classes.ChemReactionTemplate.getReactionRate")
	def testChangesInReactants(self, mockGetReactionRate):
		timeStep, rate = 3, 2
		temperature = 50
		inputReactants = [self.reactantA, self.productA]
		reactants, products = ["Mg"], ["X"]
		mockGetReactionRate.side_effect = lambda *args,**kwargs: rate

		expDict = {"Mg":-1*timeStep*rate, "X":timeStep*rate}
		actDict = self.testObjA.getChangesInReactants(timeStep, inputReactants, temperature)

		self.assertEqual(expDict, actDict)

	def testGetReactantConcFactor(self):
		reactants = ["A","B"]
		concs = [2,3]
		self.testObjA.reactants = reactants
		reactantA, reactantB = [tCode.ChemSpeciesStd(name,conc) for name,conc in zip(reactants,concs)]
		expFactor = 2*3
		actFactor = self.testObjA._getReactantConcRateFactor([reactantA,reactantB])
		self.assertEqual(expFactor,actFactor)

	@mock.patch("simple_reactions_lib.core.core_classes.unitHelp")
	def testGetk0(self,mockedUnitHelp):
		fakeBoltzFactor = 2 #Makes math less annoying
		mockedUnitHelp.BOLTZ_EV = fakeBoltzFactor
		exp_k0 = self.prefactor * math.exp( (-1*self.barrier) / (fakeBoltzFactor*self.temperature) )
		act_k0 = self.testObjA._getk0(self.temperature)
		self.assertAlmostEqual(exp_k0,act_k0)



#TODO: Do equality method for ChemSpeciesStd first
class TestReactionController(unittest.TestCase):

	def setUp(self):
		self.timeStep = 2
		self.temperature = 300
		self.pH = 4
		self.potential = 2.5
		self.startReactants = [tCode.ChemSpeciesStd(name,conc) for name,conc in zip(["Mg","H"], [1,1])]

		self.changeA = {"Mg":0.2,"H":-0.1}
		self.changeB = {"Mg":-0.3, "H":0.1}
		self.createTestObjs()

	def createTestObjs(self):
		self.reactionA = mock.Mock()
		self.reactionB = mock.Mock()

		self.reactionA.getChangesInReactants.side_effect = lambda *args, **kwargs: self.changeA
		self.reactionB.getChangesInReactants.side_effect = lambda *args, **kwargs: self.changeB

		#Create testObjA
		reactions = [self.reactionA, self.reactionB]
		args = [self.startReactants, reactions, self.timeStep]
		kwargs = {"temperature":self.temperature, "pH":self.pH, "potential":self.potential}

		self.testObjA = tCode.ReactionControllerStandard(*args, **kwargs)

	def testExpectedNextStepA(self):
		expOutReactants = [tCode.ChemSpeciesStd(name,conc) for name,conc in zip(["Mg","H"], [1+0.2-0.3, 1-0.1+0.1])]
		expStepNumber = 1
		self.testObjA._doNextStep()
		actOutReactants = self.testObjA.currentReactants
		actStepNumber = self.testObjA.step
		self.assertEqual(expOutReactants,actOutReactants)
		self.assertEqual(expStepNumber, actStepNumber)


	def testSimpleTrackStepCallbackFunct(self):
		nSteps = 3
		expStepList = [x for x in range(nSteps)]
		actStepList = list()
		def _appendStepNumberToList(instance):
			actStepList.append( instance.step )

		self.testObjA.callbackFunct = _appendStepNumberToList
		self.testObjA.doNextNSteps(nSteps)
		self.assertEqual(expStepList, actStepList)



class TestNetReactionTemplate(unittest.TestCase):

	def setUp(self):
		self.reactants = ["X","Y"] 
		self.products = ["Z"]
		self.reactantConcs = [2,2]
		self.productConcs = [2]

		self.timeStep = 2
		self.temperature = 300
		self.pH = 0 #Dud argument these days
		self.potential = 0

		self.barrier = 2
		self.prefactor = 10

		self.forwardRate = 10
		self.backwardRate = 15

		self.createTestObjs()

	def createTestObjs(self):
		#Get input reactants
		self.inpReactants = list()
		allSpecies = self.reactants + self.products
		allConcs = self.reactantConcs + self.productConcs
		for spec,conc in it.zip_longest(allSpecies,allConcs):
			self.inpReactants.append( tCode.ChemSpeciesStd(spec,conc) )


		self.forwardReaction  = tCode.ChemReactionTemplate(self.reactants, self.products, self.barrier, self.prefactor)
		self.backwardReaction = tCode.ChemReactionTemplate(self.products, self.reactants, self.barrier, self.prefactor)


		#Mock out the get-rate part
		self.forwardReaction.getReactionRate = mock.Mock()
		self.backwardReaction.getReactionRate = mock.Mock()
		self.forwardReaction.getReactionRate.side_effect = lambda *args,**kwargs: self.forwardRate
		self.backwardReaction.getReactionRate.side_effect = lambda *args,**kwargs: self.backwardRate

		#Create the overall object
		self.testObjA = tCode.NetReactionTemplate(self.forwardReaction, self.backwardReaction)

	def testGetReactionRate(self):
		# Setup
		expRate = self.forwardRate-self.backwardRate

		#run
		actRate = self.testObjA.getReactionRate(self.inpReactants, self.temperature, pH=self.pH, potential=self.potential)

		#check calls/value
		self.forwardReaction.getReactionRate.assert_called_with( self.inpReactants, self.temperature, pH=self.pH, potential=self.potential )
		self.backwardReaction.getReactionRate.assert_called_with( self.inpReactants, self.temperature, pH=self.pH, potential=self.potential )
		self.assertEqual(expRate, actRate)

	def testGetChangesInReactants(self):
		#setup
		netRate = self.forwardRate-self.backwardRate
		expChanges = {"X": -1*netRate*self.timeStep, "Y": -1*netRate*self.timeStep,
		              "Z": netRate*self.timeStep}
		
		#run. Note forward "reactants" should increase in conc, since netRate is negative
		actChanges = self.testObjA.getChangesInReactants(self.timeStep, self.inpReactants, self.temperature, pH=self.pH, potential=self.potential)

		self.assertEqual(expChanges, actChanges)

	def testGetChangesInReactants_speciesOnEachSide(self):
		#We have to do this without redoing sel.createTestObjs; since that specifically assumed no duplication of species
		self.products.append("X")

		netRate = self.forwardRate-self.backwardRate
		expChanges = {"X": 0, "Y": -1*netRate*self.timeStep, "Z": netRate*self.timeStep}
		actChanges = self.testObjA.getChangesInReactants(self.timeStep, self.inpReactants, self.temperature, pH=self.pH, potential=self.potential)

		self.assertEqual(expChanges, actChanges)

	def testGetChangesForRepeatedReactant(self):
		self.reactants.append("X")

		netRate = self.forwardRate-self.backwardRate
		expChanges = {"X":-2*netRate*self.timeStep, "Y":-1*netRate*self.timeStep, "Z":netRate*self.timeStep}
		actChanges = self.testObjA.getChangesInReactants(self.timeStep, self.inpReactants, self.temperature, pH=self.pH, potential=self.potential)

		self.assertEqual(expChanges, actChanges)


class TestBetterChemReactionTemplate(unittest.TestCase):

	def setUp(self):
		self.reactants = ["H","H"]
		self.products = ["H2"]
		self.barrier = 2
		self.prefactor = 10
		self.refPot = 2
		self.netElectronTransfer = 0
		self.temperature = 300
		self.potential = 1.5
		self.symFactor = 0.6
		self.inpReactants = mock.Mock()
		self.createTestObjs()

	def createTestObjs(self):
		args = [self.reactants, self.products, self.barrier, self.prefactor]
		kwargs = {"refPot":self.refPot, "nElecTransfer":self.netElectronTransfer, "symFactor":self.symFactor}
		self.testObjA = tCode.BetterReactionTemplate(*args, **kwargs)

	def _getReactionRate(self):
		return self.testObjA.getReactionRate( self.inpReactants, self.temperature, self.potential)

	@mock.patch("simple_reactions_lib.core.core_classes.BetterReactionTemplate._getReactantConcRateFactor")
	@mock.patch("simple_reactions_lib.core.core_classes.BetterReactionTemplate._getk0")
	def testGetReactionRate_Potential(self, mockGetk0, mockGetConcFactor):
		k0Val = 3
		concFactor = 4
		mockGetk0.side_effect = lambda *args,**kwargs: k0Val
		mockGetConcFactor.side_effect = lambda *args,**kwargs: concFactor

		expTotal = k0Val * concFactor
		actTotal = self._getReactionRate()

		mockGetk0.assert_called_with( self.temperature )
		mockGetConcFactor.assert_called_with( self.inpReactants ) 
		self.assertAlmostEqual(expTotal, actTotal)

	@mock.patch("simple_reactions_lib.core.core_classes.BetterReactionTemplate._getTafelFactor")
	@mock.patch("simple_reactions_lib.core.core_classes.BetterReactionTemplate._getReactantConcRateFactor")
	@mock.patch("simple_reactions_lib.core.core_classes.BetterReactionTemplate._getk0")
	def testGetReactionRate_mockTafelFactor(self, mockGetk0, mockGetConcFactor, mockGetTafelFactor):
		k0Val, concFactor, tafelFactor = 3, 4, 5
		mockGetk0.side_effect = lambda *args,**kwargs: k0Val
		mockGetConcFactor.side_effect = lambda *args,**kwargs: concFactor
		mockGetTafelFactor.side_effect = lambda *args,**kwargs: tafelFactor

		expTotal = k0Val*concFactor*tafelFactor
		actTotal = self._getReactionRate()
		self.assertAlmostEqual(expTotal, actTotal)

	def testGetTafelFactorCathodicA(self):
		self.netElectronTransfer = 1
		self.potential = self.refPot - 0.5
		self.createTestObjs()
		rt = unitHelp.IDEAL_GAS_R_JOULES * self.temperature
		faradOverRT = unitHelp.FARADAY_CONST / rt
		potDiff = self.potential - self.refPot
		expTafelFactor = math.exp( -1*(potDiff*self.symFactor*faradOverRT) )
		actTafelFactor = self.testObjA._getTafelFactor( self.temperature, self.potential )
		self.assertAlmostEqual(expTafelFactor, actTafelFactor)

	def testGetTafelFactorAnodicA_2Electrons(self):
		self.netElectronTransfer = -2
		self.potential = self.refPot + 0.4
		self.createTestObjs()

		rt = unitHelp.IDEAL_GAS_R_JOULES * self.temperature
		faradOverRT = unitHelp.FARADAY_CONST / rt
		potDiff = self.potential - self.refPot
		expTafelFactor = math.exp( 2*(potDiff*self.symFactor*faradOverRT) )
		actTafelFactor = self.testObjA._getTafelFactor( self.temperature, self.potential )

		self.assertAlmostEqual(expTafelFactor, actTafelFactor)





