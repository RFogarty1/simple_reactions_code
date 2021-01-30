

import copy
import itertools as it
import math
from . import core_units as unitHelp

class ChemSpeciesStd():
	""" Simple class for holding the name and concentration of a chemical species

	Attributes:
		name: (str) Name used as a label (e.g. "h_ads" for adsorbed hydrogen) 
		conc: (float) Number representing the concentration of these species (units up to the user)

	"""
	def __init__(self, name, conc, strConcFmt="{:.3g}"):
		self._eqTol = 1e-5
		self.name = name
		self.conc = conc
		self.strConcFmt = strConcFmt

	def __str__(self):
		return "ChemSpeciesStd: {} conc=".format(self.name) + self.strConcFmt.format(self.conc)

	def __eq__(self, other):
		eqTol = min(self._eqTol, other._eqTol)
		if self.name != other.name:
			return False

		concDiff = abs(self.conc-other.conc)
		if concDiff>eqTol:
			return False

		return True



class ChemReactionTemplate():

	def __init__(self, reactants, products, barrier, prefactor):
		self.reactants = reactants
		self.products = products
		self.barrier = barrier
		self.prefactor = prefactor

	def getChangesInReactants(self, timeStep, inputReactants, temperature, pH=0, potential=0):
		""" Get the changes in concentration (or amounts???) of inputReactants
		
		Args:
			timeStep: (float) Timestep to use 
			inputReactants: (iter of ChemSpeciesStd) These include the INPUT concentrations
			temperature: (float) in Kelvin

		Returns
			reactantChanges: (dict) Keys are reactant names, vals are changes in concentrations from THIS timestep
	 
		"""
		reactRate = self.getReactionRate(inputReactants, temperature ,pH, potential)
		allSpecies = set(self.reactants + self.products)
		outDict = {k:0 for k in allSpecies}

		for reactant in self.reactants:
			outDict[reactant] -= reactRate*timeStep
		for product in self.products:
			outDict[product] += reactRate*timeStep

		return outDict

	def getReactionRate(self, inputReactants, temperature, pH=0, potential=0):
		""" Gets the rate for this reaction at input conditions
		
		Args:
			inputReactants: (iter of ChemSpeciesStd) 
			temperature: (float)
			pH: (float)
			potential: (float)

		Returns
			rate: The rate of this reaction per unit-time
	 
		"""
		k0Val = self._getk0(temperature)
		tafelFactor = self._getTafelFactor(inputReactants, temperature, pH, potential)
		concFactor = self._getReactantConcRateFactor(inputReactants)
		return k0Val*tafelFactor*concFactor

	#This is actually doable
	def _getk0(self, temperature):
		return self.prefactor*math.exp(  (-1*self.barrier)/(unitHelp.BOLTZ_EV*temperature) )

	def _getTafelFactor(self, inputReactants, temperature, pH, potential):
		raise NotImplementedError("")

	def _getReactantConcRateFactor(self, inputReactants):
		""" The factor in rate which comes from the reactant concentrations, e.g. for H + H -> H2 its [H]^2, where [] denotes concentration """
		outFactor = 1
		for reactant in self.reactants:
			inpConcs = [x.conc for x in inputReactants if x.name==reactant]
			assert len(inpConcs)<2
			if len(inpConcs)==1:
				outFactor *= inpConcs[0]

		return outFactor


class ReactionControllerStandard():

	def __init__(self, startReactants, reactions, timeStep, temperature=300, pH=0, potential=0, callbackFunct=None):
		self.startReactants = startReactants
		self.reactions = reactions
		self.timeStep = timeStep
		self.temperature = temperature
		self.pH = pH
		self.potential = potential
#		self.currentReactants = copy.deepcopy(self.startReactants)
#		self.step = 0
		self.callbackFunct = callbackFunct
		self.reset()

	def reset(self):
		self.step = 0
		self.currentReactants = copy.deepcopy(self.startReactants)

	def doNextNSteps(self, n):
		for x in range(n):
			self._doNextStep()

	def _doNextStep(self):

		if self.callbackFunct is not None:
			self.callbackFunct(self)



		concChangeDict = self._getConcDiffs()
		#Update the currentReactants accordingly
		for key,val in concChangeDict.items():
			for reactant in self.currentReactants:
				if reactant.name==key:
					reactant.conc += val

		#Also update the step number
		self.step += 1

	def _getConcDiffs(self):
		outList = list()
		#Handle all the reactions
		currArgs = [self.timeStep, self.currentReactants, self.temperature]
		currKwargs = {"pH":self.pH, "potential":self.potential}
		for reaction in self.reactions:
			currDict = reaction.getChangesInReactants(*currArgs, **currKwargs)
			outList.append(currDict)

		#Figure out the total concentration diffs
		allKeys = [ [x for x in x.keys()] for x in outList]
		allKeys = set( [x for x in it.chain(*allKeys)] )

		outDict = {k:0 for k in allKeys}
		for currDict in outList:
			for key,val in currDict.items():
				outDict[key] += val
		return outDict
