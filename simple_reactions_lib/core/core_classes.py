

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


class ChemReactionBase():

	def getChangesInReactants(self, timeStep, inputReactants, temperature, potential=0, **kwargs):
		raise NotImplementedError("")

	def getReactionRate(self, inputReactants, temperature, potential=0, **kwargs):
		raise NotImplementedError("")

class ChemReactionTemplate(ChemReactionBase):

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


class BetterReactionTemplate(ChemReactionTemplate):
	""" Improved reaction template class which has a built in method for dealing with the Tafel factor """

	def __init__(self,reactants, products, barrier, prefactor, refPot=0, nElecTransfer=0, symFactor=0.5):
		""" Initializer
		
		Args:
			reactants: (iter of str) Each represents a chemical species (e.g. h_2) 
			products: (iter of str) Each represents a chemical species
			barrier: (float) The barrier for the reaction. Needs to be in electron volts
			prefactor: (float) The prefactor for the Arrhenius term ("attempt frequency") Units are inverse time; generally inverse seconds
			refPot: (float) The reference potential for this reaction. This corresponds to the potential at which the barrier is correct (generally 0)
			nElecTransfer: (int) The net number of electrons transferred in this reaction. Cathodic reactions will be positive; anodic negative. E.g. Mg->Mg2+ + 2e- would have nElecTransfer=2
			symFactor: (float) Should be between 0 and 1. Symmetry factor (may als be called charge transfer coefficient?) that appears in the Tafel term. A value of 1 means the transition state has the same charge as the product, 0 means it has the same charge as reactant while 0.5 means its midway between the two (alternatively "charge" can be substituted by "same dependence on potential as"; i.e. 1 means the transition state depends on potential in the same way as the product).
 
		"""
		self.reactants = reactants
		self.products = products
		self.barrier = barrier
		self.prefactor = prefactor
		self.refPot = refPot
		self.nElecTransfer = nElecTransfer
		self.symFactor = symFactor

	def getReactionRate(self, inputReactants, temperature, pH=0, potential=0):
		k0Val = self._getk0(temperature)
		tafelFactor = self._getTafelFactor(temperature, potential)
		concFactor = self._getReactantConcRateFactor(inputReactants)
		return k0Val*tafelFactor*concFactor

	def _getTafelFactor(self, temperature, potential):
		overPotential = potential - self.refPot
		expTerm = (-1*self.nElecTransfer*unitHelp.FARADAY_CONST*overPotential*self.symFactor) / (temperature*unitHelp.IDEAL_GAS_R_JOULES) #The -1 is because we define nElecTransfer to be +ve for cathodic processes
		return math.exp(expTerm)

class NetReactionTemplate(ChemReactionTemplate):
	""" Class representing both forward and backwards reactions simultaneously, with the interface of a single reaction """

	def __init__(self, forwardReaction, backwardReaction):
		""" Initializer
		
		Args:
			forwardReaction: (ChemReactionTemplate)  
			backwardReaction: (ChemReactionTemplate)
				 
		"""
		self.forwardReaction = forwardReaction
		self.backwardReaction = backwardReaction

		assert self.forwardReaction.reactants==self.backwardReaction.products
		assert self.forwardReaction.products==self.backwardReaction.reactants

	@property
	def reactants(self):
		return self.forwardReaction.reactants

	@property
	def products(self):
		return self.forwardReaction.products

	def getChangesInReactants(self, timeStep, inputReactants, temperature, pH=0, potential=0):
		netChange = timeStep*self.getReactionRate(inputReactants, temperature, pH=pH, potential=potential)
		allSpecies = set(self.reactants + self.products)
		outDict = {k:0 for k in allSpecies}
		
		for reactant in self.reactants:
			outDict[reactant] += -1*netChange
		for product in self.products:
			outDict[product] += netChange

		return outDict

	def getReactionRate(self, inputReactants, temperature, pH=0, potential=0):
		forwardRate = self.forwardReaction.getReactionRate(inputReactants, temperature, pH=pH, potential=potential)
		backwardRate = self.backwardReaction.getReactionRate(inputReactants, temperature, pH=pH, potential=potential)
		return forwardRate-backwardRate





#DEPRECATED; Implementation forces use of the Euler integration method in effect; this means its often very unstable
class ReactionControllerStandard():

	def __init__(self, startReactants, reactions, timeStep, temperature=300, pH=0, potential=0, callbackFunct=None, constantConcReactants=None):
		""" DEPRECATED CLASS: DO NOT USE
		
		Args:
			startReactants: (iter of ChemSpeciesStd objects)
			reactions: (iter of ChemReactionTemplate objects)
			timeStep: (float) Timestep for updating concentrations
			temperature: (float) Temperature to run at; higher temperature means faster rates (enters barriers via RT)
			pH: (DEPRECATED) For the Taylor2016 model this sets the pH for the simulation. If using BetterReactionTemplate reactions then it does nothing (pH is handled in that case using the OH- and H+ concentrations/activities)
			potential: (float) Potential to run at. Usually all input reactions have barriers defined at potential=0; in which case +ve potential means polarise anodically while -ve potential means polarise cathodically
			callbackFunct: f(instance) Function called every step; useful for printing various intermediate states or logging every n-steps etc
			constantConcReactants: (iter of ChemSpeciesStd objects) These should ALSO be present in startReactants. After every step the concentration of these species is set to their initial value. This is my way of making sure our system stays in equilibrium with a reservoir containig set concentrations of non-surface species (e.g. [H+], [Mg2+])
	 
		Raises:
			Errors
		"""
		self.startReactants = startReactants
		self.reactions = reactions
		self.timeStep = timeStep
		self.temperature = temperature
		self.pH = pH
		self.potential = potential
#		self.currentReactants = copy.deepcopy(self.startReactants)
#		self.step = 0
		self.callbackFunct = callbackFunct
		self.constantConcReactants = constantConcReactants
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

		#Reset any concentrations we need to
		if self.constantConcReactants is not None:
			for constSpecies in self.constantConcReactants:
				for reactant in self.currentReactants:
					if reactant.name==constSpecies.name:
						reactant.conc = constSpecies.conc

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



