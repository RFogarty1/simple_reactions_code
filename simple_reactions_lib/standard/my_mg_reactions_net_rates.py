

import math
from ..core import core_classes as coreObjs
from ..core import core_units as unitHelp


class StandardNetReactionTemplate(coreObjs.NetReactionTemplate):

	def __init__(self, forwardBarrier, prefactor, reactionEnergy, u0=0, symFactorForward=0.5):
		""" Initializer
		
		Args:
			forwardBarrier: (float) Forward reaction barrier in eV
			prefactor: (float) The prefactor in 
			reactionEnergy: (float) The energy to go from reactants to products. This should come from the same calculation as the barrier
			u0: (float) Potential at which barrier/reactionEnergy are calculated. Should generally be left to zero and only matters if electron transfer number != 0
			symFactorForward: (float) Value between 0 and 1 which describes how well the potential-dependence of transition state relates to reactant(0 means reacts the same)  and product (1 means the same as product)
				 
		"""
		self._createForwardReaction(forwardBarrier, prefactor, u0, symFactorForward)
		self._createBackwardReaction(forwardBarrier, prefactor, reactionEnergy, u0, symFactorForward)
	
	def _createForwardReaction(self, barrier, prefactor,u0, symFactorForward):
		args = [self.reactants, self.products, barrier, prefactor]
		kwargs = {"refPot":u0, "symFactor":symFactorForward, "nElecTransfer":self.nElecTransfer}
		self.forwardReaction = coreObjs.BetterReactionTemplate(*args, **kwargs)
	
	def _createBackwardReaction(self, forwardBarrier, prefactor, reactionEnergy, u0, symFactorForward):
		barrier = forwardBarrier - reactionEnergy
		args = [self.products, self.reactants, barrier, prefactor]
		kwargs = {"refPot":u0, "symFactor":1-symFactorForward, "nElecTransfer":-1*self.nElecTransfer}
		self.backwardReaction = coreObjs.BetterReactionTemplate(*args, **kwargs)

	@property
	def reactants(self):
		raise NotImplementedError("")

	@property
	def products(self):
		raise NotImplementedError("")

	@property
	def nElecTransfer(self):
		raise NotImplementedError("")


class TafelReactionNet(StandardNetReactionTemplate):
	
	@property
	def reactants(self):
		return ["h_ads", "h_ads"]

	@property
	def products(self):
		return ["free", "free", "h2"] 

	@property
	def nElecTransfer(self):
		return 0

class WaterAssistReaction_twoElectronXferNet(StandardNetReactionTemplate):

	@property
	def reactants(self):
		return ["free"]

	@property
	def products(self):
		return ["mg2+", "free"]

	@property
	def nElecTransfer(self):
		return 2


class OHAssistedDissolutionReaction_twoElectronXferNet(StandardNetReactionTemplate):

	@property
	def reactants(self):
		return ["oh_ads"]

	@property
	def products(self):
		return ["free", "mg2+"]

	@property
	def nElecTransfer(self):
		return 2


class VolmerReactionNet(StandardNetReactionTemplate):

	@property
	def reactants(self):
		return ["free", "free"]

	@property
	def products(self):
		return ["h_ads", "oh_ads"]

	@property
	def nElecTransfer(self):
		return 0


class Heyrovsky_waterAssistedNet(StandardNetReactionTemplate):

	@property
	def reactants(self):
		return ["h_ads"]

	@property
	def products(self):
		return ["free", "oh-", "h2"]

	@property
	def nElecTransfer(self):
		return -1


class HydrogenBulkDiffusionNet(StandardNetReactionTemplate):

	@property
	def reactants(self):
		return ["h_ads"]

	@property
	def products(self):
		return ["free", "h_diffused"]

	@property
	def nElecTransfer(self):
		return 0


class CathodicOHDesorption(StandardNetReactionTemplate):

	@property
	def reactants(self):
		return ["oh_ads"]

	@property
	def products(self):
		return ["free", "oh-"]

	@property
	def nElecTransfer(self):
		return -1



