
import math
from ..core import core_classes as coreObjs
from ..core import core_units as unitHelp


class TafelReaction(coreObjs.ChemReactionTemplate):

	def __init__(self, barrier, prefactor):
		self.barrier = barrier
		self.prefactor = prefactor

		self.reactants = ["h_ads","h_ads"]
		self.products = ["free","free","h2"]

	def _getTafelFactor(self, inputReactants, temperature, pH, potential):
		return 1


class VolmerReaction(coreObjs.ChemReactionTemplate):

	def __init__(self, barrier, prefactor):
		self.barrier = barrier
		self.prefactor = prefactor

		self.reactants = ["free", "free"]
		self.products = ["oh_ads","h_ads"]

	def _getTafelFactor(self, inputReactants, temperature, pH, potential):
		return 1


#SUSPECT i have these numbers wrong aswell and am using the wrong Nernst equation in essence
class CathodicOHDesorption(coreObjs.ChemReactionTemplate):

	def __init__(self, barrier, prefactor, symFactor=0.5, u0=-0.83):
		self.barrier = barrier
		self.prefactor = prefactor
		self.symFactor = symFactor

		self.reactants = ["oh_ads"]
		self.products = ["free", "oh_minus"]

		self.electronsLost = -1
		self.u0 = u0

	#TODO: I MAY be able to factor this out? Except tricky for non-electron transfers (return 1 if self.electronsLost=0 else calcVal)
	def _getTafelFactor(self, inputReactants, temperature, pH, potential):
		u0Val = self._getModifiedU0Value(inputReactants, temperature, pH)
		expTerm = (-1*self.symFactor*unitHelp.FARADAY_CONST*(potential-u0Val)) / (unitHelp.IDEAL_GAS_R_JOULES*temperature)
		return math.exp(expTerm)
		
	def _getModifiedU0Value(self, inputReactants, temperature, pH):
		prefactor = -1 * unitHelp.NERNST_PREFACTOR_LOG10 * (temperature/self.electronsLost)

		#Figure out the lnQ part of the Nernst equation
		protonConc = 10**(-1*pH)
		oHConc = 1e-14/protonConc
		logPart = math.log10(oHConc)

		return self.u0 - (prefactor*logPart)

#TODO: This is PROBABLY wrong, my assumption of pH=0 at the standard U0 is probably wrong
#Also - are 1 or 2 electrons lost in an electron transfer? Pretty sure its two for Nernst equation, but only 1 for the potential
#TODO: Look at the H reduction reactions here
#https://chem.libretexts.org/Bookshelves/General_Chemistry/Book%3A_Chem1_(Lower)/16%3A_Electrochemistry/16.04%3A_The_Nernst_Equation


#I'm REASONABLY sure that Taylor sees the OH as just a ligand on [MgOH]+
# This means he REALLY is assuming a 1 electron loss for the OH assisted pathway. Which is dumb AF but w/e
class OHAssistedDissolution_Taylor2016(coreObjs.ChemReactionTemplate):

	def __init__(self, barrier, prefactor, symFactor=0.5, u0=-2.38):
		self.barrier = barrier
		self.prefactor = prefactor
		self.symFactor = symFactor

		self.reactants = ["oh_ads"]
		self.products = ["free", "mg_2+", "oh_minus", "electron"]
#		self.electronsLost = 2
		self.electronsLost = 1 #TODO: THIS REQUIRES SOME ACTUAL THOUGHT

		self.u0 = u0

	def _getTafelFactor(self, inputReactants, temperature, pH, potential):
		u0Val = self._getModifiedU0Value(inputReactants, temperature, pH)
#		expTerm = (-1*self.symFactor*unitHelp.FARADAY_CONST*(potential-u0Val)) / (unitHelp.IDEAL_GAS_R_JOULES*temperature)

		expTerm = (1*self.symFactor*unitHelp.FARADAY_CONST*(potential-u0Val)) / (unitHelp.IDEAL_GAS_R_JOULES*temperature)


		return math.exp(expTerm)

	#TODO: THIS REQUIRES SERIOUS THOUGHT
	def _getModifiedU0Value(self, inputReactants, temperature, pH):
		prefactor = -1 * unitHelp.NERNST_PREFACTOR_LOG10 * (temperature/self.electronsLost)

		#We keep the constant mg2+ concentration with inputReactants, separate from "mg_2+" that gets generated in the simulation
		for reactant in inputReactants:
			if reactant.name=="fixed_mg_2+":
				mg2PlusConc = reactant.conc

		#Figure out the lnQ part of the nernst equation
		protonConc = 10**(-1*pH)
		oHConc = 1e-14/protonConc
#		oHTerm = math.log10(oHConc/(1e-14)) #At "standard" conditions [OH]=1e-14, so I THINK we need to divide [OH] by that
		oHTerm = math.log10(oHConc) #Probably also wrong and stupid
		oHTerm = 0
		mg2PlusTerm = math.log10(mg2PlusConc)

		logPart = oHTerm + mg2PlusTerm		
		outU0 = self.u0
		return self.u0 - (prefactor*logPart)





