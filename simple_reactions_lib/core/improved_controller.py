
import copy

class ReactionControllerBase():

	def moveForwardByT(self, time):
		""" Updates the reactants concentration to what they would be at t+\Delta t. How we integrate forward is implementation dependent; and generally expected to be controlled by a propagator object
		
		Args:
			time: (float) The time we want to propagate forward to
				 
		"""
		raise NotImplementedError("")

	def reset(self):
		raise NotImplementedError("")

class ConcsPropagatorBase():
	""" Job of this class is to take a list of current concentrations, and propagate it forward by a timestep; Modfying concentrations IN PLACE. See .propagate() for interface """

	def propagate(self, inputReactants, timeStep, temperature=300, potential=0):
		""" Given a dictionary of reactants/concs at t=t_init this UPDATES IN PLACE concentrations up to t=t_init + timeStep
		
		Args:
			inputReactants: (dict) Keys are labels for chemical species, values are (ChemSpeciesStd) objects. Note the key should match ChemSpeciesStd.name
			timeStep: (float) Amount of time to move forward by
			temperature: (float) The temperature of the system in Kelvin
			potential: (float) Potential of the system relative to the reference
 
		"""

		raise NotImplementedError("")

class ConcChangesFinderBase():
	""" The job of this class is to figure out the changes in concentration for a maximum timestep. Stability reasons mean its not gauranteed to move forward the requested time step (extrapolating current rates over a long timestep can lead to unphysical concentrations)

	See getConcChangesForNextTimeStep for interface

	"""
	pass

	def getConcChangesForNextTimeStep(self, inputReactants, maxTimeStep, temperature=300, potential=0):
		""" Gets changes in concentration for reactants 
		
		Args:
			inputReactants: (dict) Keys are labels for chemical species, values are (ChemSpeciesStd) objects. Note the key should match ChemSpeciesStd.name
			maxTimeStep: (float) Maximum step forward to take. May take a smaller step for stability reasons
			temperature: (float) The temperature of the system in Kelvin
			potential: (float) Potential of the system relative to the reference			
	 
		Returns
			stepTaken: (float) The actual step taken forward (<=maxTimeStep)
			concChanges: Keys are the same as inputReactants keys. Values are floats representing d[X] for that reactant (i.e. the change in that reactants concentration)
	 
		"""
		raise NotImplementedError("") 


class RateCalculatorBase():
	""" Base class for calculating rates of change in reactants from their concentrations

	"""
	def getRates(self, inputReactants, temperature=300, potential=0):
		""" Gets the rates of change (d[X]/dt) for each chemical species 
		
		Args:
			inputReactants: (dict) Keys are labels for chemical species, values are (ChemSpeciesStd) objects. Note the key should match ChemSpeciesStd.name
			temperature: (float) The temperature of the system in Kelvin
			potential: (float) Potential of the system relative to the reference

		Returns
			outRates: (dict) Keys are the same as inputReactants keys. Values are floats representing d[X]/dt for that reactant
	 
		"""
		raise NotImplementedError("")



class ReactionControllerImproved(ReactionControllerBase):

	def __init__(self, propagator, startReactants, temperature=300, potential=0):
		""" Iniitalizer
		
		Args:
			propagator: (ConcsPropagatorBase) Class specialising in moving concentrations forward by \delta T. This also implicitly contains all info on reactions
			startReactants: (iter of ChemSpeciesStd objects) All the reactants and concentrations
			temperature: (float) The temperature in Kelvin
			potential: (float) Electric potential of the system in volts. Usually the zero value is defined under whatever conditions you have barriers/reaction energies calculated at 

		"""
		self.propagator = propagator
		self.startReactants = startReactants
		self.currentReactants = copy.deepcopy(startReactants)
		self.temperature = temperature
		self.potential = potential

	def reset(self):
		self.currentReactants = copy.deepcopy(self.startReactants)

	def moveForwardByT(self, time):
		self.propagator.propagate(self.currentReactants, time, temperature=self.temperature, potential=self.potential)



#TODO: Could likely just merge this with ConcChangesFinderStandard and add a .propagte to THAT class...
#Not sure theres ever going to be much varying configuration on this class?
class ConcsPropagatorStandard(ConcsPropagatorBase):

	def __init__(self, concChangesFinder, minConc=None, maxConc=None):
		""" Initializer
		
		Args:
			concChangesFinder: (ConcChangesFinderBase) object
				 
		"""
		self.concChangesFinder = concChangesFinder
		self.relTimeTolerance = 1e-4 #
#		self.minConc = minConc
#		self.maxConc = maxConc

	#There may be some float error issues here....
	def propagate(self, inputReactants, timeStep, temperature=300, potential=0):
		doneSteps = False
		timeTaken, maxStep = 0, timeStep


		while doneSteps is False:
			currStep, concChanges = self.concChangesFinder.getConcChangesForNextTimeStep(inputReactants, maxStep, temperature=temperature, potential=potential)
			timeTaken += currStep
			maxStep -= currStep
			self._updateConcs( inputReactants, concChanges )

			if ( abs((timeTaken/timeStep)-1)< self.relTimeTolerance):
				doneSteps = True
			elif abs(timeTaken/timeStep) > 1:
				doneSteps = True
		

	def _updateConcs(self, reactants, concChanges):
		for key in concChanges.keys():
			for reactant in reactants:
				if reactant.name==key:
					reactant.conc += concChanges[key]
#					#Check conc has remained within limits if requested
#					if (self.minConc is not None):
#						if reactant.conc<self.minConc:
#							reactant.conc = self.minConc
#					if (self.maxConc is not None):
#						if reactant.conc>self.


class ConcChangesFinderStandard(ConcChangesFinderBase):

	def __init__(self, rateCalculator, variableConcSpecies, maxConcChange=0.1):
		""" Iniitalizer
		
		Args:
			rateCalculator: (RateCalculatorBase)
			variableConcSpecies: (iter of str) Names of species for which concentration is allowed to vary.
			maxConcChange: (float) The maximum change in concentrations of a "variableConcSpecies" for a timestep. A maximum value is needed for stability reasons; In the extreme if a conc change is >1 it means we'll get either negative concentration or fractional concentrations > 1 for a species.

		NOTE:
			I THINK the value of maxConcChange will basically determine how accurate a [conc] vs [time] plot is. In effect it tells us how often to sample the derivative along the true curve (as maxConcChange->0 we approach the "real" curve)

				 
		"""
		self.rateCalculator = rateCalculator
		self.variableConcSpecies = variableConcSpecies
		self.maxConcChange = maxConcChange

	def getConcChangesForNextTimeStep(self, inputReactants, maxTimeStep, temperature=300, potential=0):

		#Figure out the largest change in reactant concentrations
		rates = self.rateCalculator.getRates(inputReactants, temperature=temperature, potential=potential)
		maxChange = 0
		for key in rates.keys():
			if key in self.variableConcSpecies:
				currChange = abs(rates[key])
				if currChange > maxChange:
					maxChange = abs(rates[key])

		#Figure out the step we use
		try:
			maxPossibleTimeStep = self.maxConcChange / maxChange
		except ZeroDivisionError:
			maxPossibleTimeStep = maxTimeStep
		outStep = min(maxTimeStep, maxPossibleTimeStep)

		#Figure out the changes in input reactants accordingly
		outVals = dict()
		for key in self.variableConcSpecies:
			outVals[key] = rates[key]*outStep


		return outStep, outVals

class RateCalculatorStandard(RateCalculatorBase):

	def __init__(self, reactions):
		""" Initializer
		
		Args:
			reactions: (iter of ChemReactionBase objects) 
				 
		"""
		self.reactions = reactions


	def getRates(self, inputReactants, temperature=300, potential=0):
		outDict = dict()
		dudTimeStep = 1 #We want the rate of change rather than the amount; hence we ALWAYS use a timestep of 1
		for reaction in self.reactions:
			currChanges = reaction.getChangesInReactants(dudTimeStep, inputReactants, temperature, potential=potential)
			for key in currChanges.keys():
				if key in outDict:
					outDict[key] += currChanges[key]
				else:
					outDict[key]  = currChanges[key]
		return outDict


