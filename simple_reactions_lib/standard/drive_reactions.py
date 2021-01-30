

def runControllerUntilSteadyState(controller, steadySpeciesNames, nStepsPerCheck, tolerancePerStep):
	""" Run the reaction controller until the specified reactants reach a steady state; determined by their change per step
	
	Args:
		controller: (ReactionControllerStandard object) 
		steadySpeciesNames: (iter of str) The species we want to have reached a steady state
		nStepsPerCheck: (int) Check for convergence every N-steps. Too-small numbers probably slow things down (though i dont know)
		tolerancePerStep: (float) Amount of changes in reaction 
 
	Returns
		Nothing, but controller will have been run N times and contain concentrations/rates etc at the steady state
 
	"""
	while _isCurrentStepConverged(controller, steadySpeciesNames, tolerancePerStep) is False:
		controller.doNextNSteps(nStepsPerCheck)


def _isCurrentStepConverged(controller, steadySpeciesName, tolerancePerStep):
	allConcDiffs = controller._getConcDiffs()
	maxConcDiff = max( [abs(allConcDiffs[key]) for key in steadySpeciesName] )

	if maxConcDiff < tolerancePerStep:
		return True
	else:
		return False

	




