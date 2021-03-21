
import scipy.integrate as integrateHelp

from . import improved_controller as contrHelp


class ConcsPropagator_DOP853(contrHelp.ConcsPropagatorTemplate):

	def __init__(self, rateCalculator, variableConcSpecies, aTol=None, rTol=None):
		self.rateCalculator = rateCalculator
		self.variableConcSpecies = variableConcSpecies
		self.aTol = aTol
		self.rTol = rTol

	def _propagateVectorisedFunctionToNextTimeStep(self, startConcs, timeStep, vectorisedFunction):

		#Run the integrator until it reaches the desired timestep
		t0,tEnd = 0,timeStep
		y0 = startConcs
		solverOptions = {}
		if self.aTol is not None:
			solverOptions["atol"] = self.aTol
		if self.rTol is not None:
			solverOptions["rtol"] = self.rTol

		outObj = integrateHelp.solve_ivp(vectorisedFunction, [t0,tEnd], y0, method="DOP853", **solverOptions)
		outVals = [y[-1] for y in outObj.y]

		return outVals


class ConcsPropagator_Radau(contrHelp.ConcsPropagatorTemplate):

	def _propagateVectorisedFunctionToNextTimeStep(self, startConcs, timeStep, vectorisedFunction):
		t0,tEnd = 0,timeStep
		outObj = integrateHelp.solve_ivp(vectorisedFunction, [t0,tEnd], startConcs, method="Radau")
		outVals = [ y[-1] for y in outObj.y ]
		assert (abs(outObj.t[-1]-timeStep)/timeStep)<0.01

		return outVals

