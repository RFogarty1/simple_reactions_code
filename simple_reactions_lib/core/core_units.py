
import math

BOLTZ_EV = 8.617333262e-5

FARADAY_CONST = 96485.33212
IDEAL_GAS_R_JOULES = 8.31446261815324


NERNST_PREFACTOR_LOG10 = (math.log(10)*IDEAL_GAS_R_JOULES) / FARADAY_CONST
