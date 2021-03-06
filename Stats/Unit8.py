import math

from Stats.Unit7 import normalDistribution, findX, findBounds


def samplingDistribution(mean, std, size, side: str, *values):
    newStd = std / (math.sqrt(size))
    # solution = normalDistribution(mean, newStd, side, *values)
    # solution = findX(mean, newStd, side, values)
    # Find bounds for this
    solution = findBounds(mean, newStd, values[0])
    return solution