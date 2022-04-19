import math

from Stats.Unit7 import normalDistribution


def samplingDistribution(mean, std, size, side: str, *values):
    newStd = std / (math.sqrt(size))
    solution = normalDistribution(mean, newStd, side, *values)
    return solution


# sol = samplingDistribution(1.25, 0.05, 100, "below", 1.24, 1.26)
# print(sol)