import math

from Stats.Unit7 import findX

# from Stats.Unit8 import samplingDistribution


def confiInterval(sampleMean, n, std, threeSigma: str = "95", prob: float = 0, times: int = 1):
    rightHand = 0
    if threeSigma == "95":
        rightHand = 2 * (std / math.sqrt(n))
    elif threeSigma == "68":
        rightHand = (std / math.sqrt(n))
    elif threeSigma == "99.7":
        rightHand = 3 * (std / math.sqrt(n))
    elif threeSigma == "c":
        zStar = findCInterval(prob)
        rightHand = zStar * (std / math.sqrt(n))
    rightHand = (1 / math.sqrt(times)) * rightHand
    print(f"sampleMean ({sampleMean}) /pm {rightHand}")
    print(f"[{round(sampleMean - rightHand, 2)}, {round(sampleMean + rightHand, 2)}]")


def findCInterval(prob: float):
    area = 1 - prob
    indivArea = area / 2
    zStar = round(findX(0, 1, "below", indivArea) * -1, 3)
    print(f"sampleMean (x) /pm {zStar} * (sigma / sqrt(n))")
    return zStar


def findSampleMean(*args):
    n = len(args)
    solution = 0
    for i in args:
        solution = solution + i
    solution = solution / n
    return solution


def findMinSampleSize(prob, std, m):
    zStar = findCInterval(prob)
    n = (zStar * std / m) ** 2
    print(f"Minimum sample size (n): {round(n, 2)} approx= {round(n)}")
    return n


findMinSampleSize(0.98, 0.1, 0.04)
