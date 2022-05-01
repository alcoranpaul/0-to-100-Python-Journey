import numpy as np
import pylab as pl
import scipy.stats as stats
import matplotlib.pyplot as plt
from varname import nameof
from Unit1 import *
import seaborn as sb


def showScatterPlot(explanatory, response, showRegression=False):
    correlation = getCorrelation(explanatory, response)
    if correlation > 0:
        print("Association: Positive")
    else:
        print("Association: Negative")

    plt.scatter(explanatory, response)
    if showRegression:
        plt.plot(explanatory, showLeastSquareRegression(explanatory, response, correlation))
    plt.xlabel(nameof(explanatory))
    plt.ylabel(nameof(response))
    # plt.plot(showLeastSquareRegression(explanatory, response, correlation))
    plt.show()


def getCorrelation(explanatory, response):
    xSTD = (math.sqrt(variance(explanatory)))
    ySTD = (math.sqrt(variance(response)))
    xDeviations = calculateDeviaions(explanatory)
    yDeviations = calculateDeviaions(response)
    nProducts = 0
    for i in range(len(explanatory)):
        nProducts = nProducts + (xDeviations[i] * yDeviations[i])
    return round(nProducts / ((len(explanatory) - 1) * xSTD * ySTD), 4)


def calculateDeviaions(array):
    mean = getMean(array)
    deviations = []
    for i in array:
        deviations.append(i - mean)
    return deviations


def showLeastSquareRegression(explanatory, response, correlation):
    xSTD = (math.sqrt(variance(explanatory)))
    ySTD = (math.sqrt(variance(response)))
    b_1 = correlation * (ySTD / xSTD)
    print(b_1)
    b_0 = getMean(response) - (b_1 * getMean(explanatory))
    print(b_0)
    yReg = []
    for i in explanatory:
        yReg.append(b_0 + (b_1 * i))
    return yReg


def leastSquareRegFormula(explanatory, response, xValue):
    correlation = getCorrelation(explanatory, response)
    xSTD = (math.sqrt(variance(explanatory)))
    ySTD = (math.sqrt(variance(response)))
    b_1 = correlation * (ySTD / xSTD)
    b_0 = getMean(response) - (b_1 * getMean(explanatory))
    yReg = b_0 + (b_1 * xValue)
    return yReg


def getResidual(explanatory, response, xValue):
    yHat = leastSquareRegFormula(explanatory, response, xValue)
    y = np.interp(xValue, explanatory, response)
    return round(y - yHat, 4)


# if __name__ == '__main__':
#     xc = [5.5, 6.1, 6.7, 7.0, 7.5, 7.9, 8.6, 9.2, 9.8, 10.3]
#     yc = [0.11, 0.19, 0.24, 0.37, 0.36, 0.49, 0.59, 0.60, 0.81, 0.78]
#     print(getResidual(xc, yc, 7.0))
