import math
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt


def getMean(array, weight=None):
    sumItem = 0
    for item in array:
        sumItem = sumItem + item
    if weight is None:
        solution = sumItem / len(array)
    else:
        numerator = 0
        for item, wItem in array, weight:
            numerator = numerator + (item * wItem)
        solution = numerator / sumItem
    return solution


def quartiles(array):
    midVal = int(np.round((len(array)) / 2))
    leftHalf = array[0: midVal]
    rightHalf = array[midVal: -1]
    firstQuartile = np.median(leftHalf)
    thirdQuartile = np.median(rightHalf)
    return firstQuartile, thirdQuartile


def outlierBoxplot(array):
    q1, q3 = quartiles(array)
    IQR = q3 - q1
    lowerFence = q1 - (1.5 * IQR)
    upperFence = q3 + (1.5 * IQR)
    median = np.median(array)
    outliers = []
    for i in array:
        if i < lowerFence or i > upperFence:
            outliers.append(i)
    nonExtreme = []
    for i in array:
        if lowerFence < i < upperFence:
            nonExtreme.append(i)
    nonExtreme = np.array(nonExtreme)
    skew = stats.skew(nonExtreme)
    if skew == 0:
        print("Outlier graph skewness is normally distributed")
    elif skew > 1:
        print("Outlier graph Skewness is Right-Skewed")
    else:
        print("Outlier graph Skewness is Left-Skewed")
    print(f"Five number Summary: \nLower Fence: {lowerFence} "
          f"Q1: {q1} Median: {median} Q3: {q3} Upper Fence: {upperFence}")
    print(f"Outliers: {outliers}")
    plt.boxplot(array, vert=0)
    plt.show()


def variance(array):
    mean = getMean(array)
    nDeviations = []
    for i in array:
        nDeviations.append((i - mean) ** 2)
    sumSol = 0
    for u in nDeviations:
        sumSol = sumSol + u
    solution = sumSol / (len(array) - 1)
    print(f"Average Squared deviation: {solution}")
    print(f"Sample Deviation: {math.sqrt(solution)}")
    return solution


def fineNum(arrayOrig):
    array = np.array(arrayOrig)
    mode = stats.mode(array)
    median = np.median(array)

    mean = getMean(array)
    maxVal = array[-1]
    minVal = array[0]
    spread = maxVal - minVal
    q1, q3 = quartiles(arrayOrig)
    skew = stats.skew(array)
    if skew == 0:
        print("Original graph skewness is normally distributed")
    elif skew > 1:
        print("Original graph Skewness is Right-Skewed")
    else:
        print("Original graph Skewness is Left-Skewed")
    print(f"Mode: {mode[0]}")
    print(f"Spread (R): {spread}")
    print(f"Mean: {mean}")
    print(f"Five number Summary: \n Min:{minVal} Q1:{q1} Median:{median} Q3:{q3} Max:{maxVal}")
    plt.hist(arrayOrig, bins='auto', rwidth=0.85)
    plt.show()

