import math
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import stemgraphic


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


def quartiles(array, missing=0):
    midVal = int(np.round((len(array)+missing) / 2))
    leftHalf = array[0: midVal]
    rightHalf = array[midVal: -1]
    firstQuartile = np.median(leftHalf)
    thirdQuartile = np.median(rightHalf)
    return firstQuartile, thirdQuartile


def outlierBoxplot(array, missing=0, plot=False):
    q1, q3 = quartiles(array,missing)
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
    print(f"IQR: {IQR}")
    if plot:
        plt.boxplot(array, vert=0)
        plt.show()


def variance(array, printInfo=False):
    mean = getMean(array)
    nDeviations = []
    for i in array:
        nDeviations.append((i - mean) ** 2)
    sumSol = 0
    for u in nDeviations:
        sumSol = sumSol + u
    solution = sumSol / (len(array) - 1)
    if printInfo:
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


def getMedian(arrayOrig):
    return np.median(arrayOrig)


def specificRelativeFreq(array, lowerB, higherB, side):
    sumVal = 0
    if side == "right":
        for i in array:
            if lowerB <= i < higherB:
                sumVal = sumVal + 1
    elif side == "left":
        for i in array:
            if lowerB < i <= higherB:
                sumVal = sumVal + 1
    return sumVal / len(array)


def plotStemLeaf(array):
    prev = float("-inf")
    mySorted = True
    for i in array:
        if prev > i:
            mySorted = False
            break
        prev = i
    if not mySorted:
        mergeSort(array, 0, len(array))
    prev = 0
    for i in array:
        if prev != str(i)[0]:
            print("")
            prev = str(i)[0]
            print(f"{str(i)[0]} |", end=" ")
        print(str(i)[1:], end=" ")




def merge(arr, l, m, r):
    n1 = m - l + 1
    n2 = r - m

    # create temp arrays
    L = [0] * n1
    R = [0] * n2

    # Copy data to temp arrays L[] and R[]
    for i in range(0, n1):
        L[i] = arr[l + i]

    for j in range(0, n2):
        R[j] = arr[m + 1 + j]

    # Merge the temp arrays back into arr[l..r]
    i = 0  # Initial index of first subarray
    j = 0  # Initial index of second subarray
    k = l  # Initial index of merged subarray

    while i < n1 and j < n2:
        if L[i] <= R[j]:
            arr[k] = L[i]
            i += 1
        else:
            arr[k] = R[j]
            j += 1
        k += 1

    # Copy the remaining elements of L[], if there
    # are any
    while i < n1:
        arr[k] = L[i]
        i += 1
        k += 1

    # Copy the remaining elements of R[], if there
    # are any
    while j < n2:
        arr[k] = R[j]
        j += 1
        k += 1


# l is for left index and r is right index of the
# sub-array of arr to be sorted


def mergeSort(arr, l, r):
    if l < r:
        # Same as (l+r)//2, but avoids overflow for
        # large l and h
        m = l + (r - l) // 2

        # Sort first and second halves
        mergeSort(arr, l, m)
        mergeSort(arr, m + 1, r)
        merge(arr, l, m, r)
