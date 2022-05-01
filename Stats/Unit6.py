from decimal import Decimal
from fractions import Fraction


def pmf(func, *xValues):
    sumVal = 0
    for i in xValues:
        x = func(i)
        sumVal = sumVal + x
    c = 1 / sumVal
    return c


