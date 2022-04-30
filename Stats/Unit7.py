import math
from math import comb
from decimal import Decimal

# from scipy.stats import poisson as poisson2
# DISTRIBUTIONS
from scipy.stats import norm


def binomial(x, n, p, info=False):
    """
    Binomial

    :param info: Show information
    :param x: number of successes seeking
    :param n: total population
    :param p: probability of success
    :return:
    """
    comp = {'probability': (comb(n, x)) * (p ** x) * ((1 - p) ** (n - x)),
            'expected': n * p,
            'variance': (n * p) * (1 - p)}
    comp['STD'] = math.sqrt(comp['variance'])
    if info:
        print(f"P(X = {x}) = {comp['probability']:.4f}")
        print(f"E(X) = {n * p:.4f}")
        print(f"E(x^2) = {p:.4f}")
        print(f"Var(X) = {comp['variance']:.4f}")
        print(f"StandardD = {comp['STD']:.4f}")
    return comp


def geometric(x, p, info=False):
    """
    Geometric

    :param info: Show information
    :param x: index of first success
    :param p: probability of success
    :return:
    """
    comp = {'probability': ((1 - p) ** (x - 1)) * p,
            'cdf': 1 - ((1 - p) ** x),
            'variance': (1 - p) / (p ** 2),
            'expected': 1 / p,
            'expectedSquared': (2 - p) / (p ** 2)}
    comp['STD'] = math.sqrt(comp['variance'])
    if info:
        print(f"P(X = {x}) = {comp['probability']:.4f}")
        print(f"CDF: F(X <= {x}) = {comp['cdf']:.4f} ")
        print(f"E(X) = {comp['expected']:.4f}")
        print(f"E(x^2) = {comp['expectedSquared']:.4f}")
        print(f"Var(X) = {comp['variance']:.4f}")
        print(f"StandardD = {comp['STD']:.4f}")
    return comp


def negativeBinomial(x, r, p, info=False):
    """
    Negative Binomial

    :param info: Show information
    :param x: number of trials
    :param r: number of successes on x th trial
    :param p: probability of success
    :return:
    """
    comp = {'probability': comb(x - 1, r - 1) * (1 - p) ** (x - r) * p ** r,
            'expected': r / p,
            'variance': (r * (1 - p)) / (p ** 2)}
    comp['STD'] = math.sqrt(comp['variance'])
    if info:
        print(f"P(X = {x}) = {comp['probability']:.4f}")
        print(f"E(X) = {comp['expected']:.4f}")
        print(f"Var(X) = {comp['variance']:.4f}")
        print(f"StandardD = {comp['STD']:.4f}")
    return comp


def bernoulliTrials(x, p, info=False):
    """
    bernoulliTrials

    :param info: Show information
    :param x:
    :param p:
    :return:
    """
    comp = {
        'probability': (p ** x) * ((1 - p) ** (1 - x)),
        'expected': p,
        'expectedSquared': p,
        'expression_variance': p * (1 - p),
    }
    comp['STD'] = math.sqrt(comp['variance'])
    if info:
        print(f"P(X = {x}) = {comp['probability']:.4f}")
        print(f"E(X) = {comp['expected']:.4f}")
        print(f"E(x^2) = {comp['expectedSquared']:.4f}")
        print(f"Var(X) = {comp['variance']:.4f}")
        print(f"StandardD = {comp['STD']:.4f}")
    return comp


def hypergeometric(N, n, r, info=False, x=-1):
    """
    Hyper-geometric

    :param info: Show information
    :param x: number of successes seeking of n items
    :param N: total population size
    :param n: number of selections made from N
    :param r: total number of successes.deesired items in N
    """
    comp = {
        'expected': (n * r) / N,
        'variance': ((N - n) / (N - 1)) * n * (r / N) * (1 - (r / N)), }
    if x >= 0:
        comp['probability'] = (comb(r, x) * comb(N - r, n - x)) / comb(N, n)
    comp['STD'] = math.sqrt(comp['variance'])
    if info:
        print(f"P(X = {x}) = {comp['probability']:.4f}")
        print(f"E(X) = {comp['expected']:.4f}")
        print(f"Var(X) = {comp['variance']:.4f}")
        print(f"StandardD = {comp['STD']:.4f}")
    return comp


# Week 8


def poisson(lmbda, x=0, t=1, info=False):
    """
    poisson
    Number of events in a time interval
    Discrete: Exact value for x

    :param lmbda: Mean number of successses that occur during a specific interval
    :param x: Number of successes
    :param t: time
    :param info: Show information
    :return:
    """
    comp = {}
    lambda_value = lmbda * t
    # print(lambda_value)
    if x >= 0:
        try:
            comp['probability'] = (((lambda_value ** x) * (math.e ** (-lambda_value))) / math.factorial(x))
        except OverflowError:
            comp['probability'] = (((lambda_value ** x) * Decimal(math.e ** (-lambda_value))) / math.factorial(x))
        if info:
            print(f"P(X = {x}) = {comp['probability']:.4f}")

    comp['expected'] = lambda_value
    comp['variance'] = lambda_value
    comp['STD'] = math.sqrt(comp['variance'])
    if info:
        print(f"E(X) = {comp['expected']:.4f}")
        print(f"Var(X) = {comp['variance']:.4f}")
        print(f"StandardD = {comp['STD']:.4f}")
    return comp


def uniform(a, b, x=1, info=False):
    """
    Uniform

    :param a: Lower bound
    :param b: Higher bound
    :param x: interested value
    :param info: Show information
    :return:
    """
    comp = {'pdf': 1 / (b - a)}
    if x < a:
        comp['cdf'] = 0
    elif x > b:
        comp['cdf'] = 1
    else:
        comp['cdf'] = (x - a) / (b - a)
    comp['expected'] = (a + b) / 2
    comp['expectedSquared'] = ((a ** 2) + (a * b) + (b ** 2)) / 3
    comp['variance'] = ((b - a) ** 2) / 12
    comp['STD'] = math.sqrt(comp['variance'])
    if info:
        print(f"PDF: f(x) = {comp['pdf']} for [{a}, {b}]")
        print(f"CDF: F(x) = {comp['cdf']} for {a} <= {x} <= {b}")
        print(f"E(X) = {comp['expected']:.4f}")
        print(f"E(x^2) = {comp['expectedSquared']:.4f}")
        print(f"Var(X) = {comp['variance']:.4f}")
        print(f"StandardD = {comp['STD']:.4f}")
    return comp


def exponential(lmbda, x, proba=0, info=False):
    """
    exponential
    Time between two eevents
    Continuous on an interval

    :param lmbda: Events rate
    :param x: Interested Value
    :param proba: Probability
    :param info: Show information
    :return:
    """
    # Untested
    comp = {}
    if x == "?":
        comp['x'] = (-math.log((1 - proba))) / lmbda
        if info:
            print(f"Median or    x = {x}")
    else:
        if x < 0:
            comp['cdf'] = 0
            comp['pdf'] = 0
        else:
            comp['pdf'] = (math.e ** (-lmbda * x))
            comp['cdf'] = 1 - (math.e ** (-lmbda * x))
        if info:
            print(f"P(X > {x}) = PDF: P(X = {x}) = f(x) = {comp['pdf']}")
            print(f"P(X < {x}) = CDF: F(x) = {comp['cdf']} for 0 to {x}")

    comp['expected'] = 1 / lmbda
    comp['variance'] = 1 / (lmbda ** 2)
    comp['STD'] = math.sqrt(comp['variance'])

    if info:
        print(f"E(X) = {comp['expected']:.4f}")
        print(f"Var(X) = {comp['variance']:.4f}")
        print(f"StandardD = {comp['STD']:.4f}")
    return comp


def add(*args):
    sumValue = 0
    for i in args:
        sumValue = sumValue + i
    return sumValue


def addHyperGeometric(N, n, r, num):
    sumValue = 0
    for i in range(num, r + 1):
        sumValue = sumValue + hypergeometric(N, n, r, x=i).get('probability')
    return sumValue


def solve(func, **kwargs):
    solution = None
    components = {}
    if 'info' not in kwargs:
        kwargs['info'] = False
    if func == "binomial":
        components = binomial(kwargs['x'], kwargs['n'], kwargs['p'], kwargs['info'])
    elif func == "geometric":
        components = geometric(kwargs['x'], kwargs['p'], kwargs['info'])
    elif func == "negativeBinomial":
        components = negativeBinomial(kwargs['x'], kwargs['r'], kwargs['p'], kwargs['info'])
    elif func == "hypergeometric":
        if 'x' not in kwargs:
            kwargs['x'] = -1
        components = hypergeometric(kwargs['N'], kwargs['n'], kwargs['r'],
                                    kwargs['info'], kwargs['x'])
    elif func == "bernoulliTrials":
        components = bernoulliTrials(kwargs['x'], kwargs['p'], kwargs['info'])
    elif func == "poisson":
        if 't' not in kwargs:
            kwargs['t'] = 1
        if 'x' not in kwargs:
            kwargs['x'] = 0
        components = poisson(kwargs['lmbda'], kwargs['x'], kwargs['t'], kwargs['info'])
    elif func == "uniform":
        if 'x' not in kwargs:
            kwargs['x'] = 1
        components = uniform(kwargs['a'], kwargs['b'], kwargs['x'], kwargs['info'])
    elif func == "exponential":
        if 'proba' not in kwargs:
            kwargs['proba'] = 0
        components = exponential(kwargs['lmbda'], kwargs['x'], kwargs['proba'], kwargs['info'])
    print("=================================")
    print("Solution: ")
    if 'comp' not in kwargs:
        for i in components.keys():
            print(f"{i}: {components[i]}")
    elif kwargs['comp'] == "probability":
        if "probability" not in components:
            print(f"    There is no such thing as '{kwargs['comp']}' in the components!!")
            probs = input("    You might be interested in Pdf or Cdf, shall print them? ")
            if probs == "yes":
                print(f"P(X > {kwargs['x']}) = PDF: {round(components.get('pdf'), 4)}")
                print(f"P(X < {kwargs['x']}) = CDF: {round(components.get('cdf'), 4)}")
        else:
            solution = round(components.get('probability'), 4)
    elif kwargs['comp'] == "pdf":
        solution = round(components.get('pdf'), 4)
    elif kwargs['comp'] == "cdf":
        solution = round(components.get('cdf'), 4)
    elif kwargs['comp'] == "expected":
        solution = round(components.get('expected'), 4)
    elif kwargs['comp'] == "expectedSquared":
        solution = round(components.get('expectedSquared'), 4)
    elif kwargs['comp'] == "variance":
        solution = round(components.get('variance'), 4)
    elif kwargs['comp'] == "STD":
        solution = round(components.get('STD'), 4)
    elif kwargs['comp'] == "all":
        for i in components.keys():
            print(f"{i}: {components[i]}")
    if solution is None:
        solution = "All information is printed Above"
    return solution


def distributionSolve():
    function = (input("What kind of Distribution? "))
    if function == "binomial":
        xb = int(input("(x) number of successes seeking: "))
        nb = int(input("(n) total population: "))
        pb = input("(p) probability of success: ")
        if len(pb.split()) == 2:
            if pb.split()[0] == "sqrt":
                pb = math.sqrt(float(pb.split()[1]))
        else:
            pb = float(pb)
        infob = input("     Do you want to print all? ").lower()
        if infob == "no":
            compp = input("     What information do you need? ")
            print(" ")
            print(solve(function, x=xb, n=nb, p=pb, comp=compp))
        else:
            print(" ")
            print(solve(function, x=xb, n=nb, p=pb, info=infob))
    elif function == "geometric":
        xg = int(input("(x) index of first success: "))
        pg = input("(p) probability of success: ")
        if len(pg.split()) == 2:
            if pg.split()[0] == "sqrt":
                pg = math.sqrt(float(pg.split()[1]))
        else:
            pg = float(pg)
        infog = input("     Do you want to print all? ").lower()
        if infog == "no":
            compg = input("     What information do you need? ")
            print(" ")
            print(solve(function, x=xg, p=pg, comp=compg))
        else:
            print(" ")
            print(solve(function, x=xg, p=pg, info=infog))

    elif function == "negativeBinomial":
        xn = int(input("(x) number of successes seeking: "))
        rn = int(input("(r) number of successes on x th trial: "))
        pn = input("(p) probability of success: ")
        if len(pn.split()) == 2:
            if pn.split()[0] == "sqrt":
                pn = math.sqrt(float(pn.split()[1]))
        else:
            pn = float(pn)
        infon = input("     Do you want to print all? ").lower()
        if infon == "no":
            compn = input("     What information do you need? ")
            print(" ")
            print(solve(function, x=xn, r=rn, p=pn, comp=compn))
        else:
            print(" ")
            print(solve(function, x=xn, r=rn, p=pn, info=infon))
    elif function == "hypergeometric":
        Nh = int(input("(N) total population size: "))
        nh = int(input("(n) number of selections made: "))
        rh = int(input("(r) total number of successes: "))
        infoh = input("     Do you want to print all? ").lower()
        isx = input("   Is there an x value? ")
        if isx == "yes":
            xh = int(input("(x) Total number of successes seeking "))
            if infoh == "no":
                comph = input("     What information do you need? ")
                print(" ")
                print(solve(function, N=Nh, n=nh, r=rh, comp=comph, x=xh))
            print(" ")
            print(solve(function, N=Nh, n=nh, r=rh, info=infoh, x=xh))
        else:
            print(" ")
            print(solve(function, N=Nh, n=nh, r=rh, info=infoh))
            if infoh == "no":
                comph = input(" What information do you need? ")
                print(" ")
                print(solve(function, N=Nh, n=nh, r=rh, comp=comph))
            print(" ")
            print(solve(function, N=Nh, n=nh, r=rh, info=infoh))
    elif function == "bernoulliTrials":
        xbe = int(input("(x) index of first success: "))
        pbe = input("(p) probability of success: ")
        if len(pbe.split()) == 2:
            if pbe.split()[0] == "sqrt":
                pbe = math.sqrt(float(pbe.split()[1]))
        else:
            pbe = float(pbe)
        infobe = input("    Do you want to print all? ").lower()
        if infobe == "no":
            compbe = input("    What information do you need? ")
            print(" ")
            print(solve(function, x=xbe, p=pbe, comp=compbe))
        else:
            print(" ")
            print(solve(function, x=xbe, p=pbe, info=infobe))

    elif function == "poisson":
        lmbdap = input("(lambda) Events rate: ")
        if len(lmbdap.split()) == 2:
            if lmbdap.split()[0] == "sqrt":
                lmbdap = math.sqrt(float(lmbdap.split()[1]))
        else:
            lmbdap = float(lmbdap)
        isx = input("   Is there a (x) value? ").lower()
        ist = input("   Is there a (t) value? ").lower()
        infop = input("     Do you want to print all? ").lower()
        if isx == "yes" and ist == "yes":
            xp = int(input("(x) Number of events: "))
            tp = float(input("(t) Time: "))
            if infop == "no":
                compp = input("     What information do you need? ")
                print(" ")
                print(solve(function, lmbda=lmbdap, x=xp, t=tp, comp=compp))
            else:
                print(" ")
                print(solve(function, lmbda=lmbdap, x=xp, t=tp, info=infop))
        elif isx == "no" and ist == "yes":
            tp = int(input("(t) Time: "))
            if infop == "no":
                compp = input("     What information do you need? ")
                print(" ")
                print(solve(function, lmbda=lmbdap, t=tp, comp=compp))
            else:
                print(" ")
                print(solve(function, lmbda=lmbdap, t=tp, info=infop))
        elif isx == "yes" and ist == "no":
            xp = int(input("(x) Number of events: "))
            if infop == "no":
                compp = input("     What information do you need? ")
                print(" ")
                print(solve(function, lmbda=lmbdap, x=xp, comp=compp))
            else:
                print(" ")
                print(solve(function, lmbda=lmbdap, x=xp, info=infop))
        else:
            print(" ")
            print(solve(function, lmbda=lmbdap, info=infop))
    elif function == "uniform":
        au = int(input("(a) Lower bound: "))
        bu = int(input("(b) Higher bound: "))
        infou = input("     Do you want to print all? ").lower()
        isx = input("   Is there a (x) value? ").lower()
        if isx == "yes":
            xu = float(input("(x) Interested value: "))
            if infou == "no":
                compu = input("     What information do you need? ")
                print(" ")
                print(solve(function, a=au, b=bu, x=xu, comp=compu))
            else:
                print(" ")
                print(solve(function, a=au, b=bu, x=xu, info=infou))
        else:
            if infou == "no":
                compu = input("     What information do you need? ")
                print(" ")
                print(solve(function, a=au, b=bu, comp=compu))
            else:
                print(" ")
                print(solve(function, a=au, b=bu, info=infou))
    elif function == "exponential":
        lmdae = input("(lambda) Events rate: ")
        if len(lmdae.split()) == 2:
            if lmdae.split()[0] == "sqrt":
                lmdae = math.sqrt(float(lmdae.split()[1]))
        else:
            lmdae = float(lmdae)
        xe = int(input("(x) Interested value: "))
        infoe = input("     Do you want to print all? ").lower()
        isProba = input("   Is there a probability? ").lower()
        if isProba == "yes":
            probae = input("(Probability) probability rate ")
            if len(probae.split()) == 2:
                if probae.split()[0] == "sqrt":
                    lmdae = math.sqrt(float(lmdae.split()[1]))
            else:
                probae = float(probae)
            if infoe == "no":
                compe = input("     What information do you need? ")
                print(" ")
                print(solve(function, lmbda=lmdae, x=xe, proba=probae, comp=compe))
            else:
                print(" ")
                print(solve(function, lmbda=lmdae, x=xe, proba=probae, info=infoe))
        else:
            if infoe == "no":
                compe = input("     What information do you need? ")
                print(" ")
                print(solve(function, lmbda=lmdae, x=xe, comp=compe))
            else:
                print(" ")
                print(solve(function, lmbda=lmdae, x=xe, info=infoe))


def normalDistribution(mean, stdDev, side: str, *values):
    solution = 0
    if side == "below":
        if len(values) == 2:
            high = standardize(values[1], mean, stdDev)
            low = standardize(values[0], mean, stdDev)
            solution = round(norm.cdf(high), 4) - round(norm.cdf(low), 4)
        else:
            num = standardize(values[0], mean, stdDev)
            solution = norm.cdf(num)
    elif side == "above":
        num = standardize(values[0], mean, stdDev)
        solution = norm.cdf(num)
        solution = 1 - solution
    elif side == "equal":
        solution = 0
    return solution


def findProbabilityNormDist(side: str, *prob):
    solution = 0
    if side == "below":
        if len(prob) == 2:
            solution = norm.ppf(prob[1]) - norm.ppf(prob[0])
        elif len(prob) == 1:
            solution = norm.ppf(prob[0])
    elif side == "above":
        if len(prob) == 1:
            solution = norm.isf(prob[0])
    return solution


def standardize(value, mean, stdDev):
    newValue = (value - mean) / stdDev
    return newValue


def findX(mean, stdDev, side: str, *arg):
    solution = 0
    if side == "below":
        if len(arg) == 2:
            x = findProbabilityNormDist(side, arg[0])
            y = findProbabilityNormDist(side, arg[1])
            x = (x * stdDev) + mean
            y = (y * stdDev) + mean
            solution = y - x
        else:
            prob = findProbabilityNormDist(side, arg[0])
            solution = (prob * stdDev) + mean
    elif side == "above":
        prob = findProbabilityNormDist(side, arg[0])
        solution = (prob * stdDev) + mean
    return solution


def findStd(value, mean, prob, side: str):
    solution = 0
    if side == "above":
        prob = 1 - prob
        solution = (value - mean) / round(norm.ppf(prob), 2)
    elif side == "below":
        solution = (value - mean) / round(norm.ppf(prob), 2)
    return round(solution, 4)


def findMean(value, std, prob, side: str):
    solution = 0
    if side == "above":
        prob = 1 - prob
        solution = value - (round(norm.ppf(prob), 3) * std)
    elif side == "below":
        solution = value - (round(norm.ppf(prob), 3) * std)
    return round(solution, 4)


def linearCombMean(fMean, *mean):
    solution = fMean
    for i in mean:
        solution = solution + i
    return solution


def linearCombSTD(fStd, *std):
    solution = fStd ** 2
    for i in std:
        solution = solution + (i ** 2)
    solution = math.sqrt(solution)
    return solution


def findBounds(mean, std, probability):
    solution = 1 - probability
    bound = solution / 2
    solution = math.fabs(findX(mean, std, "below", bound))
    return solution


#
# mean1= 65
# stdd = 5
#
# solution = findX(mean1, stdd, "below", 0.025)
# print(solution)
# sol = findX(mean1, stdd, "below", 0.6)
# # # distributionSolve()
# # # sol = normalDistribution(meeaan, stdd, "below", 75)
# print(sol)
