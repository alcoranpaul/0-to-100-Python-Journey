from math import comb

# DISTRIBUTIONS
from tokenize import String


def binomial(x, n, p):
    """
    Binomial

    :param x: number of successes seeking
    :param n: total population
    :param p: probability of success
    :return:
    """
    # expression_probability = (comb(n, x)) * (p ** x) * ((1 - p) ** (n - x))
    return (comb(n, x)) * (p ** x) * ((1 - p) ** (n - x))
    # expression_variance = p * (1 - p)
    # print(f"P(X = {x}) = {expression_probability:.4f}")
    # print(f"E(X) = {p:.4f}")
    # print(f"E(x^2) = {p:.4f}")
    # print(f"Var(X) = {expression_variance:.4f}")


def geometric(x, p):
    """
    Geometric

    :param x: index of first success
    :param p: probability of success
    :return:
    """
    expression_probability = ((1 - p) ** (x - 1)) * p
    expression_variance = (1 - p) / (p ** 2)
    expression_expected = 1 / p
    expression_expectedSquared = (2 - p) / (p ** 2)
    print(f"P(X = {x}) = {expression_probability:.4f}")
    print(f"E(X) = {expression_expected:.4f}")
    print(f"E(x^2) = {expression_expectedSquared:.4f}")
    print(f"Var(X) = {expression_variance:.4f}")


def negativeBinomial(x, r, p):
    """
    Negative Binomial

    :param x: number of trials
    :param r: number of successes on x th trial
    :param p: probability of success
    :return:
    """
    expression_probability = comb(x - 1, r - 1) * (1 - p) ** (x - r) * p ** r
    expression_expected = r / p
    expression_variance = (r * (1 - p)) / (p ** 2)
    print(f"P(X = {x}) = {expression_probability:.4f}")
    print(f"E(X) = {expression_expected:.4f}")
    print(f"Var(X) = {expression_variance:.4f}")


def hypergeometric(x, N, n, r):
    """
    Hyper-geometric

    :param x: number of successes seeking
    :param N: total population size
    :param n: number of selections made
    :param r: total number of successes
    """
    expression_probability = (comb(r, x) * comb(N - r, n - x)) / comb(N, n)
    print(f"P(X = {x}) = {expression_probability:.4f}")


def bernoulliTrials(x, p):
    """
    bernoulliTrials

    :param x:
    :param p:
    :return:
    """
    expression_probability = (p ** x) * ((1 - p) ** (1 - x))
    expression_expected = p
    expression_expectedSquared = p
    expression_variance = p * (1 - p)
    print(f"P(X = {x}) = {expression_probability:.4f}")
    print(f"E(X) = {expression_expected:.4f}")
    print(f"E(x^2) = {expression_expectedSquared:.4f}")
    print(f"Var(X) = {expression_variance:.4f}")


# function = (input("What kind of Distribution? "))
# if function == "binomial":
#     x = int(input("(x) number of successes seeking: "))
#     n = int(input("(n) total population: "))
#     p = float(input("(p) probability of success: "))
#     print(" ")
#     binomial(x, n, p)
# elif function == "geometric":
#     x = int(input("(x) index of first success: "))
#     p = float(input("(p) probability of success: "))
#     print(" ")
#     geometric(x, p)
# elif function == "negativeBinomial":
#     x = int(input("(x) number of successes seeking: "))
#     r = int(input("(r) number of successes on x th trial: "))
#     p = float(input("(p) probability of success: "))
#     print(" ")
#     negativeBinomial(x, r, p)
# elif function == "hypergeometric":
#     x = int(input("(x) number of successes seeking: "))
#     N = int(input("(N) total population size: "))
#     n = int(input("(n) number of selections made: "))
#     r = int(input("(r) total number of successes: "))
#     print(" ")
#     hypergeometric(x, N, n, r)
# elif function == "bernoulliTrials":
#     x = int(input("(x) index of first success: "))
#     p = float(input("(p) probability of success: "))
#     print(" ")
#     bernoulliTrials(x, p)
