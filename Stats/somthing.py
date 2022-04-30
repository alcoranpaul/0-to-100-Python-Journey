from Tools.scripts.dutree import display
from sympy import *
init_printing()
eq = Eq(sympify('(a**2 + sqrt(b)) / log(x)'))
display(eq)