import sympy as sym
import math

x, y, v, p, lx, ly = sym.symbols('x y v p lx ly')
f = v * math.exp( (2 * math.sin(math.pi / p * abs(x - y)) ** 2) / lx ** 2 ) * math.exp( (x - y)^2 ** 0.5 / (2 * ly ** 2) )
deriv_f = f.diff(x)
deriv_f


