from sympy import *
import numpy as np

if __name__ == '__main__':
    T = symbols('T')
    J = symbols('J', cls=Function)
    startx, starty, startz = symbols('startx, starty, startz')
    startxv, startyv, startzv = symbols('startxv, startyv, startzv')
    targetx, targety, targetz = symbols('targetx, targety, targetz')

    deltapx = targetx - startxv * T - startx
    deltapy = targety - startyv * T - starty
    deltapz = targetz - startzv * T - startz

    targetxv = 0
    targetyv = 0
    targetzv = 0

    deltavx = targetxv - startxv
    deltavy = targetyv - startyv
    deltavz = targetzv - startzv

    alph1 = -12 / T**3 * deltapx + 6 / T**2 * deltavx
    alph2 = -12 / T**3 * deltapy + 6 / T**2 * deltavy
    alph3 = -12 / T**3 * deltapz + 6 / T**2 * deltavz
    beta1 = 6 / T**2 * deltapx - 2 / T * deltavx
    beta2 = 6 / T ** 2 * deltapy - 2 / T * deltavy
    beta3 = 6 / T ** 2 * deltapz - 2 / T * deltavz

    J = T + (1/3 * alph1**2*T**3 + alph1*beta1*T**2 + beta1**2*T) + \
            (1/3 * alph2**2*T**3 + alph2*beta2*T**2 + beta2**2*T) + \
            (1/3 * alph3**2*T**3 + alph3*beta3*T**2 + beta3**2*T)

    J_expand = expand(J)
    J_expr = collect(J_expand, T)
    print(J_expr)
    print('***************')

    J_diff = diff(J, T)
    J_diff_expand = expand(J_diff)
    J_diff_expr = collect(J_diff_expand, T)
    print(J_diff_expr)

    coeff = [1, -2, 1]
    print(np.roots(coeff))