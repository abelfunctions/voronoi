r"""
Kranich Homotopy
================

Methods for analytically continuing a y-fibre along a path in the complex
x-plane based on the method of Kranich.

Functions
---------

analytically_continue

"""

from sage.all import QQ, QQbar, CDF, CC, Permutation, Infinity, sqrt

import scipy
import scipy.optimize


def matching_permutation(a, b):
    r"""Determines if a is an "approximate" permuation of b.

    (Copied from Abelfunctions.)
    """
    N = len(a)
    if N != len(b):
        raise ValueError, "Lists must be of same length."

    perm = [-1]*N
    eps  = 0.5*min([abs(a[i]-a[j]) for i in range(N) for j in range(i)])
    for i in xrange(N):
        for j in xrange(N):
            dist = abs(a[i] - b[j])
            if dist < eps:
                perm[j] = i+1
                break

    # an element is not accounted for if the numerical distance is too large
    if -1 in perm:
        raise ValueError, "Could not compute matching permutation " \
                          "between %s and %s." %(a,b)

    return Permutation(perm)


def compute_epsilon(yfibre):
    r"""
    Note that there will be at least one y-fibre element. Otherwise, analytic
    continuation is trivial.
    """
    d = len(yfibre)
    epsilon = Infinity

    # compute the pairwise distances between all y-fibre elements
    for i in xrange(d):
        for j in xrange(i):
            dist = abs(yfibre[i] - yfibre[j])
            epsilon = dist if dist < epsilon else epsilon


    epsilon *= 0.5
    return epsilon


def analytically_continue(f, gamma, y0):
    r"""
    Analytically continues a y-fibre `y0` along the path `gamma`.

    Parameters
    ----------
    f : complex algebraic curve
    gamma : parameterized complex path
    y0 : y-fibre

    Returns
    -------
    y1 : orded y-fibre above end of gamma
    """
    # first verify that y0 is a y-fibre above x0
    x,y = f.parent().gens()
    x0 = gamma(0)
    f_x0 = f(x0,y).univariate_polynomial()
    y0_test = f_x0.roots(ring=CDF, multiplicities=False)
    phi = matching_permutation(y0, y0_test)  # asserts the permutation exists

    # obtain necessary tools from the curve
    fx = f.derivative(x)
    fy = f.derivative(y)
    ratio = fx/fy
    a = f.polynomial(y).coefficients(sparse=False)
    disc = f.discriminant(y).univariate_polynomial()
    p = a[0]*disc
    rts = p.roots(CDF, multiplicities=False)

    # loop
    T = 0.0
    x1 = x0
    y1 = y0
    while (T < 1.0):
        epsilon = compute_epsilon(y1)

        # determine a sufficient value of rho
        rho = min(abs(x1 - rt) for rt in rts) - 1e-8

        # obtain the parameter Y
        Y = max(ratio(x1,y1j) for y1j in y1)

        # obtian the parameter M
        bounds = [
            max(abs(coeff)*(abs(x1) + rho)**expon
                for expon,coeff in ak.dict().iteritems())
            for ak in a if ak
        ]
        M = max(bounds)

        # compute the x-step size, delta
        delta = rho*(sqrt((rho*Y-epsilon)**2 + 4*epsilon*M) - (rho*Y+epsilon))
        delta /= (2*(M-rho*Y))
        delta = abs(delta)

        # maximize Ts such that |x(T) - X(Ts)| < delta
        g = lambda Ts:  abs(gamma(T) - gamma(Ts)) - delta
        if (T < 1.0):
            T = CDF(scipy.optimize.bisect(g,T+1e-14,1))

        # update
        #
        # TODO: this part is still pretty broken
        x1 = gamma(T)
        f_x1 = f(x1,y).univariate_polynomial()
        y1_next = f_x1.roots(ring=CDF, multiplicities=False)
        phi = matching_permutation(y1, y1_next)
        y1 = phi.action(y1_next)

    return y1


if __name__ == '__main__':
    R = QQ['x,y']
    x,y = R.gens()
    f = y**3 - 2*x**3*y + x**7

    start = CDF(-1)
    end = CDF(-0.5)
    gamma = lambda s: start*(1-s) + end*s

    y0 = f(start,y).univariate_polynomial().roots(
        ring=CDF, multiplicities=False
    )

    y1 = analytically_continue(f, gamma, y0)
    print y1
