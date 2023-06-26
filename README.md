# Adaptive-BVP
Implementation of the [A fast adaptive numerical method for stiff two-point boundary value problems](https://math.nyu.edu/~greengar/lee_gr_tpbvp.pdf) by June-Yub Lee and Leslie Greengard.


## Contents
- [Introduction](#introduction)
- [Linear solver](#linear-solver)
- [Nonlinear Solver](#nonlinear-solver)
- [Examples](#examples)
- [FAQ](#faq)
- [Future work](#future-work)

## Introduction
This project contains a single Python program named `BVPsolver.py`. To use the program, run
```
import BVPsolver
```

This program contains 2 main functions, `LinearBVPSolver` for linear BVPs, and `NewtonNonlinearSolver` for nonlinear BVPs.

## Linear solver

The function `LinearBVPSolver` has the following parameters:
```
def LinearBVPSolver(p, q, f, a, c, zetal0, zetal1, gammal, zetar0, zetar1, gammar,
    C=4.0, TOL=1e-8, iters=40, eval_points=32, force_double=False):
```
This function aims to solve the following ordinary differential equation:
$$u''(x) + p(x) u'(x) + q(x) u(x) = f(x)$$
with respect to boundary conditions:
$$ \zeta_{l0} u(a) + \zeta_{l1} u'(a) = \gamma_l$$
$$ \zeta_{r0} u(c) + \zeta_{r1} u'(c) = \gamma_r$$

- `p`: A function defined on $[a,c]$ and returns a `float` number.
- `q`: A function defined on $[a,c]$ and returns a `float` number.
- `f`: A function defined on $[a,c]$ and returns a `float` number.
- `a`: A `float` number for the left boundary of the interval.
- `c`: A `float` number for the right boundary of the interval.
- `zetal0, zetal1, gammal, zetar0, zetar1, gammar`: `float` numbers determining the boundary conditions $ \zeta_{l0} u(a) + \zeta_{l1} u'(a) = \gamma_l$ and $ \zeta_{r0} u(c) + \zeta_{r1} u'(c) = \gamma_r$.
- `C`: A `float` number with the same meaning described in the article. Default to 4.0.
- `TOL`: Error tolerance, default to 1e-8.
- `iters`: Maximum allowed iteration to solve the equation to the given error tolerance. Default to 40.
- `eval_points`: This program calculates the function values at `eval_points` random points on $[a,c]$ to estimate the error. Default to 32.
- `force_double`: The algorithm in the article contains one last doubling step, splitting every subinterval node into two to ensure that error stays within `TOL`. In my experiments, this step is usually unnecessary and doubles the execution time. Default to `False`. Switch it to `True` to force the same behavior described in the article.

This function returns a class `LinearBVPSolver.Node` that is directly callable. If you have already set `u = LinearBVPSolver(...)`, then you can:
- Access function value at `x`, i.e., $u(x)$, by direclty calling `u(x)`.
- Access function derivative at `x`, i.e., $u'(x)$, by calling `u(x, deriv=True)`.

## Nonlinear solver






