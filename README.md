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
This project contains a single Python program named `BVPSolver.py`. To use the program, run
```
import BVPsolver
```

This program contains 2 main functions, `LinearBVPSolver` for linear BVPs, and `NewtonNonlinearSolver` for nonlinear BVPs.

## Linear solver

The function `LinearBVPSolver` has the following parameters:
```Python
def LinearBVPSolver(p, q, f, a, c, zetal0, zetal1, gammal, zetar0, zetar1, gammar,
    C=4.0, TOL=1e-8, iters=40, eval_points=32, force_double=False):
```
This function aims to solve the following ordinary differential equation:
$$u''(x) + p(x) u'(x) + q(x) u(x) = f(x)$$
with respect to boundary conditions:
$$\zeta_{l0} u(a) + \zeta_{l1} u'(a) = \gamma_l$$
and
$$\zeta_{r0} u(c) + \zeta_{r1} u'(c) = \gamma_r$$.

- `p`: A function defined on $[a,c]$ and returns a `float` number.
- `q`: A function defined on $[a,c]$ and returns a `float` number.
- `f`: A function defined on $[a,c]$ and returns a `float` number.
- `a`: A `float` number for the left boundary of the interval.
- `c`: A `float` number for the right boundary of the interval.
- `zetal0, zetal1, gammal, zetar0, zetar1, gammar`: `float` numbers determining the boundary conditions $\zeta_{l0} u(a) + \zeta_{l1} u'(a) = \gamma_l$
 and 
$\zeta_{r0} u(c) + \zeta_{r1} u'(c) = \gamma_r$.
- `C`: A `float` number with the same meaning described in the article. Default to 4.0.
- `TOL`: Error tolerance, default to 1e-8.
- `iters`: Number of allowed refinement iterations to solve the equation to the given error tolerance. Default to 40.
- `eval_points`: This program calculates the function values at `eval_points` random points on $[a,c]$ to estimate the error. Default to 32.
- `force_double`: The algorithm in the article contains one last doubling step, splitting every subinterval node into two to ensure that error stays within `TOL`. In my experiments, this step is usually unnecessary and doubles the execution time. Default to `False`. Switch it to `True` to force the same behavior described in the article.

This function returns a class `LinearBVPSolver.Node` that is directly callable. If you have already set `u = LinearBVPSolver(...)`, then you can:
- Access function value at `x`, i.e., $u(x)$, by direclty calling `u(x)`.
- Access function derivative at `x`, i.e., $u'(x)$, by calling `u(x, deriv=True)`.

Example usage:

```Python
from scipy.special import jv
from matplotlib import pyplot as plt
import numpy as np
from BVPSolver import LinearBVPSolver
u = LinearBVPSolver(
    lambda x: 1/x,
    lambda x: 1-(100/x)**2,
    lambda x: 0,
    0, 600,
    1,0,0,
    1,0,1,
    TOL=5e-10,
)
xs = np.linspace(0, 600, 600)
plt.plot(xs, [u(xi) for xi in xs])
plt.title("Bessel 100 on [0,600]")
plt.show()
```

## Nonlinear solver

The function `NewtonNonlinearSolver` has the following parameters:
```Python
def NewtonNonlinearSolver(
    f, f_2, f_3, a, c, zetal0, zetal1, gammal, zetar0, zetar1, gammar,
    initial=None, C=4.0, TOL=1e-6, iters=10, eval_points=32, force_double=False,
):
```
This function aims to solve the following ordinary differential equation:
$$u''(x) = f(x, u, u')$$
with respect to boundary conditions:
$$\zeta_{l0} u(a) + \zeta_{l1} u'(a) = \gamma_l$$
and
$$\zeta_{r0} u(c) + \zeta_{r1} u'(c) = \gamma_r$$.

- `f`: A function that accepts 3 `float` numbers (namely `x, u, du` where `du` is $u'(x)$) and returns a `float` number.
- `f_2`: Partial derivative of `f` with respect to the second parameter, which is $\partial_2 f$. It also accepts 3 `float` numbers (namely `x, u, du` where `du` is $u'(x)$) and returns a `float` number.
- `f_3`: Partial derivative of `f` with respect to the third parameter, which is $\partial_3 f$. It also accepts 3 `float` numbers (namely `x, u, du` where `du` is $u'(x)$) and returns a `float` number.
- `a`: Same as above
- `c`: Same as above
- `zetal0, zetal1, gammal, zetar0, zetar1, gammar`: Same as above
- `initial`: Provide an initial guess for the Newton's method. Default to `None`. If you want to provide an initial guess, set `initial=(u0,du0)` where `u0` is the initial guess function and `du0` is its derivative.
- `C`: Same as above
- `TOL`: Error tolerance, default to 1e-6.
- `iters`: Maximum number of Newton iterations to solve the equation. Default to 10.
- `eval_points`: Same as above
- `force_double`: Same as above

The formula of Newton's iteration is not explicitly specified in the article. Here is the formula:

Let $u_{n+1}(x)$ solve 
```math
u''_{n+1} -\partial_3 f(x, u_n, u'_n) u'_{n+1} -\partial_2 f(x, u_n, u'_n) u_{n+1}=f(x, u_n, u'_n) -\partial_3 f(x, u_n, u'_n) u'_n -\partial_2 f(x, u_n, u'_n) u_n
```

Try this example:
```Python
from matplotlib import pyplot as plt
import numpy as np
from BVPSolver import NewtonNonlinearSolver
r = NewtonNonlinearSolver(
    lambda x, u, du: 2*u*du,
    lambda x, u, du: 2*du,
    lambda x, u, du: 2*u,
    -np.pi/4, np.pi/4,
    1,0,-1,
    1,0,1,
    initial=(lambda x: x**2, lambda x: 2*x),
    TOL=1e-13
)
xs = np.linspace(-np.pi/4, np.pi/4, 100)
plt.plot(xs, [r(xi) for xi in xs])
plt.title("u''=2uu', u=tanx")
plt.show()
```
