from matplotlib import pyplot as plt
import numpy as np
from math import *
from BVPSolver import *
import time


if __name__ == "__main__":

    from scipy.special import jv

    t0 = time.time_ns()
    root = LinearBVPSolver(
        lambda x: 1/x,
        lambda x: 1-(100/x)**2,
        lambda x: 0,
        0, 600,
        1,0,0,
        1,0,1,
        TOL=5e-10,
        eval_points=100
    )
    print(f"time: {(time.time_ns()-t0)/1e9}")
    xs = np.linspace(0,600,600)
    plt.plot(xs, [(root(xi) - jv(100,xi)/jv(100,600)) for xi in xs])
    plt.title("Bessel 100 on [0,600]")
    plt.savefig("bessel.jpg")
    plt.show()
    

    t0 = time.time_ns()
    root = LinearBVPSolver(
        lambda x: 0,
        lambda x: 1-x**2,
        lambda x: 0,
        -6,6,
        1,0,exp(-18),
        1,0,exp(-18),
        TOL=1e-3,
    )
    print(f"time: {(time.time_ns()-t0)/1e9}")
    xs = np.linspace(-6,6,100)
    plt.plot(xs, [(root(xi) - exp(-xi**2/2)) for xi in xs])
    plt.title("u''+(1-x^2)u=0, exp(-x^2/2)")
    plt.savefig("exp-x2.jpg")
    plt.show()
    
    t0 = time.time_ns()
    root = LinearBVPSolver(
        lambda x: 0,
        lambda x: 1,
        lambda x: 0,
        0, pi/2,
        0,1,1,
        0,1,0,
    )
    print(f"time: {(time.time_ns()-t0)/1e9}")
    xs = np.linspace(0, pi/2, 100)
    plt.plot(xs, [root(xi)-sin(xi) for xi in xs])
    plt.title("u''+u=0,u'(0)=1,u'(pi/2)=0 (Pure Neumann)")
    plt.savefig("sinx.jpg")
    plt.show()
    
    t0 = time.time_ns()
    root = NewtonNonlinearSolver(
        lambda x, u, du: 2*u*du, 
        lambda x, u, du: 2*du, 
        lambda x, u, du: 2*u, 
        -pi/4, pi/4,
        1,0,-1,
        1,0,1,
        initial=(lambda x: x**2, lambda x: 2*x),
        TOL=1e-13
    )
    print(f"time: {(time.time_ns()-t0)/1e9}")
    xs = np.linspace(-pi/4, pi/4, 100)
    plt.plot(xs, [root(xi)-tan(xi) for xi in xs])
    plt.title("u''=2uu', u=tanx")
    plt.savefig("tanx.jpg")
    plt.show()
    
    t0 = time.time_ns()
    root = NewtonNonlinearSolver(
        lambda x, u, du: du**2 / (2*u), 
        lambda x, u, du: - (du / u)**2/2, 
        lambda x, u, du: du / u, 
        -1, 1,
        1,0,1,
        1,0,1,
        initial = (lambda x: x, lambda x: 1),
        TOL=1e-8,
    )
    print(f"time: {(time.time_ns()-t0)/1e9}")
    xs = np.linspace(-1, 1, 100)
    plt.plot(xs, [root(xi)-xi**2 for xi in xs])
    plt.title("u''=u'^2/(2u), u=x^2")
    plt.savefig("x^2.jpg")
    plt.show()

    t0 = time.time_ns()
    root = NewtonNonlinearSolver(
        lambda x, u, du: du * (1-du), 
        lambda x, u, du: 0, 
        lambda x, u, du: 1 - 2 * du, 
        -log(2), log(2),
        0, 1, 1/3,
        1, 0, log(3),
        initial = (lambda x: (x+1)/2, lambda x: 1/2),
        TOL=1e-8,
    )
    print(f"time: {(time.time_ns()-t0)/1e9}")
    xs = np.linspace(-log(2), log(2), 100)
    plt.plot(xs, [root(xi)-log(1+exp(xi)) for xi in xs])
    plt.title("u''=u'(1-u'), u=ln(1+exp(x))")
    plt.savefig("logistic.jpg")
    plt.show()
    
