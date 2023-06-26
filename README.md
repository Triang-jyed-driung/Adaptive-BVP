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

This program contains 2 main functions, `LinearBVPSolver()` for linear BVPs, and `NewtonNonlinearSolver()` for nonlinear BVPs.

## Linear solver




