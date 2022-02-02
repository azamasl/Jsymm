#Code for experiments in J-symm paper


solvers.jl: contains implementation for EGM, Broyden, and  various versions of J-symm methods.

sp_functions.jl contains the data structure and auxiliary functions

CC: contains  code to generate, run and plot a quadratic convex-concave example:

run settings.jl, to generate and save a  random problem

run run_all to solve and save the date with each of the methods

run plot_all to plot the data



AC: contains the code to generate, run and plot an instance of unconstrained Analytic Center (AC) problem.
See Convex Optimization Boyd and Vandenberghe, page 141 for the problem definition.

run setting_analytic_center to generate a random problem

run run_all_analytic_center to run various methods on the problem

run plot_all to plot the result


2D: contains code for nonconvex-nonconcave problem

run run_plot_2D_X where X is the method name, will run the X method from 12 initial points.
