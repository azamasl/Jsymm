#J-symm code

CC: contians the code to generate run and plot a quadratic convex-concave example

AC: to generate run and plot an instance of unconstrained Analytic Center (AC) problem.
See Convex Optimization Boyd and Vandenberghe, page 141 for the problem definition.

2D: nonconvex-nonconcave problem

solver.jl: contains implementation for EGM, Broden, and  various versions of J-symm methods.

sp_functions.jl contains the data structure and auxiliary functions

To run a Bilinear/Strongly-convex-concave example:

run settings.jl, to generate and save a  random problem
run run_all to solve and save the date with each of of the methods
run plot_all to plot the data


To run a AC example:
run setting_analytic_center to generate a random problem
run run_all_analytic_center to run various methods on the problem
run plot_all to plot the result

To run the 2D example:
run run_plot_2D_X where X is the method name will run the X method from 12 initial point on the nonconvex-nonconcave problem.
