#Code for experiments in J-symm paper


solvers.jl: contains implementation for EGM, Broyden, and  various versions of J-symm methods.

sp_functions.jl contains the data structure and auxiliary functions


CC: contains  code to generate, run and plot a quadratic convex-concave example:

run settings.jl, to generate and save a  random problem

run run_all to solve and save the data with each of the methods

run plot_all to plot the data


bilinear: contains the code to run and plot a bilinear zero-sum example:

run setting to load the instance and setup

run run_all to run various methods to solve the bilinear problem

run plot_all to plot the result


AC: contains the code to generate, run and plot an instance of unconstrained Analytic Center (AC) problem.
See Convex Optimization Boyd and Vandenberghe, page 141 for the problem definition.

run setting_analytic_center to generate a random problem

run run_all_analytic_center to run various methods on the problem

run plot_all to plot the result

run real_AC to run various methods on real instances of analytic center 

run real_plot to plot the result


2D: contains code for nonconvex-nonconcave problem

run run_plot_2D_X where X is the method name, will run the X method from 12 initial points.
