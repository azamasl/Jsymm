using CPUTime, JLD2, LinearAlgebra, Convex, Random
import FirstOrderLp
include("solver_analytic_center_specific.jl")
include("../sp_function.jl")
include("EGM.jl")

problem_name= "ran13x13" 
instance_path=joinpath("../data", "$(problem_name).mps.gz")
lp = FirstOrderLp.qps_reader_to_standard_form(instance_path)
lp.variable_upper_bound .= Inf
lp.variable_lower_bound .= 0.0

m,n = size(lp.constraint_matrix)
println([m;n;lp.num_equalities])

c, b, A=lp.objective_vector, lp.right_hand_side, lp.constraint_matrix
F_tol = 1e-4

sp  = Analytic_Center_Problem(A, b)
obj = unconst_Analytic_Center_ObjeGrad(sp)
dummy = 0;prt=0;
max_it = 2000
stepsize = 0.008
stepsize_tol = 1e-4 

e = 0.01; 
ran_seed= 1
#Random.seed!(ran_seed);
rng = MersenneTwister(ran_seed)
x0 = rand(rng,n); y0 = -A*x0+b .+ e; w0 = zeros(m,1); 
indices = findall(t -> t.<0, y0);y0[indices] .=0.1; 
@show minimum(abs.(y0))


function save_lssec()
    lsF_tol = 1e-4
    x_sol, y_sol, w_sol,iter, nfs, val, ng, ZZ  = secant_inv_analytic_center_backtrackingLS(x0,y0,w0,obj,dummy,sp,max_it,prt,lsF_tol)
    @save "lssec.jld2"  x_sol y_sol w_sol iter nfs val ng ZZ  lsF_tol
end

function save_nosec()
    stepsize = 0.04; max_it=2000;
    x_sol, y_sol, w_sol,iter, nfs, val, ng, ZZ  = secant_inv_analytic_center_fixedstep(x0,y0,w0,obj,stepsize,sp,max_it,prt,F_tol)
    @save "nosec.jld2"  x_sol y_sol w_sol iter nfs val ng stepsize ZZ max_it
end

function save_nobro()
    stepsize = 0.008;
    x_sol, y_sol, w_sol,iter, nfs, val, ng, ZZ  = broyden_analytic_center_fixedstep(x0,y0,w0,obj,stepsize,sp,max_it,prt,F_tol)
    @save "nobro.jld2"  x_sol y_sol w_sol iter nfs val ng stepsize ZZ max_it
end

function save_egm()
    F_tol = 1e-4
    x_sol, y_sol, w_sol,iter, nfs, val, ng, ZZ  = EGM(x0,y0,w0,obj,stepsize,sp,max_it,prt,F_tol)
    jldsave("egm.jld2";  x_sol, y_sol, w_sol, iter, nfs, val, ng, ZZ, F_tol)
end

tim= @elapsed save_egm();sol = load("egm.jld2");  sol["time"] = tim; @save "egm.jld2" sol
tim= @elapsed save_lssec();sol = load("lssec.jld2");  sol["time"] = tim; @save "lssec.jld2" sol
tim= @elapsed save_nosec();sol = load("nosec.jld2");  sol["time"] = tim; @save "nosec.jld2" sol
tim= @elapsed save_nobro();sol = load("nobro.jld2");  sol["time"] = tim; @save "nobro.jld2" sol

