using CPUTime, JLD2, LinearAlgebra, Convex, Random
include("solver_analytic_center_specific.jl")
include("../sp_function.jl")

#Loading the problem
filename = "UCAC-600"
#filename = "UCAC-600-rs7"
prob = load("output/problem-"*filename*".jld2")
n = prob["n"];m = prob["m"];
A = prob["A"];b =prob["b"];xstar=prob["xstar"];
F_tol = prob["F_tol"]
x0   = prob["x0"]; y0   = prob["y0"]; w0   = prob["w0"];
sp  = Analytic_Center_Problem(A, b)#initialize prob def struct
obj = unconst_Analytic_Center_ObjeGrad(sp)
dummy = 0;prt=0;
max_it = 2000

function save_newton()
    nt_F_tol = 1e-7
    x_sol, y_sol, w_sol,iter, nfs, val, ng, ZZ  = Elimination_Newton_analytic_center_backtrackingLS(x0,y0,w0,obj,dummy,sp,50,prt,nt_F_tol)
    @save "output/"*filename*"-newton.jld2"  x_sol y_sol w_sol iter nfs val ng ZZ nt_F_tol
end

function save_lssec()
    lsF_tol = 1e-4
    x_sol, y_sol, w_sol,iter, nfs, val, ng, ZZ  = secant_inv_analytic_center_backtrackingLS(x0,y0,w0,obj,dummy,sp,max_it,prt,lsF_tol)
    @save "output/"*filename*"-lssec.jld2"  x_sol y_sol w_sol iter nfs val ng ZZ  lsF_tol
end

function save_nosec()
    #stepsize = 0.01; max_it=2000;
    x_sol, y_sol, w_sol,iter, nfs, val, ng, ZZ  = secant_inv_analytic_center_fixedstep(x0,y0,w0,obj,stepsize,sp,max_it,prt,F_tol)
    @save "output/"*filename*"-nosec.jld2"  x_sol y_sol w_sol iter nfs val ng stepsize ZZ max_it
end

function save_nobro()
    #stepsize = 0.01;max_it=2000;
    x_sol, y_sol, w_sol,iter, nfs, val, ng, ZZ  = broyden_analytic_center_fixedstep(x0,y0,w0,obj,stepsize,sp,max_it,prt,F_tol)
    @save "output/"*filename*"-nobro.jld2"  x_sol y_sol w_sol iter nfs val ng stepsize ZZ max_it
end


tim= @elapsed save_newton();sol = load("output/"*filename*"-newton.jld2"); sol["time"] = tim; @save "output/"*filename*"-newton.jld2" sol
tim= @elapsed save_lssec();sol = load("output/"*filename*"-lssec.jld2");  sol["time"] = tim; @save "output/"*filename*"-lssec.jld2" sol
tim= @elapsed save_nosec();sol = load("output/"*filename*"-nosec.jld2");  sol["time"] = tim; @save "output/"*filename*"-nosec.jld2" sol
tim= @elapsed save_nobro();sol = load("output/"*filename*"-nobro.jld2");  sol["time"] = tim; @save "output/"*filename*"-nobro.jld2" sol
