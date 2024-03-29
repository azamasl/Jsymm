using CPUTime, JLD2, LinearAlgebra, Convex, Random,Plots
include("../solvers.jl")
include("../sp_function.jl")


filename = "1-1000-BC-0.0001"
prob = load("x.jld2")

n = prob["n"];m = prob["m"];
B =prob["B"];C =prob["C"];A = Matrix(prob["A"]);
xstar=prob["xstar"];ystar=prob["ystar"];
x0   = prob["x0"];     y0=prob["y0"];
sp   = Saddle_Point(B, A, C  , xstar, ystar)
obj  = saddle_point_objective(sp)
nbFc = prob["nbFc"]
prt = 0;
stepsize = prob["stepsize"];max_it=prob["max_it"];  F_tol = prob["F_tol"];
my_tr    =prob["my_tr"];    max_tr= prob["max_tr"];    eta= prob["eta"];
F_tol= 1e-4; max_it = 2000;


function save_lssec()
    x_sol, y_sol,iter, nfs, val, ng  = secant_inv(x0,y0,obj,dummy,sp,max_it,prt,1,F_tol)
    @save "lssec.jld2"  x_sol y_sol iter nfs val ng nbFc F_tol max_it
end
function save_nosec()
    stepsize = 0.08; 
    x_sol, y_sol,iter, nfs, val, ng = secant_inv(x0,y0,obj,stepsize,sp,max_it, prt, 0,F_tol)
    @save "nosec.jld2"  x_sol y_sol iter nfs val ng stepsize max_it nbFc F_tol
end

function save_nobro()
    stepsize = 0.08; 
    x_sol, y_sol,iter, nfs, val, ng = Broyden(x0,y0,obj,stepsize,sp,max_it,prt, 0,F_tol)
    @save "nobro.jld2" x_sol y_sol iter nfs val ng stepsize max_it nbFc F_tol
end

function save_egm()
    x_sol, y_sol,iter, nfs, val, ng = EGM(x0,y0,obj,dummy,sp,max_it,prt,dummy,F_tol)
    @save "egm.jld2"  x_sol y_sol iter nfs val ng max_it nbFc F_tol
end

function save_tr()
    x_sol, y_sol,iter, nfs, val, ng = tr_dogleg(x0,y0, obj,sp, max_it, prt, F_tol, my_tr, max_tr, eta)
    @save "tr.jld2"  x_sol y_sol iter nfs val ng my_tr max_tr eta nbFc F_tol max_it
end



tim= @elapsed save_lssec();sol = load("lssec.jld2");  sol["time"] = tim; @save "lssec.jld2" sol
tim= @elapsed save_tr();   sol = load("tr.jld2");  sol["time"] = tim;    @save "tr.jld2" sol
tim= @elapsed save_nosec();sol = load("nosec.jld2");  sol["time"] = tim; @save "nosec.jld2" sol
tim= @elapsed save_nobro();sol = load("nobro.jld2");  sol["time"] = tim; @save "nobro.jld2" sol
tim= @elapsed save_egm();  sol = load("egm.jld2"); sol["time"] = tim;    @save "egm.jld2" sol
