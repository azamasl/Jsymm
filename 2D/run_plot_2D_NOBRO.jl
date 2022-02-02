using CPUTime, JLD2, LinearAlgebra, Convex, Plots
include("../solvers.jl")
include("../sp_function.jl")

prob_type= "2D" # Resource Allocation problem
interaction = 100.0
bd = 150 # A = 100
#bd = 50# A =1000
#bd = 4 # for interaction A = 1,10
xmin, xmax = -bd,bd
ymin, ymax = -bd,bd

max_it = 50
prt = 0 # 0 don't print grad norm at every iterations; 1, do it.
F_tol =1e-6
reset_in_pt = 0
dummy=0

"local alg settings"
c1 = 1e-4#Armijo parameter
stepsize = 0.01#The fixed step_size, only used when do_ls = 0

"global alg settings"
tr = 0.1  #initial Δ
max_tr = 0.5 #maximum allowed Δ
eta = 0.01
TYPE="Nonconvex-nonconcave 2D function"


obj = sad_point2D_objective(interaction)
prt=0
#################################################################
# Run&plot for 12 initial points
#################################################################
x_sol, y_sol,iter, nfs, val, ng,ZZ  = Broyden([-4.0],[-2.0],obj,stepsize,dummy,max_it,prt,0,F_tol)
plot( ZZ[1,:],ZZ[2,:], marker = (:circle, 4))
F_sol = [obj.∇xL(x_sol,y_sol); -obj.∇yL(x_sol,y_sol)]
nabla_F_sol  = obj.∇F(x_sol,y_sol);
g_sol = nabla_F_sol'*F_sol;
norm_g2 = norm(g_sol)
print("||g|| = "); println(norm_g2)


x_sol, y_sol,iter, nfs, val, ng,ZZ  = Broyden([-4.0],[0.0],obj,stepsize,dummy,max_it,prt,0,F_tol)
plot!( ZZ[1,:],ZZ[2,:], marker = (:circle, 4))
F_sol = [obj.∇xL(x_sol,y_sol); -obj.∇yL(x_sol,y_sol)]
nabla_F_sol  = obj.∇F(x_sol,y_sol);
g_sol = nabla_F_sol'*F_sol;
norm_g2 = norm(g_sol)
print("||g|| = "); println(norm_g2)

x_sol, y_sol,iter, nfs, val, ng,ZZ  = Broyden([-4.0],[2.0],obj,stepsize,dummy,max_it,prt,0,F_tol)
plot!( ZZ[1,:],ZZ[2,:], marker = (:circle, 4))
F_sol = [obj.∇xL(x_sol,y_sol); -obj.∇yL(x_sol,y_sol)]
nabla_F_sol  = obj.∇F(x_sol,y_sol);
g_sol = nabla_F_sol'*F_sol;
norm_g2 = norm(g_sol)
print("||g|| = "); println(norm_g2)

x_sol, y_sol,iter, nfs, val, ng,ZZ  = Broyden([-2.0],[-4.0],obj,stepsize,dummy,max_it,prt,0,F_tol)
#@save "output/$prob_type-$interaction-lssec.jld2"  x_sol y_sol iter nfs val ng ZZ
plot!( ZZ[1,:],ZZ[2,:], marker = (:circle, 4))
F_sol = [obj.∇xL(x_sol,y_sol); -obj.∇yL(x_sol,y_sol)]
nabla_F_sol  = obj.∇F(x_sol,y_sol);
g_sol = nabla_F_sol'*F_sol;
norm_g2 = norm(g_sol)
print("||g|| = "); println(norm_g2)

x_sol, y_sol,iter, nfs, val, ng,ZZ  = Broyden([-2.0],[4.0],obj,stepsize,dummy,max_it,prt,0,F_tol)
plot!( ZZ[1,:],ZZ[2,:], marker = (:circle, 4))
F_sol = [obj.∇xL(x_sol,y_sol); -obj.∇yL(x_sol,y_sol)]
nabla_F_sol  = obj.∇F(x_sol,y_sol);
g_sol = nabla_F_sol'*F_sol;
norm_g2 = norm(g_sol)
print("||g|| = "); println(norm_g2)

x_sol, y_sol,iter, nfs, val, ng,ZZ  = Broyden([0.0],[-4.0],obj,stepsize,dummy,max_it,prt,0,F_tol)
plot!( ZZ[1,:],ZZ[2,:], marker = (:circle, 4))
F_sol = [obj.∇xL(x_sol,y_sol); -obj.∇yL(x_sol,y_sol)]
nabla_F_sol  = obj.∇F(x_sol,y_sol);
g_sol = nabla_F_sol'*F_sol;
norm_g2 = norm(g_sol)
print("||g|| = "); println(norm_g2)

x_sol, y_sol,iter, nfs, val, ng,ZZ  = Broyden([0.0],[4.0],obj,stepsize,dummy,max_it,prt,0,F_tol)
plot!( ZZ[1,:],ZZ[2,:], marker = (:circle, 4))
F_sol = [obj.∇xL(x_sol,y_sol); -obj.∇yL(x_sol,y_sol)]
nabla_F_sol  = obj.∇F(x_sol,y_sol);
g_sol = nabla_F_sol'*F_sol;
norm_g2 = norm(g_sol)
print("||g|| = "); println(norm_g2)


x_sol, y_sol,iter, nfs, val, ng,ZZ  = Broyden([2.0],[-4.0],obj,stepsize,dummy,max_it,prt,0,F_tol)
plot!( ZZ[1,:],ZZ[2,:], marker = (:circle, 4))
F_sol = [obj.∇xL(x_sol,y_sol); -obj.∇yL(x_sol,y_sol)]
nabla_F_sol  = obj.∇F(x_sol,y_sol);
g_sol = nabla_F_sol'*F_sol;
norm_g2 = norm(g_sol)
print("||g|| = "); println(norm_g2)


x_sol, y_sol,iter, nfs, val, ng,ZZ  = Broyden([2.0],[4.0],obj,stepsize,dummy,max_it,prt,0,F_tol)
plot!( ZZ[1,:],ZZ[2,:], marker = (:circle, 4))
F_sol = [obj.∇xL(x_sol,y_sol); -obj.∇yL(x_sol,y_sol)]
nabla_F_sol  = obj.∇F(x_sol,y_sol);
g_sol = nabla_F_sol'*F_sol;
norm_g2 = norm(g_sol)
print("||g|| = "); println(norm_g2)

x_sol, y_sol,iter, nfs, val, ng,ZZ  = Broyden([4.0],[-2.0],obj,stepsize,dummy,max_it,prt,0,F_tol)
plot!( ZZ[1,:],ZZ[2,:], marker = (:circle, 4))
F_sol = [obj.∇xL(x_sol,y_sol); -obj.∇yL(x_sol,y_sol)]
nabla_F_sol  = obj.∇F(x_sol,y_sol);
g_sol = nabla_F_sol'*F_sol;
norm_g2 = norm(g_sol)
print("||g|| = "); println(norm_g2)

x_sol, y_sol,iter, nfs, val, ng,ZZ  = Broyden([4.0],[0.0],obj,stepsize,dummy,max_it,prt,0,F_tol)
plot!( ZZ[1,:],ZZ[2,:], marker = (:circle, 4))
F_sol = [obj.∇xL(x_sol,y_sol); -obj.∇yL(x_sol,y_sol)]
nabla_F_sol  = obj.∇F(x_sol,y_sol);
g_sol = nabla_F_sol'*F_sol;
norm_g2 = norm(g_sol)
print("||g|| = "); println(norm_g2)

x_sol, y_sol,iter, nfs, val, ng,ZZ  = Broyden([4.0],[2.0],obj,stepsize,dummy,max_it,prt,0,F_tol)
plot!( ZZ[1,:],ZZ[2,:], marker = (:circle, 4))
F_sol = [obj.∇xL(x_sol,y_sol); -obj.∇yL(x_sol,y_sol)]
nabla_F_sol  = obj.∇F(x_sol,y_sol);
g_sol = nabla_F_sol'*F_sol;
norm_g2 = norm(g_sol)
print("||g|| = "); println(norm_g2)

#sF= @sprintf(", cond(∇F)= %2.1f",nbFc)
#sL= @sprintf(", cond(∇L)= %2.1f",nbLc)

plot!(
#title = string(funcName, ", ", TYPE,", A = $interaction, fixed stepsize=$stepsize "), titlefont=font(8),
legend = false,aspect_ratio=:equal
,xlims=(xmin,xmax), ylims=(ymin,ymax)
)

savefig("plots/$prob_type-NOLS_Broyden-$interaction.png")
