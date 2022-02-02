using CPUTime, JLD2, LinearAlgebra, Convex, Plots
include("../sp_function.jl")
ymin, ymax = -4,4
xmin, xmax = -4,4

prob_type= "2D" # Resource Allocation problem
interaction = 1000.0
funcName = "EGM w optimal stpsz"
max_it = 100
prt = 0 # 0 don't print grad norm at every iterations; 1, do it.
F_tol =1e-6
reset_in_pt = 0
dummy=0

"local alg settings"
c1 = 1e-4#Armijo parameter
stepsize = 0.05#The fixed step_size, only used when do_ls = 0

"global alg settings"
tr = 0.1  #initial Δ
max_tr = 0.5 #maximum allowed Δ
eta = 0.01
TYPE="Nonconvex-nonconcave 2D function"
obj = sad_point2D_objective(interaction)
#@show nbc =  cond(obj.∇F)#sanity check

#################################
# EGM with the optimal fixed stepsize for the 2D nonconvex-nonconcave function
# See Sean&Grimmer landscap paper page 27:
# 2D function is 172-smooth, i.e. : β = 172. The optimal fixed stepsize is defined
# eta = 0.5/(β+A)
#################################
function EGM_eta(x,y,obj,dummy0,sp,max_it,prt, use_ls,Ftol)
    println("dummy print, to make sure next line is not being eatten up")
    println("EGM with optimal fixed stepsize with interaction = $interaction")
    nablaF  = obj.∇F(x,y)
    #eta  = 1/LinearAlgebra.norm(nablaF)
    eta  = 1/(2*(interaction + 172))
    val, ngx, ngy =0,0,0
    iter=0
    F = [obj.∇xL(x,y); -obj.∇yL(x,y)]
    normF=LinearAlgebra.norm(F)
    normFAll=[]
    Zs=[x;y]
    while ( normF> Ftol ) && iter < max_it
        xk = x -eta*obj.∇xL(x,y)
        yk = y +eta*obj.∇yL(x,y)

        x= x -eta*obj.∇xL(xk,yk)
        y= y +eta*obj.∇yL(xk,yk)

        Zs = hcat(Zs,[x;y])
        F = [obj.∇xL(x,y); -obj.∇yL(x,y)]
        append!(normFAll, normF)
        normF = LinearAlgebra.norm(F)

        val =obj.L(x,y)
        ngx = LinearAlgebra.norm(obj.∇xL(x,y))
        ngy = LinearAlgebra.norm(obj.∇yL(x,y))
        iter=iter+1
        if prt==1
            println("L = $val")
            println("|∇xL| = $ngx")
            println("|∇yL| = $ngy")
            println("iter $iter")
        end
    end
    ng = sqrt(ngx^2+ngy^2)
    return x,y,iter,normFAll, val, ng,Zs
end

#################################################################
# Run&plot for 12 initial points
#################################################################

x_sol, y_sol,iter, nfs, val, ng,ZZ  = EGM_eta([-4.0],[-2.0],obj,dummy,dummy,max_it,prt,dummy,F_tol)
plot( ZZ[1,:],ZZ[2,:], marker = (:circle, 4))
#
x_sol, y_sol,iter, nfs, val, ng,ZZ  = EGM_eta([-4.0],[0.0],obj,dummy,dummy,max_it,prt,dummy,F_tol)
plot!( ZZ[1,:],ZZ[2,:], marker = (:circle, 4))

x_sol, y_sol,iter, nfs, val, ng,ZZ  = EGM_eta([-4.0],[2.0],obj,dummy,dummy,max_it,prt,dummy,F_tol)
plot!( ZZ[1,:],ZZ[2,:], marker = (:circle, 4))

x_sol, y_sol,iter, nfs, val, ng,ZZ  = EGM_eta([-2.0],[-4.0],obj,dummy,dummy,max_it,prt,dummy,F_tol)
plot!( ZZ[1,:],ZZ[2,:], marker = (:circle, 4))

x_sol, y_sol,iter, nfs, val, ng,ZZ  = EGM_eta([-2.0],[4.0],obj,dummy,dummy,max_it,prt,dummy,F_tol)
plot!( ZZ[1,:],ZZ[2,:], marker = (:circle, 4))

x_sol, y_sol,iter, nfs, val, ng,ZZ  = EGM_eta([0.0],[-4.0],obj,dummy,dummy,max_it,prt,dummy,F_tol)
plot!( ZZ[1,:],ZZ[2,:], marker = (:circle, 4))

x_sol, y_sol,iter, nfs, val, ng,ZZ  = EGM_eta([0.0],[4.0],obj,dummy,dummy,max_it,prt,dummy,F_tol)
plot!( ZZ[1,:],ZZ[2,:], marker = (:circle, 4))
x_sol, y_sol,iter, nfs, val, ng,ZZ  = EGM_eta([2.0],[-4.0],obj,dummy,dummy,max_it,prt,dummy,F_tol)
plot!( ZZ[1,:],ZZ[2,:], marker = (:circle, 4))

x_sol, y_sol,iter, nfs, val, ng,ZZ  = EGM_eta([2.0],[4.0],obj,dummy,dummy,max_it,prt,dummy,F_tol)
plot!( ZZ[1,:],ZZ[2,:], marker = (:circle, 4))

x_sol, y_sol,iter, nfs, val, ng,ZZ  = EGM_eta([4.0],[-2.0],obj,dummy,dummy,max_it,prt,dummy,F_tol)
plot!( ZZ[1,:],ZZ[2,:], marker = (:circle, 4))

x_sol, y_sol,iter, nfs, val, ng,ZZ  = EGM_eta([4.0],[0.0],obj,dummy,dummy,max_it,prt,dummy,F_tol)
plot!( ZZ[1,:],ZZ[2,:], marker = (:circle, 4))

x_sol, y_sol,iter, nfs, val, ng,ZZ  = EGM_eta([4.0],[2.0],obj,dummy,dummy,max_it,prt,dummy,F_tol)
plot!( ZZ[1,:],ZZ[2,:], marker = (:circle, 4))


#sF= @sprintf(", cond(∇F)= %2.1f",nbFc)
#sL= @sprintf(", cond(∇L)= %2.1f",nbLc)

plot!(
#title = string(funcName, ", ", TYPE,", A = $interaction "), titlefont=font(8),
legend = false, xlims=(xmin,xmax), ylims=(ymin,ymax),aspect_ratio=:equal)
#,aspect_ratio=:equal

savefig("plots/$prob_type-EGM-$interaction.png")
