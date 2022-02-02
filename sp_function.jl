struct Saddle_Point
    # L = 0.5(x-xstar)'*B*(x-xstar) +(y-ystar)'*A*(x-xstar) -0.5(y-ystar)'*C*(y-ystar)
    B::Array{Float64,2} # PSD matrix
    A::Array{Float64,2}
    C::Array{Float64,2} # PSD matrix
    xstar::Array{Float64,1}
    ystar::Array{Float64,1}
    #∇F::Array{Float64,2} # PSD matrix
end

struct fun_return #a better name for this struc whould be alg_return
    x_sol::Array{Float64}
    y_sol::Array{Float64}
    iter::Int
    nfs::Array{Float64}
    val::Float64
    ng::Float64
end


struct ObjectiveFunction
    L::Function # objective function
    ∇xL::Function # (sub)gradient of objective
    ∇yL::Function # (sub)gradient of objective
    #∇F::Array{Float64,2}
    ∇F::Function
end

########################################################
#  2D nonconvex-nonconcave problem
########################################################

struct Obje_Fun2D
    L::Function # objective function
    ∇xL::Function # (sub)gradient of objective
    ∇yL::Function # (sub)gradient of objective
    ∇F::Function
end


"saddle-point function L(x,y)= x'Bx+y'Ax-y'Cy with (sub)gradient ∇L(x,y)."
function saddle_point_objective(sp::Saddle_Point)
    B, A,C,xstar,ystar = sp.B, sp.A, sp.C,sp.xstar, sp.ystar
    function L(x,y)
        return  0.5*(x-xstar)'*B*(x-xstar) +(y-ystar)'*A*(x-xstar) -0.5*(y-ystar)'*C*(y-ystar)
    end
    function ∇xL(x,y)
        return B*(x-xstar) + A'*(y-ystar)
    end
    function ∇yL(x,y)
        return A*(x-xstar) - C'*(y-ystar)
    end
    function ∇F(x,y)
        ∇F = [sp.B    sp.A'
              -sp.A    sp.C]
        return ∇F
    end
    return ObjectiveFunction(L,∇xL,∇yL,∇F)
end

function sad_point2D_objective(interaction)
    function L(x,y)
        return  (x.^2 .- 1).*(x.^2 .- 9) .+ interaction*x'*y .- (y.^2 .- 1).*(y.^2 .- 9)
    end
    function ∇xL(x,y)
        return 4*x.^3 .- 20*x .+ interaction*y
    end
    function ∇yL(x,y)
        return -4*y.^3 .+ 20*y .+ interaction*x
    end
    function ∇F(x,y)
        return  [12*x'*x-20    interaction
                  -interaction    12*y'*y-20]
    end
    return Obje_Fun2D(L,∇xL,∇yL,∇F)
end


########################################################
# Unconstrained Analytic center see BV page 141
########################################################

struct Analytic_Center_Problem #problem definition struct
    A::Array{Float64,2}
    b::Array{Float64,1}
end

struct unconst_Analytic_Center_ObjectiveFunGrad #return struct
    L::Function #Lagrangian
    ∇xL::Function
    ∇yL::Function
    ∇wL::Function
    ∇F::Function
    f::Function #objective function
end

function unconst_Analytic_Center_ObjeGrad(sp::Analytic_Center_Problem)
    A, b= sp.A,sp.b
    #computes the analytical center of the linear inequalities {Ax <= b},
    # f = - sum_{i=1,...,m}log(b_i - a_i^Tx)
    # L = f + w*(y+Ax-b)
    function f(y)
        #y is the slack variable
        return -sum(log.(y))
    end
    function L(x,y,w)#primal, slack and dual vars
        return   f(y) + (w'*(y+ A*x-b))[1]
    end
    function ∇xL(x,y,w)
        return  A'*w
    end
    function ∇yL(x,y,w)
        return  (-1 ./ y)+w
    end
    function ∇wL(x,y,w)
        return y+A*x-b
    end
    function ∇F(x,y,w)
        Lyy = Diagonal(1 ./ (y .^2));
        n = length(x);m = length(b);
        return Matrix([ zeros(n,n)     zeros(n,m)         A';
                        zeros(m,n)      Lyy             I(m);
                        -A             -I(m)      zeros(m,m)])
    end
    return unconst_Analytic_Center_ObjectiveFunGrad(L,∇xL,∇yL, ∇wL ,∇F,f)
end


## ###########################################################
# auxiliary functions
##############################################################

#compute the condtion number as the ratio of numerically nonzero singularvalues
function myCond(M)
    U,s,V = svd(M)# s is a vector and sorted descendingly: ORDER is important
    i= length(s)
    while s[i] < 1e-8
        i=i-1
    end
    s1=s[1];sn=s[i];conM =s1/sn;
end

function random_PSD(N, scale)
    S = randn(N, N) / (N ^ 0.25)
    S = scale*(S * S')#scaling S down to avoid overflow for larger matrices
    return S
end

function random_PD(N,scale)
    @show den =  (N ^ 0.25)
    S = randn(N, N) /den
    S = scale*(S * S')
    eig_vals, eig_vecs = eigen(S)
    #adding 1e-6 to small eigenvalues to make sure S is PD
    eig_vals[findall(eig_vals .< 1E-6)] =  eig_vals[findall(eig_vals .< 1E-6)] .+ 1E-6
    S = eig_vecs*Diagonal(eig_vals)*eig_vecs'
    return S
end

function random_normalized(M,N)
    @show den =(N ^ 0.25)
    S = randn(M, N)/den
    return S
end

function random_PD2(N)
    @show smallVal = 1#1E-3 #1E-6 #1
    S = random_normalizedv2(N,N)
    S = 0.5*(S+S')
    eig_vals, eig_vecs = eigen(S)
    @show minEig = minimum(eig_vals)
    S = S + (abs(minEig) + smallVal)*I(N)
    return S
end

function random_scaled(M,N,scale)
    S = scale*randn(M, N)
    return S
end
