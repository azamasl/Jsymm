using CPUTime, JLD2, LinearAlgebra, Convex, Random,Plots
include("../sp_function.jl")

#"general settings"
ran_seed= 11
prob_type= "UCAC" # Analytic center problem
n,m= 500,100;
ts = m+n #total size
max_it = 2000#.5*1e3
prt = 1 # 0 don't print grad norm at every iterations; 1, do it.
F_tol =1e-4;
dummy=0
#"local alg settings"
stepsize = 0.02#The fixed step_size, only used when do_ls = 0
#"global alg settings"
my_tr  = 0.05  #Initial Δ, i.e. initial radius
max_tr = 1.0   #Maximum allowed Δ
eta = 0.001    #Minimum required decrease in ρ (ρ is the ratio of actual decrease vs predicated decrease.)
e = 0.01;
scale = 1;
stepsize_tol = 1e-6 #"Stepsize smaller than 1e-6 is too small"

function main()
    Random.seed!(ran_seed);
    #creating a feasible random problem instance
    A =random_scaled(m,n,scale); xstar = randn(n);b  = randn(n);
    x0 = 10*rand(n); y0 = A*x0-b .+ e; w0 = zeros(m,1); # dual var
    ##find negatives in y0 and set them to 1
    indices = findall(t -> t.<0, y0);y0[indices] .=1;
    @show minimum(abs.(y0))
    TYPE="UnConstrained Analytic Center Problem"#w domain =R++
    nbF=Matrix([ zeros(n,n)     zeros(n,m)             A';
             zeros(m,n)   Diagonal(1 ./ (y0 .^2))    I(m);
              -A             -I(m)             zeros(m,m)]);

    @show nbFc0 = cond(Array(nbF),2)
    @show nbFc = myCond(nbF)

    filename = "$prob_type-$ts-rs$ran_seed"
    display(filename)
    @save "output/problem-"*filename*".jld2"  TYPE prob_type ran_seed max_it F_tol stepsize m n scale  b A x0 y0 w0 xstar dummy  prt ts my_tr max_tr eta nbFc

end

main()
