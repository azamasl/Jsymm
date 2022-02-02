using JLD2,Random,LinearAlgebra
include("../sp_function.jl")

#"general settings"
    ran_seed= 345#233;
    #"0: random bilinear, 1: random quadratic CC,"
    prob_type= 1
    eq_spaced= 0 # 1:creates a random problem with the equi-spaced eigenvals in the final hessian, 0: not equi spaced
    n,m =  500,500;#5,5;#
    scale = 1#e-1; #scaling down the entires of the random matrices to avoid overflow
    scaleBC = 100#e-2
    if prob_type == 0
        scaleBC = 0
    end
    ts = m+n #total size
    #"Reciprocal condition number"
    condNum = 1e7;rec_cond =1/condNum;
    max_it = 4000
    prt = 1 # 0 don't print grad norm at every iterations; 1, do it.
    F_tol =1e-8
    reset_in_pt = 0;dummy=0;
#"local alg settings"
    c1 = 1e-4 #line-search parameter
    stepsize = 0.01#The fixed step_size, only used when do_ls = 0
    stepsize_tol = 1e-8 #Stepsize smaller than 1e-8 is too small
#"Trust-region alg settings"
    my_tr  = 0.05  #Initial Δ, i.e. initial radius
    max_tr = 1.0   #Maximum allowed Δ
    eta = 0.001    #Minimum required decrease in ρ (ρ is the ratio of actual decrease vs predicated decrease.)

function main()
    #"Creating the random problme instance"
    Random.seed!(ran_seed);
    xstar = randn(n);ystar = randn(m);
    println("Uncomment the next lines to see the x^* and y^* ")
    #println("\n xstar and ystar:")
    #display(xstar);display(ystar)
    x0 = randn(n);y0 = randn(m);

    TYPE="";
    B = zeros(n,n);
    C = zeros(m,m)
    A = random_normalized(m,n);
    if prob_type ==1        #CC
        println("The Problem is strongly convex-concave quadratic\n")
        B = scaleBC*random_PD2(n);
        C = scaleBC*random_PD2(m);
        @show isposdef(B)
        @show isposdef(C)
        TYPE="Strongly convex-concave,"
    elseif prob_type ==0   #bilinear
        println("The Problem is bilinear")
        TYPE = "Bilinear,"
    else
        println("Unspecified problem type")
    end

    nbF = [B  A'
         -A  C]
    @show nbFc0 = cond(Array(nbF),2)
    @show nbFc = myCond(nbF)
    filename = "$prob_type-$ts-BC-$scaleBC"
    display(filename)
    @save "output/problem-"*filename*".jld2"  TYPE prob_type ran_seed max_it F_tol stepsize m n scale scaleBC c1 A B C xstar ystar x0 y0  dummy reset_in_pt prt condNum nbFc ts my_tr max_tr eta
end

main()
