using  Plots, JLD2, Printf
ymin, ymax = 1E-10,1E5
xmin, xmax = -100, 4000

#Loading the problem
filename = "0-1000-BC-0"
filename = "1-1000-BC-0.0001"
filename = "1-1000-BC-0.01"
filename = "1-1000-BC-1"
 #filename = "1-1000-BC-100"
prob = load("output/problem-"*filename*".jld2")
TYPE = prob["TYPE"];
F_tol = prob["F_tol"]
alph = prob["scaleBC"]
lb = Dict()
lb["egm"]  ="EGM";lb["nobro"]="Broyden";lb["lssec"]="J-symm-LS";lb["nosec"]="J-symm";lb["tr"]   ="J-symm-Tr";
A = prob["A"];B = prob["B"];C = prob["C"]
nbF = [B  A'
     -A  C]
@show norm(B);@show norm(C); @show norm(A)
nbFc = myCond(nbF)
eig_valsB, eig_vecsB = eigen(B)
eig_valsC, eig_vecsC = eigen(C)
@show μB = minimum(abs.(eig_valsB))
@show μC = minimum(abs.(eig_valsC))
μ = min(μB, μC)
println("strong convexity-concavity: μ = $μ")
U,s,V = svd(nbF)# s is a vector and sorted descending
β =s[1] #@show Smax = maximum(s)
println("smothness : β = $β");ratio = μ/β; println("μ/β = $ratio") ; #compute and print  μ/β
println("conditon number = $nbFc ")


###################################################
# load and print
###################################################

method_data=load("output/"*filename*"-egm.jld2");sol = method_data["sol"]
plot(range(1,size(sol["nfs"],1),step=1), sol["nfs"], yscale = :log10, label = lb["egm"], color=:orange)
it = sol["iter"];ti = sol["time"];println("EGM iter=$it, time = $ti s")
sol_egm = sol["nfs"]

#@show sol_egm[end]
method_data=load("output/"*filename*"-nobro.jld2");sol = method_data["sol"]
plot!(range(1,size(sol["nfs"],1),step=1), sol["nfs"], yscale = :log10, label = lb["nobro"], color=:red)#, markershape=:diamond)
it = sol["iter"]; ti = sol["time"];println("Broyden iter=$it, time = $ti s")

method_data=load("output/"*filename*"-nosec.jld2");sol = method_data["sol"]
plot!(range(1,size(sol["nfs"],1),step=1), sol["nfs"], yscale = :log10, label = lb["nosec"], color=:black)#, markershape=:circle)
it = sol["iter"];ti = sol["time"];println("Jsymm iter=$it, time = $ti s")
stsz = sol["stepsize"]

method_data=load("output/"*filename*"-lssec.jld2"); sol = method_data["sol"]
plot!(range(1,size(sol["nfs"],1),step=1), sol["nfs"], yscale = :log10, label = lb["lssec"], color=:purple)
it = sol["iter"];ti = sol["time"];println("Jsymm-LS iter=$it, time = $ti s")
F_tol = sol["F_tol"]

method_data=load("output/"*filename*"-tr.jld2");sol = method_data["sol"]
plot!(range(1,size(sol["nfs"],1),step=1), sol["nfs"], yscale = :log10, label = lb["tr"], color = :blue)
it = sol["iter"];ti = sol["time"];println("Jsymm-Tr iter=$it, time = $ti s")
sol_tr = sol["nfs"]
#@show sol_tr[end]

sF= @sprintf(", cond(∇F)= %2.1f",sol["nbFc"])
plot!(xlabel = "Iteration", ylabel = "||F||",
#title = string(TYPE, " m = $m, n=$n, fixed stepsize=$stsz ",sF, ", tol=$F_tol"), titlefont=font(8),
title = string("α=$alph"), titlefont=font(12),
legend = :topright,
#legend = :bottomleft,
legendfont = font(7), ylims=(ymin,ymax), xlims=(xmin,xmax))
savefig("plots/"*filename*".png")
