using  Plots, JLD2, Printf
ymin, ymax = 1E-5,2E6
xmin, xmax = -50, 2050

problem_name= "ran13x13" 

#Loading the problem
prob = load("x.jld2")
TYPE = prob["TYPE"];
F_tol = prob["F_tol"]
alph = prob["scaleBC"]
lb = Dict()
lb["egm"]  ="EGM";lb["nobro"]="Broyden";lb["lssec"]="J-symm-LS";lb["nosec"]="J-symm";lb["tr"]   ="J-symm-Tr";
A = prob["A"];B = prob["B"];C = prob["C"]
nbF = [B  A'
     -A  C]
@show norm(B);@show norm(C); @show norm(A)


###################################################
# load and print
###################################################

method_data=load("egm.jld2");sol = method_data["sol"]
plot(range(1,size(sol["nfs"],1),step=1), sol["nfs"], yscale = :log10, label = lb["egm"], color=:orange)
it = sol["iter"];ti = sol["time"];println("EGM iter=$it, time = $ti s")
sol_egm = sol["nfs"]

#@show sol_egm[end]
method_data=load("nobro.jld2");sol = method_data["sol"]
plot!(range(1,size(sol["nfs"],1),step=1), sol["nfs"], yscale = :log10, label = lb["nobro"], color=:red)#, markershape=:diamond)
it = sol["iter"]; ti = sol["time"];println("Broyden iter=$it, time = $ti s")

method_data=load("nosec.jld2");sol = method_data["sol"]
plot!(range(1,size(sol["nfs"],1),step=1), sol["nfs"], yscale = :log10, label = lb["nosec"], color=:black)#, markershape=:circle)
it = sol["iter"];ti = sol["time"];println("Jsymm iter=$it, time = $ti s")
stsz = sol["stepsize"]

method_data=load("lssec.jld2"); sol = method_data["sol"]
plot!(range(1,size(sol["nfs"],1),step=1), sol["nfs"], yscale = :log10, label = lb["lssec"], color=:purple)
it = sol["iter"];ti = sol["time"];println("Jsymm-LS iter=$it, time = $ti s")
F_tol = sol["F_tol"]

method_data=load("tr.jld2");sol = method_data["sol"]
plot!(range(1,size(sol["nfs"],1),step=1), sol["nfs"], yscale = :log10, label = lb["tr"], color = :blue)
it = sol["iter"];ti = sol["time"];println("Jsymm-Tr iter=$it, time = $ti s")
sol_tr = sol["nfs"]
#@show sol_tr[end]

sF= @sprintf(", cond(∇F)= %2.1f",sol["nbFc"])
plot!(xlabel = "Iteration", ylabel = "||F||",
title = string(problem_name), titlefont=font(12),
legend = :topright,
#legend = :bottomleft,
legendfont = font(7), ylims=(ymin,ymax), xlims=(xmin,xmax))

fig = "bilinear"*string(problem_name)*".png"
savefig(fig)
