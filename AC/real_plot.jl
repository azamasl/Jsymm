using  Plots, JLD2, Printf
import FirstOrderLp
ymin, ymax = 1e-5,1e8

tol = 1e-4;

problem_name= "ran13x13" 
TYPE = "Analytic Center Problem"
lb = Dict()

lb["lssec"]="J-symm-LS";
lb["nosec"]="J-symm";
lb["nobro"]="Broyden";
lb["egm"]="EGM"

method_data=load("egm.jld2"); sol = method_data["sol"]
plot(range(1,size(sol["nfs"],1),step=1), sol["nfs"], yscale = :log10, label = lb["egm"], color=:yellow)
it = sol["iter"];ti = sol["time"];println("EGM iter=$it, time = $ti s")

method_data=load("nobro.jld2"); sol = method_data["sol"]
plot!(range(1,size(sol["nfs"],1),step=1), sol["nfs"], yscale = :log10, label = lb["nobro"], color=:red)
it = sol["iter"];ti = sol["time"];println("Broyden iter=$it, time = $ti s")
stsz = sol["stepsize"];

method_data=load("nosec.jld2"); sol = method_data["sol"]
plot!(range(1,size(sol["nfs"],1),step=1), sol["nfs"], yscale = :log10, label = lb["nosec"], color=:black)
it = sol["iter"];ti = sol["time"];println("Jsymm iter=$it, time = $ti s")
stsz = sol["stepsize"];

method_data=load("lssec.jld2"); sol = method_data["sol"]
plot!(range(1,size(sol["nfs"],1),step=1), sol["nfs"], yscale = :log10, label = lb["lssec"], color=:purple)
it = sol["iter"];ti = sol["time"];println("Jsymm-LS iter=$it, time = $ti s")


plot!(xlabel = "Iteration", ylabel = "||F||",
title = string(problem_name), 
#title = string(TYPE, " m = $m, n=$n, fixed stepsize=$stsz ",sF, ", tol=$tol"), titlefont=font(8),
legend = :topright,legendfont = font(7),
ylims=(ymin,ymax))

fig = "AC_plots_"*string(problem_name)*".png"
savefig(fig)


