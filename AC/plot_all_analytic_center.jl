using  Plots, JLD2, Printf
ymin, ymax = 1E-5,1E5

#Loading the problem
filename = "UCAC-600"
#filename = "UCAC-600-rs7"
prob = load("output/problem-"*filename*".jld2")
TYPE = prob["TYPE"]; tol = prob["F_tol"];
TYPE = "Analytic Center Problem"
lb = Dict()

lb["lssec"]="J-symm-LS";
lb["nosec"]="J-symm";
lb["nobro"]="Broyden";
lb["newton"]="Newton";

method_data=load("output/"*filename*"-newton.jld2"); sol = method_data["sol"]
plot(range(1,size(sol["nfs"],1),step=1), sol["nfs"], yscale = :log10, label = lb["newton"], color=:cyan)
it = sol["iter"];ti = sol["time"];println("Newton iter=$it, time = $ti s")
@show fstar =(sol["nfs"])[end]
stsz = prob["stepsize"];

#abs.(sol["nfs"] .- fstar)
method_data=load("output/"*filename*"-nobro.jld2"); sol = method_data["sol"]
plot!(range(1,size(sol["nfs"],1),step=1), sol["nfs"], yscale = :log10, label = lb["nobro"], color=:red)
it = sol["iter"];ti = sol["time"];println("Broyden iter=$it, time = $ti s")
stsz = sol["stepsize"];

method_data=load("output/"*filename*"-nosec.jld2"); sol = method_data["sol"]
plot!(range(1,size(sol["nfs"],1),step=1), sol["nfs"], yscale = :log10, label = lb["nosec"], color=:black)
it = sol["iter"];ti = sol["time"];println("Jsymm iter=$it, time = $ti s")
stsz = sol["stepsize"];

method_data=load("output/"*filename*"-lssec.jld2"); sol = method_data["sol"]
plot!(range(1,size(sol["nfs"],1),step=1), sol["nfs"], yscale = :log10, label = lb["lssec"], color=:purple)
it = sol["iter"];ti = sol["time"];println("Jsymm-LS iter=$it, time = $ti s")

sF= @sprintf(", cond(∇F)= %2.1f",prob["nbFc"])
plot!(xlabel = "Iteration", ylabel = "||F||",
#title = string(TYPE, " m = $m, n=$n, fixed stepsize=$stsz ",sF, ", tol=$tol"), titlefont=font(8),
legend = :topright,legendfont = font(7),
ylims=(ymin,ymax))

savefig("plots/"*filename*".png")
