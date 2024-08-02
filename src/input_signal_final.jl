#=
- Compute input signal u(t) for aluminum and steel 38Si7
=#

L = 0.2; # Length of 1D rod

# Aluminium
λa = 237;  # Thermal conductivity
ρa = 2700; # Density
ca = 900;  # Specific heat capacity
α_a = λa / (ρa * ca) # Diffusivity
γ_a = L^2 / α_a

# Steel 38Si7
λs = 40;   # Thermal conductivity
ρs = 7800; # Density
cs = 460;  # Specific heat capacity
α_s = λs / (ρs * cs) # Diffusivity
γ_s = L^2 / α_s

η(L,α,i) = BigFloat(L)^(2i+1) / (BigFloat(α)^(i+1) * factorial(big(2i+1)))

idx_grid = 0:40;
eta_al = zeros(BigFloat, length(idx_grid))
eta_st = zeros(BigFloat, length(idx_grid))

for (n, iter) in enumerate(idx_grid)
    eta_al[n] = η(L,α_a,iter)
    eta_st[n] = η(L,α_s,iter)
end


w =  2.0;
using FastGaussQuadrature
T = 1000;
bump(t) = exp(-1 / (t/T - (t/T)^2)^w)
t_gq, weights_gq = FastGaussQuadrature.gausslegendre(1000)
p = T/2;
Ω_int = p *FastGaussQuadrature.dot( weights_gq ,bump.(p*t_gq .+ p))

diff_ref = 100; # (y_f - y_0) = 100 Kelvin

using DelimitedFiles
T = 1000; 
path_2_file = string("results/h_results/h_results_T_", T, "_w_", round(Int,w*10), ".txt")
h_data = readdlm(path_2_file, '\t', BigFloat, '\n')

u_al = λa*(diff_ref/Ω_int) * (h_data * eta_al)
u_st = λs*(diff_ref/Ω_int) * (h_data * eta_st)

path_al = string("u_aluminum_w_", round(Int,w*10))
path_st = string("u_steel_w_", round(Int,w*10))
path_u_al = "results/u_input/"*path_al*".txt"
path_u_st = "results/u_input/"*path_st*".txt"

open(path_u_al, "w") do io
    writedlm(io,  u_al)
end

open(path_u_st, "w") do io
    writedlm(io, u_st)
end


using CairoMakie
J = 1000;   
tgrid = T/J : T/J : T-T/J;

begin
    fig1 = Figure(size=(800,600), fontsize=26)
    ax1 = Axis(fig1[1, 1], xlabel = "Time t in [s]", ylabel = L"Input u $\times 10^5$", 
        xlabelsize = 30, ylabelsize = 30,  xgridstyle = :dash, ygridstyle = :dash, 
        xtickalign = 1., xticksize = 10, 
        xminorgridvisible = true, xminorticksvisible = true, xminortickalign = 1,
        yminorgridvisible = true, yminorticksvisible = true, yminortickalign = 1,
        ytickalign = 1, yticksize = 10, xlabelpadding = 0)
    
    u_al_scaled = u_al * 1e-5;
    ax1.xticks = 0 : 100 : T;
    ax1.yticks = 0 : 1 : 5;
    lines!(tgrid, u_al_scaled;linewidth = 3)
    fig1
    save("results/figures/"*path_al*".pdf", fig1, pt_per_unit = 1)
end

begin
    fig2 = Figure(size=(800,600), fontsize=26)
    ax2 = Axis(fig2[1, 1], xlabel = "Time t in [s]", ylabel = L"Input u $\times 10^6$", 
        xlabelsize = 30, ylabelsize = 30,
        xgridstyle = :dash, ygridstyle = :dash, 
        xtickalign = 1., xticksize = 10, 
        xminorgridvisible = true, xminorticksvisible = true, xminortickalign = 1,
        yminorgridvisible = true, yminorticksvisible = true, yminortickalign = 1,
        ytickalign = 1, yticksize = 10, xlabelpadding = 0)
    

    u_st_scaled = u_st * 1e-6;
    ax2.xticks = 0 : 100 : T;
    ax2.yticks = -8 : 2 : 8;
    lines!(tgrid, u_st_scaled;linewidth = 3)
    fig2
    
    save("results/figures/"*path_st*".pdf", fig2, pt_per_unit = 1)    
end
