#=
- Compute sequence μ_i and ratio μ_i / max μ_j for j ∈ {1,...,i} 
- Compute input signal u(t) for aluminum and steel 38Si7 for variable final iteration N
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

using FastGaussQuadrature
w = 2;
T = 1000;
bump(t) = exp(-1 / (t/T - (t/T)^2)^w)
t_gq, weights_gq = FastGaussQuadrature.gausslegendre(1000)
p = T/2;
Ω_int = p *FastGaussQuadrature.dot( weights_gq ,bump.(p*t_gq .+ p))

diff_ref = 100; # (y_f - y_0) = 100 Kelvin

using DelimitedFiles
T = 1000; 
path_2_file = string("results/h_results/h_results_T_", T, ".txt")
h_data = readdlm(path_2_file, '\t', BigFloat, '\n')


u_al_all = λa*(diff_ref/Ω_int) * (h_data' .* eta_al)'
u_st_all = λs*(diff_ref/Ω_int) * (h_data' .* eta_st)'

u_al_sum = similar(u_al_all)
u_st_sum = similar(u_st_all)

u_al_norm_2 = zeros(size(u_al_all)[2])
u_st_norm_2 = zeros(size(u_st_all)[2])

using LinearAlgebra
for i in axes(u_al_all)[2]
    # Compute sum for input signal
    u_al_sum[:,i] = sum(u_al_all[:,1:i], dims=2)
    u_st_sum[:,i] = sum(u_st_all[:,1:i], dims=2)
    
    # Compute norm of input signal
    u_al_norm_2[i] = norm(u_al_all[:,i],2)
    u_st_norm_2[i] = norm(u_st_all[:,i],2)
end

u_al_norm_2_log10 = log10.(u_al_norm_2)
u_st_norm_2_log10 = log10.(u_st_norm_2)

# Ratio μ_i / max μ_j for j ∈ {1,...,i} 
# -> to find maximum iteration number
norm_rel_al = similar(u_al_norm_2)
norm_rel_st = similar(u_st_norm_2)

for (idx, el) in enumerate(u_al_norm_2)
    norm_rel_al[idx] = el / maximum(u_al_norm_2[begin:idx])
end

for (idx, el) in enumerate(u_st_norm_2)
    norm_rel_st[idx] = el / maximum(u_st_norm_2[begin:idx])
end


using CairoMakie
J = 1000;   
tgrid = T/J : T/J : T-T/J;

fig1 = Figure(fontsize=12)
ax1 = Axis(fig1[1, 1], xlabel = "Time t in [s]", ylabel = L"\text{Input } u_{N}", ylabelsize = 22,
    xlabelsize = 24, xgridstyle = :dash, ygridstyle = :dash, 
    xtickalign = 1., xticksize = 10, 
    xminorgridvisible = true, xminorticksvisible = true, xminortickalign = 1,
    yminorgridvisible = true, yminorticksvisible = true, yminortickalign = 1,
    ytickalign = 1, yticksize = 10, xlabelpadding = 0)

ax1.xticks = 0 : 100 : T;
#ax1.yticks = -40 : 10 : 40;
lines!(tgrid, u_al_sum[:,1]; linestyle = :dash,linewidth = 3, label = "N=1")
lines!(tgrid, u_al_sum[:,3];linestyle = :dashdot, linewidth = 3, label = "N=3")
lines!(tgrid, u_al_sum[:,7];linewidth = 3, label = "N=7")
axislegend(; position = :lt, bgcolor = (:grey90, 0.1));
fig1
save("results/figures/"*"u_input_aluminum_progress.pdf", fig1, pt_per_unit = 1)


fig2 = Figure(fontsize=12)
ax2 = Axis(fig2[1, 1], xlabel = "Time t in [s]", ylabel = L"\text{Input } u_{N}", ylabelsize = 22,
    xlabelsize = 24, xgridstyle = :dash, ygridstyle = :dash, 
    xtickalign = 1., xticksize = 10, 
    xminorgridvisible = true, xminorticksvisible = true, xminortickalign = 1,
    yminorgridvisible = true, yminorticksvisible = true, yminortickalign = 1,
    ytickalign = 1, yticksize = 10, xlabelpadding = 0)

ax2.xticks = 0 : 100 : T;
#ax2.yticks = -10e6 : 2e6 : 8e6;
lines!(tgrid, u_st_sum[:,5]; linestyle = :dash,linewidth = 3, label = "N=5")
lines!(tgrid, u_st_sum[:,10];linestyle = :dashdot, linewidth = 3, label = "N=10")
lines!(tgrid, u_st_sum[:,15];linewidth = 3, label = "N=15")
axislegend(; position = :lt, bgcolor = (:grey90, 0.1));
fig2
save("results/figures/"*"u_input_steel_progress.pdf", fig2, pt_per_unit = 1)


fig3 = Figure(fontsize=12)
ax3 = Axis(fig3[1, 1], xlabel = "Iteration i", ylabel = L"\log_{10}(\mu_{i})", ylabelsize = 22,
    xlabelsize = 24, xgridstyle = :dash, ygridstyle = :dash, 
    xtickalign = 1., xticksize = 10, 
    xminorgridvisible = true, xminorticksvisible = true, xminortickalign = 1,
    yminorgridvisible = true, yminorticksvisible = true, yminortickalign = 1,
    ytickalign = 1, yticksize = 10, xlabelpadding = 0);
    
idx_grid = 0:size(u_al_all)[2]-1;
ax3.xticks = idx_grid[1] : 10 : idx_grid[end];
ax3.yticks = -50 : 10 : 20;
lines!(idx_grid, u_al_norm_2_log10; linestyle = :dashdotdot,linewidth = 3, label = "Aluminum")
lines!(idx_grid, u_st_norm_2_log10; linestyle = :dashdotdot, linewidth = 3, label = "Steel 38Si7")
axislegend(; position = :rt, bgcolor = (:grey90, 0.1));
fig3
save("results/figures/"*"norm_mu_absolute.pdf", fig3, pt_per_unit = 1)


fig4 = Figure(fontsize=12)
ax4 = Axis(fig4[1, 1], xlabel = "Iteration i", ylabel = L" \mu_{i} / \max(\mu_{j}) \text{ , } j\in \{ 11,...,i \}", ylabelsize = 22,
    xlabelsize = 24, xgridstyle = :dash, ygridstyle = :dash, 
    xtickalign = 1., xticksize = 10, 
    xminorgridvisible = true, xminorticksvisible = true, xminortickalign = 1,
    yminorgridvisible = true, yminorticksvisible = true, yminortickalign = 1,
    ytickalign = 1, yticksize = 10, xlabelpadding = 0);
    
idx_grid = 0:size(u_al_all)[2]-1;
ax4.xticks = idx_grid[1] : 10 : idx_grid[end];
#ax4.yticks = 0 : 0.1 : 1;
lines!(idx_grid, norm_rel_al; linestyle = :dashdotdot,linewidth = 3, label = "Aluminum")
lines!(idx_grid, norm_rel_st; linestyle = :dashdotdot, linewidth = 3, label = "Steel 38Si7")
axislegend(; position = :rt, bgcolor = (:grey90, 0.1));
fig4
save("results/figures/"*"norm_mu_relative.pdf", fig4, pt_per_unit = 1)