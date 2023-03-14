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

log_eta_al = log10.(eta_al)
log_eta_st = log10.(eta_st)

# Ratio η_{i+1} / η_{i}
ratio_al = eta_al[begin+1:end] ./  eta_al[begin:end-1]
ratio_st = eta_st[begin+1:end] ./  eta_st[begin:end-1]

ratio_al_log = log10.(ratio_al)
ratio_st_log = log10.(ratio_st)

# η reaches its maximum
eta_max_al = η(L,α_a,9)
eta_max_st = η(L,α_s,29)

# η drops below 1
eta_0_al = η(L,α_a,28)
eta_0_st = η(L,α_s,83)

log10(eta_0_al)
log10(eta_0_st)

path2results = "results/figures/"

using CairoMakie
fig = Figure(fontsize=12)
ax1 = Axis(fig[1, 1], xlabel = "Iteration i", ylabel = L"\log_{10}(\eta_{i})", ylabelsize = 22,
    xlabelsize = 24, xgridstyle = :dash, ygridstyle = :dash, 
    xtickalign = 1., xticksize = 10, 
    xminorgridvisible = true, xminorticksvisible = true, xminortickalign = 1,
    yminorgridvisible = true, yminorticksvisible = true, yminortickalign = 1,
    ytickalign = 1, yticksize = 10, xlabelpadding = 0)

ax1.xticks = idx_grid[1] : 4 : idx_grid[end];
ax1.yticks = -15 : 5 : 30;
lines!(idx_grid, log_eta_al; linestyle = :dashdotdot,linewidth = 3, label = "Aluminum")
lines!(idx_grid, log_eta_st; linestyle = :dashdotdot, linewidth = 3, label = "Steel 38Si7")

axislegend(; position = :lb, bgcolor = (:grey90, 0.1));
fig
save(path2results*"sequence_eta_al_st.pdf", fig, pt_per_unit = 1)

fig2 = Figure(fontsize=12)
ax2 = Axis(fig2[1, 1], xlabel = "Iteration i", ylabel = L"\log_{10}(\eta_{i+1} / \eta_{i})", ylabelsize = 22,
    xlabelsize = 24, xgridstyle = :dash, ygridstyle = :dash, 
    xtickalign = 1., xticksize = 10, 
    xminorgridvisible = true, xminorticksvisible = true, xminortickalign = 1,
    yminorgridvisible = true, yminorticksvisible = true, yminortickalign = 1,
    ytickalign = 1, yticksize = 10, xlabelpadding = 0)

ax2.xticks = idx_grid[1] : 4 : idx_grid[end];
ax2.yticks = -1.5 : 0.5 : 3.0;
lines!(idx_grid[begin:end-1], ratio_al_log; linestyle = :dashdotdot,linewidth = 3, label = "Aluminum")
lines!(idx_grid[begin:end-1], ratio_st_log; linestyle = :dashdotdot, linewidth = 3, label = "Steel 38Si7")

axislegend(; position = :lb, bgcolor = (:grey90, 0.1));
fig2
save(path2results*"ratio_eta_next_al_st.pdf", fig2, pt_per_unit = 1)