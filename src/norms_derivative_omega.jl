# Create Figure: Logarithmic norms of d/dt Ω(t) with respect to T and ω
using DelimitedFiles, LinearAlgebra, FastGaussQuadrature
t_gq, weights_gq = FastGaussQuadrature.gausslegendre(10000)

function bump_fun(t; T=1, w=2)

    if (t<=0) || (t>=T)
        return 0;
    else
        return exp(-1 / (t/T - (t/T)^2)^w)
    end
end

# Fixed ω=2 and variable T
T_set = [10, 100, 1000];
h_norm_T = zeros(41,length(T_set))
w = 2.0

for (idx, T) in enumerate(T_set)
    p = T/2;
    Ω_den = p *FastGaussQuadrature.dot( weights_gq ,bump_fun.(p*t_gq .+ p,T=T, w=w))
    
    path_2_file = string("results/h_results/h_results_T_", T, ".txt")
    h_data = readdlm(path_2_file, '\t', BigFloat, '\n')
    h_data = h_data[:,1:41]/Ω_den
    h_norm_T[:,idx] = log10.( mapslices(x->norm(x,2), h_data[2:end-1,:], dims=1))
end


# Fixed T=1000 and variable ω
w_set = [11, 15, 20, 25, 30]
h_norm_w = zeros(41,length(w_set))
T = 1000;
p = T/2;

for (idx, w) in enumerate(w_set)
    Ω_den = p *FastGaussQuadrature.dot( weights_gq ,bump_fun.(p*t_gq .+ p,T=T, w=0.1*w))

    path_2_file = string("results/h_results/h_results_T_1000_w_", w, ".txt")
    h_data = readdlm(path_2_file, '\t', BigFloat, '\n')
    h_data = h_data[:,1:41]/Ω_den
    h_norm_w[:,idx] = log10.( mapslices(x->norm(x,2), h_data[2:end-1,:], dims=1))
end


using CairoMakie

begin
    fig = Figure(size=(800,600), fontsize=26)
    ax1 = Axis(fig[1, 1], xlabel = "Iteration i", ylabel = L"\log_{10}(||\Omega_{\omega,T}^{(i)}||/\hat{\Omega}_{\omega,T})", 
            xlabelsize = 30,  ylabelsize = 30,
            xgridstyle = :dash, ygridstyle = :dash, 
            xtickalign = 1., xticksize = 10, 
            xminorgridvisible = true, xminorticksvisible = true, xminortickalign = 1,
            yminorgridvisible = true, yminorticksvisible = true, yminortickalign = 1,
            ytickalign = 1, yticksize = 10, xlabelpadding = 0)

    
    idx_grid = 0:size(h_norm_T)[1]-1;    
    ax1.xticks = idx_grid[1] : 10 : idx_grid[end];
    ax1.yticks = -40 : 10 : 40;
    lines!(idx_grid, h_norm_T[:,1]; linestyle = :dot,linewidth = 3, label = "T=10")
    lines!(idx_grid, h_norm_T[:,2]; linestyle = :dash, linewidth = 3, label = "T=100")
    lines!(idx_grid, h_norm_T[:,3]; linewidth = 3, label = "T=1000")
    
    axislegend(; position = :lb, backgroundcolor = (:grey90, 0.1), labelsize=30);
    fig
    save("results/figures/"*"norm_omega_der_fix_w_2.pdf", fig, pt_per_unit = 1)
end

begin
    fig2 = Figure(size=(800,600), fontsize=26)
    ax2 = Axis(fig2[1, 1], xlabel = "Iteration i", ylabel =L"\log_{10}(||\Omega_{\omega,T}^{(i)}||/\hat{\Omega}_{\omega,T})", 
        xlabelsize = 30, ylabelsize = 30,
        xgridstyle = :dash, ygridstyle = :dash, 
        xtickalign = 1., xticksize = 10, 
        xminorgridvisible = true, xminorticksvisible = true, xminortickalign = 1,
        yminorgridvisible = true, yminorticksvisible = true, yminortickalign = 1,
        ytickalign = 1, yticksize = 10, xlabelpadding = 0)
    
    idx_grid = 0:size(h_norm_w)[1]-1;    
    ax2.xticks = idx_grid[1] : 10 : idx_grid[end];
    #ax2.yticks = -60 : 10 : 0;
    lines!(idx_grid, h_norm_w[:,1]; linestyle = :dot,linewidth = 3, label = "w=1.1")
    lines!(idx_grid, h_norm_w[:,2]; linestyle = :dash, linewidth = 3, label = "w=1.5")
    lines!(idx_grid, h_norm_w[:,3]; linewidth = 3, label = "w=2.0")
    lines!(idx_grid, h_norm_w[:,4]; linestyle = :dash, linewidth = 3, label = "w=2.5")
    lines!(idx_grid, h_norm_w[:,5]; linestyle = :dot, linewidth = 3, label = "w=3.0")
    axislegend(; position = :lb, backgroundcolor = (:grey90, 0.1), labelsize=30);
    fig2
    save("results/figures/"*"norm_omega_der_fix_T_1000.pdf", fig2, pt_per_unit = 1)
end



