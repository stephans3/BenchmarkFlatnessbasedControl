
L = 1.; # Length of 1D rod

# Material
λ = 1;  # Thermal conductivity
ρ = 1; # Density
c = 1;  # Specific heat capacity
α = λ / (ρ * c) # Diffusivity
γ = L^2 / α


η(L,α,i) = BigFloat(L)^(2i+1) / (BigFloat(α)^(i+1) * factorial(big(2i+1)))

idx_grid = 0:20;
eta = zeros(BigFloat, length(idx_grid))

for (n, iter) in enumerate(idx_grid)
    eta[n] = η(L,α,iter)
end


w =  2.0;
using FastGaussQuadrature
T = 1.0;
bump(t) = exp(-1 / (t/T - (t/T)^2)^w)
t_gq, weights_gq = FastGaussQuadrature.gausslegendre(1000)
p = T/2;
Ω_int = p *FastGaussQuadrature.dot( weights_gq ,bump.(p*t_gq .+ p))

diff_ref = 1; # (y_f - y_0) = 100 Kelvin

using DelimitedFiles
path_2_file = string("results/h_results/h_results_T_", round(Int,T), "_w_", round(Int,w*10), ".txt")
h_data = readdlm(path_2_file, '\t', BigFloat, '\n')

# Simulation
u_data = λ*(diff_ref/Ω_int) * (h_data * eta)
u_data = vcat(0, u_data, 0)

function input_signal(t,dt)
    if t <= 0
        return u_data[1]
    elseif t >= Tf
        return u_data[end]
    end
    τ = t/dt + 1
    t0 = floor(Int, τ)
    t1 = t0 + 1;

    u0 = u_data[t0]
    u1 = u_data[t1]

    a = u1-u0;
    b = u0 - a*t0

    return a*τ + b;
end


# Diffusion: x-direction
function diffusion_x!(dx,x,Nx, Ny, Δx) # in-place
    
    for iy in 1 : Ny
        for ix in 2 : Nx-1
            i = (iy-1)*Nx + ix
            dx[i] =  (x[i-1] - 2*x[i] + x[i+1])/Δx^2
        end
        i1 = (iy-1)*Nx + 1      # West
        i2 = (iy-1)*Nx + Nx     # East
        dx[i1] = (-2*x[i1] + 2*x[i1+1])/Δx^2
        dx[i2] = (2*x[i2-1] - 2*x[i2])/Δx^2
    end

    nothing 
end


# 1D heat equation
function heat_eq!(dx,x,p,t)       
    # time = t/Tf;
    #u = input_signal(time, p)
    
    u = input_signal(t,ts)
    
    diffusion_x!(dx,x,Nx,1,Δx)
  
    dx .= α * dx
    dx[1] = dx[1] + 2α/(λ * Δx) * u
end


# Discretization  
const Nx = 101;     # Number of elements x-direction
const Δx = L/(Nx-1) # Spatial sampling
const Tf = T;  # Final time
const ts = T / (length(u_data)-1) # 1.0;     # Time step width

# Simulation without optimization
using OrdinaryDiffEq

x0 = zeros(Nx) # 300 * ones(Nx) # Intial values
tspan =  (0.0, Tf)   # Time span

prob = ODEProblem(heat_eq!,x0,tspan)
#alg = KenCarp4()    # Numerical integrator
#sol = solve(prob,alg, saveat = dt)

alg = Euler()    # Numerical integrator
sol = solve(prob,alg,dt=0.2*Δx^2, saveat = ts)


using CairoMakie

tgrid = sol.t;
fig1 = Figure(fontsize=12)
ax1 = Axis(fig1[1, 1], xlabel = "Time t in [s]", ylabel = "Temperature in [K]", ylabelsize = 22,
    xlabelsize = 24, xgridstyle = :dash, ygridstyle = :dash, 
    xtickalign = 1., xticksize = 10, 
    xminorgridvisible = true, xminorticksvisible = true, xminortickalign = 1,
    yminorgridvisible = true, yminorticksvisible = true, yminortickalign = 1,
    ytickalign = 1, yticksize = 10, xlabelpadding = 0)

#ax1.xticks = 0 : 100 : Tf;
#ax1.yticks = -40 : 10 : 40;
lines!(tgrid, sol[26,:]; linestyle = :dash,linewidth = 3, label = "x=0.05 m")
lines!(tgrid, sol[51,:]; linestyle = :dashdotdot,linewidth = 3, label = "x=0.1 m (center)")
lines!(tgrid, sol[end,:];linewidth = 3, label = "x=0.2 m (right side)")
axislegend(; position = :lt, bgcolor = (:grey90, 0.1));
fig1




























path_u_in = string("results/u_input/u_normalized_w_", round(Int,w*10),".txt")
path_plot_u_in = string("results/figures/u_normalized_w_", round(Int,w*10),".pdf")

#=
open(path_u_in, "w") do io
    writedlm(io,  u_al)
end
=#

using CairoMakie
J = 1000;   
tgrid = T/J : T/J : T-T/J;

fig1 = Figure(fontsize=12)
ax1 = Axis(fig1[1, 1], xlabel = "Time t in [s]", ylabel = "Input u", ylabelsize = 22,
    xlabelsize = 24, xgridstyle = :dash, ygridstyle = :dash, 
    xtickalign = 1., xticksize = 10, 
    xminorgridvisible = true, xminorticksvisible = true, xminortickalign = 1,
    yminorgridvisible = true, yminorticksvisible = true, yminortickalign = 1,
    ytickalign = 1, yticksize = 10, xlabelpadding = 0)

ax1.xticks = 0 : T/10 : T;
#ax1.yticks = -40 : 10 : 40;
lines!(tgrid, u_data;linewidth = 3)
fig1

# save(path_plot_u_in, fig1, pt_per_unit = 1)
