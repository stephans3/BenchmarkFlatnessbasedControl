#=
    Simulation of a 1D heat condition model with flatness-based control input
    Physical properties: aluminum
=#

L = 0.2; # Length of 1D rod

# Aluminium
λ = 237;  # Thermal conductivity
ρ = 2700; # Density
c = 900;  # Specific heat capacity
α = λ / (ρ * c) # Diffusivity


using DelimitedFiles
w = 2.0
path_2_file = string("results/u_input/"*"u_aluminum_w_" * string(round(Int64, w*10)) * ".txt")
u_data = readdlm(path_2_file, '\t', BigFloat, '\n')[:]

# path_2_file = string("u_aluminum.txt")
# u_data = readdlm(path_2_file, '\t', BigFloat, '\n')[:]
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
const Tf = 1000.0;  # Final time
const ts = 1.0;     # Time step width

# Simulation without optimization
using OrdinaryDiffEq

x0 = 300 * ones(Nx) # Intial values
tspan = (0.0, Tf)   # Time span
prob = ODEProblem(heat_eq!,x0,tspan)

alg = KenCarp4()    # Numerical integrator
sol = solve(prob,alg, saveat = ts)

# alg = Euler()    # Numerical integrator
# sol = solve(prob,alg,dt=0.2*(Δx^2/α), saveat = ts)

using CairoMakie
tgrid = sol.t;

begin
    fig1 = Figure(size=(800,600), fontsize=26)
    ax1 = Axis(fig1[1, 1], xlabel = "Time t in [s]", ylabel = "Temperature in [K]", 
        xlabelsize = 30,  ylabelsize = 30,
        xgridstyle = :dash, ygridstyle = :dash, 
        xtickalign = 1., xticksize = 10, 
        xminorgridvisible = true, xminorticksvisible = true, xminortickalign = 1,
        yminorgridvisible = true, yminorticksvisible = true, yminortickalign = 1,
        ytickalign = 1, yticksize = 10, xlabelpadding = 0)
    
    ax1.xticks = 0 : 100 : Tf;
    ax1.yticks = 300 : 20 : 460;
    lines!(tgrid, sol[26,:]; linestyle = :dot,linewidth = 3, label = "x=0.05 m")
    lines!(tgrid, sol[51,:]; linestyle = :dash,linewidth = 3, label = "x=0.1 m (center)")
    lines!(tgrid, sol[end,:];linewidth = 3, label = "x=0.2 m (right side)")
    axislegend(; position = :lt, backgroundcolor = (:grey90, 0.1), labelsize=30);
    fig1
    save("results/figures/"*"he_aluminum.pdf", fig1, pt_per_unit = 1)    
end

