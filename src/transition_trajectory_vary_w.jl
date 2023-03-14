using FastGaussQuadrature

t_gq, weights_gq = FastGaussQuadrature.gausslegendre(10000)

function bump_fun(t; T=1, w=2)

    if (t<=0) || (t>=T)
        return 0;
    else
        return exp(-1 / (t/T - (t/T)^2)^w)
    end
end

function transition_fun(t; T=1, w=2)
    if t<= 0
        return 0;
    elseif t >= T
        return 1;
    else
        q = t/2;
        p = T/2;
        Ω_num = q *FastGaussQuadrature.dot( weights_gq ,bump_fun.(q*t_gq .+ q,T=T, w=w))
        Ω_den = p *FastGaussQuadrature.dot( weights_gq ,bump_fun.(p*t_gq .+ p,T=T, w=w))
        return Ω_num / Ω_den 
    end

end

Tf = 100;
J = 1000;
tgrid = 0.0 : Tf/J : Tf;
p = Tf/2;
w_vec = [11,15,20,25,30]
# Transition
Φ = zeros(length(tgrid), length(w_vec))
d1Φ = zeros(length(tgrid), length(w_vec))

for (idx, w) in enumerate(w_vec)
    # Original function
    Φ[:,idx] = transition_fun.(tgrid,T=Tf,w=0.1*w) 

    # First derivative
    Ω_den = p *FastGaussQuadrature.dot( weights_gq ,bump_fun.(p*t_gq .+ p,T=Tf, w=0.1*w))
    Ω_num = bump_fun.(tgrid,T=Tf, w=0.1*w)
    d1Φ[:,idx] = Ω_num / Ω_den
end

using CairoMakie
fig1 = Figure(fontsize=12)
ax1 = Axis(fig1[1, 1], xlabel = "Time t in [s]", ylabel= L"\text{Transition } \Phi_{\omega,T}", ylabelsize = 22,
    xlabelsize = 24, xgridstyle = :dash, ygridstyle = :dash, 
    xtickalign = 1., xticksize = 10, 
    xminorgridvisible = true, xminorticksvisible = true, xminortickalign = 1,
    yminorgridvisible = true, yminorticksvisible = true, yminortickalign = 1,
    ytickalign = 1, yticksize = 10, xlabelpadding = 0)

ax1.xticks = 0 : 10 : Tf;
ax1.yticks = 0. : 0.25 : 1;
lines!(tgrid, Φ[:,1];linestyle = :dot, linewidth = 3, label = L"\omega=1.1")
lines!(tgrid, Φ[:,2];linestyle = :dash, linewidth = 3, label=L"\omega=1.5")
lines!(tgrid, Φ[:,3];linewidth = 3, label=L"\omega=2.0")
lines!(tgrid, Φ[:,4];linestyle = :dash, linewidth = 3, label=L"\omega=2.5")
lines!(tgrid, Φ[:,5];linestyle = :dot, linewidth = 3, label=L"\omega=3.0")
axislegend(; position = :lt, bgcolor = (:grey90, 0.1));
fig1
save("results/figures/"*"transition_trajectory_orig.pdf", fig1, pt_per_unit = 1)

fig2 = Figure(fontsize=12)
ax1 = Axis(fig2[1, 1], xlabel = "Time t in [s]", ylabel=L"\text{Derivative  } \frac{d}{dt} \Phi_{\omega,T}", ylabelsize = 22,
    xlabelsize = 24, xgridstyle = :dash, ygridstyle = :dash, 
    xtickalign = 1., xticksize = 10, 
    xminorgridvisible = true, xminorticksvisible = true, xminortickalign = 1,
    yminorgridvisible = true, yminorticksvisible = true, yminortickalign = 1,
    ytickalign = 1, yticksize = 10, xlabelpadding = 0)

ax1.xticks = 0 : 10 : Tf;
#ax1.yticks = 0. : 0.25 : 1;
lines!(tgrid, d1Φ[:,1];linestyle = :dot, linewidth = 3, label = L"\omega=1.1")
lines!(tgrid, d1Φ[:,2];linestyle = :dash, linewidth = 3, label=L"\omega=1.5")
lines!(tgrid, d1Φ[:,3];linewidth = 3, label=L"\omega=2.0")
lines!(tgrid, d1Φ[:,4];linestyle = :dash, linewidth = 3, label=L"\omega=2.5")
lines!(tgrid, d1Φ[:,5];linestyle = :dot, linewidth = 3, label=L"\omega=3.0")
axislegend(; position = :lt, bgcolor = (:grey90, 0.1));
fig2
save("results/figures/"*"transition_trajectory_deriv.pdf", fig2, pt_per_unit = 1)



# 20. derivative
p = Tf/2;
Ω_den = p *FastGaussQuadrature.dot( weights_gq ,bump_fun.(p*t_gq .+ p,T=Tf, w=2.0))
Ω_num = bump_fun.(tgrid,T=Tf, w=2.0)

d1_Φ = Ω_num / Ω_den


using DelimitedFiles

trans = zeros(J-1,0)

for (idx, w) in enumerate(w_vec)
    path_2_file = string("results/h_results/h_results_T_1000_w_", w, ".txt")
    h_data = readdlm(path_2_file, '\t', BigFloat, '\n')
    Ω_den = p *FastGaussQuadrature.dot( weights_gq ,bump_fun.(p*t_gq .+ p,T=Tf, w=0.1*w))
    trans = hcat(trans, h_data[:,1]/Ω_den)
end

trans = vcat(zeros(1, length(w_vec)), trans,zeros(1, length(w_vec)))

using CairoMakie
fig1 = Figure(fontsize=12)
ax1 = Axis(fig1[1, 1], xlabel = "Time t in [s]", ylabelsize = 22,
    xlabelsize = 24, xgridstyle = :dash, ygridstyle = :dash, 
    xtickalign = 1., xticksize = 10, 
    xminorgridvisible = true, xminorticksvisible = true, xminortickalign = 1,
    yminorgridvisible = true, yminorticksvisible = true, yminortickalign = 1,
    ytickalign = 1, yticksize = 10, xlabelpadding = 0)

ax1.xticks = 0 : 10 : Tf;
#ax1.yticks = 0. : 0.25 : 1;
lines!(tgrid, trans[:,1];linewidth = 3, label = L"\omega=1.1")
lines!(tgrid, trans[:,2];linestyle = :dash, linewidth = 3, label=L"\omega=1.5")
lines!(tgrid, trans[:,3];linestyle = :dash, linewidth = 3, label=L"\omega=2.0")
lines!(tgrid, trans[:,4];linestyle = :dash, linewidth = 3, label=L"\omega=2.5")
lines!(tgrid, trans[:,5];linestyle = :dash, linewidth = 3, label=L"\omega=3.0")
axislegend(; position = :lt, bgcolor = (:grey90, 0.1));
fig1


path_2_file = string("results/h_results/h_results_T_", Tf, ".txt")
h_data = readdlm(path_2_file, '\t', BigFloat, '\n')
dΩ_19 = h_data[:,20]
d20_Φ = dΩ_19 / Ω_den


# using Plots
# plot(tgrid, Φ )


# Transition
using CairoMakie
fig1 = Figure(fontsize=12)
ax1 = Axis(fig1[1, 1], xlabel = "Time t in [s]", ylabelsize = 22,
    xlabelsize = 24, xgridstyle = :dash, ygridstyle = :dash, 
    xtickalign = 1., xticksize = 10, 
    xminorgridvisible = true, xminorticksvisible = true, xminortickalign = 1,
    yminorgridvisible = true, yminorticksvisible = true, yminortickalign = 1,
    ytickalign = 1, yticksize = 10, xlabelpadding = 0)

ax1.xticks = 0 : 10 : Tf;
ax1.yticks = 0. : 0.25 : 1;
lines!(tgrid, Φ;linewidth = 3, label = L"\text{Transition } \Phi_{\omega,T}")
lines!(tgrid, d1_Φ;linestyle = :dash, linewidth = 3, label=L"\text{1. Derivative } \frac{d}{dt} \Phi_{\omega,T}")
axislegend(; position = :lt, bgcolor = (:grey90, 0.1));
fig1
save("transition_trajectory.pdf", fig1, pt_per_unit = 1)

# 20. derivative
fig2 = Figure(fontsize=12)
ax2 = Axis(fig2[1, 1], xlabel = "Time t in [s]", ylabelsize = 22,
    xlabelsize = 24, xgridstyle = :dash, ygridstyle = :dash, 
    xtickalign = 1., xticksize = 10, 
    xminorgridvisible = true, xminorticksvisible = true, xminortickalign = 1,
    yminorgridvisible = true, yminorticksvisible = true, yminortickalign = 1,
    ytickalign = 1, yticksize = 10, xlabelpadding = 0)

tgrid2 = T/J : T/J : T-T/J;
ax2.xticks = 0 : 10 : T;
#ax2.yticks = -1e-10 : 1.0e-12 : 1.e-10;
lines!(tgrid2, d20_Φ;linewidth = 3, label=L"\text{Derivative } \Phi_{\omega,T}^{(20)}")
axislegend(; position = :lt, bgcolor = (:grey90, 0.1));
fig2
save("transition_derivative_20.pdf", fig2, pt_per_unit = 1)