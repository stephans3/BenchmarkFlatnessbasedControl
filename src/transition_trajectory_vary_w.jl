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
begin
    fig1 = Figure(size=(800,600), fontsize=26)
    ax1 = Axis(fig1[1, 1], xlabel = "Time t in [s]", ylabel= L"\text{Transition } \Phi_{\omega,T}", 
            xlabelsize = 30,  ylabelsize = 30,
            xgridstyle = :dash, ygridstyle = :dash, 
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
    axislegend(; position = :lt, backgroundcolor = (:grey90, 0.1), labelsize=30);
    fig1
    save("results/figures/"*"transition_trajectory_orig.pdf", fig1, pt_per_unit = 1)    
end


begin
    fig2 = Figure(size=(800,600), fontsize=26)
    ax1 = Axis(fig2[1, 1], xlabel = "Time t in [s]", ylabel=L"\text{Derivative  } \frac{d}{dt} \Phi_{\omega,T}", 
            xlabelsize = 30,  ylabelsize = 30,
            xgridstyle = :dash, ygridstyle = :dash, 
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
    axislegend(; position = :lt, backgroundcolor = (:grey90, 0.1), labelsize=30);
    fig2
    save("results/figures/"*"transition_trajectory_deriv.pdf", fig2, pt_per_unit = 1)        
end

