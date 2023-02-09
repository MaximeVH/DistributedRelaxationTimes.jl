function plot_DRT(taus,gammas;lab="DRT",title_ = "Distribution of relaxation times",color = :auto)
    fig = Plots.plot(taus,gammas,xaxis=:log,label = lab,ylabel = "γ (Ω)",xlabel = "τ (s)", title = title_,linecolor=color)
    return fig
end

function plot_DRT!(taus,gammas;lab="DRT",color = :auto)
    fig = Plots.plot!(taus,gammas,xaxis=:log,label = lab,ylabel = "γ (Ω)",xlabel = "τ (s)",linecolor=color)
    return fig
end

# methods with dashed lines at the peaks

function plot_DRT(peak_times,peak_amps,taus,gammas;lab="DRT",title_ = "Distribution of relaxation times",color = :auto)
    fig = Plots.plot(taus,gammas,xaxis=:log,label = lab,ylabel = "γ (Ω)",xlabel = "τ (s)", title = title_,linecolor=color)
    for (time,amp) in zip(peak_times,peak_amps)
        Plots.plot!([time,time],[0,amp],linestyle=:dash,linewidth = 0.5, linecolor=:black, label=nothing)
    end
    return fig
end

function plot_DRT!(peak_times,peak_amps,taus,gammas;lab="DRT",color = :auto)
    fig = Plots.plot!(taus,gammas,xaxis=:log,label = lab,ylabel = "γ (Ω)",xlabel = "τ (s)",linecolor=color)
    for (time,amp) in zip(peak_times,peak_amps)
        Plots.plot!([time,time],[0,amp],linestyle=:dash,linewidth = 0.5, linecolor=:black, label=nothing)
    end
    return fig
end