function plot_DRT(taus,gammas,lab="DRT",title_ = "Distribution of relaxation times")
    fig = Plots.plot(taus,gammas,xaxis=:log,label = lab,ylabel = "γ",xlabel = "τ", title = title_)
    return fig
end

function plot_DRT!(taus,gammas,lab="DRT")
    fig = Plots.plot!(taus,gammas,xaxis=:log,label = lab,ylabel = "γ",xlabel = "τ")
    return fig
end