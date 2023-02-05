function plot_DRT(taus,gammas)
    fig = Plots.plot(taus,gammas,xaxis=:log,label = "DRT",ylabel = "γ",xlabel = "τ", title = "Distribution of relaxation times")
    return fig
end