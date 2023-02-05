function plot_DRT(taus,gammas)
    fig = plot(taus,gammas,xlim=[0,taus[index]],xaxis=:log,label = "DRT",ylabel = "γ",xlabel = "τ", title = "Distribution of relaxation times")
    return fig
end