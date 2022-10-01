function plot_DRT(taus,gammas)
    factor = 2 #zooming in on interesting part of the function.
    index = findlast(gammas .> factor*mean(gammas))
    fig = plot(taus,gammas,xlim=[0,taus[index]],label = "DRT",ylabel = "γ",xlabel = "τ", title = "Distribution of relaxation times")
    return fig
end