# methods with dashed lines at the peaks
"""
    plot_DRT(peak_times, peak_amps, τs, γs)

Visualises the output of the `compute_DRT` function.

Keyword arguments:
    - `lab::String="DRT"`: label for the plot.
    - `title::String="Distribution of relaxation times"`: The title of the plot.
    - `color::Symbol=:auto`: The color of the DRT plot.
    - `style::Symbol=:solid` : The linestyle of the DRT plot.
"""
function plot_DRT(peak_times, peak_amps, τs, γs; lab="DRT",
        title_ = "Distribution of relaxation times",color = :auto, style= :solid)
    fig = Plots.plot(τs,γs,xaxis=:log,label = lab,ylabel = "γ (Ω)",xlabel = "τ (s)", title = title_,linecolor=color,linestyle= style)
    for (time,amp) in zip(peak_times,peak_amps)
        Plots.plot!([time,time],[0,amp],linestyle=:dash,linewidth = 0.5, linecolor=:black, label=nothing)
    end
    return fig
end

function plot_DRT!(peak_times,peak_amps,τs,γs;lab="DRT",color = :auto, style= :solid)
    fig = Plots.plot!(τs,γs,xaxis=:log,label = lab,ylabel = "γ (Ω)",xlabel = "τ (s)",linecolor=color,linestyle= style)
    for (time,amp) in zip(peak_times,peak_amps)
        Plots.plot!([time,time],[0,amp],linestyle=:dash,linewidth = 0.5, linecolor=:black, label=nothing)
    end
    return fig
end

function plot_DRT(τs,γs;lab="DRT",title_ = "Distribution of relaxation times",color = :auto, style= :solid)
    fig = Plots.plot(τs,γs,xaxis=:log,label = lab,ylabel = "γ (Ω)",xlabel = "τ (s)", title = title_,linecolor=color,linestyle= style)
    return fig
end

function plot_DRT!(τs,γs;lab="DRT",color = :auto,style = :solid)
    fig = Plots.plot!(τs,γs,xaxis=:log,label = lab,ylabel = "γ (Ω)",xlabel = "τ (s)",linecolor=color,linestyle= style)
    return fig
end
