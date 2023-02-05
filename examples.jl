using PlotlyJS, DataFrames, CSV

Z_data = CSV.read("example_measurements.csv",DataFrame)
frequencies = Z_data[:,3]
measurements = Z_data[:,1] .+ (Z_data[:,2] .* im)

# Nyquist plot so see what the original data looks like.
Plots.plot(real(measurements),-imag(measurements))

# DRT calculation.
@time relaxation_times, taus_out, drt = compute_DRT(frequencies, measurements)

# Visualisation
plot_DRT(taus_out,drt)
vline!(relaxation_times, label = "Peaks", linestyle = :dot)

