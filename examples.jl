using PlotlyJS, DataFrames, CSV

Z_data = CSV.read("example_measurements.csv",DataFrame)
freqs = Z_data[:,3]
measurements = Z_data[:,1] .+ (Z_data[:,2] .* im)

# Nyquist plot so see what the original data looks like.
Plots.plot(real(measurements),-imag(measurements))

# Find the optimal lambda parameter using Saccoccio et al.'s Re-Imcross-validation test functions, this may take several minutes.

@time λ_opt =  find_optimal_lambda(freqs, measurements)

# DRT calculation.
relaxation_times, peak_amplitudes, taus_out, drt = compute_DRT(freqs, measurements, λ = λ_opt)

# Visualisation
plot_DRT(relaxation_times, peak_amplitudes, taus_out, drt)
vline!(relaxation_times, label = "Peaks", linestyle = :dot)

