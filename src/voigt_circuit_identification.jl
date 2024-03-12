
"""
    simulate_voigt_eis(relaxation_time, R_N,f)      
    a function to simulate the EIS of a single Voigt element.
    the relaxation time, R_N, and the frequency are the inputs.
    the function returns the impedance spectrum for the given frequencies.
"""
function simulate_voigt_eis(relaxation_time, R_N,f)
    return @.  R_N / (1 + 2*pi *im *f *relaxation_time)
end


"""
    joint_simulation(relaxation_times, R_Ns,f)      
    a function to simulate the EIS of a series of voigt elements with relaxation times.
    the relaxation times, R_Ns, and the frequencies are the inputs.
    the function returns the impedance spectrum for the given frequencies.
"""
function joint_simulation(relaxation_times, R_Ns,f)
    eis_subspectras = Array{Any,1}(undef, length(relaxation_times)) 
    for i in eachindex(relaxation_times)
        eis_subspectras[i] = simulate_voigt_eis(relaxation_times[i], R_Ns[i], f)
    end
    return sum(eis_subspectras)
end

# function objective_function(relaxation_times,x,frequencies,impedances)
#     model_output = joint_simulation(relaxation_times, x,frequencies)
#     return  mean((abs.(impedances - model_output).^2)./(abs.(impedances).^2 .+ abs.(model_output).^2))
# end

"""
obtain_RNs(frequencies, impedances, relaxation_times)
a function to obtain the R_N values corresponding to the peaks in the DRT.
the inputs are the frequencies, impedances, and relaxation times.
the function returns the R_N values.

"""
function obtain_RNs(frequencies, impedances, relaxation_times) #,R_0
    n = length(relaxation_times)
    impedances = impedances .- R_0
    obj = x -> sum(abs.(joint_simulation(relaxation_times, x,frequencies) .- impedances)) #.+ R_0_
    x_init = fill(10.0, n)
    upper = 1e8ones(n) # Arbitrarily large upper bound.
    lower = zeros(n)
    results = optimize(obj, lower, upper, x_init, Fminbox(BFGS()), Optim.Options(g_tol=1e-8), autodiff=:forward)
    return results.minimizer
end

"""
voigtparametertuple(R0,RNs,relaxation_times)
a function to create a NamedTuple of the Voigt parameters.
the inputs are the R_0, R_Ns, and relaxation times.
the function returns the NamedTuple of the Voigt circuit parameters.
"""
function voigtparametertuple(R0,RNs,relaxation_times)
    name_strings = ["R_0"]
    param_values = [R0]
    number = 1
    for (rn,rt) in zip(RNs,relaxation_times)
        push!(name_strings, "R_$number")
        push!(param_values, rn)
        push!(name_strings, "C_$number")
        push!(param_values, rt/rn)
        number += 1
    end
    name_symbols = Tuple(map(Symbol, name_strings))
    return NamedTuple{name_symbols}(Tuple(param_values))
end

"""
DRT_circuit_identify(impedances,frequencies)
a function to identify the Voigt circuit parameters from the impedance measurements.
the inputs are the impedances and frequencies.
the function returns the NamedTuple of the Voigt circuit parameters.
"""
function DRT_circuit_identify(impedances,frequencies)
    # step one: calculate the DRT
        R_0 = minimum(real(impedances))
        relaxation_times, peak_amplitudes, taus_out, drt = compute_DRT(frequencies,impedances)
    # step two: determine the number of peaks
        n = length(relaxation_times)
    # step there: get R_N values corresponding to the peaks
        R_Ns = obtain_RNs(frequencies, impedances, relaxation_times)       
        return voigtparametertuple(R_0,RNs,relaxation_times)
end

"""
DRT_circuit_identify_nyquist(impedances,frequencies)
A diagnostic tool to evaluate the performance of the DRT voigt circuit identification.
the inputs are the impedances and frequencies.
the function returns a comparison of the voigt-model and the original measurements.
"""
function DRT_circuit_identify_nyquist(impedances,frequencies)
    # step one: calculate the DRT
        R_0 = minimum(real(impedances))
        relaxation_times, peak_amplitudes, taus_out, drt = compute_DRT(frequencies,impedances)
    # step two: determine the number of peaks
        n = length(relaxation_times)
    # step there: get R_N values corresponding to the peaks
        R_Ns = obtain_RNs(frequencies, impedances, relaxation_times)
        # step four: simulate the EIS using the obtained R_N values
        spectrum = joint_simulation(relaxation_times, R_Ns,frequencies) .+ R_0
        fig_ = nyquist(impedances)
        return nyquist!(spectrum)
end