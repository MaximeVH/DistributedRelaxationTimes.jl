function compute_DRT(frequencies,impedances,λ)

    N_freqs = length(frequencies)
    τ = 1 ./ frequencies
    
    #Real part of Z is used for the DRT calculation.
    b_re = real(impedances)
    ϵ = calculate_shape_factor(frequencies) 

    #Assemble matrices.
    A_real_temp = construct_A_real(frequencies,ϵ)
    M_temp = construct_M_1(frequencies,ϵ)

    #Add column for resistance
    A_real = zeros(N_freqs,N_freqs+1)
    A_real[:,2:end] = A_real_temp
    A_real[:,1] .= 1
    M = zeros(N_freqs+1,N_freqs+1)
    M[2:end,2:end] = M_temp

    #Put the expression in quadratic form.
    H_re, c_re = quad_format(A_real, b_re, M, λ)

    #Define objective function to be minimized.
    objective_fnct(x) = 0.5*(transpose(x)*H_re*x) + (c_re*x)[1]
    
    # Optimisation, also try other optimisation approaches.
    results = optimize(objective_fnct, zeros(length(c_re)), ones(length(c_re))*1000, ones(length(c_re)))

    τ_out,γ_out = calculate_gamma(results.minimizer[2:end], τ, ϵ)
    return τ_out, γ_out
end