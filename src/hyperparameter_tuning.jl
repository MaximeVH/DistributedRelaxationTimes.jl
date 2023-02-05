function find_optimal_lambda(frequencies, measurements, width_coeff = 0.1) #Alternative method: choose the lambda value for which the peaks of real and imag methods are as similar as possible.
    ϵ = calculate_shape_factor(frequencies,width_coeff,rbf_kernel)
    Z_exp_real = real(measurements)
    Z_exp_imag = imag(measurements)
    lambda_values = [10.0^i for i in  range(-6, stop = 6, step=1)]
    CV_objectives = Array{Float64}(undef,length(lambda_values))
    discrepancy_values = Array{Float64}(undef,length(lambda_values))
    Z_drt_imag = construct_Z_imag(frequencies, ϵ)
    Z_drt_real = construct_Z_real(frequencies, ϵ)
    n = length(frequencies) + 2
    θ = fill(0.0, n)

    for i in eachindex(lambda_values)
        λ = lambda_values[i]
        results_imag = optimize(x -> objective_fun(Z_drt_imag, -Z_exp_imag, x, λ), lower, upper, θ,
        Fminbox(BFGS()),autodiff=:forward);
        results_real = optimize(x -> objective_fun(Z_drt_real, Z_exp_real, x, λ), lower, upper, θ,
        Fminbox(BFGS()),autodiff=:forward);
        x_primeprime = results_imag.minimizer
        x_prime = results_real.minimizer
        re_im_cv_prime = loss_fun(Z_drt_real, Z_exp_real, x_primeprime)
        re_im_cv_primeprime = loss_fun(Z_drt_imag, -Z_exp_imag, x_prime)
        CV_objectives[i] = re_im_cv_prime + re_im_cv_primeprime
        discrepancy_values[i] = norm(x_primeprime - x_prime)^2
    end

    min_idx_cv = findmin(CV_objectives)[2]
    min_idx_disc = findmin(discrepancy_values)[2]
    println("Minimum for CV: $(lambda_values[min_idx_cv])")
    println("Minimum for discrepancy: $(lambda_values[min_idx_disc])")
    return CV_objectives, discrepancy_values
end