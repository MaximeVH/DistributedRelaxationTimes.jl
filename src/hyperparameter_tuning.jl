"""
find_optimal_lambda(frequencies::Array{Float64,1}, measurements::Array{Complex{Float64},1}; method="re_im_cv", width_coeff = 0.1,rbf_kernel = SqExponentialKernel())

Suggests values for the hyperparameter λ using Saccoccio et al.'s discrepancy or Re-Im-crossvalidation methods.

The function's inputs and keyword arguments are similar to the those of the `compute_DRT` function,
with the exception of the `method` keyword argument, which allows users to choose between

- `"re_im_cv"`: Re-Im-crossvalidation
- `"discrepancy"`: minimisation of the discrepancy between the weights θ calculated using the real and imaginary parts of the impedance spectrum.

"""
function find_optimal_lambda(frequencies, measurements; method="re_im_cv", width_coeff = 0.1,rbf_kernel = SqExponentialKernel()) #Do to: method based on found peaks
    ϵ = calculate_shape_factor(frequencies,width_coeff,rbf_kernel)
    Z_exp_real = real(measurements)
    Z_exp_imag = imag(measurements)
    lambda_values = [10.0^i for i in  range(-6, stop = 1, step=1)]
    hyperparam_objective_values = Array{Float64}(undef,length(lambda_values))
    Z_drt_imag = construct_Z_imag(frequencies, ϵ, rbf_kernel)
    Z_drt_real = construct_Z_real(frequencies, ϵ, rbf_kernel)
    n = length(frequencies) + 2
    θ = fill(0.05, n)
    upper = 1e8 .*ones(n) 
    lower = zeros(n)

    for i in eachindex(lambda_values)
        λ = lambda_values[i]
        results_imag = optimize(x -> objective(Z_drt_imag, -Z_exp_imag, x, λ), lower, upper, θ,
        Fminbox(BFGS()),autodiff=:forward);
        results_real = optimize(x -> objective(Z_drt_real, Z_exp_real, x, λ), lower, upper, θ,
        Fminbox(BFGS()),autodiff=:forward);
        x_primeprime = results_imag.minimizer
        x_prime = results_real.minimizer
        re_im_cv_prime = loss(Z_drt_real, Z_exp_real, x_primeprime)
        re_im_cv_primeprime = loss(Z_drt_imag, -Z_exp_imag, x_prime)
        if method == "re_im_cv"
            hyperparam_objective_values[i] = re_im_cv_prime + re_im_cv_primeprime
        elseif method == "discrepancy"
            hyperparam_objective_values[i] = norm(x_primeprime - x_prime)^2
        end
    end

    min_idx = findmin(hyperparam_objective_values)[2]

    return lambda_values[min_idx]
end