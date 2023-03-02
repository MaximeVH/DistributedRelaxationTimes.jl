"""
    find_optimal_lambda(frequencies, measurements;
            method="re_im_cv",
            width_coeff=0.1,
            rbf_kernel=SqExponentialKernel())

Suggests values for the hyperparameter `λ`` using Saccoccio et al.'s discrepancy or Re-Im-crossvalidation methods.

The function's inputs and keyword arguments are similar to the those of the `compute_DRT` function,
with the exception of the `method` keyword argument, which allows users to choose between

- `re_im_cv` : Re-Im-crossvalidation
- `discrepancy`: minimisation of the discrepancy between the weights `θ` calculated
                    using the real and imaginary parts of the impedance spectrum.

"""
function find_optimal_lambda(frequencies, measurements;
                    method="re_im_cv",
                    width_coeff=0.1,
                    rbf_kernel=SqExponentialKernel()) #TODO: method based on found peaks
    ϵ = calculate_shape_factor(frequencies, width_coeff, rbf_kernel)
    Z_exp_real = real(measurements)
    Z_exp_imag = imag(measurements)

    #Select a range of λ values to evaluate.
    lambda_values = [10.0^i for i in  -6:1]

    #initialize Vector of hyperparameter objective values.
    hyperparam_objective_values = Vector{Float64}(undef, length(lambda_values))

    Z_drt_imag = construct_Z_imag(frequencies, ϵ, rbf_kernel)
    Z_drt_real = construct_Z_real(frequencies, ϵ, rbf_kernel)

    #initialize thetas.
    n = length(frequencies) + 1
    θ = fill(0.05, n)
    upper = 1e8ones(n) 
    lower = zeros(n)

    for i in eachindex(lambda_values)
        λ = lambda_values[i]
        # Get the estimated θ vectors using the real and imaginary parts of the impedance.
        results_imag = optimize(x -> objective(Z_drt_imag, -Z_exp_imag, x, λ), lower, upper, θ,
                                Fminbox(BFGS()),autodiff=:forward);
        results_real = optimize(x -> objective(Z_drt_real, Z_exp_real, x, λ), lower, upper, θ,
                                Fminbox(BFGS()),autodiff=:forward);
        # Select the relevant parts.
        θ_imag = results_imag.minimizer
        θ_real = results_real.minimizer

        #reconstructed impedances where real part is calculated with thetas from imaginary part, and vice versa.
        re_im_cv_Re = norm(Z_drt_real*θ_imag - Z_exp_real)^2
        re_im_cv_Im = norm(Z_drt_imag*θ_real + Z_exp_imag)^2

        #The differences between θ_imag and θ_real, as well as their reconstructed Z's should be minimized.
        if method == "re_im_cv"
            hyperparam_objective_values[i] = re_im_cv_Re + re_im_cv_Im
        elseif method == "discrepancy"
            hyperparam_objective_values[i] = norm(θ_imag - θ_real)^2
        end
    end

    min_idx = findmin(hyperparam_objective_values)[2]

    return lambda_values[min_idx]
end