"""
    calculate_shape_factor(frequencies,coefficient, rbf_kernel)

Calculates the shape factor of the RBF using the Full Width at Half Maximum (FWHM).

The inputs are the `frequencies`, the width coefficient hyperparameter `coefficient` and
the used RBF function `rbf_kernel`.
"""
function calculate_shape_factor(frequencies,coefficient, rbf_kernel)
    rbf_fwhm(x) = rbf_kernel(x, 0) - 0.5
    fwhm_coeff = 2find_zero(rbf_fwhm, 1);
    D_f = mean(diff(log.(1 ./ frequencies))) 
    ϵ  = fwhm_coeff * coefficient / D_f 
    return ϵ
end