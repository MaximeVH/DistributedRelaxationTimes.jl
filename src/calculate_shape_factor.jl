function calculate_shape_factor(frequencies,coefficient,rbf_kernel)
    rbf_fwhm(x) = (rbf_kernel(x,0))-0.5
    FWHM_coeff = 2*find_zero(rbf_fwhm, 1);
    D_f = mean(diff(log.(1 ./ frequencies))) 
    epsilon  = FWHM_coeff*coefficient/D_f 
    return epsilon
end