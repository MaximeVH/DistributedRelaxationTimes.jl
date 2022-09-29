function calculate_shape_factor(frequencies)
    rbf_inverse_quadric_4_FWHM(x) =  1 ./sqrt(1+(x).^2)-1/2;
    FWHM_coeff = 2*find_zero(rbf_inverse_quadric_4_FWHM, 1);
    delta = mean(diff(log.(1 ./frequencies)));
    coeff = 0.5
    out_val  = coeff*FWHM_coeff/delta;
    return out_val
end