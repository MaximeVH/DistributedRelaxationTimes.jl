"""
integration_Z_imag(ϵ,f_r,f_c, rbf_kernel)

Performs the numerical integration required to calculate the values of the reconstructed imaginary impedance values in the DRT equation.
"""
function integration_Z_imag(ϵ,f_r,f_c, rbf_kernel) 
    tmp = 2*π*(f_r/f_c)
    integrand(x) = tmp/(1 ./exp(x)+(tmp^2)*exp(x))*rbf_kernel(x*ϵ,0)
    out_integral = quadgk(x -> integrand(x), -Inf, Inf, rtol=1e-9)
    return out_integral[1]
end 

"""
integration_Z_real(ϵ,f_r,f_c, rbf_kernel)

Performs the numerical integration required to calculate the values of the reconstructed real impedance values in the DRT equation.
"""
function integration_Z_real(ϵ,f_r,f_c, rbf_kernel) 
    tmp = 2*π*(f_r/f_c)
    integrand(x) = (1 ./(1+(tmp^2)*exp(2*x)))*rbf_kernel(x*ϵ,0)
    out_integral = quadgk(x -> integrand(x), -Inf, Inf, rtol=1e-9)
    return out_integral[1]
end 

"""
construct_Z_imag(freq, ϵ, rbf_kernel)

Calculates the matrix to be multiplied with the weights Θ to obtain the imaginary part of the reconstructed impedance values in the DRT equation.
"""
function construct_Z_imag(freq, ϵ, rbf_kernel)
    R = zeros((1,length(freq)))
    C = zeros((length(freq),1))

    for i in eachindex(freq)
        freq_R =  freq[i]
        freq_C = freq[1]
        C[i,1] = integration_Z_imag(ϵ,freq_R,freq_C,rbf_kernel)
    end

    for j in eachindex(freq)
        freq_R =  freq[1]
        freq_C = freq[j]
        R[1,j] = integration_Z_imag(ϵ,freq_R,freq_C,rbf_kernel)
    end

    temp_Z_im = Toeplitz(vec(C),vec(R))
    out_Z_im = hcat(zeros((length(freq), 2)), temp_Z_im)

    return out_Z_im
end 

"""
construct_Z_real(freq, ϵ, rbf_kernel)

Calculates the matrix to be multiplied with the weights Θ to obtain the real part of the reconstructed impedance values in the DRT equation.
"""
function construct_Z_real(freq, ϵ, rbf_kernel)
    R = zeros((1,length(freq)))
    C = zeros((length(freq),1))
    
    for i in eachindex(freq)
        freq_R =  freq[i]
        freq_C = freq[1]
        C[i,1] = integration_Z_real(ϵ,freq_R,freq_C, rbf_kernel)
    end
        
    for j in eachindex(freq)
        freq_R =  freq[1]
        freq_C = freq[j]
        R[1,j] = integration_Z_real(ϵ,freq_R,freq_C, rbf_kernel)
    end
        
        temp_Z_re = Toeplitz(vec(C),vec(R))
    
        temp_Z_re_2 = hcat(ones(length(freq)), temp_Z_re)
        out_Z_re = hcat(zeros(length(freq)),temp_Z_re_2)
    return out_Z_re
 end 