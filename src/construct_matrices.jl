"""
integration_Z_imag(ϵ,f_r,f_c, rbf_kernel)

Performs the numerical integration required to calculate the values of the reconstructed imaginary impedance values in the DRT equation.
"""
function integration_Z_imag(ϵ,f_r,f_c,rbf_kernel) 
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
As this matrix can be Toeplitz factorized, it can be efficiently constructed using the first column and first rows.

"""
function construct_Z_imag(freq, ϵ, rbf_kernel)
     #Initialize first row and column of out_Z_im
    R = zeros((1,length(freq)))
    C = zeros((length(freq),1))

     # Calculate the values of the first column.
    for i in eachindex(freq)
        freq_R =  freq[i]
        freq_C = freq[1]
        C[i,1] = integration_Z_imag(ϵ,freq_R,freq_C,rbf_kernel)
    end

    # Calculate the values of the first row.
    for j in eachindex(freq)
        freq_R =  freq[1]
        freq_C = freq[j]
        R[1,j] = integration_Z_imag(ϵ,freq_R,freq_C,rbf_kernel)
    end

    temp_Z_im = Toeplitz(vec(C),vec(R))
    out_Z_im = hcat(zeros((length(freq), 1)), temp_Z_im)

    return out_Z_im
end 

"""
construct_Z_real(freq, ϵ, rbf_kernel)

Calculates the matrix to be multiplied with the weights Θ to obtain the real part of the reconstructed impedance values in the DRT equation.
As this matrix can be Toeplitz factorized, it can be efficiently constructed using the first column and first rows.
"""
function construct_Z_real(freq, ϵ, rbf_kernel)
    #Initialize first row and column of out_Z_re
    R = zeros((1,length(freq)))
    C = zeros((length(freq),1))
    
    # Calculate the values of the first column.
    for i in eachindex(freq)
        freq_R =  freq[i]
        freq_C = freq[1]
        C[i,1] = integration_Z_real(ϵ,freq_R,freq_C, rbf_kernel)
    end
        
    # Calculate the values for the first row.
    for j in eachindex(freq)
        freq_R =  freq[1]
        freq_C = freq[j]
        R[1,j] = integration_Z_real(ϵ,freq_R,freq_C, rbf_kernel)
    end
        
    # Construct the Toeplitz matrix using the first column and row.
        temp_Z_re = Toeplitz(vec(C),vec(R))
        out_Z_re = hcat(ones(length(freq)),temp_Z_re)

    return out_Z_re
 end 

# Slow functions that don't make use of the Toeplitz factorization.
 function construct_Z_real_full(freq, ϵ, rbf_kernel)
    out_Z_re = Matrix{Float64}(undef,length(freq),length(freq))
    for i in eachindex(freq)
        for j in eachindex(freq)
            out_Z_re[i,j] = integration_Z_real(ϵ,freq[i],freq[j], rbf_kernel)
        end
    end
    return hcat(ones(length(freq)),out_Z_re)
end

function construct_Z_im_full(freq, ϵ, rbf_kernel)
    out_Z_re = Matrix{Float64}(undef,length(freq),length(freq))
    for i in eachindex(freq)
        for j in eachindex(freq)
            out_Z_re[i,j] = integration_Z_imag(ϵ,freq[i],freq[j], rbf_kernel)
        end
    end
    return hcat(zeros(length(freq)),out_Z_re)
end