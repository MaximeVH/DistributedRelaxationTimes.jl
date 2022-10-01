function construct_A_real(frequencies,ϵ)
    
    τ = 1 ./ frequencies

    N_freqs = length(frequencies)
    out_A_re = zeros(N_freqs, N_freqs)

    #check for log-spacedness of frequencies.
    std_diff_freq = std(diff(log.(1 ./frequencies)))
    mean_diff_freq = mean(diff(log.(1 ./ frequencies)))
    toeplitz_trick = std_diff_freq/mean_diff_freq < 0.01
    R = zeros(N_freqs)
    C = zeros(N_freqs)
    if toeplitz_trick
        for p in 1:N_freqs
            C[p] = A_real_element(frequencies[p], τ[1], ϵ) #can later be expanded to include other RBFs.
        end
        for q in 1:N_freqs
            R[q]  = A_real_element(frequencies[1], τ[q], ϵ) 
        end
        out_A_re =  Toeplitz(C,R)
    else
        for p in 1:N_freqs
            for q in 1:N_freqs
                out_A_re[p, q]= A_real_element(frequencies[p], τ[q], ϵ)
            end
        end
    end
    return out_A_re
end

function A_real_element(freq_n,τ_m,ϵ)
    α = 2*pi*freq_n*τ_m  
    rbf_inverse_quadric_4_FWHM(x) =   1 / sqrt(1 + (ϵ*x)^2);
    integrand(x) = 1/(1+(α^2)*exp(2*x))*rbf_inverse_quadric_4_FWHM(x)
    return quadgk(x -> integrand(x), -50, 50, rtol=1e-9)[1]
end 

function A_imag_element(freq_n,τ_m,ϵ)
    α = 2*pi*freq_n*τ_m  
    rbf_inverse_quadric_4_FWHM(x) =   1 / sqrt(1 + (ϵ*x)^2);
    integrand(x) = α/(1 ./exp(x)+(α^2)*exp(x))*rbf_inverse_quadric_4_FWHM(x)
    return quadgk(x -> integrand(x), -50, 50, rtol=1e-9)[1]
end 


function construct_M_1(frequencies,ϵ)
    τ = 1 ./ frequencies
    N_freqs = length(frequencies)

    out_M = zeros(N_freqs, N_freqs)

    std_diff_freq = std(diff(log.(1 ./frequencies)))

    #check for log-spacedness of frequencies.
    mean_diff_freq = mean(diff(log.(1 ./ frequencies)))
    toeplitz_trick = std_diff_freq/mean_diff_freq < 0.01
    R = zeros(N_freqs)
    C = zeros(N_freqs)
    if toeplitz_trick
        for n in 1:N_freqs
            C[n] = inner_prod_rbf_1(frequencies[1], frequencies[n], ϵ)
        end
        for m in 1:N_freqs
            R[m] = inner_prod_rbf_1(frequencies[m], frequencies[1], ϵ)
        end
        return Toeplitz(C,R) 

    else
        for n in 0:N_freqs
            for m in 0:N_freqs            
                out_M[n,m] = inner_prod_rbf_1(frequencies[n], frequencies[m], ϵ)
            end
        end
end
end

function inner_prod_rbf_1(freq_n, freq_m, ϵ) 
    y_n = -log(freq_n)
    y_m = -log(freq_m)
    
    rbf_n(y) =  1/sqrt(1+(ϵ*(y-y_n))^2)
    rbf_m(y) =  1/sqrt(1+(ϵ*(y-y_m))^2)
    # compute derivative, may be improved with Zygote.
    delta = 1E-8
    sqr_drbf_dy(y) =  1/(2*delta)*(rbf_n(y+delta)-rbf_n(y-delta))*1/(2*delta)*(rbf_m(y+delta)-rbf_m(y-delta)) 
    out_val = quadgk(y -> sqr_drbf_dy(y), -50, 50, rtol=1e-9)[1]
    return out_val
end

function construct_A_imag(frequencies,ϵ)
    
    τ = 1 ./ frequencies

    N_freqs = length(frequencies)
    out_A_im = zeros(N_freqs, N_freqs)

    #check for log-spacedness of frequencies.
    std_diff_freq = std(diff(log.(1 ./frequencies)))
    mean_diff_freq = mean(diff(log.(1 ./ frequencies)))
    toeplitz_trick = std_diff_freq/mean_diff_freq < 0.01
    R = zeros(N_freqs)
    C = zeros(N_freqs)

    if toeplitz_trick
        for p in 1:N_freqs
            C[p] = A_imag_element(frequencies[p], τ[1], ϵ) #can later be expanded to include other RBFs.
        end
        for q in 1:N_freqs
            R[q]  = A_imag_element(frequencies[1], τ[q], ϵ) 
        end
        out_A_im =  Toeplitz(C,R)
    else
        for p in 1:N_freqs
            for q in 1:N_freqs
                out_A_re[p, q]= A_imag_element(frequencies[p], τ[q], ϵ)
            end
        end
    end
    return out_A_im
end