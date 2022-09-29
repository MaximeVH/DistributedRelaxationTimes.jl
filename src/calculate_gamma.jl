function calculate_gamma(x, τ, ϵ) 
    τ_map = exp10.(range(log10(minimum(τ))-0.5, log10(maximum(τ))+0.5, length=10*length(τ)))   
    rbf(t) =   1 / sqrt(1 + (ϵ*t)^2);
        
    N_taus = length(τ)
    N_tau_map = length(τ_map)
    γ = zeros(N_tau_map, 1)
    B = zeros(N_tau_map, N_taus)
    delta_log_tau = 0
        for p in 1:N_tau_map
            for q in 1:N_taus
                delta_log_tau = log(τ_map[p])-log(τ[q])
                B[p,q] = rbf(delta_log_tau)       
            end
        end       
    
        γ  = B * x
        out_τ = τ_map 
        
    return out_τ,γ 
end
