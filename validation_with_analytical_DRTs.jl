# In this script we validate the DRT deconvolution method through comparison with two circuits for which the
# expression for the analytical DRT is available.
using DRT, Plots

# I) Zarc model

# function to calculate the impedance of the Zarc model.
function zarc_impedance(R_inf,R_ct,τ_0,ϕ,f)
    Z = R_inf + (R_ct/(1 + (2*π*im*f*τ_0)^ϕ))
    return Z 
end

# function for the analytical DRT calculation for the Zarc model.
function zarc_analytical_DRT(R_ct,τ_0,ϕ,f) #R_inf unused
    τ = 1/f
    gamma = (R_ct/(2*π))*(sin((1-ϕ)*π)/(cosh(ϕ*log(τ/τ_0)) - cos((1-ϕ)*π)))
    return gamma
end

# define frequencies and parameters of the Zarc model.
frequencies_test = [10.0^i for i in LinRange(-2, 6, 100)]
R_inf, R_ct, τ_0, ϕ = [20,100,0.015,0.65]

# simulate and visualize the EIS measurements of the Zarc model.
measurements_test = [zarc_impedance(R_inf,R_ct,τ_0,ϕ,f_) for f_ in frequencies_test]
plot(real(measurements_test),-imag(measurements_test),label= "Zarc Nyquist")

# simulate and visualize the analytical DRT of the Zarc model.
analyrical_DRT_test = [zarc_analytical_DRT(R_ct,τ_0,ϕ,f_) for f_ in frequencies_test]
taus = 1 ./ frequencies_test
scatter(taus,analyrical_DRT_test,xaxis=:log, label = "Zarc analytical DRT")

# Calculate the DRT from the impedance measurements and compare with the analytical solution.
tau_relax_test  , amps_test  , taus_test  , drt_test  = compute_DRT(frequencies_test,measurements_test,λ =  10^-6, width_coeff = 0.1 ,method= "re_im")
plot!(taus_test  , drt_test,xaxis = :log, label = "Zarc numerical DRT")
xlabel!("τ (s)")
ylabel!("γ (Ω)")

# II) Double Zarc model

# function to calculate the impedance of the double Zarc model.
function double_zarc_imepdance(R_inf,R_ct,τ_0,ϕ,τ_02,f)
    Z = R_inf + (R_ct/(1 + (2*π*im*f*τ_0)^ϕ)) + (R_ct/(1 + (2*π*im*f*τ_02)^ϕ))
    return Z 
end

# function for the analytical DRT calculation for the double Zarc model.
function double_zarc_analytical_DRT(R_ct,τ_0,τ_02,ϕ,f) #R_inf unused
    τ = 1/f
    gamma = (R_ct/(2*π))*sin((1-ϕ)*π) *(1/(cosh(ϕ*log(τ/τ_0)) - cos((1-ϕ)*π)) + 1/(cosh(ϕ*log(τ/τ_02)) - cos((1-ϕ)*π)))
    return gamma
end

# Define additional parameters for double Zarc model.
τ_0,τ_02 = 0.02, 0.0008

# simulate and visualize the EIS measurements of the double Zarc model.
measurements_test2 = [double_zarc_imepdance(R_inf,R_ct,τ_0,ϕ,τ_02,f_) for f_ in frequencies_test]
plot(real(measurements_test2),-imag(measurements_test2),label= "Zarc Nyquist") #Note that the relaxation processes are badly resolved in the Nyquist plot.

# simulate and visualize the analytical DRT of the double Zarc model.
analyrical_DRT_test = [double_zarc_analytical_DRT(R_ct,τ_0,τ_02,ϕ,f_) for f_ in frequencies_test]
taus = 1 ./ frequencies_test
scatter(taus,analyrical_DRT_test,xaxis=:log, label="analytical DRT") # With two clear peaks, the relaxation processes are well resolved here.

# Calculate the DRT from the impedance measurements and compare with the analytical solution.
tau_relax_test_double  , amps_test_double  , taus_test_double  , drt_test_double  = compute_DRT(frequencies_test,measurements_test2,λ =  10^-6, width_coeff = 0.1 ,method= "re_im")
plot!(taus_test_double , drt_test_double, xaxis = :log, label="numerical DRT")
xlabel!("τ (s)")
ylabel!("γ (Ω)")
