"""
regulariser(λ::Float64,Θ::Vector{Float64})

Calculates the regularisation term using the L2-norm of the weight vector `Θ` and the tuning hyperparameter `λ`.

"""
function regulariser(λ,θ)
    return λ *norm(θ)^2
end 

"""
objective(X::Matrix{Float64}, Y::Vector{Float64}, θ::Vector{Float64}, λ::Float64)

Objective function for the Tikhonov regularisation, where `X` is the matrix for the reconstruction of the real or imaginary impedance values,
`Y` is the real or imaginary part of the impedance measurements, `θ` is a vector of weights to be optimised, and `λ` is the regularisation hyperparameter.
"""
function objective(X, Y, θ, λ)
    return norm(X*θ - Y)^2 + regulariser(λ,θ)
end 

"""
joint_objective(X1::Matrix{Float64}, Y1::Vector{Float64}, X2::Matrix{Float64}, Y2::Vector{Float64}, θ::Vector{Float64}, λ, weights::Vector{Int64} = [1,1])

Objective function for the DRT calculation using both the real (`Y1`) and imaginary (`Y2`) parts of the impedance measurements.
- `X1` and `X2` are the matrices for the reconstruction of the real and imaginary impedance values, respectively.
- `θ` is a vector of weights to be optimised.
- `λ` is the regularisation hyperparameter.
- `weights` provides the option to provide more weight to the real or imaginary parts during the optimisation.
"""
function joint_objective(X1, Y1, X2, Y2, θ, λ, weights = [1,1])
    return weights[1]*objective(X1, Y1, θ, λ) + weights[2]*objective(X2, Y2, θ, λ)
end