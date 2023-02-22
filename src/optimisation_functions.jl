"""
    regulariser(λ ,Θ)

Calculates the regularisation term using the L2-norm of
the weight vector `Θ` and the tuning hyperparameter `λ`.

"""
function regulariser(λ, θ)
    return λ * norm(θ)^2
end 

"""
    objective(X, Y, θ, λ)

Objective function for the Tikhonov regularisation, where `X` is the matrix
for the reconstruction of the real or imaginary impedance values,
`Y` is the real or imaginary part of the impedance measurements, `θ` is
a vector of weights to be optimised, and `λ` is the regularisation hyperparameter.
"""
function objective(X, Y, θ, λ)
    return norm(X * θ - Y)^2 + regulariser(λ, θ)
end 

"""
    oint_objective(X₁, Y₁, X₂, Y₂, θ, λ, weights)

Objective function for the DRT calculation using both the real (`Y₁`) and imaginary (`Y₂`)
        parts of the impedance measurements.
    - `X₁` and `X₂` are the matrices for the reconstruction of the real and imaginary impedance values, respectively.
    - `θ` is a vector of weights to be optimised.
    - `λ` is the regularisation hyperparameter.
    - `weights` provides the option to provide more weight to the real or imaginary parts during the optimisation.
"""
function joint_objective(X₁, Y₁, X₂, Y₂, θ, λ, weights = [1,1])
    return weights[1] * objective(X₁, Y₁, θ, λ) + weights[2] * objective(X₂, Y₂, θ, λ)
end