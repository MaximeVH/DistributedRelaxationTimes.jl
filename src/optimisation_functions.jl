function loss(X, Y, θ)
    l2_norm = norm(X*θ - Y)^2
    return l2_norm
end

function regulariser(θ)
    return norm(θ)^2
end 

function objective(X, Y, θ, λ)
    return loss(X, Y, θ) + λ * regulariser(θ)
end 

function joint_objective(X1, Y1, X2, Y2, θ, λ, weights = [1,1])
    return weights[1]*objective(X1, Y1, θ, λ) + weights[2]*objective(X2, Y2, θ, λ)
end