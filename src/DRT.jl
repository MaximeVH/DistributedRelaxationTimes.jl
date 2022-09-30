module DRT

    export compute_DRT
    
    using Optim, QuadGK, Roots, Statistics, ToeplitzMatrices, LinearAlgebra, Plots

    include("construct_matrices.jl")
    include("quad_format.jl")
    include("calculate_gamma.jl")
    include("calculate_shape_factor.jl")
    include("compute_DRT.jl") 

end