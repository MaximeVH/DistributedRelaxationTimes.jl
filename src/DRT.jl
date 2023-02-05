module DRT

    export compute_DRT, find_optimal_lambda, plot_DRT

    using Optim, QuadGK, Roots, FindPeaks1D
    using Statistics, ToeplitzMatrices, LinearAlgebra, KernelFunctions
    using Plots

    include("construct_matrices.jl")
    include("optimisation_functions.jl")
    include("hyperparameter_tuning.jl")
    include("calculate_shape_factor.jl")
    include("compute_DRT.jl") 

end