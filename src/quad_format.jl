function quad_format(A,b,M,λ)

    H = 2 .*(transpose(A)*A + λ .* M)
    H = (transpose(H)+H) ./2
    c = -2 .* transpose(b)*A
    
    return H,c
end 