export equiv_upto_phase, subspace_indices, different_bits

function equiv_upto_phase(A::AbstractMatrix{T}, B::AbstractMatrix{G}) where {T,G}
    size(A) != size(B) && error("A and B must have the same size")

    # Find the first non-zero element ratio
    first_non_zero_ratio = nothing
    for i in eachindex(A)
        if abs(A[i]) > 1e-8 && abs(B[i]) > 1e-8
            first_non_zero_ratio = A[i] / B[i]
            break
        end
    end

    if first_non_zero_ratio === nothing
        throw(ArgumentError("At least one non-zero element should exist in both matrices."))
    end

    # Check if all other elements have the same ratio up to a tolerance
    return all(
        isapprox(A[i] / B[i], first_non_zero_ratio, atol = 1e-6) for
        i in eachindex(A) if abs(A[i]) > 1e-8 && abs(B[i]) > 1e-8
    )
end

function subspace_indices(n::Int, l::Int)
    bases = Int[]
    for loc in Combinatorics.combinations(1:n, l)
        b = 0
        for l in loc
            b += 1 << (l - 1)
        end
        push!(bases, b + 1)
    end
    return bases
end

function different_bits(n::Int, a::Int, b::Int)
    pos = []
    diff_bits = a ⊻ b
    for i = 1:n
        if diff_bits & (1 << (i - 1)) != 0
            push!(pos, i)
        end
    end
    return pos
end
