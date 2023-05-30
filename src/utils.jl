export equiv_upto_phase

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
    return all(isapprox(A[i] / B[i], first_non_zero_ratio, atol=1e-6) for i in eachindex(A) if abs(A[i]) > 1e-8 && abs(B[i]) > 1e-8)
end
