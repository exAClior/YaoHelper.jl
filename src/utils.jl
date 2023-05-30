export equiv_upto_phase

function equiv_upto_phase(A::AbstractMatrix{T}, B::AbstractMatrix{G}) where {T,G}
    size(A) != size(B) && error("A and B must have the same size")
    # test if each element of A is equal to the corresponding element of B upto a global phase
    target_phase = A[1,1] / B[1,1]
    return all(isapprox(target_phase,x,atol=1e-5) for x in (A ./ B))
end
