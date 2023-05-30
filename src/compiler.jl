export decomp_4x4_diag

# implement decomposition for arbitrary diagonal 4x4 matrices
function decomp_4x4_diag(A::AbstractMatrix{T}) where {T}
    !isdiag(A) && error("A must be diagonal")

    # columns are
    # rotation angle factors for Rz on 1st qubit
    # rotation angle factors for Rz on 2nd qubit
    # rotation angle factors for cRz on 2nd qubit, control on 1st qubit
    # global phase
    rotation_factors = [-1 -1 0 1; 1 -1 -1 1; -1 1 0 1; 1 1 1 1]

    # get diagonal part of A into a vector
    return rotation_factors \ (imag(log.(diag(A))) .* 2)
end
