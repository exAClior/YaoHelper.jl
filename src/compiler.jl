export decomp_subspace_diag, euler_angles_2x2, decomp_dxd_unitary

# need to change this to decompose for arbit weight sector of 2^n x 2^n diagonal
# implement decomposition for arbitrary diagonal 4x4 matrices
function decomp_subspace_diag(A::AbstractVector{T}, n::Int, charge::Int) where {T}
    length(A) == binomial(n, charge) || error("A must be of size $binomial(n,charge)")

    bases = subspace_indices(n, charge)
    sort!(bases)

    rotation_factors = zeros(length(bases), n + binomial(n, 2))

    for i = 1:n
        rotation_factors[:, i] = diag(mat(-kron(n, i => Z))[bases, bases])
    end
    idx = n + 1
    for (i, j) in Combinatorics.combinations(1:n, 2)
        rotation_factors[:, idx] = diag(mat(-kron(n, i => Z, j => Z))[bases, bases])
        idx += 1
    end

    # get diagonal part of A into a vector
    return rotation_factors \ (imag.(log.(A)) .* 2.0)
end


function euler_angles_2x2(A::AbstractMatrix{T}) where {T}
    size(A) != (2, 2) && error("A must be 2x2")

    # perform euler angle decomposition on a 2x2 matrix

    phase = det(A)^(-1.0 / 2.0)
    A = phase .* A
    theta = atan(abs(A[1, 2]), abs(A[1, 1])) * 2.0
    alpha = angle(A[1, 1]) * 2.0
    beta = angle(A[1, 2]) * 2.0

    psi = (alpha + beta) / 2.0
    delta = (alpha - beta) / 2.0

    return phase, theta, psi, delta
end


# use givens rotation to compile a unitary to series of 2-level unitary gates
function decomp_dxd_unitary(A::AbstractMatrix{T}) where {T}
    two_level_unitaries = SparseMatrixCSC{ComplexF64,Int64}[]
    subspaces = Vector{Int64}[]
    isunitary(A) || error("A must be unitary")
    # either follow  https://arxiv.org/pdf/1210.7366.pdf
    for i = 1:size(A, 1)
        for j = (i+1):size(A, 1)
            G, r = LinearAlgebra.givens(A, i, j, i)
            G_phase = sparse(Matrix{ComplexF64}(I, size(A)))
            if j != size(A, 1)
                r = 1.0
            end
            G_phase[i, i] = G.c / r
            G_phase[j, j] = conj(G.c) / r
            G_phase[i, j] = G.s / r
            G_phase[j, i] = -conj(G.s) / r
            pushfirst!(two_level_unitaries, G_phase)
            pushfirst!(subspaces, [i, j])
            A = G_phase * A
        end
    end
    # diagonal adjusting for phase
    rot_angle = angle(A[end, end]) / float(size(A, 1))
    for i = 1:(size(A, 1)-1)
        D_phase = sparse(Matrix{ComplexF64}(I, size(A)))
        D_phase[i, i] = exp(im * rot_angle)
        D_phase[end, end] = exp(-im * rot_angle)
        pushfirst!(two_level_unitaries, D_phase)
        pushfirst!(subspaces, [i, size(A, 1)])
        A = D_phase * A
    end
    return two_level_unitaries, subspaces
end
