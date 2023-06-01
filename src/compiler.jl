export decomp_4x4_diag, euler_angles_2x2, arbi_2x2unitary_gate, decomp_dxd_unitary

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
    return rotation_factors \ (imag(log.(diag(A))) .* 2.0)
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

function arbi_2x2unitary_gate(A::AbstractMatrix{T}, n::G, i::G, j::G) where {T,G}
    phase, theta, psi, delta = euler_angles_2x2(A)
    return put(n, (i, j) => diag_2x2_gate(Diagonal([exp(im * psi / 2.0), exp(-im * psi / 2.0)]))) *
           Rxy(n, i, j, theta, -π / 2.0) *
           put(3, (1, 2) => diag_2x2_gate(Diagonal([exp(im * delta / 2), exp(-im * delta / 2)])))
end

# use givens rotation to compile a unitary to series of 2-level unitary gates
function decomp_dxd_unitary(A::AbstractMatrix{T}) where {T}
    two_level_unitaries = SparseMatrixCSC{ComplexF64,Int64}[]
    isunitary(A) || error("A must be unitary")
    # either follow  https://arxiv.org/pdf/1210.7366.pdf
    for i in 1:size(A, 1)
        for j in (i+1):size(A, 1)
            G, r = LinearAlgebra.givens(A, i, j, i)
            G_phase = sparse(Matrix{ComplexF64}(I, size(A)))
            if j != size(A, 1)
                r = 1.0
            end
            G_phase[i, i] = G.c / r
            G_phase[j, j] = G.c / r
            G_phase[i, j] = G.s / r
            G_phase[j, i] = -conj(G.s) / r
            push!(two_level_unitaries, G_phase)
            # display(G * Matrix(I,size(A)))
            # display(G_phase)
            # display(G * A)
            # display(G_phase * A)
            # display(isunitary(G * A))
            # display(isunitary(G_phase * A))
            A = G_phase * A
        end
    end
    # diagonal adjusting for phase
    rot_angle = angle(A[end, end]) / float(size(A, 1))
    for i in 1:(size(A, 1)-1)
        D_phase = sparse(Matrix{ComplexF64}(I, size(A)))
        D_phase[i, i] = exp(im * rot_angle)
        D_phase[end, end] = exp(-im * rot_angle)
        display(D_phase)
        push!(two_level_unitaries, D_phase)
        A = D_phase * A
    end
    # display(isunitary(A))
    # display(diag(A))
    # display(equiv_upto_phase(A, Matrix(I, size(A))))
    return two_level_unitaries
end
