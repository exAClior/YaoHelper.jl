export Rxx, Ryy, Rxy, Rzz, cRz, diag_subspace_gate, arbi_2x2unitary_gate

function Rxx(n::Int, i::Int, j::Int, θ::T) where {T}
    return rot(kron(n, i => X, j => X), θ)
end

function Ryy(n::Int, i::Int, j::Int, θ::T) where {T}
    return rot(kron(n, i => Y, j => Y), θ)
end

function Rxy(n::Int, i::Int, j::Int, θ::T, β::D = 0.0) where {T,D}
    # gave two types in case θ is an irrational value
    # because XX and YY commute, we can apply BCH formula
    return kron(n, i => Rz(-β)) *
           Rxx(n, i, j, θ / 2) *
           Ryy(n, i, j, θ / 2) *
           kron(n, i => Rz(β))
end

function cRz(n::Int, ctrl::Int, tgt::Int, θ::T) where {T}
    return control(n, ctrl, tgt => Rz(θ))
end

function Rzz(n::Int, i::Int, j::Int, θ::T) where {T}
    return rot(kron(n, i => Z, j => Z), θ)
end

function diag_subspace_gate(A::AbstractVector{T}, n::Int, charge::Int) where {T}

    # calculate the diagonal angles here in the subspace, reduce to subspace
    # need to complexify
    angles = decomp_subspace_diag(A, n, charge)

    # need to modify here
    gate = kron(n, [i => Rz(angles[i]) for i = 1:n]...)
    idx = n + 1
    for (i, j) in Combinatorics.combinations(1:n, 2)
        gate = chain(n, gate, Rzz(n, i, j, angles[idx]))
        idx += 1
    end

    return Yao.Optimise.simplify(gate)
end

function arbi_2x2unitary_gate(
    A::AbstractMatrix{T},
    n::G,
    charge::G,
    subspace_idx::Vector{G},
) where {T,G}
    bases = subspace_indices(n, charge)
    i, j = different_bits(n, (bases[subspace_idx] .- 1)...)
    phase, theta, psi, delta = euler_angles_2x2(A)
    # just chain the diagonal gates + Rxy + diagonal
    psi_diag = ones(ComplexF64, binomial(n, charge))
    delta_diag = ones(ComplexF64, binomial(n, charge))
    psi_diag[subspace_idx[1]] = exp(im * (psi / 2.0)) / phase
    psi_diag[subspace_idx[2]] = exp(-im * (psi / 2.0)) / phase
    delta_diag[subspace_idx[1]] = exp(im * delta / 2.0)
    delta_diag[subspace_idx[2]] = exp(-im * delta / 2.0)
    return diag_subspace_gate(psi_diag, n, charge) *
           Rxy(n, i, j, theta, -π / 2.0) *
           diag_subspace_gate(delta_diag, n, charge)
end
