export Rxx, Ryy, Rxy, cRz, diag_2x2_gate

function Rxx(n::Int, i::Int, j::Int, θ::T) where {T}
    return rot(kron(n, i => X, j => X), θ)
end

function Ryy(n::Int, i::Int, j::Int, θ::T) where {T}
    return rot(kron(n, i => Y, j => Y), θ)
end

function Rxy(n::Int, i::Int, j::Int, θ::T, β::D=0.0) where {T,D}
    # gave two types in case θ is an irrational value
    # because XX and YY commute, we can apply BCH formula
    return kron(n, i => Rz(-β)) * Rxx(n, i, j, θ / 2) * Ryy(n, i, j, θ / 2) * kron(n, i => Rz(β))
end

function cRz(n::Int, ctrl::Int, tgt::Int, θ::T) where {T}
    return control(n, ctrl, tgt => Rz(θ))
end

function diag_2x2_gate(A::AbstractMatrix{T}) where {T}
    !isdiag(A) && error("A must be diagonal")

    angles = decomp_4x4_diag(Diagonal([0 A[1, 1] A[2, 2] 0]))

    return chain(2, put(2, 1 => Rz(angles[1])), put(2, 2 => Rz(angles[2])), cRz(2, 1, 2, angles[3]))
end
