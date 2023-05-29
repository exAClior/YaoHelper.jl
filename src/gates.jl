export Rxx, Ryy, Rxy

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
