
@testset "Gates" begin
    θ = rand()
    β = rand()
    @test mat(Rxx(2, 1, 2, θ)) ≈ [
        cos(θ / 2) 0 0 -im*sin(θ / 2)
        0 cos(θ / 2) -im*sin(θ / 2) 0
        0 -im*sin(θ / 2) cos(θ / 2) 0
        -im*sin(θ / 2) 0 0 cos(θ / 2)
    ]
    @test mat(Ryy(2, 1, 2, θ)) ≈ [
        cos(θ / 2) 0 0 im*sin(θ / 2)
        0 cos(θ / 2) -im*sin(θ / 2) 0
        0 -im*sin(θ / 2) cos(θ / 2) 0
        im*sin(θ / 2) 0 0 cos(θ / 2)
    ]
    @test mat(Rxy(2, 1, 2, θ, β)) ≈ [
        1 0 0 0
        0 cos(θ / 2) -im*sin(θ / 2)*exp(-im * β) 0
        0 -im*sin(θ / 2)*exp(im * β) cos(θ / 2) 0
        0 0 0 1
    ]
    @test mat(Rxy(2, 1, 2, π)) ≈
          [1 0 0 0; 0 cos(π / 2) -im*sin(π / 2) 0; 0 -im*sin(π / 2) cos(π / 2) 0; 0 0 0 1]

    @test mat(Rzz(2, 1, 2, θ)) ≈ [
        exp(-im * θ / 2) 0 0 0
        0 exp(im * θ / 2) 0 0
        0 0 exp(im * θ / 2) 0
        0 0 0 exp(-im * θ / 2)
    ]

    bases = subspace_indices(3, 2)
    diag_vec = [exp(im * rand()), exp(im * rand()), exp(im * rand())]
    @test mat(diag_subspace_gate(diag_vec, 3, 2))[bases, bases] = Diagonal(diag_vec)

    rdu1 = rand_unitary(2)
    @test equiv_upto_phase(
        mat(arbi_2x2unitary_gate(rdu1, 3, 2, [1, 2]))[bases[1:2], bases[1:2]],
        rdu1,
    )

end
