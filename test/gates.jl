
@testset "Gates" begin
    θ = rand()
    β = rand()
    @test mat(Rxx(2, 1, 2, θ)) ≈ [cos(θ / 2) 0 0 -im*sin(θ / 2); 0 cos(θ / 2) -im*sin(θ / 2) 0; 0 -im*sin(θ / 2) cos(θ / 2) 0; -im*sin(θ / 2) 0 0 cos(θ / 2)]
    @test mat(Ryy(2, 1, 2, θ)) ≈ [cos(θ / 2) 0 0 im*sin(θ / 2); 0 cos(θ / 2) -im*sin(θ / 2) 0; 0 -im*sin(θ / 2) cos(θ / 2) 0; im*sin(θ / 2) 0 0 cos(θ / 2)]
    @test mat(Rxy(2, 1, 2, θ, β)) ≈ [1 0 0 0; 0 cos(θ / 2) -im*sin(θ / 2)*exp(-im * β) 0; 0 -im*sin(θ / 2)*exp(im * β) cos(θ / 2) 0; 0 0 0 1]
    @test mat(Rxy(2, 1, 2, π)) ≈ [1 0 0 0; 0 cos(π / 2) -im*sin(π / 2) 0; 0 -im*sin(π / 2) cos(π / 2) 0; 0 0 0 1]
    @test equiv_upto_phase(mat(diag_2x2_gate(Diagonal([1.0, 1.0]))), Matrix{ComplexF64}(I, 4, 4))
end
