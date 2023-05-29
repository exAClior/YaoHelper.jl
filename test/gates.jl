
@testset "Gates" begin
    θ = rand()
    β = rand()
    @test mat(Rxx(2, 1, 2, θ)) ≈ [cos(θ / 2) 0 0 -im*sin(θ / 2); 0 cos(θ / 2) -im*sin(θ / 2) 0; 0 -im*sin(θ / 2) cos(θ / 2) 0; -im*sin(θ / 2) 0 0 cos(θ / 2)]
    @test mat(Ryy(2, 1, 2, θ)) ≈ [cos(θ / 2) 0 0 im*sin(θ / 2); 0 cos(θ / 2) -im*sin(θ / 2) 0; 0 -im*sin(θ / 2) cos(θ / 2) 0; im*sin(θ / 2) 0 0 cos(θ / 2)]
    @test mat(Rxy(2, 1, 2, θ, β)) ≈ [1 0 0 0; 0 cos(θ / 2) -im*sin(θ / 2)*exp(-im * β) 0; 0 -im*sin(θ / 2)*exp(im * β) cos(θ / 2) 0; 0 0 0 1]
    @test mat(Rxy(2, 1, 2, π)) ≈ [1 0 0 0; 0 cos(π / 2) -im*sin(π / 2) 0; 0 -im*sin(π / 2) cos(π / 2) 0; 0 0 0 1]

end
