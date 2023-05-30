@testset "Compiler" begin
    mtx = Diagonal(exp.(im .* rand(4) .* 2 * π))

    angles = decomp_4x4_diag(mtx)
    @test exp.(im * angles[4] / 2 .* ones(4)) .* mat(chain(2, put(2, 1 => Rz(angles[1])), put(2, 2 => Rz(angles[2])), cRz(2, 1, 2, angles[3]))) ≈ mtx

    rdu = rand_unitary(2)
    phase, theta, psi, delta = euler_angles_2x2(rdu)
    @test phase .* rdu ≈ [exp(im * psi / 2.0) 0.0; 0.0 exp(-im * psi / 2.0)] * [cos(theta / 2.0) sin(theta / 2.0); -sin(theta / 2.0) cos(theta / 2.0)] * [exp(im * delta / 2.0) 0.0; 0.0 exp(-im * delta / 2.0)]
end
