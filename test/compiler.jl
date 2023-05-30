@testset "Compiler" begin
    mtx = Diagonal(exp.(im .* rand(4) .* 2 * π))

    angles = decomp_4x4_diag(mtx)
    @test exp.(im * angles[4] / 2 .* ones(4)) .* mat(chain(2, put(2, 1 => Rz(angles[1])), put(2, 2 => Rz(angles[2])), cRz(2, 1, 2, angles[3]))) ≈ mtx
end
