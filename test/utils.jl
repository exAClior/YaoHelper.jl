
@testset "Utils" begin

    A = rand(4,4)
    phase = exp(im * 2 * π * rand())
    B = phase .* A

    C = deepcopy(B)
    C[rand([1 2 3 4]),rand([1 2 3 4])] *= exp(im * 2 * π * rand())

    @test equiv_upto_phase(A, B)
    @test !equiv_upto_phase(A, C)
end
