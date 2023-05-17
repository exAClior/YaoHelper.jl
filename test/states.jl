

@testset "States" begin
    n = 2
    @test one_state(n).state ≈ [0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 1.0 + 0.0im]'
    @test plus_state(n).state ≈ [0.5 0.5 0.5 0.5]'
    @test minus_state(n).state ≈ [0.5 -0.5 -0.5 0.5]'

end
