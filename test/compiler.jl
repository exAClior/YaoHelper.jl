@testset "Compiler" begin
    # mtx = Diagonal(exp.(im .* rand(4) .* 2 * π))

    # angles = decomp_4x4_diag(mtx)
    # @test exp.(im * angles[4] / 2 .* ones(4)) .* mat(chain(2, put(2, 1 => Rz(angles[1])), put(2, 2 => Rz(angles[2])), cRz(2, 1, 2, angles[3]))) ≈ mtx

    rdu = rand_unitary(2)
    phase, theta, psi, delta = euler_angles_2x2(rdu)
    @test phase .* rdu ≈
          [exp(im * psi / 2.0) 0.0; 0.0 exp(-im * psi / 2.0)] *
          [cos(theta / 2.0) sin(theta / 2.0); -sin(theta / 2.0) cos(theta / 2.0)] *
          [exp(im * delta / 2.0) 0.0; 0.0 exp(-im * delta / 2.0)]


    qubits = 3
    charge = 2
    bases = subspace_indices(qubits, charge)
    rdu = Matrix{ComplexF64}(I, 2^qubits, 2^qubits)
    rdu[bases, bases] = rand_unitary(length(bases))
    two_level_unitaries = decomp_dxd_unitary(rdu[bases, bases])
    reconstruct = Matrix{ComplexF64}(I, length(bases), length(bases))
    for u4 in two_level_unitaries
        reconstruct = u4' * reconstruct
    end
    @test equiv_upto_phase(reconstruct, rdu[bases, bases])
end
