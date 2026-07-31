using ClusterAlgebras
using AbstractAlgebra
using Test
using Aqua

@testset "ClusterAlgebras.jl" begin
    @testset "Aqua" begin
        Aqua.test_all(ClusterAlgebras)
    end

    include("test_quiver.jl")
    include("test_canonical_form.jl")
    include("test_random_quiver.jl")
    include("test_seed.jl")
    include("test_mutation.jl")
    include("test_named_quivers.jl")
    include("test_root_system.jl")
    include("test_folding.jl")
    include("test_greedy.jl")
    include("test_bounds.jl")
    include("test_denominator_vector.jl")
    include("test_frieze.jl")
    include("test_mutation_class.jl")
    include("test_vectors.jl")
    include("test_fpolynomials.jl")
    include("test_enumerative.jl")
    include("test_coefficients.jl")
    include("test_green_sequences.jl")
    include("test_periodicity.jl")
    include("test_dt_transformation.jl")
    include("test_ks_dilog.jl")
    include("test_grassmannian.jl")
    include("test_symbol_alphabet.jl")
    include("test_mutation_type.jl")
    include("test_interop.jl")
    include("test_block_decomposition.jl")
    include("test_makie_ext.jl")
end
