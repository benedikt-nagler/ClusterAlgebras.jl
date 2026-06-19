using ClusterAlgebras
using AbstractAlgebra
using Test
using Aqua

@testset "ClusterAlgebras.jl" begin
    @testset "Aqua" begin
        Aqua.test_all(ClusterAlgebras)
    end

    include("test_quiver.jl")
    include("test_seed.jl")
    include("test_mutation.jl")
    include("test_named_quivers.jl")
    include("test_root_system.jl")
    include("test_denominator_vector.jl")
    include("test_frieze.jl")
    include("test_mutation_class.jl")
    include("test_vectors.jl")
    include("test_fpolynomials.jl")
    include("test_enumerative.jl")
    include("test_coefficients.jl")
    include("test_green_sequences.jl")
end
