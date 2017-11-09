using Microbiome
using Distances
using Base.Test

@testset "Abundances" begin
    abund = AbundanceTable(
        rand(100, 10), ["sample_$x" for x in 1:10],
        ["feature_$x" for x in 1:100])

    relab_fract = relativeabundance(abund)
    @test typeof(relab_fract) <: AbundanceTable

    relab_perc = relativeabundance(abund, kind=:percent)
    @test typeof(relab_perc) <: AbundanceTable

    filt = filterabund(relab_fract, 5)
    @test typeof(filt) <: AbundanceTable

    @test size(abund) == (100, 10)
    @test size(relab_fract) == (100, 10)
    @test size(relab_perc) == (100, 10)
    @test size(filt) == (6, 10)

    for i in 1:10
        @test sum(relab_fract[:, i]) ≈ 1
        @test sum(relab_perc[:, i]) ≈ 100
        @test sum(filt[:, i]) ≈ 1
    end

    @test filt.features[end] == "other"
end

@testset "Distances" begin
    abund = AbundanceTable(
        rand(100, 10), ["sample_$x" for x in 1:10],
        ["feature_$x" for x in 1:100])

    dm = getdm(abund, Jaccard()) # TODO: change to BrayCurtis after PR merges in Distances.jl
    p = pcoa(dm, correct_neg=true)
    rowdm = getrowdm(abund, Jaccard()) # TODO: change to BrayCurtis after PR merges in Distances.jl

    @test size(dm) == (10, 10)
    @test size(rowdm) == (100, 100)
    for i in 1:10; @test dm[i,i] == 0; end
    for i in 1:100; @test rowdm[i,i] == 0; end

    @test sum(p.variance_explained) ≈ 1
    for i in p.eigenvalues; @test i > 0; end
end
