using Microbiome
using Distances
using DataFrames
using Base.Test

@testset "Abundances" begin
    # Constructors
    M = rand(100, 10)
    df = hcat(DataFrame(x=collect(1:100)), DataFrame(M))

    a1 = AbundanceTable(M)
    a2 = AbundanceTable(df)

    abund = AbundanceTable(
        M, ["sample_$x" for x in 1:10],
        ["feature_$x" for x in 1:100])

    @test a1.table == a2.table == abund.table

    # Normalization functions
    relab_fract = relativeabundance(abund)
    @test typeof(relab_fract) <: AbundanceTable

    relab_perc = relativeabundance(abund, kind=:percent)
    @test typeof(relab_perc) <: AbundanceTable

    rnorm = rownormalize(abund)
    cnorm = colnormalize(abund)

    @test size(abund) == (100, 10)
    @test size(relab_fract) == (100, 10)
    @test size(relab_perc) == (100, 10)
    @test size(rnorm) == (100, 10)
    @test size(cnorm) == (100, 10)

    for j in 1:10
        @test sum(relab_fract[:, j]) ≈ 1
        @test sum(relab_perc[:, j]) ≈ 100
    end

    for i in 1:size(rnorm,1)
        @test maximum(rnorm[i,:]) ≈ 1
    end

    for j in 1:size(cnorm,2)
        @test maximum(cnorm[:,j]) ≈ 1
    end

    # Filtering Functions
    filt = filterabund(relab_fract, 5)
    @test typeof(filt) <: AbundanceTable
    @test typeof(filterabund(df, 5)) <: AbundanceTable

    @test size(filt) == (6, 10)
    for i in 1:10
        @test sum(filt[:, i]) ≈ 1
    end

    @test filt.features[end] == "other"
end

@testset "Distances" begin
    # Constructors
    M = rand(100, 10)
    df = hcat(DataFrame(x=collect(1:100)), DataFrame(M))
    abund = AbundanceTable(
        M, ["sample_$x" for x in 1:10],
        ["feature_$x" for x in 1:100])

    dm1 = getdm(M, Jaccard()) # TODO: switch to BrayCurtis when Distances.jl adds release
    dm2 = getdm(df, Jaccard()) # TODO: switch to BrayCurtis when Distances.jl adds release
    dm = getdm(abund, Jaccard()) # TODO: switch to BrayCurtis when Distances.jl adds release

    @test dm.dm == dm1.dm == dm2.dm

    rowdm1 = getrowdm(M, Jaccard()) # TODO: switch to BrayCurtis when Distances.jl adds release
    rowdm2 = getrowdm(df, Jaccard()) # TODO: switch to BrayCurtis when Distances.jl adds release
    rowdm = getrowdm(abund, Jaccard()) # TODO: switch to BrayCurtis when Distances.jl adds release

    @test rowdm.dm == rowdm1.dm == rowdm2.dm

    @test size(dm) == (10, 10)
    @test size(rowdm) == (100, 100)
    for i in 1:10; @test dm[i,i] == 0; end
    for i in 1:100; @test rowdm[i,i] == 0; end

    # PCoA
    p = pcoa(dm, correct_neg=true)
    @test sum(p.variance_explained) ≈ 1
    for i in p.eigenvalues; @test i > 0; end

    @test length(principalcoord(p, 1)) == size(dm, 1)
    for i in 1:9; @test typeof(eigenvalue(p, i)) <: Real; end
end
