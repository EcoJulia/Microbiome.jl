using Microbiome
using Test
using SparseArrays
using Tables
using Dictionaries

@testset "Samples and Features" begin
    @testset "MicriobiomeSamples and metadata" begin
        ms = MicrobiomeSample("sample")
        @test name(ms) == "sample"
        @test isempty(metadata(ms))
        @test metadata(ms) isa Dictionary
        @test_throws Dictionaries.IndexError ms.thing = "metadata"
        @test_throws Dictionaries.IndexError ms[:thing] = "metadata"
        
        @test insert!(ms, :thing, "metadata") isa MicrobiomeSample
        @test ms.thing == ms[:thing] == "metadata"
        @test let
            ms.thing = "other metadata"
            @test ms.thing == ms[:thing] == "other metadata"
            ms[:thing] = "metadata"
            ms.thing == ms[:thing] == "metadata"
        end
        @test_throws MethodError ms["thing"] = "metadata"
        @test_throws Dictionaries.IndexError insert!(ms, :thing, "still other metadata")
        @test_throws Dictionaries.IndexError delete!(ms, :thing2)
        @test let
            delete!(ms, :thing)
            unset!(ms, :thing2)
            true
        end
        @test let
            set!(ms, :thing2, "metadata2")
            ms.thing2 == "metadata2"
        end
    end
    
    @testset "Taxa" begin
        _clades = (:domain, :kingdom, :phylum, :class, :order, :family, :genus, :species, :subspecies, :strain)
        txm = Taxon("taxon", missing)
        @test txm === Taxon("taxon")
        @test ismissing(clade(txm))
        @test !hasclade(txm)

        for (i, c) in enumerate(_clades)
            tx = Taxon("taxon", c)
            @test clade(tx) == c
            @test tx === Taxon("taxon", i-1)
            @test tx !== txm
            @test tx == txm
        end
        
        @test_throws ErrorException Taxon("taxon", :invalid)
        @test_throws ErrorException Taxon("taxon", 10)
        @test let tx = Taxon("taxon", :kingdom)
            hasclade(tx)
        end
    end

    @testset "Gene Functions" begin
        gfm = GeneFunction("gene", missing)
        @test name(gfm) == "gene"
        @test gfm === GeneFunction("gene")
        @test ismissing(taxon(gfm))
        @test !hastaxon(gfm)

        gf1 = GeneFunction("gene", Taxon("sp1", :species))
        gf2 = GeneFunction("gene", Taxon("sp1"))
        @test name(gf1) == "gene"
        
        @test gf1 == gf2 == gfm
        @test gf1 !== gf2
        @test gf1 !== gfm
        @test hastaxon(gf1)
        @test !ismissing(taxon(gf1))
        @test taxon(gf1) == taxon(gf2)
    end
end

@testset "Profiles" begin
    _clades = (:domain, :kingdom, :phylum, :class, :order, :family, :genus, :species, :subspecies, :strain)
    mss = [MicrobiomeSample("sample$i") for i in 1:5]
    txs = [Taxon("taxon$i", _clades[i]) for i in 1:9]
    push!(txs, Taxon("taxon10", missing))
    
    mat = spzeros(10,5)
    for i in 1:5; mat[i,i] = 1.; end

    comm = CommunityProfile(mat, txs, mss)
    @test nsamples(comm) == 5
    @test nfeatures(comm) == 10
    for (i, col) in enumerate(Tables.columns(comm))
        if i == 1
            @test col == name.(txs)
        else
            @test col ==  mat[:, [i-1]] 
        end
    end
    for (i, row) in enumerate(Tables.rows(comm))
        @test row == (; :features => name(txs[i]), (Symbol("sample$(j)") => mat[i, j] for j in 1:5)...)
    end
    @test features(comm) == txs
    for i in 1:5
        @test abundances(comm[:, "sample$i"]) == mat[:, [i]]
        @test abundances(comm["taxon$i", :]) == mat[[i], :]
    end

    tbl = Tables.columntable(comm)
    @test tbl.features == featurenames(comm)
    @test keys(tbl) == (:features, Symbol.(samplenames(comm))...)
    @test tbl.sample1 == abundances(comm[:, "sample1"])
    

end 

# Abundance Tables

# @testset "Abundances" begin
#     # Constructors
#     M = rand(100, 10)
#     df = hcat(DataFrame(x=collect(1:100)), DataFrame(M))

#     a1 = abundancetable(M)
#     a2 = abundancetable(df)
#     @test typeof(a1) <:AbstractComMatrix
#     @test typeof(a2) <:AbstractComMatrix

#     abund = abundancetable(
#         M, ["sample_$x" for x in 1:10],
#         ["feature_$x" for x in 1:100])

#     @test a1.occurrences == a2.occurrences == abund.occurrences

#     # Normalization functions
#     relab_fract = relativeabundance(abund)
#     @test typeof(relab_fract) <: AbstractComMatrix

#     relab_perc = relativeabundance(abund, kind=:percent)
#     @test typeof(relab_perc) <: AbstractComMatrix

#     rnorm = rownormalize(abund)
#     cnorm = colnormalize(abund)

#     @test size(abund) == (nfeatures(abund), nsamples(abund))
#     @test size(relab_fract) == (nfeatures(relab_fract), nsamples(relab_fract))
#     @test size(relab_perc) == (nfeatures(relab_perc), nsamples(relab_perc))
#     @test size(rnorm) == (nfeatures(rnorm), nsamples(rnorm))
#     @test size(cnorm) == (nfeatures(cnorm), nsamples(cnorm))

#     for j in 1:10
#         @test sum(getsample(relab_fract, j)) ≈ 1
#         @test sum(getsample(relab_perc, j)) ≈ 100
#     end

#     for i in 1:nfeatures(rnorm)
#         @test maximum(getfeature(rnorm,i)) ≈ 1
#     end

#     for j in 1:nsamples(cnorm)
#         @test maximum(getsample(cnorm, j)) ≈ 1
#     end

#     # Filtering Functions
#     filt = filterabund(relab_fract, 5)
#     @test typeof(filt) <: AbstractComMatrix
#     @test typeof(filterabund(df, 5)) <: AbstractComMatrix

#     @test size(filt) == (6, 10)
#     for i in 1:10
#         @test sum(getsample(filt, i)) ≈ 1
#     end

#     @test featurenames(filt)[end] == "other"

#     @test present(0.1, 0.001)
#     @test !present(0.001, 0.1)
#     @test present(rand(), 0.)

#     a = zeros(100)
#     a[randperm(100)[1:10]] .= rand(10)

#     @test prevalence(a, 0.) == 0.1
# end

# @testset "Distances" begin
#     # Constructors
#     Random.seed!(1)
#     M = rand(100, 10)
#     df = hcat(DataFrame(x=collect(1:100)), DataFrame(M))
#     abund = abundancetable(
#         M, ["sample_$x" for x in 1:10],
#         ["feature_$x" for x in 1:100])

#     # Diversity indicies
#     R = 100
#     s1 = rand(R) # high diversity
#     s2 = [i % 10 == 0 ? s1[i] : 0 for i in 1:R] # low diversity
#     s3 = ones(R) # uniform
#     s4 = [1., zeros(R-1)...] # no diversity

#     @test shannon(s1) > shannon(s2)
#     @test shannon(s3) ≈ log(R)
#     @test shannon(s4) ≈ 0.

#     @test ginisimpson(s1) > ginisimpson(s2)
#     @test ginisimpson(s3) ≈ 1. - 1/R
#     @test ginisimpson(s4) ≈ 0.

#     @test length(shannon(abund)) == 10
#     @test length(ginisimpson(abund)) == 10
# end
