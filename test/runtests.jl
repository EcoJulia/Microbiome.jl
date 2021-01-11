using Microbiome
using Test
using Microbiome.SparseArrays
using Microbiome.Tables
using Microbiome.Dictionaries
import Microbiome.MultivariateStats: MDS

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
    _clades = (:domain, :kingdom, :phylum, :class, :order, :family, :genus, :species, :subspecies)
    mss = [MicrobiomeSample("sample$i") for i in 1:5]
    txs = [Taxon("taxon$i", _clades[i]) for i in 1:9]
    push!(txs, Taxon("taxon10", missing))
    
    mat = spzeros(10,5)
    for i in 1:5; mat[i,i] = 1.; end
    for i in 1:5; mat[i+5,i] = 0.6; end

    comm = CommunityProfile(mat, txs, mss)
    
    @testset "Profile operations" begin
        @test CommunityProfile{Float64, Taxon, MicrobiomeSample}(mat, txs, mss) isa CommunityProfile
        @test nsamples(comm) == 5
        @test nfeatures(comm) == 10
        @test size(comm) == (10, 5)
        @test profiletype(comm) == Taxon
        @test clades(comm)[1:9] == [_clades...]
        @test sampletotals(comm) == [1.6 1.6 1.6 1.6 1.6]
        @test featuretotals(comm) == reshape([1.0, 1.0, 1.0, 1.0, 1.0, 0.6, 0.6, 0.6, 0.6, 0.6], 10, 1)
        @test features(comm) == txs
        @test samples(comm) == mss

        @test present(0.1)
        @test !present(0.1, 0.2)
        @test_throws DomainError present(-0.1)
        @test_throws DomainError present(0., -0.1)

        @test prevalence([0.0, 0.1, 0.2, 0.3]) ≈ 0.75
        @test prevalence([0.0, 0.1, 0.2, 0.3], 0.15) ≈ 0.5

        @test all(≈(0.2), prevalence(comm))
        @test all(prevalence(comm, 0.7) .≈ [0.2, 0.2, 0.2, 0.2, 0.2, 0, 0, 0, 0, 0])

        @test let c2 = deepcopy(comm)
            relativeabundance!(c2)
            @test all(≈(0.625), featuretotals(c2)[1:5])
            all(≈(0.375), featuretotals(c2)[6:end])
        end
        @test let c2 = deepcopy(comm)
            relativeabundance!(c2, kind=:percent)
            @test all(≈(62.5), featuretotals(c2)[1:5])
            all(≈(37.5), featuretotals(c2)[6:end])
        end

        @test_throws ArgumentError relativeabundance!(comm, kind=:invalid)
    end
    
    @testset "Indexing and Tables integration" begin
        for i in 1:5
            @test abundances(comm[:, "sample$i"]) == mat[:, [i]]
            @test abundances(comm["taxon$i", :]) == mat[[i], :]
        end

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

        tbl = Tables.columntable(comm)
        @test tbl.features == featurenames(comm)
        @test keys(tbl) == (:features, Symbol.(samplenames(comm))...)
        @test tbl.sample1 == abundances(comm[:, "sample1"])
    end

    @testset "Diversity" begin
        R = 10
        s1 = collect(0:R-1) # high diversity
        s2 = [i % 2 == 0 ? s1[i] : 0 for i in eachindex(s1)] # low diversity
        s3 = ones(R) # uniform
        s4 = [10, zeros(R-1)...] # no diversity

        @test shannon(s1) > shannon(s2)
        @test shannon(s3) ≈ log(R)
        @test shannon(s4) ≈ 0.0
        for s in (s1, s2, s3, s4)
            shannon(s ./ sum(s)) ≈ shannon(s)
        end

        @test ginisimpson(s1) > ginisimpson(s2)
        @test ginisimpson(s3) ≈ 1.0 - 1/R
        @test ginisimpson(s4) ≈ 0.0
        for s in (s1, s2, s3, s4)
            ginisimpson(s ./ sum(s)) ≈ ginisimpson(s)
        end


        @test size(shannon(comm)) == (1, 5)
        @test size(ginisimpson(comm)) == (1, 5)

        @test let c2 = deepcopy(comm)
            shannon!(c2)
            @test all(s-> haskey(s, :shannon), samples(c2))
            @test_throws Dictionaries.IndexError shannon!(c2)
            shannon!(c2, overwrite=true)
            true
        end
        @test let c2 = deepcopy(comm)
            ginisimpson!(c2)
            @test all(s-> haskey(s, :ginisimpson), samples(c2))
            @test_throws Dictionaries.IndexError ginisimpson!(c2)
            ginisimpson!(c2, overwrite=true)
            true
        end

        @test all(braycurtis(comm) .== [0 1 1 1 1; 1 0 1 1 1; 1 1 0 1 1; 1 1 1 0 1; 1 1 1 1 0])
        @test all(jaccard(comm) .== [0 1 1 1 1; 1 0 1 1 1; 1 1 0 1 1; 1 1 1 0 1; 1 1 1 1 0])
        @test all(hellinger(comm) .== [0 1 1 1 1; 1 0 1 1 1; 1 1 0 1 1; 1 1 1 0 1; 1 1 1 1 0])
        @test pcoa(comm) isa MDS
    end

end
