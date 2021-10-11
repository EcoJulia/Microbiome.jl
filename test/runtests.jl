using Microbiome
using Test
using Microbiome.SparseArrays
using Microbiome.Tables
using Microbiome.Dictionaries
import Microbiome.MultivariateStats: MDS
using Documenter

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

        ms2 = MicrobiomeSample("sample2"; age=10, birthtype="vaginal", allergies=true)
        @test ms2.age == 10
        @test ms2.birthtype == "vaginal"
        @test ms2.allergies
    end
    
    @testset "Taxa" begin
        _ranks = (:domain, :kingdom, :phylum, :class, :order, :family, :genus, :species, :subspecies, :strain)
        txm = Taxon("taxon", missing)
        @test txm === Taxon("taxon")
        @test ismissing(taxrank(txm))
        @test !hasrank(txm)

        for (i, c) in enumerate(_ranks)
            tx = Taxon("taxon", c)
            @test taxrank(tx) == c
            @test tx === Taxon("taxon", i-1)
            @test tx !== txm
        end
        
        @test_throws ErrorException Taxon("taxon", :invalid)
        @test_throws ErrorException Taxon("taxon", 10)
        @test let tx = Taxon("taxon", :kingdom)
            hasrank(tx)
        end
    end

    @testset "Gene Functions" begin
        gfm = GeneFunction("gene", missing)
        @test name(gfm) == "gene"
        @test gfm === GeneFunction("gene")
        @test ismissing(taxon(gfm))
        @test !hastaxon(gfm)
        @test !hasrank(gfm)

        gf1 = GeneFunction("gene", Taxon("sp1", :species))
        gf2 = GeneFunction("gene", Taxon("sp1"))
        @test name(gf1) == "gene"
        
        @test gf1 != gfm
        @test gf1 === GeneFunction("gene", Taxon("sp1", :species))
        @test gf2 === GeneFunction("gene", Taxon("sp1"))
        @test gf1 !== gf2
        @test hastaxon(gf1)
        @test !ismissing(taxon(gf1))
        @test taxon(gf1) == Taxon("sp1", :species)
        @test taxon(gf1) != taxon(gf2)
        @test hasrank(gf1)
        @test taxrank(gf1) == :species
    end
end

@testset "Profiles" begin
    _ranks = Tuple(keys(Microbiome._ranks))[1:9]
    mss = [MicrobiomeSample("sample$i") for i in 1:5]
    txs = [Taxon("taxon$i", _ranks[i]) for i in 1:9]
    push!(txs, Taxon("taxon10", missing))
    
    mat = spzeros(10,5)
    dmat = zeros(10,5)
    for i in 1:5; 
        mat[i,i] = 1.
        dmat[i,i] = 1.
    end
    for i in 1:5
        mat[i+5,i] = 0.6
        dmat[i+5,i] = 0.6
    end

    comm = CommunityProfile(mat, txs, mss)

    @testset "Profile operations" begin
        @test repr(comm) == "CommunityProfile{Float64, Taxon, MicrobiomeSample} with 10 features in 5 samples\n\nFeature names:\ntaxon1, taxon2, taxon3...taxon9, taxon10\n\nSample names:\nsample1, sample2, sample3, sample4, sample5\n\n"
        @test CommunityProfile{Float64, Taxon, MicrobiomeSample}(mat, txs, mss) isa CommunityProfile
        @test comm == CommunityProfile(dmat, txs, mss)
        
        @test nsamples(comm) == 5
        @test nfeatures(comm) == 10
        @test size(comm) == (10, 5)
        @test profiletype(comm) == Taxon
        @test ranks(comm)[1:9] == [_ranks...]
        @test sampletotals(comm) == [1.6 1.6 1.6 1.6 1.6]
        @test featuretotals(comm) == reshape([1.0, 1.0, 1.0, 1.0, 1.0, 0.6, 0.6, 0.6, 0.6, 0.6], 10, 1)
        @test features(comm) == txs
        @test samples(comm) == mss

        @test size(rankfilter(comm, :species), 1) == 1
        @test size(rankfilter(comm, :genus; keepempty=true), 1) == 2

        @test size(rankfilter(comm, 5), 1) == 1
        @test size(rankfilter(comm, 6; keepempty=true), 1) == 2
        @test_throws ErrorException rankfilter(comm, :foo)
        @test_throws ErrorException rankfilter(comm, 10)
        @test_throws ErrorException rankfilter(rankfilter(comm, :species), :genus) # will be empty

        @test featurenames(filter(f-> hasrank(f) && taxrank(f) == :species, comm)) == featurenames(rankfilter(comm, :species))

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

        @test_throws ErrorException commjoin(comm, comm)
        let c3 = commjoin(comm[:,1:2], comm[:, 3:4], comm[:, 5])
            @test abundances(c3) == abundances(comm)
            @test samples(c3) == samples(comm)
            @test features(c3) == features(comm)
        end

        filtertest = CommunityProfile(sparse(Float64[3 2 1 # 0.66, assuming minabundance 2
                                              2 2 2 # 1.0
                                              0 0 1 # 0.0
                                              2 0 0 # 0.33
                                       ]),
                               [Taxon(string(i)) for i in 1:4],
                               [MicrobiomeSample(string(i)) for i in 1:3]); 
        @test size(prevalence_filter(filtertest)) == (4,3)
        @test size(prevalence_filter(filtertest, minabundance=2)) == (3,3)
        @test size(prevalence_filter(filtertest, minabundance=2, minprevalence=0.4)) == (2,3)
        @test all(<=(1.0), abundances(prevalence_filter(filtertest, renorm=true)))
        @test all(x-> isapprox(x, 1.0, atol=1e-8), sum(abundances(prevalence_filter(filtertest, renorm=true)), dims=1))
    end 

    @testset "Stratified gene functions operations" begin
        strat = CommunityProfile(sparse(Float64[3 2 1 
                                                2 2 2
                                                0 0 1
                                                0 1 0
                                        ]),
                                        [GeneFunction("gene1"),
                                         GeneFunction("gene1", "species1"),
                                         GeneFunction("gene1", "species2"),
                                         GeneFunction("gene2")
                                         ],
                                        [MicrobiomeSample(string(i)) for i in 1:3]);
        
        @test filter(hastaxon, strat)         |> nfeatures == 2
        @test filter(!hastaxon, strat)        |> nfeatures == 2
        @test strat["gene1", :]               |> nfeatures == 3
        @test strat[["gene1", "gene2"], :]    |> nfeatures == 4
        @test strat[GeneFunction("gene1"), :] |> nfeatures == 1
    end


    @testset "Profile Metadata" begin
        s1 = MicrobiomeSample("sample1", Dictionary(Dict(:age=> 37, :name=>"kevin", :something=>1.0)))
        s2 = MicrobiomeSample("sample2", Dictionary(Dict(:age=> 37, :name=>"kevin", :something_else=>2.0, :still_other=>"boo")))

        @testset "Single sample" begin
            c4 = CommunityProfile(sparse([1 1; 2 2; 3 3]), [Taxon(string(i)) for i in 1:3], [s1, s2])
            md1, md2 = metadata(c4)

            @test_throws Dictionaries.IndexError insert!(c4, "sample1", :something, 3.0)
            @test_throws Dictionaries.IndexError delete!(c4, "sample1", :something_else)

            @test delete!(c4, "sample1", :something) isa MicrobiomeSample
            @test !haskey(c4, "sample1", :something)
            @test get(c4, "sample1", :something_else, 42) == 42
            @test insert!(c4, "sample1", :something, 3.0) isa MicrobiomeSample
            @test get(c4, "sample1", :something, 42) == 3.0
            set!(c4, "sample1", :something, 1.0)
            @test first(metadata(c4))[:something] == 1.0
            @test haskey(c4, "sample1", :something)
            @test unset!(c4, "sample1", :something) isa MicrobiomeSample
            @test !haskey(c4, "sample1", :something)
            set!(c4, "sample1", :something, 1.0)

            @test collect(keys(c4, "sample1")) == [:age, :name, :something]
        end

        @testset "Whole community" begin
            c5 = CommunityProfile(sparse([1 1; 2 2; 3 3]), [Taxon(string(i)) for i in 1:3], [s1, s2])
            md1, md2 = metadata(c5)
            
            @test all(row-> row[:age] == 37, [md1, md2])
            @test all(row-> row[:name] == "kevin", [md1, md2])
            @test md1[:something] == 1.0
            @test ismissing(md2[:something])
            @test md2[:something_else] == 2.0
            @test ismissing(md1[:something_else])

            @test_throws IndexError add_metadata!(c5, "sample1", Dict(:something=>3.0))
            add_metadata!(c5, "sample1", Dict(:something=>3.0), overwrite=true)
            @test first(metadata(c5))[:something] == 3
            
            @test_throws IndexError add_metadata!(c5, "sample1", (; something=4.0))
            add_metadata!(c5, "sample1", (; something=4.0), overwrite=true)
            @test first(metadata(c5))[:something] == 4

            tbl = [(sample="sample1", something=5,  newthing="bar"),
                   (sample="sample2", something=10, newthing="baz"),
                   (sample="sample3", something=42, newthing="fuz")]
            @test_throws IndexError add_metadata!(c5, :sample, tbl, overwrite=true) # for sample that doesn't exist
            @test_throws IndexError add_metadata!(c5, :sample, tbl[1:2])            # for metadata that already exists
            add_metadata!(c5, :sample, tbl[1:2]; overwrite = true)
            @test first(metadata(c5))[:something] == 5
            @test first(metadata(c5))[:newthing] == "bar"
        end
    end
    
    @testset "Indexing and Tables integration" begin
        for i in 1:5
            @test abundances(comm[:, "sample$i"]) == mat[:, [i]]
            @test abundances(comm["taxon$i", :]) == mat[[i], :]
        end

        @test abundances(comm[r"taxon1", :]) == abundances(comm[["taxon1", "taxon10"], :]) == abundances(comm[[1,10], :])
        @test abundances(comm[:, r"sample[13]"]) == abundances(comm[:,["sample1", "sample3"]]) == abundances(comm[:, [1,3]])

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

doctest(Microbiome)