function _dense_join_reference(profiles)
    feats = unique(vcat(features.(profiles)...))
    smpls = vcat(samples.(profiles)...)
    matrix = zeros(Float64, length(feats), length(smpls))
    offset = 0
    for profile in profiles
        for (row, feature) in enumerate(features(profile))
            target = findfirst(isequal(feature), feats)
            for col in axes(abundances(profile), 2)
                matrix[target, offset + col] = abundances(profile)[row, col]
            end
        end
        offset += nsamples(profile)
    end
    return feats, smpls, matrix
end

_random_join_values(rng, ::Type{Bool}, m, n) = rand(rng, Bool, m, n)
_random_join_values(rng, ::Type{T}, m, n) where {T} = T.(rand(rng, -3:3, m, n))

@testset "Sparse community joins" begin
    a, b, c, d, e = Taxon.(["a", "b", "c", "d", "e"])
    s1 = MicrobiomeSample("s1"; cohort = "ROC")
    s2 = MicrobiomeSample("s2"; cohort = "OOM")
    s3 = MicrobiomeSample("s3"; visit = 2)
    firstprofile = CommunityProfile([1 0; 0 2; 3 0; 0 0], [a, c, b, e], [s2, s1])
    secondprofile = CommunityProfile(Float32[4, 5, 6], [d, b, a], s3)
    before = copy(abundances(firstprofile))
    joined = commjoin(firstprofile, secondprofile)
    @test features(joined) == [a, c, b, e, d]
    @test samplenames(joined) == ["s2", "s1", "s3"]
    @test Matrix(abundances(joined)) == [1 0 6; 0 2 0; 3 0 5; 0 0 0; 0 0 4]
    @test eltype(abundances(joined)) == Float64
    @test nnz(abundances(joined)) == 6
    @test abundances(firstprofile) == before
    @test get(samples(joined)[1]) === get(s2)
    @test get(samples(joined)[2]) === get(s1)
    @test get(samples(joined)[3]) === get(s3)

    single = commjoin(firstprofile)
    @test single == firstprofile
    @test abundances(single) !== abundances(firstprofile)
    @test eltype(abundances(single)) == Float64
    @test commjoin([firstprofile, secondprofile]) == joined
    @test commjoin(view([firstprofile, secondprofile], :)) == joined
    @test commjoin(CommunityProfile[firstprofile, secondprofile]) == joined
    @test commjoin([firstprofile]) == single
    @test_throws ArgumentError commjoin(CommunityProfile[])
    @test_throws ErrorException commjoin([firstprofile, firstprofile])
    @test_throws ErrorException commjoin(firstprofile, firstprofile)
    thirdprofile = CommunityProfile([1], [a], MicrobiomeSample("s4"))
    @test_throws ErrorException commjoin(firstprofile, firstprofile, thirdprofile)
    @test_throws ErrorException commjoin(thirdprofile, firstprofile, firstprofile)

    nofeatures = CommunityProfile(spzeros(0, 1), Taxon[], [MicrobiomeSample("empty")])
    nosamples = CommunityProfile(spzeros(1, 0), [d], MicrobiomeSample[])
    @test size(commjoin(nofeatures, nosamples)) == (1, 1)
    @test features(commjoin(nofeatures, nosamples)) == [d]
    @test iszero(only(abundances(commjoin(nofeatures, nosamples))))
    @test size(commjoin(nosamples)) == (1, 0)
    @test size(commjoin(nofeatures)) == (0, 1)

    special = CommunityProfile(sparse([1, 2], [1, 1], [-0.0, NaN], 3, 1), [a, b, c], [s1])
    infinities = CommunityProfile([Inf, -Inf], [a, d], s2)
    combined = commjoin(special, infinities)
    @test isequal(Matrix(abundances(combined)), [0.0 Inf; NaN 0.0; 0.0 0.0; 0.0 -Inf])
    @test nnz(abundances(combined)) == 3
    @test nnz(abundances(special)) == 2

    indexed32 = SparseMatrixCSC{Float32, Int32}(sparse(Float32[0 2; 1 0]))
    viewed = CommunityProfile(indexed32, view([a, b, c], [2, 1]), [s1, s2])
    @test commjoin(viewed, secondprofile) == commjoin(CommunityProfile(Matrix(indexed32), [b, a], [s1, s2]), secondprofile)

    metabolite = Metabolite("metabolite", "annotation", 1.5, 2.5)
    metaboliteprofile = CommunityProfile([7.0], [metabolite], s3)
    mixed = commjoin(firstprofile, metaboliteprofile)
    @test features(mixed) == AbstractFeature[a, c, b, e, metabolite]
    @test Matrix(abundances(mixed)) == [1 0 0; 0 2 0; 3 0 0; 0 0 0; 0 0 7]
    @test last(features(mixed)) === metabolite

    @testset "Randomized alignment and numeric types" begin
        rng = MersenneTwister(714)
        for T in (Bool, Int16, Float32, Float64, BigFloat), repetition in 1:12
            pool = isodd(repetition) ? Taxon.(string.(1:24)) : GeneFunction.(string.(1:24))
            profiles = map(1:rand(rng, 2:5)) do index
                nf, ns = rand(rng, 0:20), rand(rng, 0:3)
                feats = pool[randperm(rng, length(pool))[1:nf]]
                smpls = [MicrobiomeSample("$(index)_$(col)"; group = index) for col in 1:ns]
                matrix = _random_join_values(rng, T, nf, ns)
                matrix[rand(rng, size(matrix)...) .< 0.7] .= zero(T)
                CommunityProfile(matrix, feats, smpls)
            end
            expected_features, expected_samples, expected = _dense_join_reference(profiles)
            output = commjoin(profiles...)
            @test commjoin(profiles) == output
            @test features(output) == expected_features
            @test samples(output) == expected_samples
            @test Matrix(abundances(output)) == expected
            @test eltype(abundances(output)) == Float64
            @test nnz(abundances(output)) == count(!iszero, expected)
            @test all(get(result) === get(original) for (result, original) in zip(samples(output), expected_samples))
        end
    end
end
