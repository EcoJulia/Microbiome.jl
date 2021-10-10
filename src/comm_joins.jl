"""
    commjoin(c1::CommunityProfile, comms::CommunityProfile...)

Join multiple `CommunityProfile`s, creating a new `CommunityProfile`.
For now, sample names cannot overlap in any of the input profiles.

```jldoctest
julia> mss = [MicrobiomeSample(string("sample",i)) for i in 1:15];

julia> txs = [Taxon(string("taxon",i)) for i in 1:20];

julia> cm1 = CommunityProfile(spzeros(10,5), txs[1:10], mss[1:5])
CommunityProfile{Float64, Taxon, MicrobiomeSample} with 10 features in 5 samples

Feature names:
taxon1, taxon2, taxon3...taxon9, taxon10

Sample names:
sample1, sample2, sample3, sample4, sample5



julia> cm2 = CommunityProfile(spzeros(10,5), txs[6:15], mss[6:10])
CommunityProfile{Float64, Taxon, MicrobiomeSample} with 10 features in 5 samples

Feature names:
taxon6, taxon7, taxon8...taxon14, taxon15

Sample names:
sample6, sample7, sample8, sample9, sample10



julia> cm3 = CommunityProfile(spzeros(10,5), txs[11:20], mss[11:15])
CommunityProfile{Float64, Taxon, MicrobiomeSample} with 10 features in 5 samples

Feature names:
taxon11, taxon12, taxon13...taxon19, taxon20

Sample names:
sample11, sample12, sample13, sample14, sample15



julia> commjoin(cm1, cm2, cm3)
CommunityProfile{Float64, Taxon, MicrobiomeSample} with 20 features in 15 samples

Feature names:
taxon1, taxon2, taxon3...taxon19, taxon20

Sample names:
sample1, sample2, sample3...sample14, sample15
```
"""
function commjoin(c1::CommunityProfile, comms::CommunityProfile...)
    length(intersect(samplenames(c1), samplenames.(comms)...)) == 0 || error("Duplicate sample names detected: $(intersect(samplenames(c1), samplenames.(comms)...))")

    all_samples = vcat(samples(c1), samples.(comms)...)
    sample_dict = dictionary(zip(all_samples, eachindex(all_samples)))
    all_features = unique(vcat(features(c1), features.(comms)...))
    feature_dict = dictionary(zip(all_features, eachindex(all_features)))

    mat = spzeros(length(all_features), length(all_samples))
    for comm in (c1, comms...)
        for sample in samples(comm)
            for feature in features(comm)
                mat[feature_dict[feature], sample_dict[sample]] = comm[feature, sample]
            end
        end
    end

    return CommunityProfile(mat, all_features, all_samples)
end
