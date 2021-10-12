"""
    commjoin(c1::CommunityProfile, comms::CommunityProfile...)

Join multiple `CommunityProfile`s, creating a new `CommunityProfile`.
For now, sample names cannot overlap in any of the input profiles.
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
