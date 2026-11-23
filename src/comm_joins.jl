"""
    commjoin(c1::CommunityProfile, comms::CommunityProfile...)

Join multiple `CommunityProfile`s, creating a new `CommunityProfile`.
Samples are concatenated in argument order; features retain their first-occurrence
order, including features with only zero abundances. Sample names must be unique
across all inputs. The result has `Float64` abundances and shares feature and sample
objects with the inputs, but owns its abundance matrix. A single input also
produces an independent abundance matrix.
"""
function commjoin(c1::CommunityProfile, comms::CommunityProfile...)
    all_samples = vcat(samples(c1), samples.(comms)...)
    seen_samples = Set{String}()
    for sample in all_samples
        sample_name = name(sample)
        sample_name in seen_samples && error("Duplicate sample names detected: $sample_name")
        push!(seen_samples, sample_name)
    end
    all_features = unique(vcat(features(c1), features.(comms)...))
    feature_dict = dictionary(zip(all_features, eachindex(all_features)))

    profiles = (c1, comms...)
    stored_entries = sum(comm -> nnz(abundances(comm)), profiles)
    rows, cols, values = Int[], Int[], Float64[]
    sizehint!(rows, stored_entries)
    sizehint!(cols, stored_entries)
    sizehint!(values, stored_entries)
    coloffset = 0
    for comm in profiles
        rowmap = [feature_dict[feature] for feature in features(comm)]
        matrix = abundances(comm)
        localrows = rowvals(matrix)
        localvalues = nonzeros(matrix)
        for col in axes(matrix, 2)
            for entry in nzrange(matrix, col)
                value = Float64(localvalues[entry])
                if !iszero(value)
                    push!(rows, rowmap[localrows[entry]])
                    push!(cols, coloffset + col)
                    push!(values, value)
                end
            end
        end
        coloffset += nsamples(comm)
    end

    mat = sparse(rows, cols, values, length(all_features), length(all_samples))
    return CommunityProfile(mat, all_features, all_samples)
end
