"""
    commjoin(c1::CommunityProfile, comms::CommunityProfile...)
    commjoin(profiles::AbstractVector{<:CommunityProfile})

Join multiple `CommunityProfile`s, creating a new `CommunityProfile`.
Samples are concatenated in argument order; features retain their first-occurrence
order, including features with only zero abundances. Sample names must be unique
across all inputs. The result has `Float64` abundances and shares feature and sample
objects with the inputs, but owns its abundance matrix. A single input also
produces an independent abundance matrix.
"""
function commjoin(c1::CommunityProfile, comms::CommunityProfile...)
    return _commjoin((c1, comms...))
end

function commjoin(profiles::AbstractVector{<:CommunityProfile})
    isempty(profiles) && throw(ArgumentError("At least one profile is required."))
    return _commjoin(profiles)
end

_join_samples(profiles::Tuple) = vcat(samples.(profiles)...)
_join_samples(profiles::AbstractVector) = reduce(vcat, samples.(profiles))
_join_features(profiles::Tuple) = vcat(features.(profiles)...)
_join_features(profiles::AbstractVector) = reduce(vcat, features.(profiles))

function _commjoin(profiles)
    all_samples = _join_samples(profiles)
    seen_samples = Set{String}()
    for sample in all_samples
        sample_name = name(sample)
        sample_name in seen_samples && error("Duplicate sample names detected: $sample_name")
        push!(seen_samples, sample_name)
    end
    all_features = unique(_join_features(profiles))
    feature_dict = dictionary(zip(all_features, eachindex(all_features)))

    stored_entries = sum(comm -> nnz(abundances(comm)), profiles)
    rows, cols, values = Int[], Int[], Float64[]
    sizehint!(rows, stored_entries)
    sizehint!(cols, stored_entries)
    sizehint!(values, stored_entries)
    coloffset = 0
    for comm in profiles
        _append_profile_entries!(rows, cols, values, feature_dict, comm, coloffset)
        coloffset += nsamples(comm)
    end

    mat = sparse(rows, cols, values, length(all_features), length(all_samples))
    return CommunityProfile(mat, all_features, all_samples)
end

function _append_profile_entries!(rows, cols, values, feature_dict, comm::CommunityProfile, coloffset)
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
    return nothing
end
