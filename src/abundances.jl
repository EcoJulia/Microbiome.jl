# Methods for absolute and relative abundances

# """
#     filterabund(abun::AbstractComMatrix, n::Int=minimum(10, nfeatures(abun)))

# Filter an abundance table to the top `n` features accross all samples

# This function also adds a row for "other", which sums the abundances of the
# remaining features.
# """
# function filterabund(abun::AbstractAbundanceTable, n::Int=minimum(10, nfeatures(abun)))
#     # TODO: add prevalence filter

#     totals = featuretotals(abun)

#     srt = sortperm(totals, rev=true)

#     newabun = getfeature(abun, srt[1:n])

#     remainder = [sum(occurrences(abun)[srt[n+1:end], i]) for i in 1:size(abun, 2)]'
#     newabun = vcat(newabun, remainder)
#     newrows = cat(featurenames(abun)[srt[1:n]], ["other"], dims=1)

#     return abundancetable(newabun, samplenames(abun), newrows)
# end

"""
    relativeabundance!(a::AbstractAbundanceTable; kind::Symbol=:fraction)

Normalize each sample in AbstractAbundanceTable to the sum of the sample.

By default, columns sum to 1.0.
Use `kind=:percent` for columns to sum to 100.
"""
function relativeabundance!(at::AbstractAbundanceTable; kind::Symbol=:fraction)
    in(kind, [:percent, :fraction]) || throw(ArgumentError("Invalid kind: $kind"))
    eltype(abundances(at)) <: AbstractFloat || throw(ArgumentError("relativeabundance! requires profile to have AbstractFloat eltype. Try relativeabundance instead"))
    abund = abundances(at)
    abund ./= sampletotals(at)
    kind == :percent && (abund .*= 100)
    dropzeros!(abund) # shouldn't need dropzeros - https://github.com/JuliaLang/julia/issues/39018
    return at
end

"""
    relativeabundance(at::AbstractAbundanceTable, kind::Symbol=:fraction)

Like [`relativeabundance!`](@ref), but does not mutate original.
"""
function relativeabundance(at::T, kind::Symbol=:fraction) where T <: AbstractAbundanceTable
    comm = T(float.(abundances(at)), deepcopy(features(at)), deepcopy(samples(at)))
    relativeabundance!(comm)
end

"""
    present(t::Union{Real, Missing}, minabundance::Real=0.0)
    present(at::AbstractAbundanceTable, minabundance::Real=0.0)

Check if a given (non-zero) value is greater than or equal to a minimum value.
If the minimum abundance is 0, just checks if value is non-zero.

If used on an `AbstractAbundanceTable`, returns a sparse boolean matrix of the same size.
"""
function present(t::Real, minabundance::Real=0.0)
    (minabundance >= 0 && t >= 0) || throw(DomainError("Only defined for positive values"))
    t == 0 ? false : t >= minabundance
end

present(::Missing, m) = missing

function present(at::AbstractAbundanceTable, minabundance::Real=0.0)
    mat = spzeros(Bool, Ssize(at)...)
    for i in eachindex(mat)
        mat[i] = present(at[i], minabundance)
    end
    return mat
end


"""
    prevalence(a::AbstractArray{<:Real}, minabundance::Real=0.0)
    prevalence(at::AbstractAbundanceTable, minabundance::Real=0.0)

Return the fraction of values that are greater than or equal to a minimum.
If the minimum abundance is 0, returns the fraction of non-zero values.

If used on an `AbstractAbundanceTable`,
returns a prevalence value for each `feature` accross the `sample`s.
"""
prevalence(a::AbstractArray{<:Real}, minabundance::Real=0.0) = mean(x-> present(x, minabundance), a)

# makes it work for any iterable
prevalence(a, minabundance::Real=0.0) = mean(x-> present(x, minabundance), (y for y in a))

function prevalence(at::AbstractAbundanceTable, minabundance::Real=0.0)
    mean(x-> present(x, minabundance), abundances(at), dims=2)
end
