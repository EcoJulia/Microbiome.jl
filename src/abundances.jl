# Methods for absolute and relative abundances

abundancetable(df::DataFrame) = ComMatrix(
    convert(Matrix{Float64}, df[2:end]),
        string.(df[1]),
        string.(names(df[2:end]))
    )

abundancetable(table::AbstractArray{T,2},
    site = ["sample_$x" for x in axes(table, 2)],
    species = ["feature_$x" for x in axes(table, 1)]
    ) where T<:Real = ComMatrix(Float64.(table), species, site)


"""
Filter an abundance table to the top `n` features accross all samples

This function also adds a row for "other", which sums the abundances of the
remaining features.
"""
function filterabund(abun::AbstractComMatrix, n::Int=minimum(10, nfeatures(abun)))
    # TODO: add prevalence filter

    totals = featuretotals(abun)

    srt = sortperm(totals, rev=true)

    newabun = getfeature(abun, srt[1:n])

    remainder = [sum(occurrences(abun)[srt[n+1:end], i]) for i in 1:size(abun, 2)]'
    newabun = vcat(newabun, remainder)
    newrows = cat(featurenames(abun)[srt[1:n]], ["other"], dims=1)

    return abundancetable(newabun, samplenames(abun), newrows)
end

filterabund(df::DataFrame, n::Int=10) = filterabund(abundancetable(df), n)

function rownormalize!(abt::AbstractComMatrix)
    for i in 1:size(abt, 1)
        rowmax = maximum(occurrences(abt)[i,:])
        for j in 1:size(abt, 2)
            occurrences(abt)[i,j] /= rowmax
        end
    end
end

function rownormalize(abt::AbstractComMatrix)
    abt = deepcopy(abt)
    rownormalize!(abt)
    return abt
end

function colnormalize!(abt::AbstractComMatrix)
    for j in 1:size(abt, 2)
        colmax = maximum(occurrences(abt)[:,j])
        for i in 1:size(abt, 1)
            occurrences(abt)[i,j] /= colmax
        end
    end
end

function colnormalize(abt::AbstractComMatrix)
    abt = deepcopy(abt)
    rownormalize!(abt)
    return abt
end

function relativeabundance!(a::AbstractComMatrix; kind::Symbol=:fraction)
    in(kind, [:percent, :fraction]) || error("Invalid kind: $kind")
    if eltype(a.occurrences) != Float64; a.occurrences = Float64.(a.occurrences) end
    for i in 1:nsamples(a)
        s = sum(getsample(a,i))
        s == 0 && continue
        for x in 1:nfeatures(a)
            kind == :fraction ? occurrences(a)[x,i] /= s : occurrences(a)[x,i] /= (s / 100.)
        end
    end
end


function relativeabundance(a::AbstractComMatrix; kind::Symbol=:fraction)
    relab = deepcopy(a)
    relativeabundance!(relab, kind=kind)
    return relab
end

"""
Check if a given (non-zero) value is greater than a minimum value.
If the minimum abundance is 0, just checks if value is non-zero
"""
present(t::Float64, minabundance::Float64=0.0001) = t == 0 ? false : t >= minabundance

"""
Return the fraction of values that are greater than a minimum
"""
prevalence(a::AbstractArray{<:Real}, minabundance::Float64=0.0001) = mean(x-> present(x, minabundance), a)
