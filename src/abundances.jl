# Methods for absolute and relative abundances

abundancetable(df::DataFrame) = ComMatrix(Matrix(df[2:end]), string.(df[1]), String.(names(df[2:end])))
abundancetable(table::AbstractArray{T,2}, site = ["sample_$x" for x in indices(table, 2)],
    species = ["feature_$x" for x in indices(table, 1)]) where T<:Real =
    ComMatrix(table, species, site)

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

    remainder = [sum(abun.occurrences[srt[n+1:end], i]) for i in 1:size(abun, 2)]'
    newabun = vcat(newabun, remainder)
    newrows = cat(1, featurenames(abun)[srt[1:n]], ["other"])

    return abundancetable(newabun, samplenames(abun), newrows)
end

filterabund(df::DataFrame, n::Int=10) = filterabund(abundancetable(df), n)

function rownormalize!(abt::AbstractComMatrix)
    for i in 1:size(abt, 1)
        rowmax = maximum(abt.occurrences[i,:])
        for j in 1:size(abt, 2)
            abt.occurrences[i,j] /= rowmax
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
        colmax = maximum(abt.occurrences[:,j])
        for i in 1:size(abt, 1)
            abt.occurrences[i,j] /= colmax
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

    for i in 1:nsamples(a)
        s = sum(getsample(a,i))
        for x in 1:nfeatures(a)
            kind == :fraction ? a.occurrences[x,i] /= s : a.occurrences[x,i] /= (s / 100.)
        end
    end
end


function relativeabundance(a::AbstractComMatrix; kind::Symbol=:fraction)
    relab = deepcopy(a)
    relativeabundance!(relab, kind=kind)
    return relab
end
