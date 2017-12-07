# Methods for absolute and relative abundances

struct AbundanceTable{T<:Real} <: AbstractArray{T,2}
    table::Array{T,2}
    samples::Vector{S} where S
    features::Vector{R} where R
end

function AbundanceTable(df::DataFrame)
    return AbundanceTable(
            Matrix(df[:,2:end]),
            names(df[2:end]),
            Vector(df[1]))
end
function AbundanceTable(table::Array{T,2}) where T<:Real
    return AbundanceTable(
             table,
             Vector{Int64}(1:size(table,2)),
             Vector{Int64}(1:size(table,1)))
end

@forward_func AbundanceTable.table Base.getindex, Base.setindex, Base.length, Base.size


"""
Filter an abundance table to the top `n` species accross all samples

This function also adds a row for "other", which sums the
"""
function filterabund(abun::AbundanceTable, n::Int=10)
    # TODO: add prevalence filter

    totals = [sum(abun[i,:]) for i in 1:size(abun, 1)]

    srt = sortperm(totals, rev=true)

    newabun = abun[srt[1:n], :]

    remainder = [sum(abun[srt[n+1:end], i]) for i in 1:size(abun, 2)]'
    newabun = vcat(newabun, remainder)
    newrows = cat(1, abun.features[srt[1:n]], ["other"])

    return AbundanceTable(newabun, abun.samples, newrows)
end

filterabund(df::DataFrame, n::Int=10) = filterabund(AbundanceTable(df), n)


function relativeabundance(a::AbundanceTable; kind::Symbol=:fraction)
    in(kind, [:percent, :fraction]) || error("Invalid kind: $kind")

    relab = a.table
    for j in 1:size(relab, 2)
        s = sum(relab[:,j])
        for i in 1:size(relab,1)
            @inbounds relab[i,j] = a[i,j] / s
        end
    end
    kind == :percent ? relab = relab .* 100. : true

    return AbundanceTable(relab, a.samples, a.features)
end
