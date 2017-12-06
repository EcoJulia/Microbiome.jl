# Methods for absolute and relative abundances

abundancetable(df::DataFrame) = ComMatrix(df, sitecolumns = true)
abundancetable(table::Array{T,2}) where T<:Real = ComMatrix(table')

@forward_func AbundanceTable.table Base.getindex, Base.setindex, Base.length, Base.size
#do something about these forward funcions

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
