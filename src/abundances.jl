# Methods for absolute and relative abundances

struct AbundanceTable{T<:Real} end #TODO placeholder

abundancetable(df::DataFrame) = ComMatrix(df, sitecolumns = true)
abundancetable(table::Array{T,2}) where T<:Real = ComMatrix(table')

@forward_func AbundanceTable.table Base.getindex, Base.setindex, Base.length, Base.size
#do something about these forward funcions

SpatialEcology.show(io::IO, com::AbstractComMatrix) = show(io, full(com.occurrences))

"""
Filter an abundance table to the top `n` species accross all samples

This function also adds a row for "other", which sums the
"""
function filterabund(abun::AbstractComMatrix, n::Int=10)
    srt = sortperm(total_abundance_species(abun), rev=true)
    keep, compact = srt[1:n], srt[n+1:end]
    # TODO: add prevalence filter

    newabun = view(abun, species = keep)
    remainder = total_abundance_site(view(abun, species = compact))

    ComMatrix(hcat(newabun.occurrences, remainder), [specnames(newabun); "other"] , sitenames(abun))
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
