function bysample(df::DataFrame, samples::Array{String, 1}, rows::Array{Int,1}=Int[])
    s = [in(String(n), samples) for n in names(df[2:end])]

    length(rows) > 0 ? df = df[rows, [true,s...]] : df = df[:, [true, s...]]

    newrows = DataFrame(samples=String.(names(df[2:end])))
    newcols = Symbol.(df[1])

    df = DataFrame(Matrix(df[2:end])', newcols)
    names!(df, newcols)
    return hcat(newrows, df)
end


function filter_rows(df::DataFrame, quant::Real; kind::Symbol=:percolumn, calc=:gt)
    in(kind, [:percolumn, :perrow]) || error("kind must be :percolumn or :perrow")
    in(calc, [:gt, :lt, :sum, :mean, :max, :min]) || error("calc must be :gt, :lt, :sum, :mean, :max, or :min")

    if kind==:perrow && calc==:sum
        keep = [sum(Matrix(df[i,2:end])) > quant ? true : false for i in 1:size(df, 1)]
    else
        error("that hasn't been implemented yet")
    end

    return df[keep, :]
end

#==============
MetaPhlAn Utils
==============#

"""
taxfilter!(df::DataFrame, level::Int=7; shortnames::Bool=true)

Filter a MetaPhlAn table (df) to a particular taxon level.
1 = Kingdom
2 = Phylum
3 = Class
4 = Order
5 = Family
6 = Genus
7 = Species
8 = Strain

If shortnames is true (default), also changes names in the first column to
remove higher order taxa
"""
function taxfilter(taxonomic_profile::DataFrame, level::Int=7; shortnames::Bool=true)
    filt = taxonomic_profile[length.(split.(taxonomic_profile[1], '|')) .== level, :]
    if shortnames
        matches = collect.(eachmatch.(r"[kpcofgs]__(\w+)", filt[1]))
        filt[1] = String.([m[level].captures[1] for m in matches])
    end
    return filt
end


#==============
PanPhlAn Utils
==============#

function panphlan_calcs(df::DataFrame)
    abun = AbundanceTable(df)
    dm = getdm(df, Jaccard())
    rowdm = getrowdm(df, Jaccard())
    clust_h = hclust(dm.dm, :single)
    clust_v = hclust(rowdm.dm, :single)
    pco = pcoa(dm)

    return abun, dm, clust_h, clust_v, pco
end
