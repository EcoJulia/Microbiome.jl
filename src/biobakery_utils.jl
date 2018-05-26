"""
Functions for working with files from Biobakery tools.
Note: these functions may not be well tested
"""

function fliptable(df::DataFrame, samples::Array{String,1})
    newrows = DataFrame(samples=samples)
    newcols = Symbol.(df[1])

    df = DataFrame(Matrix(df[2:end])', newcols)
    names!(df, newcols)
    return hcat(newrows, df)
end

function bysample(df::DataFrame, samples::Array{Symbol, 1}, rows::Array{Int,1}=Int[])
    s = [in(n, samples) for n in names(df[2:end])]
    length(rows) > 0 ? df = df[rows, [true,s...]] : df = df[:, [true, s...]]

    fliptable(df, String.(samples))
end

function bysample(df::DataFrame, samples::Array{String, 1}, rows::Array{Int,1}=Int[])
    s = [in(n, samples) for n in String.(names(df[2:end]))]
    length(rows) > 0 ? df = df[rows, [true,s...]] : df = df[:, [true, s...]]

    fliptable(df, samples)
end

bysample(df::DataFrame) = bysample(df, df[1])


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
const taxlevels = Dict([
    :kingom     => 1,
    :phylum     => 2,
    :class      => 3,
    :order      => 4,
    :family     => 5,
    :genus      => 6,
    :species    => 7,
    :subspecies => 8])

function metaphlan_import(path::String; level=0, shortnames::Bool=true)
    df = FileIO.load(path) |> DataFrame
    for n in names(df)
        df[n] = coalesce.(df[n], 0)
    end

    if typeof(level) <: Symbol
        in(level, keys(taxlevels)) || error("$level not a valid taxonomic level")
        level = taxlevels[level]
    end

    level > 0 && taxfilter!(df, level, shortnames=shortnames)
    return abundancetable(df)
end

function metaphlan_import(paths::Array{String,1}; level=0, shortnames::Bool=true)
    tax = DataFrame(SampleID=String[])
    for f in paths
        df = load(f) |> DataFrame
        rename!(df, Symbol("#SampleID"), :SampleID)
        tax = join(tax, df, on=:SampleID, kind=:outer)
    end

    for n in names(tax)
        tax[n] = coalesce.(tax[n], 0)
    end

    if typeof(level) <: Symbol
        in(level, keys(taxlevels)) || error("$level not a valid taxonomic level")
        level = taxlevels[level]
    end

    level > 0 && taxfilter!(tax, level, shortnames=shortnames)
    return abundancetable(tax)
end

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
8 = Subspecies

If shortnames is true (default), also changes names in the first column to
remove higher order taxa
"""
function taxfilter!(taxonomic_profile::DataFrames.DataFrame, level::Int=7; shortnames::Bool=true)
    filter!(row->length(split(row[1], '|')) == level, taxonomic_profile)
    if shortnames
        matches = collect.(eachmatch.(r"[kpcofgs]__(\w+)", taxonomic_profile[1]))
        taxonomic_profile[1] = String.([m[level].captures[1] for m in matches])
    end
    return taxonomic_profile
end

function taxfilter!(taxonomic_profile::DataFrames.DataFrame, level::Symbol; shortnames::Bool=true)
    in(level, keys(taxlevels)) || error("$level not a valid taxonomic level")
    taxfilter!(taxonomic_profile, taxlevels[level], shortnames=shortnames)
end


function taxfilter(taxonomic_profile::DataFrames.DataFrame, level::Int=7; shortnames::Bool=true)
    filt = deepcopy(taxonomic_profile)
    taxfilter!(filt, level, shortnames=shortnames)
    return filt
end

function taxfilter(taxonomic_profile::DataFrames.DataFrame, level::Symbol; shortnames::Bool=true)
    filt = deepcopy(taxonomic_profile)
    taxfilter!(filt, level, shortnames=shortnames)
    return filt
end

"""
Given a dataframe with a column that has a pvalue column, perform
Benjamini Hochberch correction to generate q value column with given Q.
"""
function qvalue!(df::DataFrame, q::Float64=0.2; pcol::Symbol=:p_value, qcol::Symbol=:q_value)
    if eltype(df[pcol]) <:StatsBase.PValue
        ranks = invperm(sortperm(map(x->x.v,df[pcol])))
    else
        ranks = invperm(sortperm(map(x->x,df[pcol])))
    end
    m = length(ranks)
    df[qcol] = [i / m * q for i in eachindex(df[pcol])]
end

#==============
PanPhlAn Utils
==============#

function panphlan_calcs(df::DataFrame)
    abun = abundancetable(df)
    dm = getdm(df, Jaccard())
    rowdm = getrowdm(df, Jaccard())
    col_clust = hclust(dm.dm, :single)
    row_clust = hclust(rowdm.dm, :single)
    optimalorder!(col_clust, dm.dm)
    optimalorder!(row_clust, rowdm.dm)

    pco = pcoa(dm)

    return abun, dm, col_clust, row_clust, pco
end
