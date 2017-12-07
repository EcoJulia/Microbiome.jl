# Methods for absolute and relative abundances

abundancetable(df::DataFrame) = ComMatrix(df, sitecolumns = true)
abundancetable(table::Array{T,2}, site = ["sample_$x" for x in indices(table, 2)],
    species = ["feature_$x" for x in indices(table, 1)]) where T<:Real =
    ComMatrix(table', species, site)

@forward_func ComMatrix.occurrences Base.length

"""
Filter an abundance table to the top `n` species accross all samples

This function also adds a row for "other", which sums the abundances of
the remaining species
"""
function filterabund(abun::AbstractComMatrix, n::Int=minimum(10, nspecies(abun)))
    srt = sortperm(total_abundance_species(abun), rev=true)
    keep, compact = srt[1:n], srt[n+1:end]
    # TODO: add prevalence filter

    newabun = view(abun, species = keep)
    remaining_abundances = total_abundance_site(view(abun, species = compact))

    ComMatrix(hcat(newabun.occurrences, remaining_abundances), [specnames(newabun); "other"] , sitenames(abun))
end

filterabund(df::DataFrame, n::Int=10) = filterabund(abundancetable(df), n)

"""
Express abundances as relative abundances, summing to 1 at each site (or to 100
if `kind` is `:percent`)
"""
function relativeabundance!(abun::AbstractComMatrix{Float64}; kind::Symbol=:fraction)
    in(kind, [:percent, :fraction]) || error("Invalid kind: $kind")

    for j in 1:nspecies(abun)
       s = sum(getspecies(abun, j))
       for i in 1:nsites(abun)
           abun.occurrences[i,j] = abun.occurrences[i,j] / s
       end
    end

    kind == :percent && (abun.occurrences .*= 100)
    abun
end

relativeabundance(abun::AbstractComMatrix) =
    relativeabundance!(ComMatrix(convert(Array{Float}, abun.occurrences),
                                 specnames(abun), sitenames(abun)))
