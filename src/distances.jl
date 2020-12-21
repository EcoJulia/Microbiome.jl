## Diversity Measures

function shannon(v::AbstractVector{T}) where T<:Real
    total = sum(v)
    relab = map(x-> x/total, v)
    return -sum([log(x^x) for x in relab])
end

function shannon(abt::AbstractAbundanceTable)
    return (; (s=> shannon(abt[:, s]) for s in Symbol.(samplenames(abt)))...)
end

function ginisimpson(v::AbstractVector{T}) where T<:Real
    total = sum(v)
    relab = map(x-> x/total, v)
    return 1 - sum([x^2 for x in relab])
end

function ginisimpson(abt::AbstractAbundanceTable)
    return (; (s=> ginisimpson(abt[:, s]) for s in Symbol.(samplenames(abt)))...)
end

braycurtis(abt::AbstractAbundanceTable) = pairwise(BrayCurtis(), collect(abundances(abt)), dims=2)
pcoa(abt::AbstractAbundanceTable) = fit(MDS, braycurtis(abt), distances = true)