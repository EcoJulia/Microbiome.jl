## Diversity Measures

function shannon(v::AbstractVector{T}) where T<:Real
    total = sum(v)
    relab = map(x-> x/total, v)
    return -sum([log(x^x) for x in relab])
end

function shannon(abt::AbstractComMatrix)
    return shannon.(eachcol(occurrences(abt)))
end

function ginisimpson(v::AbstractVector{T}) where T<:Real
    total = sum(v)
    relab = map(x-> x/total, v)
    return 1 - sum([x^2 for x in relab])
end

function ginisimpson(abt::AbstractComMatrix)
    return ginisimpson.(eachcol(occurrences(abt)))
end
