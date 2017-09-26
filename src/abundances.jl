struct AbundanceTable{T<:Real} <: AbstractArray{T,2}
    t::Array{T,2}
    samples::Vector{S} where S
    rows::Vector{R} where R
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

@forward_func AbundanceTable.t Base.getindex, Base.setindex, Base.length, Base.size

cat(1, [3],[1,2])

function filterabund(abun::AbundanceTable, n::Int=10)
    2 < n < 12 || error("n must be between 2 and 12")

    totals = [sum(abun[i,1:end]) for i in 1:size(abun, 1)]
    remainder = [100-t for t in totals]

    srt = sortperm(totals, rev=true)

    newabun = abun[srt[1:10], :]
    remainder = [100-sum(newabun[:, i]) for i in 1:size(newabun, 2)]'
    newabun = vcat(newabun, remainder)
    newrows = cat(1, abun.rows[srt[1:10]], ["total"])

    return AbundanceTable(newabun, newrows, abun.samples)
end

filterabund(df::DataFrame, n::Int=10) = filterabund(AbundanceTable(df), n)

## For use with abundance tables generated from Humman2
