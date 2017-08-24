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

## For use with abundance tables generated from Humman2
