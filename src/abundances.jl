struct AbundanceTable{T<:Real} <: AbstractArray{T,2}
    t::AbstractArray{T,2}
    samples::AbstractVector
    rows::AbstractVector

    function AbundanceTable(t, samples, rows)
        length(samples) == size(t, 2) || error("# of samples ≠ # of columns")
        length(rows) == size(t, 1) || error("# of row names ≠ # of rows")
        new(t, samples, rows)
    end
end

AbundanceTable(df::DataFrame) = AbundanceTable(
                                    Matrix(df[:,2:end]),
                                    names(df[2:end]),
                                    vector(df[1])
                                    )

getindex(t::AbundanceTable) = getindex(AbundanceTable.t)
setindex(t::AbundanceTable) = setindex(AbundanceTable.t)
length(t::AbundanceTable) = length(AbundanceTable.t)

## For use with abundance tables generated from Humman2
