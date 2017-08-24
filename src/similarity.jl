struct DistanceMatrix{T<:Real} <: AbstractArray{T,2}
    dm::AbstractArray{T,2}
    labels::AbstractVector{S} where S
    distance::PreMetric

    function DistanceMatrix(dm, labels, dtype)
        size(dm, 1) == size(dm, 2) || error("Should be a symetrical distance matrix")
        for i in 1:size(dm,1)
            dm[i,i] == 0 || error("Distance between a sample and itself should be 0")
        end
        length(labels) == size(dm, 1) || error("labels must have same dimension as matrix")
        new(dm, labels, dtype)
    end
end

DistanceMatrix(dm::AbstractArray, distance) = DistanceMatrix(dm, Vector(1:size(dm,1)), distance)

getindex(t::DistanceMatrix) = getindex(DistanceMatrix.dm)
setindex(t::DistanceMatrix) = setindex(DistanceMatrix.dm)
length(t::DistanceMatrix) = length(DistanceMatrix.dm)

function getdm(t::AbstractArray, distance::PreMetric)
    return pairwise(distance, t)
end
