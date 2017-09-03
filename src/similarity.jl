struct DistanceMatrix{T<:Real} <: AbstractArray{T,2}
    dm::AbstractArray{T,2}
    labels::AbstractVector{S} where S
    distance::PreMetric
end

DistanceMatrix(dm::AbstractArray, distance) = DistanceMatrix(dm, Vector(1:size(dm,1)), distance)

@forward_func DistanceMatrix.dm Base.getindex, Base.setindex, Base.length, Base.size

function getdm(t::AbstractArray, distance::PreMetric)
    return DistanceMatrix(
            pairwise(distance, t),
            Vector(1:size(t,2)),
            distance)
end

function getdm(df::DataFrame, distance::PreMetric)
    return DistanceMatrix(
            pairwise(distance, Matrix(df[2:end])),
            Vector(names(df[2:end])),
            distance)
end

function getrowdm(df::DataFrame, distance::PreMetric)
    m = Matrix(df[2:end])'
    return DistanceMatrix(
            pairwise(distance, m),
            Vector(df[:,1]),
            distance)
end
