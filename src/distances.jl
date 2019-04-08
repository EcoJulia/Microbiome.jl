struct DistanceMatrix{T<:Real} <: AbstractArray{T,2}
    dm::AbstractArray{T,2}
    samples::AbstractVector{S} where S
    distance::PreMetric
end

struct PCoA{T<:AbstractFloat} <: AbstractArray{T,2}
    eigenvectors::Array{T,2}
    eigenvalues::Array{T,1}
    variance_explained::Array{T,1}
end

DistanceMatrix(dm::AbstractArray{T,2}, distance) where {T<:Real} = DistanceMatrix(dm, Vector(1:size(dm,1)), distance)

@forward_func DistanceMatrix.dm Base.getindex, Base.setindex, Base.length, Base.size
@forward_func PCoA.eigenvectors Base.getindex, Base.setindex, Base.length, Base.size

function getdm(t::AbstractComMatrix, distance::PreMetric)
    dm = pairwise(distance, t.occurrences, dims=2)
    for i in eachindex(dm); if isnan(dm[i]); dm[i] = 1; end; end
    return DistanceMatrix(
            dm,
            samplenames(t),
            distance)
end

function getdm(t::AbstractArray, distance::PreMetric)
    dm = pairwise(distance, t, dims=2)
    for i in eachindex(dm); if isnan(dm[i]); dm[i] = 1; end; end
    return DistanceMatrix(
            dm,
            Vector(1:size(t,2)),
            distance)
end

function getdm(df::DataFrame, distance::PreMetric)
    dm = pairwise(distance, convert(Matrix, df[2:end]), dims=2)
    for i in eachindex(dm); if isnan(dm[i]); dm[i] = 1; end; end
    return DistanceMatrix(
            dm,
            Vector(names(df[2:end])),
            distance)
end

function getrowdm(abt::AbstractComMatrix, distance::PreMetric)
    m = abt.occurrences'
    return DistanceMatrix(
            pairwise(distance, m, dims=2),
            samplenames(abt),
            distance)
end

function getrowdm(arr::AbstractArray, distance::PreMetric)
    m = arr'
    return DistanceMatrix(
            pairwise(distance, m, dims=2),
            Vector(1:size(arr,1)),
            distance)
end

function getrowdm(df::DataFrame, distance::PreMetric)
    m = convert(Matrix, df[2:end])'
    return DistanceMatrix(
            pairwise(distance, m, dims=2),
            Vector(df[:,1]),
            distance)
end


function pcoa(D::DistanceMatrix; correct_neg::Bool=false)
    n = size(D,1)
    A = -1/2 * D.^2
    Δ1 = getdelta(A)

    f = sortedeig(Δ1)

    if correct_neg && f.values[end] < 0
        c = abs(f.values[end])
        for (i,h) in [(i, h) for i in 1:n for h in 1:n if i != h]
            A[h,i] = -1/2 * D[h,i]^2 - c
        end
        Δ1 = getdelta(A)
        f = sortedeig(Δ1)
    end

    vals = f.values[1:n-1]
    return PCoA(f.vectors[:,1:n-1], vals, [v/sum(vals) for v in vals])
end

function sortedeig(M::Array{Float64,2})
    f = eigen(M, scale=true, permute=true)
    v = real.(f.values)
    p = sortperm(v, rev = true)
    return LinearAlgebra.Eigen(v[p], real.(f.vectors[:,p]))
end


function getdelta(A::AbstractArray{T,2}) where T <: AbstractFloat
    n = size(A,1)
    return reshape(
        [A[h,i] - mean(A[h,:]) - mean(A[:,i]) + mean(A) for i in 1:n for h in 1:n],n,n)
end


@inline eigenvalue(p::PCoA, inds...) = p.eigenvalues[inds...]
@inline variance(p::PCoA, inds...) = p.variance_explained[inds...]
@inline principalcoord(p::PCoA, inds...) = p[:,inds...]


## Diversity Measures

function shannon(v::AbstractVector{T}) where T<:Real
    total = sum(v)
    relab = map(x-> x/total, v)
    return -sum([log(x^x) for x in relab])
end


function ginisimpson(v::AbstractVector{T}) where T<:Real
    total = sum(v)
    relab = map(x-> x/total, v)
    return 1 - sum([x^2 for x in relab])
end
