abstract type AbstractAbundanceTable{T <: Real, 
                                     F <: AbstractFeature, 
                                     S <: AbstractSample} <: EcoBase.AbstractAssemblage{T, F, S}
end

"""
    CommunityProfile{T, F, S} <: AbstractAbundanceTable{T, F, S}

An `AbstractAssemblage` from [EcoBase.jl](https://github.com/EcoJulia/EcoBase.jl)
that uses a `SparseMatrixCSC` under the hood.

`CommunityProfile`s are tables with `AbstractFeature`-indexed rows and
`AbstractSample`-indexed columns.
Note - we can use the `name` of samples and features to index.
"""
mutable struct CommunityProfile{T, F, S} <: AbstractAbundanceTable{T, F, S}
    aa::AbstractSparseMatrix{T}
    features::AbstractVector{F}
    samples::AbstractVector{S}
    fidx::Dictionary{String, Int}
    sidx::Dictionary{String, Int}

    function CommunityProfile(aa::AbstractSparseMatrix,
                              feats::AbstractVector{<:AbstractFeature},
                              samps::AbstractVector{<:AbstractSample})
        fidx = Dictionary(name.(feats), eachindex(feats))
        sidx = Dictionary(name.(samps), eachindex(samps))

        T = eltype(aa)
        F = eltype(feats)
        S = eltype(samps)
        return new{T, F, S}(aa, feats, samps, fidx, sidx)
    end
end

function CommunityProfile(tab::AbstractMatrix,
                          feats::AbstractVector{<:AbstractFeature},
                          samps::AbstractVector{<:AbstractSample})
    return CommunityProfile(sparse(tab), feats, samps)
end
## -- Convienience functions -- ##

function ==(p1::CommunityProfile, p2::CommunityProfile)
    return abundances(p1) == abundances(p2) && 
           samples(p1)    == samples(p2) &&
           features(p1)   == features(p2)
end

"""
    taxonomicprofile(mat, features, samples)
"""
function taxonomicprofile(mat, features::AbstractVector{<:AbstractString}, samples::AbstractVector{<:AbstractString})
    CommunityProfile(mat, Taxon.(features), MicrobiomeSample.(samples))
end

"""
    functionalprofile(mat, features, samples)
"""
function functionalprofile(mat, features::AbstractVector{<:AbstractString}, samples::AbstractVector{<:AbstractString})
    CommunityProfile(mat, GeneFunction.(features), MicrobiomeSample.(samples))
end

"""
    metabolicprofile(mat, features, samples)
"""
function metabolicprofile(mat, features::AbstractVector{<:AbstractString}, samples::AbstractVector{<:AbstractString})
    CommunityProfile(mat, Metabolite.(features), MicrobiomeSample.(samples))
end


"""
    features(at::AbstractAbundanceTable)

Returns features in `at`. To get featurenames instead, use [`featurenames`](@ref).
"""
features(at::AbstractAbundanceTable) = at.features

"""
    samples(at::AbstractAbundanceTable)

Returns samples in `at`. To get samplenames instead, use [`samplenames`](@ref).
"""
samples(at::AbstractAbundanceTable) = at.samples

"""
    samples(at::AbstractAbundanceTable, name::AbstractString)

Returns sample in `at` with name `name`.
"""
function samples(at::AbstractAbundanceTable, name::AbstractString)
    idx = findall(==(name), samplenames(at))
    length(idx) == 0 && throw(IndexError("No samples called $name"))
    length(idx) > 1 && throw(IndexError("More than one sample matches name $name"))
    return samples(at)[axes(at.aa, 2)][first(idx)]
end

profiletype(at::AbstractAbundanceTable) = eltype(features(at))
ranks(at::AbstractAbundanceTable) = taxrank.(features(at))

Base.size(at::AbstractAbundanceTable, dims...) = size(at.aa, dims...)

Base.copy(at::AbstractAbundanceTable) = CommunityProfile(copy(abundances(at)), copy(features(at)), deepcopy(samples(at)))

# -- Indexing -- #

function _toinds(d, inds::AbstractVector{Regex})
    return findall(a-> any(ind-> contains(a, ind), inds), keys(d))
end

function _toinds(d, inds::AbstractVector{<: Union{AbstractString}})
    return findall(a-> any(==(a), inds), keys(d))
end

function _toinds(d, inds::AbstractVector{<: Union{AbstractFeature, AbstractSample}})
    return findall(a-> any(i-> name(i) == a, inds), keys(d))
end

# fall back ↑
_toinds(arr, ind::Union{AbstractSample, AbstractFeature, AbstractString, Regex}) = _toinds(arr, [ind])

# if inds are integers, just return them
_toinds(_, ind::Int) = ind
_toinds(_, inds::AbstractVector{Int}) = inds

function Base.getindex(at::CommunityProfile, inds...)
    mat = at.aa[inds...]
    
    CommunityProfile(mat, features(at)[inds[1]], samples(at)[inds[2]])
end

Base.getindex(at::CommunityProfile, rowind::Int, colind::Int) = at.aa[rowind, colind]
    
function Base.getindex(at::CommunityProfile, rowind::Union{T, AbstractVector{<:T}} where T<:Union{AbstractString,Regex}, colind)
    rows = _toinds(featurenames(at), rowind)
    mat = at.aa[rows, colind]

    CommunityProfile(mat, features(aa)[rows], samples(aa)[colind])
end

function Base.getindex(at::CommunityProfile, rowind, colind::Union{T, AbstractVector{<:T}} where T<:Union{AbstractString,Regex})
    cols = _toinds(samplenames(at), colind)
    mat = at.aa[rowind, cols]

    CommunityProfile(mat, features(aa)[rowind], samples(aa)[cols])
end

function Base.getindex(at::CommunityProfile, rowind::Union{T, AbstractVector{<:T}} where T<:Union{AbstractString,Regex},
                                             colind::Union{S, AbstractVector{<:S}} where S<:Union{AbstractString,Regex})
    rows = _toinds(featurenames(at), rowind)
    at[rows, colind]
end

## -- EcoBase Translations -- ##
# see src/ecobase.jl for Microbiome function names
# thing => feature
# place => sample
# occurrences => abundances (or totals)

EcoBase.thingnames(at::AbstractAbundanceTable) = name.(features(at))
EcoBase.placenames(at::AbstractAbundanceTable) = name.(samples(at))
EcoBase.occurrences(at::AbstractAbundanceTable) = at.aa
EcoBase.nthings(at::AbstractAbundanceTable) = size(at, 1)
EcoBase.nplaces(at::AbstractAbundanceTable) = size(at, 2)
## todo
# EcoBase.thingoccurrences(at::AbstractAbundanceTable, things) = nothing
# EcoBase.placeoccurrences(at::AbstractAbundanceTable, places) = nothing

# for custom printing
EcoBase.thingkind(asm::AbstractAbundanceTable) = "feature"
EcoBase.placekind(asm::AbstractAbundanceTable) = "sample"
## not needed for now
# EcoBase.thingkindplural(asm::AbstractAbundanceTable) = "$(thingkind(asm))s"
# EcoBase.placekindplural(asm::AbstractAbundanceTable) = "$(placekind(asm))s"

"""
    featuretotals(at::AbstractAbundanceTable)

Returns sum of each row (feature) in `at`.
Note, return value is a nfeatures x 1 `Matrix`, not a `Vector`.
If you need 1D `Vector`, use `vec(featuretotals(at))`.
"""
featuretotals(at::AbstractAbundanceTable) = sum(abundances(at), dims=2)

"""
    sampletotals(at::AbstractAbundanceTable)

Returns sum of each row (feature) in `at`.
Note, return value is a 1 x nsamples `Matrix`, not a `Vector`.
If you need 1D `Vector`, use `vec(sampletotals(at))`.
"""
sampletotals(at::AbstractAbundanceTable) = sum(abundances(at), dims=1)

## -- Tables Interface -- ##

Tables.istable(::AbstractAbundanceTable) = true
Tables.columnaccess(::AbstractAbundanceTable) = true
Tables.rowaccess(::AbstractAbundanceTable) = true

Tables.getcolumn(at::AbstractAbundanceTable, i::Int) = i == 1 ? featurenames(at) : abundances(at[:, i-1])
Tables.getcolumn(at::AbstractAbundanceTable, i::AbstractString) = i == "features" ? features(at) : abundances(at[:, i])
Tables.getcolumn(at::AbstractAbundanceTable, i::Symbol) = Tables.getcolumn(at, string(i))

Tables.columnnames(at::AbstractAbundanceTable) = [:features, Symbol.(samplenames(at))...]

function Tables.schema(at::AbstractAbundanceTable)
    elt = eltype(abundances(at))
    coltypes = [eltype(features(at)), (elt for _ in 1:nsamples(at))...]
    return Tables.Schema(Tables.columnnames(at), coltypes)
end


Tables.columns(at::AbstractAbundanceTable) = (; (col => Tables.getcolumn(at, col) for col in Tables.columnnames(at))...)

function _makerow(row::AbstractAbundanceTable)
    size(row, 1) == 1 || error("Can't make row from table of size $(size(row))")
    NamedTuple{(:features, Symbol.(samplenames(row))...)}((first(features(row)), abundances(row)...))
end

Tables.rows(at::AbstractAbundanceTable) = (_makerow(at[i, :]) for i in 1:nfeatures(at))

## -- Methods for absolute and relative abundances -- ##

# """
#     filterabund(abun::AbstractComMatrix, n::Int=minimum(10, nfeatures(abun)))

# Filter an abundance table to the top `n` features accross all samples

# This function also adds a row for "other", which sums the abundances of the
# remaining features.
# """
# function filterabund(abun::AbstractAbundanceTable, n::Int=minimum(10, nfeatures(abun)))
#     # TODO: add prevalence filter

#     totals = featuretotals(abun)

#     srt = sortperm(totals, rev=true)

#     newabun = getfeature(abun, srt[1:n])

#     remainder = [sum(occurrences(abun)[srt[n+1:end], i]) for i in 1:size(abun, 2)]'
#     newabun = vcat(newabun, remainder)
#     newrows = cat(featurenames(abun)[srt[1:n]], ["other"], dims=1)

#     return abundancetable(newabun, samplenames(abun), newrows)
# end

"""
    relativeabundance!(a::AbstractAbundanceTable; kind::Symbol=:fraction)

Normalize each sample in AbstractAbundanceTable to the sum of the sample.

By default, columns sum to 1.0.
Use `kind=:percent` for columns to sum to 100.
"""
function relativeabundance!(at::AbstractAbundanceTable; kind::Symbol=:fraction)
    in(kind, [:percent, :fraction]) || throw(ArgumentError("Invalid kind: $kind"))
    eltype(abundances(at)) <: AbstractFloat || throw(ArgumentError("relativeabundance! requires profile to have AbstractFloat eltype. Try relativeabundance instead"))
    abund = abundances(at)
    abund ./= sampletotals(at)
    kind == :percent && (abund .*= 100)
    dropzeros!(abund) # shouldn't need dropzeros - https://github.com/JuliaLang/julia/issues/39018
    return at
end

"""
    relativeabundance(at::AbstractAbundanceTable, kind::Symbol=:fraction)

Like [`relativeabundance!`](@ref), but does not mutate original.
"""
function relativeabundance(at::AbstractAbundanceTable, kind::Symbol=:fraction)
    comm = typeof(at)(float.(abundances(at)), deepcopy(features(at)), deepcopy(samples(at)))
    relativeabundance!(comm)
end

"""
    present(t::Union{Real, Missing}, minabundance::Real=0.0)
    present(at::AbstractAbundanceTable, minabundance::Real=0.0)

Check if a given (non-zero) value is greater than or equal to a minimum value.
If the minimum abundance is 0, just checks if value is non-zero.

If used on an `AbstractAbundanceTable`, returns a sparse boolean matrix of the same size.
"""
function present(t::Real, minabundance::Real=0.0)
    (minabundance >= 0 && t >= 0) || throw(DomainError("Only defined for positive values"))
    t == 0 ? false : t >= minabundance
end

present(::Missing, m::Real=0.0) = missing

function present(at::AbstractAbundanceTable, minabundance::Real=0.0)
    mat = spzeros(Bool, size(at)...)
    for i in eachindex(mat)
        mat[i] = present(at[Tuple(i)...], minabundance)
    end
    return mat
end


"""
    prevalence(a::AbstractArray{<:Real}, minabundance::Real=0.0)
    prevalence(at::AbstractAbundanceTable, minabundance::Real=0.0)

Return the fraction of values that are greater than or equal to a minimum.
If the minimum abundance is 0, returns the fraction of non-zero values.

If used on an `AbstractAbundanceTable`,
returns a prevalence value for each `feature` accross the `sample`s.
"""
prevalence(a::AbstractArray{<:Real}, minabundance::Real=0.0) = mean(x-> present(x, minabundance), a)

# makes it work for any iterable
prevalence(a, minabundance::Real=0.0) = mean(x-> present(x, minabundance), (y for y in a))

function prevalence(at::AbstractAbundanceTable, minabundance::Real=0.0)
    mean(x-> present(x, minabundance), abundances(at), dims=2)
end

"""
    prevalence_filter(comm::AbstractAbundanceTable; minabundance=0.0; minprevalence=0.05, renorm=false)

Return a filtered `CommunityProfile` where features with prevalence lower than `minprevalence` are removed.
By default, a feature is considered "present" if > 0, but this can be changed by setting `minabundance`.

Optionally, set `renorm = true` to calculate relative abundances after low prevalence features are removed.
"""
function prevalence_filter(comm::AbstractAbundanceTable; minabundance=0.0, minprevalence=0.05, renorm=false)
    comm = comm[vec(prevalence(comm, minabundance) .>= minprevalence), :]
    return renorm ? relativeabundance(comm) : comm
end

"""
    rankfilter(comm::AbstractAbundanceTable, cl::Union{Symbol, Int}; keepempty=false)

Return a copy of `comm`, where only rows that have `taxrank(feature) == cl` are kept.
Use `keepempty = true` to also keep features that don't have a `rank` (eg "UNIDENTIFIED").
"""
function rankfilter(comm::AbstractAbundanceTable, cl::Symbol; keepempty=false)
    in(cl, keys(_ranks)) ||  error("Invalid rank $cl, must be one of $(keys(_ranks))")
    if keepempty
        return filter(f-> !hasrank(f) || taxrank(f) == cl, comm)
    else
        return filter(f-> hasrank(f) && taxrank(f) == cl, comm)
    end
end

function rankfilter(comm::AbstractAbundanceTable, rank::Int; keepempty=false)
    0 <= rank <= 9 ||  error("Invalid rank $rank, must be one of $_ranks")
    return rankfilter(comm, keys(_ranks)[rank+1]; keepempty)
end


## Metadata

"""
    metadata(commp::CommunityProfile)

Returns iterator of `NamedTuple` per sample, where keys are `:sample`
and each metadata key found in `commp`.
Samples without given metadata are filled with `missing`.

Returned values can be passed to any Tables.rowtable - compliant type,
eg `DataFrame`.
"""
function metadata(commp::CommunityProfile)
    ss = samples(commp)
    cols = unique(reduce(vcat, collect.(keys.(metadata.(samples(commp))))))
    return Tables.rowtable(merge((; sample=name(s)), 
                     NamedTuple(c => get(s, c, missing) for c in cols)
                    ) for s in ss)
end


"""
    set!(commp::CommunityProfile, sample::AbstractString, prop::Symbol, val)
    set!(commp::CommunityProfile, sample::AbstractString, md::Union{AbstractDict, NamedTuple})

Update or insert a value `val` to the metadata of `sample` in the CommunityProfile `commp` using a Symbol `prop`. 
If you want an error to be thrown if the value already exists, use [`insert!`](@ref).

Can also pass a Dictionary or NamedTuple containing key=> value pairs,
all of which will be `set!`.
"""
function set!(commp::CommunityProfile, sample::AbstractString, prop::Symbol, val)
    sample = samples(commp, sample)
    prop in _restricted_fields(sample) && error("Cannot set! $prop for $(typeof(sample)).")
    set!(sample.metadata, prop, val)
    return sample
end

function set!(commp::CommunityProfile, sample::AbstractString, md::Dictionary)
    for (key, value) in pairs(md)
        set!(commp, sample, key, value)
    end
    return nothing
end

function set!(commp::CommunityProfile, sample::AbstractString, md::Union{AbstractDict, NamedTuple})
    md = Dictionary(md)
    set!(commp, sample, md)
    return nothing
end


"""
    set!(cp::CommunityProfile, md; namecol=:sample)

Add metadata (in the form of a `Tables.jl` table) a `CommunityProfile`.
One column (`namecol`) should contain sample names that exist in `commp`,
and other columns should contain metadata that will be added to the metadata of each sample.
"""
function set!(commp::CommunityProfile, md; namecol=:sample)
    Tables.istable(md) || throw(ArgumentError("Metadata must be a Tables.table"))
    sns = Set(samplenames(commp))
    md = filter(row-> row[namecol] in sns, md)
    for row in Tables.rows(md)
        sample = samples(commp, row[namecol])
        ks = filter(!=(namecol), keys(first(md)))
        for k in ks
            set!(sample, k, row[k])
        end
    end
    return nothing
end


"""
    insert!(commp::CommunityProfile, sample::AbstractString, prop::Symbol, val)

Insert a value `val` to the metadata of `sample` in a CommunityProfile `commp` using a Symbol `prop`, 
and it will throw an error if `prop` exists. 
If you don't want an error to be thrown if the value exists, use [`set!`](@ref).
"""
function insert!(commp::CommunityProfile, sample::AbstractString, prop::Symbol, val)
    sample = samples(commp, sample)
    prop in _restricted_fields(sample) && error("Cannot insert! $prop for $(typeof(sample)).")
    insert!(sample.metadata, prop, val)
    return sample
end

function insert!(commp::CommunityProfile, sample::AbstractString, md::Dictionary)
    isempty(Set(keys(commp, sample)) ∩ Set(keys(md))) || throw(IndexError("Duplicate keys found. Use `set!` to overwrite"))
    for (key, value) in pairs(md)
        insert!(commp, sample, key, value)
    end
    return nothing
end

function insert!(commp::CommunityProfile, sample::AbstractString, md::Union{<:AbstractDict, NamedTuple})
    md = Dictionary(md)
    insert!(commp, sample, md)
    return nothing
end


"""
    insert!(cp::CommunityProfile, md; namecol=:sample)

Add metadata (in the form of a `Tables.jl` table) a `CommunityProfile`.
One column (`namecol`) should contain sample names that exist in `commp`,
and other columns should contain metadata that will be added to the metadata of each sample.

Before starting, this will check that every value in every row is `insert!`able,
and will throw an error if not.
This requires iterating over the metadata table twice, which may be slow.
If performance matters, you can use `set!` instead, 
though this will overwrite existing data.
"""
function insert!(commp::CommunityProfile, md; namecol=:sample, careful=true)
    Tables.istable(md) || throw(ArgumentError("Metadata must be a Tables.table"))
    sns = Set(samplenames(commp))
    md = filter(row-> row[namecol] in sns, md)
    for row in Tables.rows(md)
        sample = samples(commp, row[namecol])
        ks = filter(!=(namecol), keys(first(md)))
        for k in ks
            haskey(sample, k) && throw(IndexError("Duplicate metadata detected. Use `set!` to force overwrite."))
        end
    end
    for row in Tables.rows(md)
        sample = samples(commp, row[namecol])
        ks = filter(!=(namecol), keys(first(md)))
        for k in ks
            set!(sample, k, row[k])
        end
    end
    return nothing
end


"""
unset!(commp::CommunityProfile, sample::AbstractString, prop::Symbol)

Delete a metadata entry in `sample` from CommunityProfile `commp` using the Symbol `prop`. 
If you want an error to be thrown if the value does not exist, use [`delete!`](@ref).
"""
function unset!(commp::CommunityProfile, sample::AbstractString, prop::Symbol)
    sample = samples(commp, sample)
    prop in _restricted_fields(sample) && error("Cannot unset! $prop for $(typeof(sample)).")
    unset!(sample.metadata, prop)
    return sample
end


"""
    delete!(commp::CommunityProfile, sample::AbstractString, prop::Symbol)

Delete a metadata entry in `sample` from CommunityProfile `commp` using the Symbol `prop` if it exists, or throw an error otherwise.
If you don't want an error to be thrown if the value does not exist, use [`unset!`](@ref).
"""
function delete!(commp::CommunityProfile, sample::AbstractString, prop::Symbol)
    sample = samples(commp, sample)
    prop in _restricted_fields(sample) && error("Cannot delete! $prop for $(typeof(sample)).")
    delete!(sample.metadata, prop)
    return sample
end

"""
    keys(commp::CommunityProfile, sample::AbstractString)

Return an iterator over all keys of the metadata attached to `sample` in a CommunityProfile `commp`. 
`collect(keys(commp, sample))` returns an array of keys. 
"""
Base.keys(commp::CommunityProfile, sample::AbstractString) = keys(metadata(samples(commp, sample)))

"""
    haskey(commp::CommunityProfile, sample::AbstractString, key::Symbol)

Determine whether the metadata of `sample` in a CommunityProfile `commp` has a mapping for a given `key`. 
Use `!haskey` to determine whether a `sample` in a CommunityProfile doesn't have a mapping for a given `key`
"""
Base.haskey(commp::CommunityProfile, sample::AbstractString, key::Symbol) = in(key, keys(samples(commp, sample)))

"""
    get(commp::CommunityProfile, sample::AbstractString, key::Symbol, default)

Return the value of the metadata in a `sample` stored for the given `key`, or the given `default` value if no mapping for the key is present.
"""
Base.get(commp::CommunityProfile, sample::AbstractString, key::Symbol, default=missing) = get(metadata(samples(commp, sample)), key, default)

"""
    get(commp::CommunityProfile, key::Symbol, default)

Return the value of the metadata in a `sample` stored for the given `key`, or the given `default` value if no mapping for the key is present.
"""
Base.get(commp::CommunityProfile, key::Symbol, default=missing) = [get(commp, sample, key, default) for sample in samplenames(commp)]


"""
    filter(f, comm::CommunityProfile)

Apply `f` to the features of `comm`,
and return a copy where `f(feature)` is `true`.
"""
function Base.filter(f::Function, commp::CommunityProfile)
    ridx = findall(f, features(commp))
    isempty(ridx) && error("Can't return empty profile")
    return copy(commp[ridx, :])
end