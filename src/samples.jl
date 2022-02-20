# can define alternate method to change which types are blocked for setproperty! etc
_restricted_fields(::AbstractSample) = (:name, :metadata)

# Samples should at minimum have name and metadata

"""
    name(t::Union{AbstractSample, AbstractFeature})

Get the `name` field from an `AbstractSample` or `AbstractFeature`.
"""
name(as::AbstractSample) = as.name

"""
    metadata(t::AbstractSample)

Get the `metadata` field from an `AbstractSample`.
Note that this is not a copy, so modifications to the returned value
will update the parent `AbstractSample` as well.
"""
metadata(as::AbstractSample) = as.metadata

Base.:(==)(m1::AbstractSample, m2::AbstractSample) = name(m1) == name(m2)

Base.String(as::AbstractSample) = name(as)
Base.String(af::AbstractFeature) = name(af)

"""
    getindex(as::AbstractSample, prop::Symbol)

Return the `prop` value in the metadata dictionary of `as`.
This enables using bracket syntax for access, eg `as[prop]`.
"""
function Base.getindex(as::AbstractSample, prop::Symbol)
    prop in _restricted_fields(as) && error("Do not use getindex to access $prop of $(typeof(as)). Use accessor function or getfield instead.")
    getindex(as.metadata, prop)
end


function Base.setindex!(as::AbstractSample, val, prop::Symbol)
    prop in _restricted_fields(as) && error("Do not use setindex! to change $prop of $(typeof(as)).")
    setindex!(as.metadata, val, prop)
    return as
end

"""
    getproperty(as::AbstractSample, prop::Symbol)

Return the `prop` value in the metadata dictionary of `as`.
This enables using dot syntax for access, eg `as.prop`.
"""
function Base.getproperty(as::AbstractSample, prop::Symbol)
    prop in _restricted_fields(as) ? getfield(as, prop) : as.metadata[prop]
end

function Base.setproperty!(as::AbstractSample, prop::Symbol, val)
    prop in _restricted_fields(as) && error("Cannot change $prop of $(typeof(as)) using setproperty!")
    setindex!(as.metadata, val, prop)
    return as
end

"""
    set!(as::AbstractSample, prop::Symbol, val)

Update or insert a value `val` to the metadata of sample `as` using a Symbol `prop`. 
If you want an error to be thrown if the value already exists, use [`insert!`](@ref).
"""
function set!(as::AbstractSample, prop::Symbol, val)
    prop in _restricted_fields(as) && error("Cannot set! $prop for $(typeof(as)).")
    set!(as.metadata, prop, val)
    return as
end

"""
    unset!(as::AbstractSample, prop::Symbol)

Delete a metadata entry of sample `as` using the Symbol `prop`. 
If you want an error to be thrown if the value does not exist, use [`delete!`](@ref).
"""
function unset!(as::AbstractSample, prop::Symbol)
    prop in _restricted_fields(as) && error("Cannot unset! $prop for $(typeof(as)).")
    unset!(as.metadata, prop)
    return as
end

"""
    insert!(as::AbstractSample, prop::Symbol, val)

Insert a value `val` to the metadata of sample `as` using a Symbol `prop`, 
and it will throw an error if `prop` exists. 
If you don't want an error to be thrown if the value exists, use [`set!`](@ref).
"""
function insert!(as::AbstractSample, prop::Symbol, val)
    prop in _restricted_fields(as) && error("Cannot insert! $prop for $(typeof(as)).")
    insert!(as.metadata, prop, val)
    return as
end

"""
    delete!(as::AbstractSample, prop::Symbol)

Delete a metadata entry of sample `as` using the Symbol `prop` if it exists, or throw an error otherwise.
If you don't want an error to be thrown if the value does not exist, use [`unset!`](@ref).
"""
function delete!(as::AbstractSample, prop::Symbol)
    prop in _restricted_fields(as) && error("Cannot delete! $prop for $(typeof(as)).")
    delete!(as.metadata, prop)
    return as
end

"""
    keys(as::AbstractSample)

Return an iterator over all keys of the metadata attached to sample `as`. 
`collect(keys(as))` returns an array of keys. 
"""
Base.keys(as::AbstractSample) = keys(metadata(as))

"""
    haskey(as::AbstractSample, key::Symbol)

Determine whether the metadata of sample `as` has a mapping for a given `key`. 
Use `!haskey` to determine whether a sample `as` in a CommunityProfile doesn't have a mapping for a given `key`
"""
Base.haskey(as::AbstractSample, key::Symbol) = in(key, keys(as))

"""
    get(as::AbstractSample, key::Symbol, default)

Return the value of the metadata in the sample `as` stored for the given `key`, or the given `default` value if no mapping for the key is present.
"""
Base.get(as::AbstractSample, key::Symbol, default) = get(metadata(as), key, default)


function set!(as::AbstractSample, d::Union{NamedTuple, Dictionary{Symbol, <:Any}})
    for (key, value) in pairs(d)
        set!(as, key, value)
    end
    return as
end

function insert!(as::AbstractSample, d::Union{NamedTuple, Dictionary{Symbol, <:Any}})
    isempty(Set(keys(as)) ∩ Set(keys(d))) || throw(ArgumentError("Duplicate keys found. Use `set!` to overwrite"))
    for (key, value) in pairs(d)
        insert!(as, key, value)
    end
    return as
end

"""
    MicrobiomeSample(name::String, metadata::Dictionary{Symbol, T}) <: AbstractSample
    MicrobiomeSample(name::String; kwargs...)
    MicrobiomeSample(name::String)

Microbiome sample type that includes a name and a [`Dictionary`](https://github.com/andyferris/Dictionaries.jl)
of arbitrary metadata using `Symbol`s (other than `:name` or `:metadata`) as keys.

Metadata can be accessed using `getproperty` or `getindex` on the sample itself.

Samples can be instantiated with only a name, leaving the `metadata` `Dictionary` blank

Adding or changing metadata follows [the same rules](https://github.com/andyferris/Dictionaries.jl#accessing-dictionaries) as for the normal `Dictionary`.

"""
struct MicrobiomeSample <: AbstractSample
    name::String
    metadata::Dictionary{Symbol, T} where {T <: Any} # currently non-functional
end

MicrobiomeSample(n::AbstractString; kwargs...) = isempty(kwargs) ? MicrobiomeSample(n, Dictionary{Symbol, Any}()) : MicrobiomeSample(n, dictionary(kwargs))
MicrobiomeSample(n::AbstractString, d::Union{AbstractDict,NamedTuple}) = MicrobiomeSample(n; pairs(d)...)