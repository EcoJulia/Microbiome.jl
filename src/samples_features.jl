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

name(as::AbstractFeature) = as.name

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

Examples
≡≡≡≡≡≡≡≡≡≡

```jldoctest
julia> set!(ms, :thing1, "metadata1")

julia> ms.thing1
"metadata1"
```
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

Examples
≡≡≡≡≡≡≡≡≡≡

```jldoctest
julia> unset!(ms, :thing1)

julia> !haskey(ms, :thing1)
true
```
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


Examples
≡≡≡≡≡≡≡≡≡≡

```jldoctest
julia> insert!(ms, :thing, "metadata")

julia> ms.thing
"metadata"
```
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

Examples
≡≡≡≡≡≡≡≡≡≡

```jldoctest
julia> delete!(ms, :thing) 

julia> !haskey(ms, :thing)
true
 ```
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

Examples
≡≡≡≡≡≡≡≡≡≡

```jldoctest
julia> ms = MicrobiomeSample("sample1", Dictionary([:thing1, :thing2], ["metadata1", "metadata2"]))

julia> collect(keys(ms))
2-element Vector{Symbol}:
 :thing1
 :thing2
```
"""
Base.keys(as::AbstractSample) = keys(metadata(as))

"""
    haskey(as::AbstractSample, key::Symbol)

Determine whether the metadata of sample `as` has a mapping for a given `key`. 
Use `!haskey` to determine whether a sample `as` in a CommunityProfile doesn't have a mapping for a given `key`

Examples
≡≡≡≡≡≡≡≡≡≡

```jldoctest
julia> set!(ms, :thing1, "metadata1")

julia> haskey(ms, :thing1)
true

julia> delete!(ms, :thing1, "metadata1")

julia> !haskey(ms, :thing1)
true
```
"""
Base.haskey(as::AbstractSample, key::Symbol) = in(key, keys(as))

"""
    get(as::AbstractSample, key::Symbol, default)

Return the value of the metadata in the sample `as` stored for the given `key`, or the given `default` value if no mapping for the key is present.

Examples
≡≡≡≡≡≡≡≡≡≡

```jldoctest
julia> get(ms, :thing1, 42)
42 

julia> insert!(ms, :thing1, 3.0) 

julia> get(ms, :thing1, 42)
3.0
```
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

```jldoctest MicrobiomeSample
julia> ms = MicrobiomeSample("sample1", Dictionary([:gender, :age], ["female", 180]))
MicrobiomeSample("sample1", {:gender │ "female", :age │ 180})

julia> name(ms)
"sample1"

julia> ms.name
"sample1"

julia> ms.gender
"female"

julia> ms.age
180
```

Samples can be instantiated with only a name, leaving the `metadata` `Dictionary` blank

```jldoctest MicrobiomeSample
julia> ms2 = MicrobiomeSample("sample2")
MicrobiomeSample("sample2", {})
```

or using keyword arguments.

```jldoctest MicrobiomeSample
julia> ms3 = MicrobiomeSample("sample3"; age=20)
MicrobiomeSample("sample2", {:age | 20})
```

Adding or changing metadata follows [the same rules](https://github.com/andyferris/Dictionaries.jl#accessing-dictionaries) as for the normal `Dictionary`.

"""
struct MicrobiomeSample <: AbstractSample
    name::String
    metadata::Dictionary{Symbol, T} where {T <: Any} # currently non-functional
end

MicrobiomeSample(n::AbstractString; kwargs...) = isempty(kwargs) ? MicrobiomeSample(n, Dictionary{Symbol, Any}()) : MicrobiomeSample(n, dictionary(kwargs))

const _ranks = (
    domain     = 0,
    kingdom    = 1,
    phylum     = 2,
    class      = 3,
    order      = 4,
    family     = 5,
    genus      = 6,
    species    = 7,
    subspecies = 8,
    strain     = 9
)

const _shortranks = (
    d = :domain,
    k = :kingdom,
    p = :phylum,
    c = :class,
    o = :order,
    f = :family,
    g = :genus,
    s = :species,
    t = :subspecies,
    u = missing
)

"""
    Taxon(name::String, rank::Union{Missing, Symbol, Int}) <: AbstractFeature
    Taxon(name::String)

Microbial taxon with a name and a rank that can be one of 

0. `:domain`
1. `:kingom`
2. `:phylum`
3. `:class`
4. `:order`
5. `:faamily`
6. `:genus`
7. `:species`
8. `:subspecies`
9. `:strain`

or `missing`. Contructors can also use numbers 0-9, or pass a string alone
(in which case the `taxon` will be stored as `missing`)
"""
struct Taxon <: AbstractFeature
    name::String
    rank::Union{Missing, Symbol}
    
    Taxon(s::AbstractString, ::Missing) = new(s, missing)
    Taxon(s::AbstractString, rank::Symbol) = in(rank, keys(_ranks)) ? 
                                                        new(s, rank)  :
                                                        error("Invalid rank $rank, must be one of $(keys(_ranks))")
end

Taxon(n::AbstractString, rank::Int) = 0 <= rank <= 9 ?
                                            Taxon(n, keys(_ranks)[rank+1]) :
                                            error("Invalid rank $rank, must be one of $_ranks")
Taxon(n::AbstractString) = Taxon(n, missing)

function Base.String(t::Taxon)
    if hasrank(t)
        return taxrank(t) == :strain ? string("t__", name(t)) : string(first(string(taxrank(t))), "__", name(t))
    else
        return string("u__", name(t))
    end
end

function taxon(n::AbstractString)
    m = match(r"^([dkpcofgst])__(.+)", n)
    isnothing(m) && return Taxon(n)
    return Taxon(string(m.captures[2]), _shortranks[Symbol(m.captures[1])])
end

"""
    taxrank(t::Union{Taxon, missing})

Get the `rank` field from an `Taxon`.
Returns `missing` if the rank is not set.
"""
taxrank(t::Taxon) = t.rank
taxrank(::Missing) = missing

"""
    hasrank(t::Taxon)::Bool

Pretty self-explanatory.
"""
hasrank(t::Taxon) = !ismissing(taxrank(t))


"""
    GeneFunction(name::String, taxon::Union{Taxon, String, Missing}) <: AbstractFeature
    GeneFunction(name::String)

Microbial gene function object with optional stratification (`taxon`).
"""
struct GeneFunction <: AbstractFeature
    name::String
    taxon::Union{Missing, Taxon}
end

GeneFunction(n::AbstractString) = GeneFunction(n, missing)
GeneFunction(n::AbstractString, t::AbstractString) = GeneFunction(n, Taxon(t))

"""
    taxon(t::GeneFunction)

Get the `taxon` field from a `GeneFunction`.
Returns `missing` if the taxon is not set.
"""
taxon(gf::GeneFunction) = gf.taxon

"""
    hastaxon(t::GeneFunction)::Bool

Pretty self-explanatory.
"""
hastaxon(gf::GeneFunction) = !ismissing(taxon(gf))

"""
    taxrank(gf::GeneFunction)

Get the `rank` field from the Taxon, if `gf` has one.
Returns `missing` if the taxon or rank is not set.
"""
taxrank(gf::GeneFunction) = taxrank(taxon(gf))

"""
    hasrank(t::GeneFunction)::Bool

Pretty self-explanatory.
"""
hasrank(gf::GeneFunction) = hastaxon(gf) && !ismissing(taxrank(gf))

Base.String(gf::GeneFunction) = hastaxon(gf) ? string(name(gf), '|', String(taxon(gf))) : name(gf)

"""
    genefunction(n::AbstractString)

Make a gene function from a string,
Converting anything after an initial `|` as a [`Taxon`](@ref).
"""
function genefunction(n::AbstractString)
    if contains(n, '|')
        spl = split(n, '|')
        return GeneFunction(string(spl[1]), taxon(join(spl[2:end], '|')))
    else
        return GeneFunction(n)
    end
end