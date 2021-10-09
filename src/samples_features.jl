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

## getters and setters for metadata ##
function Base.getindex(as::AbstractSample, prop::Symbol)
    prop in _restricted_fields(as) && error("Do not use getindex to access $prop of $(typeof(as)). Use accessor function or getfield instead.")
    getindex(as.metadata, prop)
end

function Base.setindex!(as::AbstractSample, val, prop::Symbol)
    prop in _restricted_fields(as) && error("Do not use setindex! to change $prop of $(typeof(as)).")
    setindex!(as.metadata, val, prop)
    return as
end

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

"""
    MicrobiomeSample(name::String, metadata::Dictionary{Symbol, T}) <: AbstractSample
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

Adding or changing metadata follows the same rules as for the normal `Dictionary` type.

- to change a value use `setproperty` or `setindex`.
  Note that this will fail if the key does not already exist.
- to add or remove a value with validation that it does or does not exist already,
  use `insert!` and `delete!` respectively.
- to add or remove a value without validation ("upsert"),
  use `set!` and `unset!` respectively.

```jldoctest MicrobiomeSample
julia> ms
MicrobiomeSample("sample1", {:gender │ "female", :age │ 180})

julia> ms.age = 16 * 365 # or `ms[:age] = 16 * 365`
5840

julia> set!(ms, :gender, "nonbinary")
MicrobiomeSample("sample1", {:gender │ "nonbinary", :age │ 5840})

julia> insert!(ms, :occupation, "clerk")
MicrobiomeSample("sample1", {:gender │ "nonbinary", :age │ 5840, :occupation │ "clerk"})

julia> insert!(ms, :occupation, "bagger")
ERROR: IndexError("Dictionary already contains index: occupation")

julia> delete!(ms, :occupation)
MicrobiomeSample("sample1", {:gender │ "nonbinary", :age │ 5840})

julia> delete!(ms, :occupation)
ERROR: IndexError("Index doesn't exist: occupation")

julia> unset!(ms, :occupation)
MicrobiomeSample("sample1", {:gender │ "nonbinary", :age │ 5840})
```
"""
struct MicrobiomeSample <: AbstractSample
    name::String
    metadata::Dictionary{Symbol, T} where {T <: Any} # currently non-functional
end

MicrobiomeSample(n::AbstractString) = MicrobiomeSample(n, Dictionary{Symbol, Any}())
# Base.convert(::Type{MicrobiomeSample}, s::AbstractString) = MicrobiomeSample(s)

const _clades = (
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

"""
    Taxon(name::String, clade::Union{Missing, Symbol, Int}) <: AbstractFeature
    Taxon(name::String)

Microbial taxon with a name and a clade that can be one of 

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
    clade::Union{Missing, Symbol}
    
    Taxon(s::AbstractString, ::Missing) = new(s, missing)
    Taxon(s::AbstractString, clade::Symbol) = in(clade, keys(_clades)) ? 
                                                        new(s, clade)  :
                                                        error("Invalid clade $clade, must be one of $(keys(_clades))")
end

Taxon(n::AbstractString, clade::Int) = 0 <= clade <= 9 ?
                                            Taxon(n, keys(_clades)[clade+1]) :
                                            error("Invalid clade $clade, must be one of $_clades")
Taxon(n::AbstractString) = Taxon(n, missing)

Base.String(t::Taxon) = hasclade(t) ? string(first(string(clade(t))), "__", name(t)) : name(t)


"""
    clade(t::Union{Taxon, missing})

Get the `clade` field from an `Taxon`.
Returns `missing` if the clade is not set.
"""
clade(t::Taxon) = t.clade
clade(::Missing) = missing

"""
    hasclade(t::Taxon)::Bool

Pretty self-explanatory.
"""
hasclade(t::Taxon) = !ismissing(clade(t))


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
    taxon(t::Union{GeneFunction, missing})

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
    clade(gf::GeneFunction)

Get the `clade` field from the Taxon, if `gf` has one.
Returns `missing` if the taxon or clade is not set.
"""
clade(gf::GeneFunction) = clade(taxon(gf))

"""
    hasclade(t::GeneFunction)::Bool

Pretty self-explanatory.
"""
hasclade(gf::GeneFunction) = hastaxon(gf) && !ismissing(clade(gf))
