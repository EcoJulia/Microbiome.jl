# can define alternate method to change which types are blocked for setproperty! etc
_restricted_fields(::AbstractSample) = (:name, :metadata)

# Samples should at minimum have name and metadata
name(as::AbstractSample) = as.name
metadata(as::AbstractSample) = as.metadata
name(as::AbstractFeature) = as.name

# Allows Strings to be used to index (since `Thing("name", missing)` will be == `Thing("name", something)`)
Base.String(as::AbstractSample) = name(as)
Base.String(af::AbstractFeature) = name(af)

Base.:(==)(s1::T, s2::T) where {T <: Union{AbstractSample, AbstractFeature}} = name(s1) == name(s2)

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

function set!(as::AbstractSample, prop::Symbol, val)
    prop in _restricted_fields(as) && error("Cannot set! $prop for $(typeof(as)).")
    set!(as.metadata, prop, val)
    return as
end

function unset!(as::AbstractSample, prop::Symbol)
    prop in _restricted_fields(as) && error("Cannot unset! $prop for $(typeof(as)).")
    unset!(as.metadata, prop)
    return as
end

function insert!(as::AbstractSample, prop::Symbol, val)
    prop in _restricted_fields(as) && error("Cannot insert! $prop for $(typeof(as)).")
    insert!(as.metadata, prop, val)
    return as
end

function delete!(as::AbstractSample, prop::Symbol)
    prop in _restricted_fields(as) && error("Cannot delete! $prop for $(typeof(as)).")
    delete!(as.metadata, prop)
    return as
end

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

clade(t::Taxon) = t.clade
clade(::Missing) = missing
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

taxon(gf::GeneFunction) = gf.taxon
hastaxon(gf::GeneFunction) = !ismissing(taxon(gf))
