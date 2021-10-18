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

"""
    taxon(::AbstractString)

Return a [`Taxon`](@ref) from a string representation.
If the string contains taxonomic rank information in the form
`"x__Thename"` where `x` is the first letter of the rank,
this information will be used.

## Examples

```julia-repl
julia> taxon("Unknown")
Taxon("Unknown", missing)

julia> taxon("s__Prevotella_copri")
Taxon("Prevotella_copri", :species)
```
"""
function taxon(n::AbstractString)
    m = match(r"^([dkpcofgstu])__(.+)", n)
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