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
(in which case the `taxon` will be stored as `missing`).

See also [`taxon`](@ref Microbiome.taxon).
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

Get the `rank` field from a [`Taxon`](@ref) `t`.
Returns `missing` if the rank is not set.
"""
taxrank(t::Taxon) = t.rank
taxrank(::Missing) = missing

"""
    hasrank(t::Taxon)::Bool

Boolean function that returns `true` if the `rank`
field in [`Taxon`](@ref) `t` is not `missing`,
or `false` if it is `missing`
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
GeneFunction(n::AbstractString, t::AbstractString) = GeneFunction(n, taxon(t))

"""
    taxon(gf::GeneFunction)

Get the `taxon` field from a [`GeneFunction`](@ref), `gf`.
Returns `missing` if the taxon is not set.
"""
taxon(gf::GeneFunction) = gf.taxon

"""
    hastaxon(gf::GeneFunction)::Bool

Boolean function that returns `true` if the `taxon`
field in a [`GeneFunction`](@ref) `gf` is not `missing`,
or `false` if it is `missing`
"""
hastaxon(gf::GeneFunction) = !ismissing(taxon(gf))

"""
    taxrank(gf::GeneFunction)

Get the `rank` field from the `taxon` field of a [`GeneFunction`](@ref) `gf`
if it has one.
Returns `missing` if the `taxon` or `rank` is not set.
"""
taxrank(gf::GeneFunction) = taxrank(taxon(gf))

"""
    hasrank(gf::GeneFunction)::Bool

Boolean function that returns:

- `true` if `gf` has a [`Taxon`](@ref) with a non-missing `rank` field,
- `false` if there's no `Taxon`, or 
- `false` if the `Taxon` has no `rank`
"""
hasrank(gf::GeneFunction) = hastaxon(gf) && !ismissing(taxrank(gf))

Base.String(gf::GeneFunction) = hastaxon(gf) ? string(name(gf), '|', String(taxon(gf))) : name(gf)

"""
    genefunction(n::AbstractString)

Make a [`GeneFunction`](@ref) from a string,
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

Base.string(f::AbstractFeature) = String(f)

"""
    Metabolite(name::String, commonname::Union{Missing, String}, mz::Union{Missing, Float64}, rt::Union{Missing, Float64}) <: AbstractFeature
    Metabolite(name::String)

Represents a small-molecule metabolite coming from an LCMS.
The fields are

- `name`: required, this should be a unique identifier
- `commonname`: might refer to a chemical name like "proprionate"
- `mz`: The mass/charge ratio
- `rt`: The retention time
"""
struct Metabolite <: AbstractFeature
    name::String
    commonname::Union{Missing, String}
    mz::Union{Missing, Float64}
    rt::Union{Missing, Float64}
end

Metabolite(n::AbstractString) = Metabolite(n, missing, missing, missing)

name(m::Metabolite) = m.name

"""
    commonname(m::Metabolite)

Accessor function for the `commonname` field of a [`Metabolite`](@ref).
"""
commonname(m::Metabolite) = m.commonname


"""
    masscharge(m::Metabolite)

Accessor function for the `mz` field of a [`Metabolite`](@ref).
"""
masscharge(m::Metabolite) = m.mz

"""
    retentiontime(m::Metabolite)

Accessor function for the `rt` field of a [`Metabolite`](@ref).
"""
retentiontime(m::Metabolite) = m.rt

@testset "Metabolites" begin
    m1 = Metabolite("name", "common", 1., 1.)
    @test name(m1) == "name"
    @test commonname(m1) == "common"
    @test masscharge(m1) == 1
    @test retentiontime(m1) == 1
    m2 = Metabolite("name")
    @test name(m2) == "name"
    @test ismissing(commonname(m2))
    @test ismissing(masscharge(m2))
    @test ismissing(masscharge(m2))

end