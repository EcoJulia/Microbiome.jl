#=
Based on:
Bar-Joseph et. al. "Fast optimal leaf ordering for hierarchical clustering." _Bioinformatics_. (2001)
=#
using Distances


m = rand(5, 15)
dm = pairwise(BrayCurtis(), m)

hcl = hclust(dm, :average)
dump(hcl)

function optimalorder(hc::Hclust, dm::Array{Float64,2})
    ord = hc.order
    extents = Tuple[]

    for (vl, vr) in zip(hc.merge[:,1], hc.merge[:,2])
        @show vl, vr
        if vl < 0 && vr < 0
            m = findfirst(abs(vl) .== ord)
            k = findfirst(abs(vr) .== ord)
            push!(extents, (m, k))
        elseif vl < 0
            m = findfirst(abs(vl) .== ord)
            kidx = extents[vr][1]
            widx = extents[vr][2]

            k = ord[kidx]
            w = ord[widx]
            flp = flip1(m, k, w, dm)
            flp == 2 && reverse!(ord, kidx, widx)
            push!(extents, (m, widx))
        elseif vl > 0 && vr > 0
            uidx = extents[vl][1]
            midx = extents[vl][2]
            kidx = extents[vr][1]
            widx = extents[vr][2]

            u = ord[uidx]
            m = ord[midx]
            k = ord[kidx]
            w = ord[widx]

            flp = flip2(u, m, k, w, dm)
            @show flp
            (flp == 2 || flp == 4) && reverse!(ord, uidx, midx)
            (flp == 3 || flp == 4) && reverse!(ord, kidx, widx)
            push!(extents, (uidx, widx))

        else
            error("ooops")
        end
    end
    return ord
end

hcl.order = optimalorder(hcl, dm)

"""
For 1 multi-leaf branch and a leaf, determine if flipping branch is required
1 = do not flip
2 = flip right
"""
function flip1(m::Int, k::Int, w::Int, dm::Array{Float64,2})
    @inbounds dm[m,k] <= dm[m,w] ? 1 : 2
end

"""
For 2 multi-leaf branches, determine if one or two flips is required
1 = do not flip
2 = flip left
3 = flip right
4 = flip both
"""
function flip2(u::Int, m::Int, k::Int, w::Int, dm::Array{Float64,2})
    s = sortperm([dm[m,k], dm[u,k], dm[m,w], dm[u,w]])
    return s[1]
end
