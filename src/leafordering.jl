"""
optimalorder(hc::Hclust, dm::Array{Float64,2})

Given a hierarchical cluster, use fast algorithm to determine optimal leaf order
minimizing the distance between adjacent leaves. This is done using a heuristic
where, when combining multi-leaf sub branches, only the outermost leaves are
compared (a maximum of 4 comparisons per intersection).

Sub branches are flipped if necessary to minimize the distance between adjacent
nodes, and then the combined branches are treated as a block for future
comparisons.

Based on:
Bar-Joseph et. al. "Fast optimal leaf ordering for hierarchical clustering." _Bioinformatics_. (2001)
"""
function optimalorder(hc::Hclust, dm::Array{Float64,2})
    ord = copy(hc.order)
    ext = Tuple[]
    orderleaves!(ord, ext, hc, dm)
    return ord
end


function optimalorder!(hc::Hclust, dm::Array{Float64,2})
    ord = hc.order
    ext = Tuple[]
    orderleaves!(ord, ext, hc, dm)
end


function orderleaves!(order::Vector{Int}, extents::Vector{Tuple}, hcl::Hclust, dm::Array{Float64,2})
    for (vl, vr) in zip(hcl.merge[:,1], hcl.merge[:,2])
        (u, m, uidx, midx) = leaflocs(vl, order, extents)
        (w, k, widx, kidx) = leaflocs(vr, order, extents)
        if vl < 0 && vr < 0
            push!(extents, (m, k))
        elseif vl < 0
            flp = flip1(m, k, w, dm)
            flp == 2 && reverse!(order, widx, kidx)
            push!(extents, (m, widx))
        elseif vl > 0 && vr > 0
            flp = flip2(u, m, k, w, dm)
            (flp == 2 || flp == 4) && reverse!(order, uidx, midx)
            (flp == 3 || flp == 4) && reverse!(order, widx, uidx)
            push!(extents, (uidx, widx))
        else
            error("invalid 'merge' order in Hclust: ($vl, $vr) ")
        end
    end
end


function leaflocs(v::Int, order::Vector{Int}, extents::Vector{Tuple})
    if v < 0
        outer = 0
        inner = findfirst(abs(v) .== order)
        branchleft = 0
        branchright = 0
    elseif v > 0
        branchleft = extents[abs(v)][1]
        branchright = extents[abs(v)][2]
        outer = order[branchleft]
        inner = order[branchright]
    else
        error("leaf postition cannot be zero")
    end
    return outer, inner, branchleft, branchright
end



"""
For 1 multi-leaf branch and a leaf, determine if flipping branch is required
1 = do not flip
2 = flip right
"""
function flip1(m::Int, k::Int, w::Int, dm::Array{Float64,2})
    dm[m,k] <= dm[m,w] ? 1 : 2
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
