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
[Bar-Joseph et. al. "Fast optimal leaf ordering for hierarchical clustering." _Bioinformatics_. (2001)](https://doi.org/10.1093/bioinformatics/17.suppl_1.S22)
"""
function optimalorder(hc::Hclust, dm::Array{Float64,2})
    ord = deepcopy(hc)
    optimalorder!(ord, dm)
    return ord
end


function optimalorder!(hc::Hclust, dm::Array{Float64,2})
    ord = hc.order
    orderleaves!(ord, hc, dm)
end


function orderleaves!(order::Vector{Int}, hcl::Hclust, dm::Array{Float64,2})
    extents = Tuple{Int,Int}[]
    for (vl, vr) in zip(hcl.merges[:,1], hcl.merges[:,2])
        (u, m, uidx, midx) = leaflocs(vl, order, extents)
        (k, w, kidx, widx) = leaflocs(vr, order, extents)
        if vl < 0 && vr < 0
            # Nothing needs to be done
        elseif vl < 0
            flp = flip1(m, k, w, dm)
            flp == 2 && reverse!(order, kidx, widx)
        elseif vr < 0
            flp = flip1(k, m, u, dm)
            flp == 2 && reverse!(order, uidx, midx)
        elseif vl > 0 && vr > 0
            flp = flip2(u, m, k, w, dm)
            (flp == 2 || flp == 4) && reverse!(order, uidx, midx)
            (flp == 3 || flp == 4) && reverse!(order, kidx, widx)
        else
            error("invalid 'merge' order in Hclust: ($vl, $vr) ")
        end
        push!(extents, (uidx, widx))
    end
end


function leaflocs(v::Int, order::Vector{Int}, extents::Vector{Tuple{Int,Int}})
    ##################
    if v < 0
        leftextent = findfirst(abs(v) .== order)
        rightextent = leftextent
    elseif v > 0
        leftextent = extents[v][1]
        rightextent = extents[v][2]
    else
        error("leaf position cannot be zero")
    end
    left = order[leftextent]
    right = order[rightextent]
    return left, right, leftextent, rightextent
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
    argmin([dm[m,k], dm[u,k], dm[m,w], dm[u,w]])
end
