#=
Based on:
Bar-Joseph et. al. "Fast optimal leaf ordering for hierarchical clustering." _Bioinformatics_. (2001)
=#

leaves = [1,2,3,4,5]

dm = [0. .1 .4 .5 .6;
      .1 0. .45 .3 .55;
      .4 .45 0. .2 .7;
      .5 .3 .2 0. .8;
      .6 .55 .7 .8 0.]

current_max = -Inf64
hcl = hclust(dm, :average)


orders = Tuple[]
connects = Tuple[]

merge = [(hcl.merge[i,1], hcl.merge[i,2]) for i in 1:size(hcl.merge, 1)]
for (vl, vr) in merge
    if vl < 0 && vr < 0
        push!(orders, (abs(vl), abs(vr)))
    elseif vl < 0
        push!(orders, (abs(vl),))
    elseif vr < 0
        push!(orders, (abs(vr),))
    else
        push!(connects, (vl, vr))
    end
end

"""
For 2 multi-leaf branches, determine if one or two flips is required
1 = do not flip
2 = flip left
3 = flip right
4 = flip both
"""
function flip2(order, dm::Array{Float64,2}, idx1::Int64, idx2::Int64)
    u = order[idx1][1]
    m = order[idx1][end]
    k = order[idx2][1]
    w = order[idx2][end]

    ds = [dm[m,k], dm[u,k], dm[m,w], dm[u,w]]
    return sortperm(ds)[1]
end

flip2(orders, dm, 1, 2)

connects
orders
