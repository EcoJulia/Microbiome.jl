@recipe function plot(dm::DistanceMatrix)
    plot(pcoa(dm.dm))
end

@recipe function plot(pc::PCoA)
    xticks := false
    yticks := false
    xlabel -> "PCo1 ($(round(pc.variance_explained[1], 1) * 100)%)"
    ylabel -> "PCo2 ($(round(pc.variance_explained[2], 1) * 100)%)"
    seriestype := :scatter
    principalcoord(pc, 1), principalcoord(pc,2)
end

# function plot_test(abund::AbundanceTable, dm::DistanceMatrix, annotations::Array{T,1}) where T <: AbstractArray
#     clust = hclust(dm.dm, :single)
#
#
#
#     plot(
#         plot(hclustplot(clust, false), seriestype=:path, color=:black, xlims=(0.5,length(c)+1), framestyle=:none, legend=false, margin=-5px),
#         plotannotations(diags, colors),
#         heatmap(abund.t[:,clust.order], color=:PuBu, legend=false, ticks=false, margin=-5px), layout=grid(3,1,heights=[0.1,0.05,.85]))
# end




function hclustplot(hc::Hclust, useheight::Bool=false)
    o = indexmap(hc.order)
    n = [x for x in 1:length(o)]

    pos = treepositions(hc, useheight)

    xs = []
    ys = []
    for i in 1: size(hc.merge, 1)
        x1 = pos[hc.merge[i,1]][1]
        x2 = pos[hc.merge[i,2]][1]
        append!(xs, [x1,x1,x2,x2])

        y1 = pos[hc.merge[i,1]][2]
        y2 = pos[hc.merge[i,2]][2]
        useheight ? h = hc.height[i] : h = 1
        newy = maximum([y1,y2]) + h
        append!(ys, [y1,newy,newy,y2])
    end
    return (reshape(xs, 4, size(hc.merge, 1)), reshape(ys, 4, size(hc.merge, 1)))
end

function treepositions(hc::Hclust, useheight::Bool=false)
    order = indexmap(hc.order)
    positions = Dict{}()
    for (k,v) in order
        positions[-k] = (v, 0)
    end
    for i in 1:size(hc.merge,1)
        xpos = mean([positions[hc.merge[i,1]][1], positions[hc.merge[i,2]][1]])
        if hc.merge[i,1] < 0 && hc.merge[i,2] < 0
            useheight ? ypos = hc.height[i] : ypos = 1
        else
            currenttop = maximum([positions[hc.merge[i,1]][2], positions[hc.merge[i,2]][2]])
            useheight ? h = hc.height[i] - currenttop : h = 1
            ypos = maximum([positions[hc.merge[i,1]][2], positions[hc.merge[i,2]][2]]) + h
        end

        positions[i] = (xpos, ypos)
    end
    return positions
end

@recipe function annotations(annotations::Array{T,1}, colormap::Dict{T, Symbol}) where T
    xs = Int[]
    for i in 1:length(annotations)
        append!(xs, [0,0,1,1,0] .+ (i-1))
    end
    xs = reshape(xs, 5, length(samples))
    ys = hcat([[0,1,1,0,0] for _ in samples]...)

    fillcolor := reshape([colormap[a] for a in annotations], 1, length(annotations))
    seriestype := path
    legend := false

    color -> :black

    xs, ys
end

function plotannotations(colors::Array{T,1}) where T
    xs = Int[]
    for i in 1:length(colors)
        append!(xs, [0,0,1,1,0] .+ (i-1))
    end
    xs = reshape(xs, 5, length(colors))
    ys = hcat([[0,1,1,0,0] for _ in colors]...)

    fc = reshape(colors, 1,length(colors))

    plot(xs, ys, seriestype=:path, fill=(0,1),fillcolor=fc,
        legend=:false, color=:black, ticks=false, framestyle=false,
        top_margin=-5mm, bottom_margin=-1mm)
end


function plotabund(abun::AbundanceTable, n::Int=10; sorton::Symbol=:top)
    in(sorton, [:top, :hclust, abun.samples...]) || error("invalid sorton option")
    2 < n < 12 || error("n must be between 2 and 12")


    topabund = filterabund(abun, n)

    c = [color("#a6cee3") color("#1f78b4") color("#b2df8a") color("#33a02c") color("#fb9a99") color("#e31a1c") color("#fdbf6f") color("#ff7f00") color("#cab2d6") color("#6a3d9a") color("#ffff99") color("#b15928")]

    rows = replace.(string(topabund.rows), r"^[\w+\|]+?s__", "")
    foo = topabund.t'

    if sorton == :top
        srt = sortperm([topabund[n+1,i] for i in 1:size(topabund, 2)], rev=true)
    elseif sorton == :hclust
        DM = getdm(topabund, BrayCurtis())
        hc = hclust(DM, :single)
        srt = hc.order
    end
    groupedbar(foo[srt,:], bar_position=:stack, color=c, label=Vector(topabund.samples), legend=false)
end
