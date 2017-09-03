@recipe function plot(dm::DistanceMatrix)
    dist = dm.dm
    M = classical_mds(dist, 2)
    labels -> string.(dm.labels)
    seriestype := :scatter
    M[1,:], M[2,:]
end


function hclustplot(hc::Hclust, useheight::Bool)
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

function treepositions(hc::Hclust, useheight::Bool)
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

function plotannotations(ann::Array{T,1}, colormap::Dict{T, Symbol}) where T
    xs = Int[]
    for i in 1:length(ann)
        append!(xs, [0,0,1,1,0] .+ (i-1))
    end
    xs = reshape(xs, 5, length(ann))
    ys = hcat([[0,1,1,0,0] for _ in ann]...)

    fc = reshape([colormap[a] for a in ann], 1,length(ann))

    plot(xs, ys, seriestype=:path, fill=(0,1),fillcolor=fc,
        legend=:false, color=:black, ticks=false, framestyle=false,
        top_margin=-5mm, bottom_margin=-1mm)
end
