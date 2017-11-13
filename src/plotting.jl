@recipe function f(pc::PCoA)
    xticks := false
    yticks := false
    xlabel --> "PCo1 ($(round(pc.variance_explained[1] * 100, 2))%)"
    ylabel --> "PCo2 ($(round(pc.variance_explained[2] * 100, 2))%)"
    seriestype := :scatter
    principalcoord(pc, 1), principalcoord(pc,2)
end

@recipe function f(abun::AbundanceTable; topabund::Int=10, sorton::Symbol=:top)
    in(sorton, [:top, :hclust, abun.samples...]) || error("invalid sorton option")
    2 < topabund < 12 || error("n must be between 2 and 12")

    top = filterabund(abun, topabund)

    c = [color("#a6cee3") color("#1f78b4") color("#b2df8a") color("#33a02c") color("#fb9a99") color("#e31a1c") color("#fdbf6f") color("#ff7f00") color("#cab2d6") color("#6a3d9a") color("#ffff99") color("#b15928")]

    rows = replace.(string.(top.features), r"^[\w+\|]+?s__", "")
    foo = top.table'

    if sorton == :top
        srt = sortperm([top[topabund+1,i] for i in 1:size(top,2)], rev=true)
    elseif sorton == :hclust
        DM = getdm(top, BrayCurtis())
        hc = hclust(DM, :single)
        srt = hc.order
    else
        error("invalid sorton option")
    end


    @series begin
        bar_position := :stack
        color := c
        label := top.features
        StatPlots.GroupedBar((1:size(foo,1), foo[srt,:]))
    end
end

function treepositions(hc::Hclust; useheight::Bool=false)
    order = StatsBase.indexmap(hc.order)
    positions = Dict{}()
    for (k,v) in order
        positions[-k] = (v, 0)
    end
    for i in 1:size(hc.merge,1)
        xpos = mean([positions[hc.merge[i,1]][1], positions[hc.merge[i,2]][1]])
        if hc.merge[i,1] < 0 && hc.merge[i,2] < 0
            useheight ? ypos = hc.height[i] : ypos = 1
        else
            useheight ? h = hc.height[i] : h = 1
            ypos = maximum([positions[hc.merge[i,1]][2], positions[hc.merge[i,2]][2]]) + h
        end

        positions[i] = (xpos, ypos)
    end
    return positions
end


function annotationbar(colors::Array{T,1}) where T
    xs = Int[]
    for i in 1:length(colors)
        append!(xs, [0,0,1,1,0] .+ (i-1))
    end
    xs = reshape(xs, 5, length(colors))
    ys = hcat([[0,1,1,0,0] for _ in colors]...)

    fc = reshape(colors, 1,length(colors))

    plot(xs, ys,
        seriestype=:path,
        fill=(0,1),
        fillcolor=fc,
        legend=false,
        color=:black,
        ticks=false,
        framestyle=false)
end

@recipe function f(hc::Hclust; useheight::Bool=false)
    useheight ? yticks = true : yticks = false

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
    @series begin
        xlims := (0.5, length(hc.order) + 0.5)
        yticks := yticks
        legend := false
        color := :black
        xticks := false
        plot((reshape(xs, 4, size(hc.merge, 1)), reshape(ys, 4, size(hc.merge, 1))))
    end
end
