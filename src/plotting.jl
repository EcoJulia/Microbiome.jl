@recipe function f(pc::PCoA)
    xticks := false
    yticks := false
    xlabel --> "PCo1 ($(round(variance(pc, 1) * 100, 2))%)"
    ylabel --> "PCo2 ($(round(variance(pc, 2) * 100, 2))%)"
    seriestype := :scatter
    principalcoord(pc, 1), principalcoord(pc,2)
end

@userplot AbundancePlot
@recipe function f(plt::AbundancePlot; topabund=10, sorton=:top)
    abun = plt.args[1]
    typeof(abun) <: AbstractComMatrix || error("AbundancePlot not defined for $(typeof(abun))")

    topabund = min(topabund, nfeatures(abun))
    in(sorton, [:top, :hclust, Symbol.(samplenames(abun))...]) || error("invalid sorton option")
    2 <= topabund < 12 || error("n must be between 2 and 12")

    top = filterabund(abun, topabund)

    rows = featurenames(top)

    if sorton == :top
        srt = sortperm(getfeature(top, topabund + 1), rev=true)
    elseif sorton == :hclust
        DM = getdm(abun, BrayCurtis())
        hc = hclust(DM, :single)
        srt = hc.order
    else
        error("invalid sorton option") # TODO: add feature-specific sorting
    end

    bar_position := :stack
    label := rows
    StatPlots.GroupedBar((1:nsamples(top), occurrences(top)[:,srt]'))
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


@userplot HClustPlot
@recipe function f(plt::HClustPlot; useheight=true)
    typeof(useheight) <: Bool || error("'useheight' argument must be true or false")

    hc = plt.args[1]
    useheight ? yt = true : yt = false

    pos = treepositions(hc, useheight=useheight)
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
    xs = reshape(xs, 4, size(hc.merge, 1))
    ys = reshape(ys, 4, size(hc.merge, 1))

    xlims := (0.5, length(hc.order) + 0.5)
    legend := false
    color := :black
    yticks --> yt
    xticks --> (1:length(hc.labels), hc.labels[hc.order])
    (xs, ys)
end

# From Michael K. Borregaard (posted to slack 2018-05-18)
# Usage:
#
# x = randn(100)
# y = randn(100) + 2
# y[y.<0] = 0
#
# zeroyplot(x,y, yscale = :log10, color = :red, markersize = 8)
@userplot ZeroYPlot
@recipe function f(h::ZeroYPlot)
           length(h.args) != 2 && error("zeroyplot only defined for x and y input arguments")
           x, y = h.args
           val0 = y .== 0

           layout := @layout([a
                              b{0.05h}])

           markeralpha --> 0.4
           seriestype --> :scatter

           @series begin
               primary := false
               subplot := 2
               yscale := :identity
               ylim := (-1,1)
               yticks := ([0],[0])
               grid := false
               x[val0], y[val0]
           end

           subplot := 1
           x[.!val0], y[.!val0]
end
