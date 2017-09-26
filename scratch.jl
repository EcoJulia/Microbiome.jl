using StatPlots
using DataFrames
using Clustering

datadir = "/Users/kev/computation/science/hutlab/HMPs/data"
plotsdir = "$datadir/plots"
distdir = "$datadir/distmatrices"
isdir(plotsdir) || mkdir(plotsdir)
isdir(distdir) || mkdir(distdir)

meta = readtable("$datadir/allmeta.csv")

meta = meta[.!isna.(meta[:Project]), [:Project, :External_ID, :site_sub_coll, :data_type, :diagnosis, :site_name, :Age_at_diagnosis]]

meta = meta[meta[:data_type] .== "metagenomics", :]
meta[:External_ID] = replace.(meta[:External_ID], "_P", "")
diags = Dict(meta[i, :External_ID] => meta[i, :diagnosis] for i in 1:size(meta, 1))

function dcolor(sample)
    sample = String(sample)
    startswith(sample, "REF") && return colorant"white"

    s = replace(sample, "_panphlan_map", "")

    !in(s, keys(diags)) && error("Diagnosis dict missing $s")

    diags[s] == "CD" && return colorant"darkred"
    diags[s] == "UC" && return colorant"orange"
    diags[s] == "nonIBD" && return colorant"cornflowerblue"

    error("diagnosis not CD, UC or nonIBD")
end

function plots(profile::String)
    m = match(r"\/?([a-z0-9]+)_pa\.tsv", profile)
    name = m.captures[1]
    df = readtable(profile)
    (abun, dm, clust_h, clust_v, pco) = panphlan_calcs(df)

    c = dcolor.(dm.labels)

    length(pco.variance_explained) < 2 && return

    plot(pco, color=c, primary=false, marker=(:circle, 5, stroke(1)),
        xlabel="PCo1 ($(round(pco.variance_explained[1], 4) * 100)%)",
        ylabel="PCo2 ($(round(pco.variance_explained[2], 4) * 100)%)",
        title="$name strains (PanPhlAn)")
    savefig("$plotsdir/pcoa_$name.png")

    corder = c[clust_h.order]

    plot(
        plot(hclustplot(clust_h, false), seriestype=:path, color=:black,
            xlims=(0.5,length(clust_h.order)+0.5), framestyle=:none,
            bottom_margin=-3mm, title="$name genes (PanPhlAn)"), legend=false,
        plotannotations(corder),
        heatmap(abun.t[clust_v.order,clust_h.order], color=:YlGnBu, legend=false, ticks=false, top_margin=-4mm, bottom_margin=-1mm), layout=grid(3,1,heights=[0.1,0.05,.85]))

    savefig("$plotsdir/presabs_$name.png")
end

for f in readdir("$datadir/profiles")
    m = match(r"\/?([a-z0-9]+)_pa\.tsv", f)
    name = m.captures[1]

    m = getdm(readtable("$datadir/profiles/$f"), Jaccard())

    size(m.dm, 1) < 3 && continue
    writedlm("$distdir/dist_$name.tsv", m.dm,"\t")
end


#--------------------JDRF

using StatPlots
using DataFrames
using Clustering
using Distances

datadir = "/Users/kev/computation/science/hutlab/jdrf/data"
plotsdir = "$datadir/plots"
distdir = "$datadir/distmatrices"
isdir(plotsdir) || mkdir(plotsdir)
isdir(distdir) || mkdir(distdir)

clibrary(:colorbrewer)
function plots(profile::String)
    df = readtable(profile)
    df = filter_rows(df, 2, kind=:perrow, calc=:sum)

    (abun, dm, clust_h, clust_v, pco) = panphlan_calcs(df)

    m = match(r"\/?([a-z0-9]+)_pa\.tsv", profile)
    name = m.captures[1]

    c = [startswith(string(x), "REF") ? :white : :blue for x in dm.labels]

    length(pco.variance_explained) < 2 && return

    plot(pco, color=c, primary=false, marker=(:circle, 5, stroke(1)),
        xlabel="PCo1 ($(round(pco.variance_explained[1], 4) * 100)%)",
        ylabel="PCo2 ($(round(pco.variance_explained[2], 4) * 100)%)",
        title="$name strains (PanPhlAn)")
    savefig("$plotsdir/pcoa_$name.svg")

    corder = c[clust_h.order]

    plot(
        plot(hclustplot(clust_h, false), seriestype=:path, color=:black,
            xlims=(0.5,length(clust_h.order)+0.5), framestyle=:none,
            bottom_margin=-3mm, title="$name genes (PanPhlAn)"), legend=false,
            size=(1200,1200),
        plotannotations(corder),
        heatmap(abun.t[clust_v.order,clust_h.order], color=:OrRd, legend=false, ticks=false, top_margin=-4mm, bottom_margin=-1mm), layout=grid(3,1,heights=[0.1,0.05,.85]))

    savefig("$plotsdir/presabs_$name.svg")
end


for f in readdir("$datadir/tmp")
    m = match(r"\/?([a-z0-9]+)_pa\.tsv", f)
    name = m.captures[1]
    info("Doing $name")
    plots("$datadir/tmp/$f")
    # size(m.dm, 1) < 3 && continue
    # writedlm("$distdir/dist_$name.tsv", m.dm,"\t")
end


for f in readdir("$datadir/tmp")
    m = match(r"\/?([a-z0-9]+)_pa\.tsv", f)
    name = m.captures[1]

    info("Doing $name")
    df = readtable("$datadir/tmp/$f")
    df = filter_rows(df, 2, kind=:perrow, calc=:sum)
    (abun, dm, clust_h, clust_v, pco) = panphlan_calcs(df)

    df = DataFrame(abun.t[clust_v.order, clust_h.order])
    df[:rownames] = [randstring() for _ in 1:size(df, 1)]
    df = df[:, [:rownames, names(df[1:end-1])...]]
    writetable("/Users/kev/Desktop/jdrf/clustered_$name.csv", df)

    # size(m.dm, 1) < 3 && continue
    # writedlm("$distdir/dist_$name.tsv", m.dm,"\t")
end

tax = readtable("$datadir/metaphlan_table.tsv")
# get rows broken down to species (k,p,c,o,f,g,s)
sp = tax[length.(split.(tax[1], "|")) .== 7, :]

bact = sp[startswith.(sp[:ID], "k__Bacteria"), :]
abt = AbundanceTable(bact)

totals = [sum(abt.t[i, :]) for i in 1:size(abt.t, 1)]
srt = sortperm(totals, rev=true)

plotabund(abt, sorton=:hclust)

function plotabund(df::DataFrame, n::Int=10, sorton::Symbol=:top)
    in(sorton, [:top, :hclust, df[:ID]...]) || error("invalid sorton option")
    2 < n < 12 || error("n must be between 2 and 12")


    topabund = filterabund(df, n)

    d = Dict{Symbol,Any}(x => 100-sum(topabund[x]) for x in names(topabund[2:end]))
    d[:ID] = "other"
    df2=DataFrame()
    for n in names(topabund)
        df2[n] = [d[n]]
    end

    goodcolors = [color("a6cee3",color("1f78b4",color("b2df8a",color("33a02c",color("fb9a99",color("e31a1c",color("fdbf6f",color("ff7f00",color("cab2d6",color("6a3d9a",color("ffff99",color("b15928"]

    append!(topabund, df2)
    topabund[:ID] = replace.(topabund[:ID], r"^[\w+\|]+?s__", "")
    foo = Matrix(topabund[2:end])'

    if sorton == :top
        srt = sortperm([topabund[11,i] for i in 2:size(topabund, 2)], rev=true)
    elseif sorton == :hclust
        DM = getdm(topabund, BrayCurtis())
        hc = hclust(DM, :single)
        srt = hc.order
    end
    groupedbar(foo[srt,:], bar_position=:stack, color=c, label=Vector(topabund[:ID]), legend=false)
end

savefig("/Users/kev/Desktop/jdrf/abund_hclust.svg")

groupedbar(foo[srt,:], bar_position=:stack, color=c, label=Vector(topabund[:ID]))
savefig("/Users/kev/Desktop/jdrf/abund_legend.svg")
