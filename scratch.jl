using StatPlots
using DataFrames
using Clustering

datadir = "/Users/kev/computation/science/hutlab/iamc/data"
plotsdir = "/Users/kev/Desktop/plots"
isdir(plotsdir) || mkdir(plotsdir)
colors = Dict(
    "IA" => RGBA(color("pink")),
    "RA" => RGBA(color("red")),
    "PsA" => RGBA(color("purple")),
    "REA" => RGBA(color("orange")),
    "AS" => RGBA(color("yellow")),
    "Unc" => RGBA(color("darkgrey")),
    "REF" => RGBA(color("black"), 0.2)
    )
d = readtable("$datadir/metadata/diagnoses.csv")
diagnoses = Dict(d[x,1] => colors[d[x,2]] for x in 1:size(d,1))

function get_color(sample_id::String)
    startswith(sample_id, "REF") && return RGBA(color("black"), 0.2)

    m = match(r"(^[A-Z]+)([0-9]+)_panphlan", sample_id)
    if m == nothing
        return :black
    elseif m.captures[1] in keys(colors)
        return colors[m.captures[1]]
    elseif "$(m.captures[1])$(m.captures[2])" in keys(diagnoses)
        return diagnoses["$(m.captures[1])$(m.captures[2])"]
    else
        return RGBA(color("darkgrey"))
    end
end


function panphlan_plots(name::String)
    df = readtable("$datadir/panphlan/$(name)_pa_refs.tsv")
    df = filter_rows(df, 1; kind=:perrow, calc=:sum)

    (abun, dm, clust_h, clust_v, pco) = panphlan_calcs(df)

    c = get_color.(string.(dm.labels))

    plot(pco, color=c, primary=false, marker=(:circle, 7, stroke(0)),
        xlabel="PCo1 ($(round(pco.variance_explained[1], 4) * 100)%)",
        ylabel="PCo2 ($(round(pco.variance_explained[2], 4) * 100)%)",)
        # title="$name strains (PanPhlAn)")
    savefig("$plotsdir/iamc_$(name)_pcoa.svg")

    corder = c[clust_h.order]

    plot(
            plot(hclustplot(clust_h, false), seriestype=:path, color=:black,
                xlims=(0.5,length(clust_h.order)+0.5), framestyle=:none,
                bottom_margin=-3mm), legend=false,
            plotannotations(corder),
            heatmap(abun.t[clust_v.order,clust_h.order], color=:YlGnBu, legend=false, ticks=false, top_margin=-4mm, bottom_margin=-1mm), layout=grid(3,1,heights=[0.1,0.05,.85]))

    savefig("$plotsdir/iamc_$(name)_heatmap.svg")
end


df = readtable("$datadir/panphlan/ecoli16_pa_refs.tsv")

df = filter_rows(df, 1; kind=:perrow, calc=:sum)

sum(!startswith.(string.(names(df)), "REF"))
sum(startswith.(string.(names(df)), "REF"))

abun = AbundanceTable(df)

dm = getdm(df, Jaccard())

rowdm = getrowdm(df, Jaccard())
clust_h = hclust(dm.dm, :single)
clust_v = hclust(rowdm.dm, :single)

pco = pcoa(dm)
pco.eigenvalues[1] / sum(pco.eigenvalues)
pco.variance_explained

plot(pco, primary=false, marker=(:circle, 7, stroke(0)),
    xlabel="PCo1 ($(round(pco.variance_explained[1], 4) * 100)%)",
    ylabel="PCo2 ($(round(pco.variance_explained[2], 4) * 100)%)",)


outtable = DataFrame(rowdm.dm)
writetable("/Users/kev/Desktop/rowdm.csv", outtable)
