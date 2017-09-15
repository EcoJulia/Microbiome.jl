using StatPlots
using MultivariateStats
using DataFrames
using Clustering

datadir = "/Users/kev/computation/science/hutlab/jdrf/data"
plotsdir = "$datadir/plots"
isdir(plotsdir) || mkdir(plotsdir)

rename!(meta, :Stool_kit_id_number, :sample_id)
mapping = by(mapping, :indexing_id) do df
    DataFrame(sample_id = df[1,:sample_id])
end
# meta = meta[[in(x, mapping[:sample_id]) for x in meta[:sample_id]], :]
meta[:Diagnosis_at_baseline] = lowercase.(meta[:Diagnosis_at_baseline])
meta = join(meta, mapping, on=:sample_id, kind=:outer)

diag_map = Dict("inflammatory arthralgia" => "IA",
                "rheumatoid arthritis" => "RA",
                "psoriatic arthritis" => "PsA",
                "psa" => "PsA",
                "unknown" => "Unc",
                "not available" => "Unc",
                "palindromic inflammatory arthritis" => "IA",
                "drug-induced inflammatory arthritis" => "IA",
                "unclassified arthritis" => "Unc")

meta[:diagnosis] = [isna(x) ? "Unc" : diag_map[x] for x in meta[:Diagnosis_at_baseline]]
ids = Dict(sample => index for (index,sample) in enumerate(meta[:indexing_id]))


diags = ["IA","RA","PsA","REA","AS","Unc","REF"]
colors = Dict(
    diags[i] => x for (i,x) in enumerate(distinguishable_colors(length(diags), colorant"darkred"))
    )


function get_diag(s::Symbol)
    name = String(s)

    if ismatch(r"^(REA)|(RA)|(AS)", name)
        diag = match(r"^[A-Z]+", name).match
    elseif ismatch(r"^x(\d+)", name)
        name = match(r"x(\d+)_\w+", name).captures[1]
        diag = meta[ids[parse(Int, name)], :diagnosis]
    else
        diag = "REF"
    end
    return String(diag)
end

function plotit(filename)
    df = readtable("$datadir/profiles/$filename")
    dm = getdm(df, Jaccard())
    abun = AbundanceTable(df)

    rowdm = getrowdm(df, Jaccard())
    size(dm.dm, 1) > 1 || return false
    clust_h = hclust(dm.dm, :single)
    clust_v = hclust(rowdm.dm, :single)

    diags = get_diag.(abun.samples[clust_h.order])

    name = basename(splitext(filename)[1])
    plot(
        plot(hclustplot(clust_h, false), seriestype=:path, color=:black,
            xlims=(0.5,length(clust_h.order)+0.5), framestyle=:none,title=name,
            legend=false, bottom_margin=-3mm),
        plotannotations(diags, colors),
        heatmap(abun.t[clust_v.order,clust_h.order], color=:YlGnBu, legend=false, ticks=false, top_margin=-4mm, bottom_margin=-1mm), layout=grid(3,1,heights=[0.1,0.05,.85]))
end

D = [0 3.16228 3.16228 7.07107 7.07107;
      3.16228 0 4.47214 4.47214 6.32456;
      3.16228 4.47214 0 6.32456 4.47214;
      7.07107 4.47214 6.32456 0 4.47214;
      7.07107 6.32456 4.47214 4.47214 0]

A = -1/2 * D.^2

Δ = getdelta(A)

f = eigfact(Δ)
v = f.values
vec = f.vectors
p = sortperm(v)

dm = DistanceMatrix(D, [i for i in 1:5], Jaccard())
p = pcoa(dm)
p.eigenvalues
p.eigenvectors
p.variance_explained

plot(p)

plot(principalcoord(p, 1), principalcoord(p,2), seriestype=:scatter, xticks = false, yticks = false,
    xlabel = "PCo1 ($(round(p.variance_explained[1], 1) * 100)%)",
    ylabel = "PCo2 ($(round(p.variance_explained[2], 1) * 100)%)")


classical_mds(D, 4)

1.342/0.223607
tax = readtable("/Users/kev/computation/science/hutlab/iamc/data/uk_analysis/taxonomic_profiles.tsv")


dm = getdm(tax, BrayCurtis())
p = pcoa(dm, correct_neg=true)
classical_mds(dm,2)
p[:,1:2]

plotly()

p.eigenvalues
plot(p)
