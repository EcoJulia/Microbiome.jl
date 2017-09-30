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

meta = readtable("$datadir/diabimmune_metadata.txt", separator='\t')
gid = Dict(meta[i, :gid_wgs]=>meta[i, :subjectID] for i in 1:size(meta,1) if !isna(meta[i, :gid_wgs]))
diag = Dict(meta[i, :gid_wgs]=>meta[i, :t1d_diagnosed] for i in 1:size(meta,1) if !isna(meta[i, :gid_wgs]))

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

df = readtable("$datadir/tmp/bfragilis16_pa.tsv")
(abun, dm, clust_h, clust_v, pco) = panphlan_calcs(df)

heatmap(abun.t[clust_v.order,clust_h.order],
    color=:OrRd, legend=false, xticks=false)

genes = sort(df[clust_v.order, :x][14000:17000])





tax = readtable("$datadir/metaphlan_table.tsv")
# get rows broken down to species (k,p,c,o,f,g,s)
sp = tax[length.(split.(tax[1], "|")) .== 7, :]

bact = sp[startswith.(sp[:ID], "k__Bacteria"), :]
bact[:ID] = replace.(bact[:ID], r"^[\w\|]+s__", "")
rgnavus = bact[bact[:ID] .== "Ruminococcus_gnavus", :]
rgnavus = Dict(string(i) => rgnavus[1, i] for i in names(rgnavus)[2:end])

Set([gid[s] for s in keys(rgnavus) if rgnavus[s] > 0.1])

abt = AbundanceTable(bact)
ids = string.(names(bact)[2:end])
dcolor = [diag[i] ? colorant"orange" : colorant"cornflowerblue" for i in ids]

totals = [sum(abt.t[i, :]) for i in 1:size(abt.t, 1)]
srt = sortperm(totals, rev=true)

plotabund(abt, sorton=:hclust)

savefig("/Users/kev/Desktop/jdrf/abund_hclust.svg")

(p, srt) = plotabund(abt, sorton=:top)
plot(
    p,
    plotannotations(dcolor[srt]), layout=grid(2,1,heights=[0.9, 0.1]))
savefig("/Users/kev/Desktop/jdrf/abund_annotation.svg")



"REF_G000297735"
"REF_G000598845"



# -------------- IAMC

using DataFrames
using StatPlots

datadir = "/Users/kev/computation/science/hutlab/iamc/data"
plotsdir = "$datadir/plots"
isdir(plotsdir) || mkdir(plotsdir)

meta = readtable("$datadir/metadata/all_metadata.csv")
meta = meta[[:sample_id, :Classified_Diagnosis]]
metadict = Dict(meta[i,1]=>meta[i,2] for i in 1:size(meta,1))



tax = readtable("$datadir/taxonomic_profiles.tsv")
# get rows broken down to species (k,p,c,o,f,g,s)
sp = tax[length.(split.(tax[1], "|")) .== 7, :]
bact = sp[startswith.(sp[:_SampleID], "k__Bacteria"), :]
bact[:_SampleID] = replace.(bact[:_SampleID], r"^[\w\|]+s__", "")

delete!(bact,:x219165_taxonomic_profile)
delete!(bact,:RA1001_TCCGGAGA_taxonomic_profile)
delete!(bact,:x209190_taxonomic_profile)


abun = AbundanceTable(bact)
dm = getdm(bact, Jaccard())

cs = []
for s in dm.labels
    m = match(r"(\w+)_taxonomic_profile", string(s))
    name = string(m.captures[1])
    if !in(name, keys(metadict))
        push!(cs, color("#1b9e77"))
        continue
    end
    @show metadict[name]
    if isna(metadict[name])
        push!(cs, color("#1b9e77"))
    elseif metadict[name] == "RA"
        push!(cs, color("#d95f02"))
    elseif metadict[name] == "AS"
        push!(cs, color("#7570b3"))
    elseif metadict[name] == "PsA"
        push!(cs, color("#e7298a"))
    elseif metadict[name] == "IA"
        push!(cs, color("#66a61e"))
    else
        push!(cs, color("#1b9e77"))
    end
end
pco = pcoa(dm)
# cs = reshape(cs, 1, length(cs))
plot(pco, color=cs, primary=false, marker=(:circle, 8),
    xlabel="PCo1 ($(round(pco.variance_explained[1], 4) * 100)%)",
    ylabel="PCo2 ($(round(pco.variance_explained[2], 4) * 100)%)", size=(1200,1200))
scatter!(Float64[ ], Float64[ ], label = "RA", marker=(:circle,color("#d95f02")))
scatter!(Float64[ ], Float64[ ], label = "AS", marker=(:circle,color("#7570b3")))
scatter!(Float64[ ], Float64[ ], label = "PsA", marker=(:circle,color("#e7298a")))
scatter!(Float64[ ], Float64[ ], label = "IA", marker=(:circle,color("#66a61e")))
scatter!(Float64[ ], Float64[ ], label = "Other/unknown", marker=(:circle,color("#1b9e77")))
savefig("/Users/kev/Desktop/aruk_pco.svg")

(p, srt) = plotabund(abun, sorton=:hclust)

plot(p, legend=true)
savefig("/Users/kev/Desktop/aruk_abun_legend.svg")
plot(p, legend=false)
savefig("/Users/kev/Desktop/aruk_abun.svg")

plotannotations(cs[srt])
savefig("/Users/kev/Desktop/aruk_abun_annotations.svg")


function getcolors(dm::DistanceMatrix)
    cs = []
    for s in dm.labels
        if startswith(string(s), "REF")
            push!(cs, color("#FFFFFF"))
            continue
        end

        m = match(r"^(\w+)_panphlan_map", string(s))
        name = string(m.captures[1])
        if !in(name, keys(metadict))
            push!(cs, color("#CCCCCC"))
            continue
        end
        @show metadict[name]
        if isna(metadict[name])
            push!(cs, color("#CCCCCC"))
        elseif metadict[name] == "RA"
            push!(cs, color("#ED1C24"))
        elseif metadict[name] == "AS"
            push!(cs, color("#D95F02"))
        elseif metadict[name] == "PsA"
            push!(cs, color("#FCEE21"))
        elseif metadict[name] == "IA"
            push!(cs, color("#FF7BAC"))
        else
            push!(cs, color("#CCCCCC"))
        end
    end
    return cs
end


profile="$datadir/panphlan/prevotella_pa_refs.tsv"
df = readtable(profile)
df = filter_rows(df, 2, kind=:perrow, calc=:sum)

(abun, dm, clust_h, clust_v, pco) = panphlan_calcs(df)

m = match(r"\/?([a-z0-9]+)_pa_refs\.tsv", profile)
name = m.captures[1]

c = getcolors(dm)

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
