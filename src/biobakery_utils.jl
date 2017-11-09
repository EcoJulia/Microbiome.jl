function filter_rows(df::DataFrame, quant::Real; kind::Symbol=:percolumn, calc::Symbol=:gt)
    in(kind, [:percolumn, :perrow]) || error("kind must be :percolumn or :perrow")
    in(calc, [:gt, :lt, :sum, :mean, :max, :min]) || error("calc must be :gt, :lt, :sum, :mean, :max, or :min")

    if kind==:perrow && calc==:sum
        keep = [sum(Matrix(df[i,2:end])) > quant ? true : false for i in 1:size(df, 1)]
    else
        error("that hasn't been implemented yet")
    end

    return df[keep, :]
end

"""
PanPhlAn Utils
"""

function panphlan_calcs(df::DataFrame)
    abun = AbundanceTable(df)
    dm = getdm(df, Jaccard())
    rowdm = getrowdm(df, Jaccard())
    clust_h = hclust(dm.dm, :single)
    clust_v = hclust(rowdm.dm, :single)
    pco = pcoa(dm)

    return abun, dm, clust_h, clust_v, pco
end
