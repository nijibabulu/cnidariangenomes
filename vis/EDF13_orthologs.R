
# compute or retrieve cache of computed ALGs
infer_cached_algs(
  "((Em,Aq),(Ta,((((((Nv,Ep),Am),Xs),(Aa,Re)),(Ch,Hv)),((Tc,Cr),(Py,(Bf,Lo))))));",
  critical_nodes = list(B="Tc,Cr,Py,Bf,Lo",
                        C="Nv,Ep,Am,Xs,Aa,Re,Ch,Hv",
                        M="Em,Aq,Ta,Nv,Ep,Am,Xs,Aa,Re,Ch,Hv,Tc,Cr,Py,Bf,Lo"),
  ortholog_source = "data/orthologs/OrthologousGroups.fixed.txt")

p <- plot_orthologs_summary(ogs)
save_fig(p, "EDF13_orthologs", height=6, width=6)

