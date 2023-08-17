# compute or retrieve cache of computed ALGs
infer_cached_algs(
  "((Em,Aq),(Ta,((((((Nv,Ep),Am),Xs),(Aa,Re)),(Ch,Hv)),((Tc,Cr),(Py,(Bf,Lo))))));",
  critical_nodes = list(B="Tc,Cr,Py,Bf,Lo",
                        C="Nv,Ep,Am,Xs,Aa,Re,Ch,Hv",
                        M="Em,Aq,Ta,Nv,Ep,Am,Xs,Aa,Re,Ch,Hv,Tc,Cr,Py,Bf,Lo"),
  ortholog_source = "data/orthologs/OrthologousGroups.fixed.txt")

p <- patchwork::wrap_plots(list(
  plot_cluster_internals(clusters$M, "M", c("a","d","g")),
  plot_cluster_internals(clusters$C, "C", c("b","e","h")),
  plot_cluster_internals(clusters$B, "B", c("c","f","i"))
  ) , nrow=1)
save_fig(p, "EDF14_cluster_internals", height=6, width=12)

