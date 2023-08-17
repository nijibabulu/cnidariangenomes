#! /bin/bash

pkg_dir=$(dirname $(dirname $0))
oma_file=$(realpath $pkg_dir/data/orthologs/OrthologousGroups.fixed.txt)
Rscript -e "
setwd(\"$pkg_dir\")
devtools::load_all(\".\")
infer_cached_algs(
  \"((Em,Aq),(Ta,((((((Nv,Ep),Am),Xs),(Aa,Re)),(Ch,Hv)),((Tc,Cr),(Py,(Bf,Lo))))));\",
  critical_nodes = list(B=\"Tc,Cr,Py,Bf,Lo\",
                        C=\"Nv,Ep,Am,Xs,Aa,Re,Ch,Hv\",
                        M=\"Em,Aq,Ta,Nv,Ep,Am,Xs,Aa,Re,Ch,Hv,Tc,Cr,Py,Bf,Lo\"),
                        ortholog_source = \"$oma_file\")
"

