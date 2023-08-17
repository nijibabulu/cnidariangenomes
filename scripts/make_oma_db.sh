#! /bin/bash

# This script links the files in the data/genomes directory to the data/orthologs/db directory 
# for use with OMA.
OMA_DB="data/orthologs/db"
ln -s $(realpath data/genomes/aurelia/Aurelia.Genome_v1.2_Protein_Models_12-28-18.fasta ) $OMA_DB/aa.fa
ln -s $(realpath data/genomes/digitifera/adi_aug101220_pasa_prot.fa ) $OMA_DB/ad.fa
ln -s $(realpath data/genomes/millepora/amilgenomev1.1/amil_1.1.maker_006.proteins.fasta ) $OMA_DB/am.fa
ln -s $(realpath data/genomes/amphimedon/Amphimedon_queenslandica.Aqu1.pep.all.fa ) $OMA_DB/aq.fa
ln -s $(realpath data/genomes/Bfl.fastp ) $OMA_DB/bf.fa
ln -s $(realpath data/genomes/celegans/c_elegans.PRJNA13758.WS276.protein.fa ) $OMA_DB/ce.fa
ln -s $(realpath data/genomes/clytia/full_nr_align.fasta ) $OMA_DB/ch.fa
ln -s $(realpath data/genomes/horseshoe/Annotation/MHSC_MAKER_proteins.fasta ) $OMA_DB/cr.fa
ln -s $(realpath data/genomes/drosophila/GCF_000001215.4_Release_6_plus_ISO1_MT_protein.faa ) $OMA_DB/dm.fa
ln -s $(realpath data/genomes/ephydatia/Emu_v1_prots.fasta) $OMA_DB/em.fa
ln -s $(realpath data/genomes/exaiptasia/GCF_001417965.1_Aiptasia_genome_1.1_protein.faa ) $OMA_DB/ep.fa
ln -s $(realpath data/genomes/hydra/hydra2.0_genemodels.aa ) $OMA_DB/hm.fa
ln -s $(realpath data/genomes/hydrav/GCF_022113875.1_Hydra_105_v3_protein.faa) $OMA_DB/hv.fa
ln -s $(realpath data/genomes/human/Homo_sapiens.GRCh38.pep.all.flt.fa ) $OMA_DB/hs.fa
ln -s $(realpath data/genomes/spottedgar/Lepisosteus_oculatus.LepOcu1.pep.all.fa ) $OMA_DB/lo.fa
ln -s $(realpath data/genomes/nematostella/nve.gene_models.vie130208/nveGenes.vienna130208.protein.fasta) $OMA_DB/nv.fa
ln -s $(realpath data/genomes/Scal100/NY_Scal100_v1.protein.fasta) $OMA_DB/sc.fa
sed -e 's/\(>PY_[^_]\+\)_\S\+/\1/' ../../../data/genomes/yesso/PY-genome/PYgene-strucanno.gffread.cds.transeq.prot.filt.fa > $OMA_DB/py.fa
ln -s $(realpath data/genomes/streblspio/Sbenedicti_v2.proteins.fasta) $OMA_DB/sb.fa
ln -s $(realpath data/genomes/trichoplax/Triad1_best_proteins.Tadh_P.ids.fasta) $OMA_DB/ta.fa
ln -s $(realpath data/genomes/millipede/Trigoniulus_corallinus_hic.proteins.fa) $OMA_DB/tg.fa
ln -s $(realpath data/genomes/rhopilema/pep.fasta) $OMA_DB/re.fa
ln -s $(realpath data/genomes/xenia/xenSp1.proteins.fa) $OMA_DB/xs.fa
