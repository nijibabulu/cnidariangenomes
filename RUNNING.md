# Running the code

<!-- TOC -->

- [Running the code](#running-the-code)
    - [TL;DR](#tldr)
    - [Details](#details)
        - [Ancestral Linkgage Groups](#ancestral-linkgage-groups)
        - [Genomes](#genomes)
        - [Orthologs](#orthologs)
        - [Molecular Clock](#molecular-clock)
        - [Genome-genome alignment](#genome-genome-alignment)
        - [Size estimates](#size-estimates)
        - [QC](#qc)
        - [Contig and Scaffold Assemblies](#contig-and-scaffold-assemblies)
        - [Superscaffold/Chromosome Assembly with Hi-C](#superscaffoldchromosome-assembly-with-hi-c)
        - [Visualizing the Assembly](#visualizing-the-assembly)

<!-- /TOC -->

The figures from the paper can be readily reproduced with the code below. We follow a top-down approach here by showing how to make the figures first and proceed to generating ancestral linkage groups, orthologs and genome assemblies last.

## TL;DR

The workflow for generating the figures from the data is based on [devtools](https://devtools.r-lib.org/), working with a local clone of this repository. It will also be necessary to download the data from the supplemental data repository (not yet established). The following steps should be sufficient to generate the figures:

1. Clone this reponsitory to your local machine.
2. Install the `devtools` package using `install.packages("devtools")`.
3. Install the dependencies using `devtools::install_deps()`.
4. Load the package using `devtools::load_all()`.
5. Run the scripts in the `vis/` directory. It should be possible to start with a fresh R session for each of the scripts. Run them line-by-line with the current working directory as the top-level directory of this repository.
6. The figures will be saved in a new `figures/` directory.

## Details

Below we describe the main steps that were used as a rough guide to the pipeline used in this paper. It won't be possible to run the code exactly as it is written here, since the file names and paths are not included. However, the code should be sufficient to guide the user to the appropriate functions and scripts.

### Ancestral Linkgage Groups

For many scripts, the `infer_cached_algs()` function is called, which either returns a cached version of the ALGs if they exist in the `cache/` directory, or runs a very compute- and memory-intensive algorithm to infer them. For this purpose, we provide the script `scripts/infer_algs.sh` which can be readily run from a command line alone on an HPC. This script will generate the ALGs for all of the assemblies in the `data/` directory, and place them in the `cache/` directory. The script can be run with no arguments. The main driving function is in `R/alg_inference.R`. It runs on a single core and takes under 64GB at its peak.

### Genomes

The download locations of the non-Edwardsiid genomes are listed in the paper, but for convenience we provide the genome annotation and amino acid sequence files in the `data/genomes` directory. For inferring ALGs, the function `load_cached_gene_positions()`. The function can serve as a guide to the directory.

### Orthologs

Orthologs were generated using the OMA algorithm. The parameters are provided in the `data/orthologs/parameters.drw` file. The genome database can be made with the `script/make_oma_db.sh` command.

### Molecular Clock

We took an automated approach by finding orthologs to BUSCOs, as described in the paper, aligning the found genes across species, concatenating the alignments and finally inferring the molecular clock using `r8s`. For each of the proteomes, we use the [BUSCO](https://busco.ezlab.org/) algorithm to find orthologs to the metazoan BUSCOs (version 9). The commands used are:

```bash
for f in data/genomes/orthologs/db/*.fa; do
    run_BUSCO.py -l metazoa_odb9 -m prot -i $f -o tmp/seq/$(basename $f).buscos
done
```

The resulting BUSCOs are then concatenated and aligned using [mafft](https://mafft.cbrc.jp/alignment/software/), followed by [trimal](http://trimal.cgenomics.org/). The tree is generated using [IQ-TREE](http://www.iqtree.org/) The commands used are, e.g.:

```bash
for f in tmp/seq/*.buscos; do mafft --maxiterate 1000 --genafpair $f > tmp/aln/$(basename $f) done
for f in tmp/aln/*; do trimal -gappyout -in $f -out tmp/trm/$(basename $f); done
iqtree -bb 1000 -alrt 1000 -msub nuclear -p tmp/trm
```

At this point, the tree should be edited to be rooted as needed, then, referring to the `r8s` [manual](https://naturalis.github.io/mebioda/doc/week1/w1d5/r8s1.7.manual.pdf), fill in
the appropriate MCRA constraints. In our case, we took 2 constraints and used them to determine the
bounds of the distances:

```bash
r8s  -b -f high.nex > high.out
r8s  -b -f low.nex > low.out
```

From there, extract the newick tree from the file to use with the `EDF1_clock.R` script.

### Genome-genome alignment

To compare the dovetail scaffold-based assembly and the contig-based assembly, we used the `nucmer` and `delta-filter` programs from the [MUMmer](http://mummer.sourceforge.net/) package. The commands used are:

```bash
nucmer --prefix=tmp/nv_dovetail-nv_rb data/reviewed_assemblies/nv_bwa_lowmasked_gapless_chroms/nv_bwa_lowmasked_gapless_chroms.final.fasta data/reviewed_assemblies/nv_ph_bwa_lowmasked_gapless_chroms/nv_ph_bwa_lowmasked_gapless_chroms.final.fasta
delta-filter -r -q tmp/nv_dovetail-nv_rm.delta > nv_dovetail-nv_rm.delta.filter
```

We used a customized version of the `mummerplot` script to generate the dotplots found in `scripts/mummerplot`.

```bash
perl scripts/mummerplot -t png --medium -R reviewed_assemblies/nv_bwa_lowmasked_gapless_chroms/nv_bwa_lowmasked_gapless_chroms.final.fasta -Q reviewed_assemblies/nv_ph_bwa_lowmasked_gapless_chroms/nv_ph_bwa_lowmasked_gapless_chroms.final.fasta -p assembly/17_comparisons/nv/nv_dovetail-nv_rb assembly/17_comparisons/nv/nv_dovetail-nv_rb.delta.filter
```

### Size estimates

We made k-mer curves based on size estimates using the [GenomeScope](http://qb.cshl.edu/genomescope/) algorithm. These were performed using a separate sample of _Nematostella vectnesis_ from the one used for the assemblies. The first step is to generate the k-mer counts using [jellyfish](https://github.com/gmarcais/Jellyfish), which can done with something like the following:

```bash
cat nv*.fastq | jellyfish count  --quality-start=33 --min-quality=20 -m 56 -s 1G -o size/01_jellyfish/nvm$m.jf /dev/stdin
cat sc*.fastq | jellyfish count  --quality-start=33 --min-quality=20 -m 18 -s 1G -o size/01_jellyfish/nvm$m.jf /dev/stdin
```

Note we tried several different k-mer sizes, and the ones that worked best were 56 for _Nematostella vectensis_ and 18 for _Scolanthus callimorphus_, because we had initially underestimated the size of the genome. The output of this is can be used with the code contained in this package (which contains modified code from GenomeScope), in `vis/EDF2a_size_estimates.R`.

### QC

Assembly size, length distribution and contiguity curves are calculated directly in functions called by the `vis/EDF2cdefgh_assembly_summary.R` script. They rely on the fasta index files generated by `samtools faidx`. If the index file is present, it is not necessary to have the fasta file.

[REAPR](https://mybiosoftware.com/reapr-1-0-16-genome-assembly-evaluation.html) breaks the genome where there is conflict between the assembly and the reads, and the relative degree to which the genome is broken can serve as a proxy for the quality of the assembly. The REAPR requires a reference genome and short reads. We reuse the above reads for the assembly size assessment. It may be useful to concatinate all the reads into single files or to later concatinate the bam files produced by REAPR. In general this pipeline takes a lot of time. Example commands for a single assembly are:

```bash
R1=nv_r1.fastq
R2=nv_r2.fastq
ASSY=nv_bwa_lowmasked_gapless_chroms.final.fasta
reapr facheck $ASSY # check the assembly for formatting issues that may cause problems
reapr perfectmap $ASSY $R1 $R2 400 tmp/reapr/${ASSY}perfectmap.bam
reapr smaltmap $ASSY $R1 $R2 tmp/reapr/${ASSY}smaltmap.bam
reapr pipeline $ASSY tmp/reapr/${ASSY}perfectmap.bam tmp/reapr/${ASSY}out tmp/reapr/${ASSY}perfectmap
```

See the [manual](https://ftp.sanger.ac.uk/resources/software/reapr/Reapr_1.0.18.manual.pdf) for more details. The final `05.summary.report.tsv` is used in the `vis/EDF2cdefgh_assembly_summary.R` script.

[BUSCO](https://busco.ezlab.org/) is run in genomic mode in order to assess the relative completeness of the genomes. The `-sp` parameter is used by augustus when searching for genes in the assembly. In our case, we trained augustus on our _Nematostella vectensis_ genes, however, one can also use the built-in model. Note we used version 3.3 and the current version is 5+. Results will vary here, but the important point is relative completeness.

```bash
ASSY=nv_bwa_lowmasked_gapless_chroms.final.fasta
run_BUSCO.py -l metazoa_odb9 -m prot -i $f -o tmp/busco/$(basename $ASSY) -sp nematostella_vectensis
```

### Contig and Scaffold Assemblies

[CANU](https://github.com/marbl/canu) was used for contig assembly. Anecdotally, the newer [flye](https://github.com/fenderglass/Flye) assembler may produce better results for a similar data set. We used the following parameters for both genomes. `$READS` is a set of reads in uncompressed fastq format

```bash
canu useGrid=false genomeSize=400m maxMemory=1000 rawErrorRate=0.3 correctedErrorRate=0.045 -p nv -d tmp/assembly/nv -pacbio-raw $READS"
```

The _Nematostella_ and _Scolanthus_ genomes followed separate paths. The _Nematostella_ genome received a separate seuqencing of the same original sample of DNA by Dovetail genomics. They performed HiRise sequencing in a proprietary pipeline. The data are available via the SRA data set.

The _Nematostella_ dovetail scaffolds and _Scolanthus_ contigs had separate issues: the _Scolanthus_ contigs derived from under-covered sequence due to the unexpected large size of the genome, leading to a large number of redundant contigs due to misassemblies. We did not achieve good separation between what appears to be redundant contigs due to coverage and single-copy contigs. Therefore, for _Scolanthus_, unlike _Nematostella_ we used the [Redundans](https://github.com/Gabaldonlab/redundans) pipeline:

```bash
redundans.py -f $SC_CONTIGS -o tmp/redundans_out -t 8 --noscaffolding --nogapclosing --overlap 0.66
```

The output file from `contigs.reduced.fa` was then used as input to Hi-C assembly.

We used [Purge Haplotigs](https://bitbucket.org/mroachawri/purge_haplotigs/src/master/) in order to remove haplotigs from the _Nematostella_ assembly. This pipeline operates on an understanding of over-coverage of certain contigs over others which indicates that they are likely haplotigs. Therefore one must first map external reads or the source reads to the assembly. For the following we illustrate how to run with a single sample of reads, but in practice several reads need to be sorted and merged together. The commands are:

```bash
IND=tmp/ph/$(basename $NV_CONTIGS)
cp $NV_CONTIGS tmp/ph
bwa index -p $IND
bwa mem $IND reads1.fastq reads2.fastq | samtools view -bF $(samtools flags "UNMAP,MUNMAP,SECONDARY,SUPPLEMENTARY" | cut -f1) > $IND.bam
samtools sort $IND.bam > $IND.sorted.bam
purge_haplotigs readhist -b $IND.sorted.bam -g $IND
purge_haplotigs cov -i $IND.gencov -l 5 -m 20 -h 200 -o $IND.cov_stats.csv
purge_haplotigs purge -g $IND -c $IND.cov_stats.csv -b $IND.sorted.bam -p $IND.purged
```

The pipeline may also be run using PacBio reads.

### Superscaffold/Chromosome Assembly with Hi-C

Our Hi-C data were generated using the Phase genomics kit, and we followed the [guideline](https://phasegenomics.github.io/2019/09/19/hic-alignment-and-qc.html) they provided for aligning and filtering the data. The provided bamfiles created allwed us to do an intial assembly using [Lachesis](https://shendurelab.github.io/LACHESIS/). The commands used are:

```bash
LACHESIS_CONTIGS=tmp/$(basename $NV_CONTIGS)
cp $NV_CONTIGS tmp/lachesis
perl scripts/CountMotifsInFasta.pl $LACHESIS_CONTIGS
grep '^>' $LACHESIS_CONTIGS > tmp/lachesis/names.txt
```

From there modify the [example INI file](https://github.com/shendurelab/LACHESIS/blob/master/src/bin/INIs/test_case.ini) from the github repository to include the `DRAFT_ASSEMBLY_FASTA`, `SAM_FILES` (in `SAM_DIR`), `RE_SITE_SEQ`, and `REF_ASSEMBLY_FASTA` variables to run against this INI file (`Lachesis INI_FILE`). Lachesis does a good job of providing an initial assembly, but refinement is best performed with the [Juicebox](https://aidenlab.gitbook.io/juicebox/desktop) tool. The recommended workflow is to start with the 3d-DNA pipeline which produces a `.hic` file. We used the [jcvi library](https://pypi.org/project/jcvi/) to generate `agp` files:

```bash
python3 -mjcvi.assembly.hic agp -o $LACHESIS_CONTIGS.agp tmp/lachesis/main_results $LACHESIS_CONTIGS
```

The [Juicebox Scripts](https://github.com/phasegenomics/juicebox_scripts) provided by Phase genomics will create an assembly file needed for the 3d-DNA pipeline:

```bash
python agp2assembly.py $LACHESIS_CONTIGS.agp > $LACHESIS_CONTIGS.assembly
matlock bam2 juicer $LACHESIS_CONTIGS.bam $LACHESIS_CONTIGS.links.txt # the BAM should be generated in the phase genomics pipeline
sort -k2,2 -k6,6 $LACHESIS_CONTIGS.links.txt $LACHESIS_CONTIGS.links.srt.txt
```

The [3d-DNA pipeline](https://github.com/aidenlab/3d-dna) can produce a `.hic` file:

```bash
bash run-assembly-visualizer.sh -p false $LACHESIS_CONTIGS.assembly $LACHESIS_CONTIGS.links.srt.txt
```

This generates the `.hic` file which can be used with Juicebox. The tool will dump a new .assembly file that describes the new arrangement of contigs. The new assembly can be exported using Juicebox Scripts again:

```bash
juicebox_assembly_converter.py -g 0 -s -a $LACHESIS_CONTIGS.reviewed.assembly -f $LACHESIS_CONTIGS -p $LACHESIS_CONTIGS.reviewed
```

This process can be repeated until the assembly is satisfactory.

### Visualizing the Assembly

The assembly as visualized in `vis/F1de_contact_maps.R` is based on the exported matrices from the above. From a `.hic` file, one can dump the matrix from the command line using a script from the `3d-DNA` pipeline. For example if one wishes to dump the 100k matrix:

```bash
RESOLUTION=100
bash juicebox_tools.sh dump observed KR $LACHESIS_CONTIGS.reviewed.hic assembly assembly BP ${RESOLUTION}000 $LACHESIS_CONTIGS.reviewed.${RESOLUTION}k.mat
```

The matrix would be part of the folder which should also contain the `superscaffold_track`, `scaffold_track` emitted by Juicebox.

