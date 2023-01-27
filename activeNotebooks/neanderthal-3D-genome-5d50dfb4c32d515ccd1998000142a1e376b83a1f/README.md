# README

Created: October 28, 2021 4:54 PM
Created By: Evonne McArthur
Last Edited: February 5, 2022 3:06 PM

## About

- Code and data for McArthur et al. “Reconstructing the 3D genome organization of Neanderthals reveals that chromatin folding shaped phenotypic and sequence divergence”
- Manuscript preprint: [TBD]
- Most of the analysis and figure-generation was done in the `bin/Neanderthal_3D_Genome.ipynb` with python 3.6.10 (see dependencies in `environment.yml`). It is not recommended to run this notebook from start-to-finish, but to run individual sections for the analysis you want to replicate. All custom files used are in this repo but there many be some additional files used downloaded from other resources (e.g. eQTLs from GTEx). These publicly available resources are all linked and cited in the methods section of the paper (section “Data availability”). To run Akita on your own sequences, please see the Akita repo: [https://github.com/calico/basenji/tree/master/manuscripts/akita](https://github.com/calico/basenji/tree/master/manuscripts/akita).

## Manifest

### bin/

- `Neanderthal_3D_Genome.ipynb` - Ipython/Jupyter notebook in python for analysis and figure generating code. It is organized by sub-section and figure to align with the presentation of results from the paper. Other python scripts are also launched from this notebook so if you want example usage for the other python scripts (below), check this out. Sections include:
    
    <img src="https://github.com/emcarthur/neanderthal-3D-genome/blob/d594cbe1f01b3303b2bdb6eccfc4c72899cfa12e/toc.png" width="600" >
    
- `environment.yml`- List and version of packages/dependencies used in Anaconda environment for all python scripts in this project
- `compareCellTypes.py` - Runs Akita on positions across the genome to compare HFF and ESC results, used for a supplemental figure. See `Neanderthal_3D_Genome.ipynb` for how this was parallelized.
- `inSilicoMutagenesis.py` - Conducts in silico mutagenesis over a file of 1 Mb windows across the genome (replace archaic-specific alleles one by one and run akita and comparisons). See `Neanderthal_3D_Genome.ipynb` for how this was parallelized.
- `goEnrichment.py`- Generates empirical distribution of HPO/GWAS annotation enrichment for the null hypothesis so we can calculate p values and enrichment statistics for observed counts
- `goEnrichmentEmpiricalFDR.py` - Used to generate an expected number of false positives for ontology enrichment in order to calculate FDR-corrected P values (q)
- generate3Dgenome.py - in progress, runs Akita across genome
- comparisons3Dgenome.py - in progress, compares Akita outputs across genome
- comparisonsSequence.py - in progress, compares sequence divergence across the genome
- expected3DgenomeDivergence.py - in progress, generates empiric distribution of 3D differences based on shuffling sequence differences
- locateSeqDifferences.py - in progress, finds sequence differences between two individuals
- `akita_model/` & `basenji/` From Akita CNN model, from here: [https://github.com/calico/basenji/tree/master/manuscripts/akita](https://github.com/calico/basenji/tree/master/manuscripts/akita)

### data/

- `regionsWithFullCoverage.csv`: 1 Mb Regions (chr and start coordinate) with full coverage in 1000G modern humans (to calculate the end coordinate as a bedfile add 2^20 to the windowStartPos column `regionsWithFullCoverage.bed`). Further data analysis considers only these regions.
- `ideogram_hg19.txt`: cytoBandIdeo table from UCSC Table browser in hg19 (to color ideogram) (inspired by [https://www.biostars.org/p/9922/](https://www.biostars.org/p/9922/))
- `archaicMissingness/` (`vindija/` & `chagyrskaya/` & `altai/` & `denisova/` & `all/`) - Places  (bed coordinates) where archaic genomes have missingness. The `all/` folder is where missingness is observed in at least one archaic. Coordinates in where at least one archaic has missingness get filled in with reference to mask for appropriate comparisons
- `sequenceComparisons/` - pairwise quantification of amount of sequence differences between two individuals
    - `withArchaics/` - masked regions of archaic missingness to facilitate comparisons (fill with hg19 ref)
    - `withoutArchaics/` - not masked because we are only comparing between modern humans in regions with no missingness
- ctcf/
    - `HFF_allCTCF_hg19.bed` - CTCF-bound open chromatin candidate cis-regulatory elements (cCREs) in the HFF cell type ( [https://screen.encodeproject.org/](https://screen.encodeproject.org/) > Downloads > by cell type > HFF-Myc male newborn originated from foreskin fibroblast, lifted-over to hg19)
    - `hg19_archetype_motifs_CTCF.bed.gz` - CTCF motifs are from genome-wide motif scans v1.0 from Viestra et al: [https://resources.altius.org/~jvierstra/projects/motif-clustering/releases/v1.0/](https://resources.altius.org/~jvierstra/projects/motif-clustering/releases/v1.0/) (all models in the CTCF archetype motif cluster) lifted over to hg19
    - `AFR_ESN_female_HG03105_vindija_differences.bed.gz`- coordinates of sequence differences between HG03105 and Vindija (from `locateSeqDifferences.py`)
    - `4DN_HFFboundariesMicroC_hg19.bed` - TAD boundaries in the HFF cell type are from processed MicroC data available at the 4D nucleome data portal ( [https://data.4dnucleome.org/experiment-set-replicates/4DNES9X112GZ/](https://data.4dnucleome.org/experiment-set-replicates/4DNES9X112GZ/) , lifted-over to hg19)
- `examples/` - Data for example figures in Figs 3 & 6
    - `chr[2/7]/`
        - `chr[2/7]_ccre_hg19.bed:` ENCODE open chromatin candidate cis-regulatory elements (cCREs) to highlight promoters (promoter-like signature,) and enhancers (proximal and distal enhancer-like signature) combined from all cell types downloaded from the UCSC table browser (lifted over to hg19)
        - `encode_chr[2/7].tsv:` Transcription Factor (TF) ChIP-seq Clusters (130 cell types) from ENCODE 3
        - `akitaPredsPerIndivAcrossChr7exampleRegion.tsv` - More fine-grained/resolution predictions of 3D divergence across the chr7 example (MHs-vs-Vindija), generated in `Neanderthal_3D_Genome.ipynb`
- `ontologyEnrichment/`
    - `GWAS_Catalog_2019.txt` - ontology enrichment links for GWAS from Enrichr [https://maayanlab.cloud/Enrichr/#libraries](https://maayanlab.cloud/Enrichr/#libraries)
    - `Human_Phenotype_Ontology.txt` - ontology enrichment links for HPO from Enrichr [https://maayanlab.cloud/Enrichr/#libraries](https://maayanlab.cloud/Enrichr/#libraries)
    - `4DN_HFF_MicroC_hg19_processedTADs.bed` - TADs created from regions between TAD boundaries (created in `Neanderthal_3D_Genome.ipynb` notebook from `4DN_HFFboundariesMicroC_hg19.bed`(above))
    - `refSeqGeneCoordinatesSimplifiedProteinCoding.bed` - longest transcript from protein-coding genes (NM prefix) from NCBI RefSeq downloaded from the UCSC Table Browser (hg19)
- `introgression/`
    - `[INTERSECT/UNIQUE]_[sprime_segments/sprime_segments_neanMatchingFilter/vernot].bed.gz` - Sprime Segments from Browning et al. or S* segments from Vernot et al. given if they are observed in all 3 1kgp superpopulations (EAS, EUR, SAS) versus if they are unique to only one super-population
    - `individualIntrogressionCalls/chen20_neanderthal_seq_in_1kG_50kb.txt.gz` - 1000 Genome individual data downloaded from [https://drive.google.com/drive/folders/1mDQaDFS-j22Eim5_y7LAsTTNt5GWsoow](https://drive.google.com/drive/folders/1mDQaDFS-j22Eim5_y7LAsTTNt5GWsoow) , renamed

## results/

- `3DgenomeComparisons/` - pairwise quantification of 3D divergence between two individuals
    - `withArchaics/` - masked regions of archaic missingness to facilitate comparisons (fill with hg19 ref)
    - `withoutArchaics/` - not masked because we are only comparing between modern humans in regions with no missingness
    - `comparisonAcrossCellTypes/3dcompAcrossCellTypes_AFR_ESN_female_HG03105_vs_vindija.tsv.gz` - 3d comparisons across the genome using different metrics to compare between HFF and h1ESC contexts, generated with `compareCellTypes.py`, used for a supplemental fig
- `empiricalExpected3dDivergence.tsv` - 100 shuffles per 1 Mb window generated by `expected3DgenomeDivergence.py`. Used to determine the average expected 3D genome divergence based on shuffling observed sequence differences between two individuals (HG03105 & vindija)
- `ontologyEnrichment/` - intermediate and results files from functional annotation enrichment analyses, everything was done with 3 different sets of AH-MH divergent windows (”intersectNeanderthal” = windows divergent in all 3 Neanderthals, “unionNeanderthal” = windows divergent in any of the 3 Neanderthals, “denisova” = windows divergent in Denisovan) and with 2 sets of annotations (”hpo” = human phenotype ontology, “gwas” = genome-wide association study)
    - `observedCounts/[gwas/hpo]_ontology_counts_[unionNeanderthal/intersectNeanderthal/denisova].tsv` - Per divergent window set and annotation-type, specifies the observed links to each ontology annotation (gene-phenotype links)
    - `empiricCounts/empiric_[gwas/hpo]_ontology_counts_[unionNeanderthal/intersectNeanderthal/denisova].tsv.gz` - Per divergent window set and annotation-type, empirical counts of HPO/GWAS annotation links for the null hypothesis (n=500,000) so we can calculate p values and enrichment statistics for observed counts, generated by `goEnrichment.py`
    - `minPval_fdrCorrection: [gwas/hpo]_fdr_go_enr.csv.gz` - p-values based on empiric shuffles for annotation enrichment, used to generate an expected number of false positives for ontology enrichment in order to calculate FDR-corrected P values (q), generated by `goEnrichmentEmpiricalFDR.py`
    - `results/ontology_[hpo/gwas]_pvals_[unionNeanderthal/intersectNeanderthal/denisova].tsv` - ontology enrichment results and p-values from comparing the observedCounts (above) (per set and annotation-type) with the empiricCounts (above)
        - `results/ontology_[hpo/gwas]_pvals_intersectNeanderthal_bySystem_manual.tsv` - manually annotated version of `ontology_[hpo/gwas]_pvals_intersectNeanderthal.tsv` with phenotype trait categories
- `inSilicoMutagenesis/`
    - `rawFiles/inSilicoMut_[archaic]_chr[#]_[#].tsv.gz` - in silico mutagenesis for each 1 Mb window of interest. Each archaic variant is inserted into a MH (african) background one-by-one and the resulting 3D genome comparison is measured, generated with `inSilicoMutagenesis.py`
    - `intermediateFiles/` - a variety of interemediate files used to filter the insilico mutagenesis variants by those that cause 3D structure change and are either archaic specific or introgressed (all generated by `Neanderthal_3D_Genome.ipynb`)
- `divergentWindows/` - files with the `_lessStrict` suffix are using less stringent criteria to identify AH-MH divergent windows than those in the main text
    - `AH-MH_divergedWindows_filtered[_lessStrict].tsv` - AH-MH divergent windows called from the 1 Mb sliding window
    - `AH-MH_divergedWindows_filtered_mergeOverlapping[_lessStrict].tsv` - AH-MH divergent windows where overlapping windows are merged (these are supplemental tables)
    - `AH-MH_divergedWindows_[lessStrict_]inSilicoMut.tsv` - coordinates of 3D modifying variants in each AH-MH divergent window
    - `AH-MH_divergedWindows_inSilicoMut_withAnnots.tsv` - annotations of the 3D modifying variants in `AH-MH_divergedWindows_inSilicoMut.tsv` (eg. overlap with eQTL, TADs, genes, 1kgp AF, ctcf motifs, chen adaptive haplotypes) (supplemental table)
