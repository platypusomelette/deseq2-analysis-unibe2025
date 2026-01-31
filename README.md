This repository contains code for DESeq2 interaction analysis of bulk RNAseq data from Singhania et al. 2019. The data is subsetted from blood samples of _Ifnar−/− × Ifngr−/−_ vs control mice, either healthy or infected with _Toxoplasma gondii_.



## Part 1: BASH Scripts

Pipeline:

1. Trim (Calliefastp)
2. Quality Control (Calliefastqc \& Calliemultiqc)
3. Map reads to reference genome (Calliehisat)
4. Quality Control (Calliefastqc \& Calliemultiqc)
5. Count number of reads per gene via FeatureCounts (Calliesubread)

## Part 2: R Scripts

All analyses are carried out using a 2-factor design.

### /R\_Scripts

### /raw\_inputfiles

Raw data input files are used by deseq\_v3.R and customGO\_v1:

1. **counts.txt**: contains read counts, output from FeatureCounts
2. **samplenames.txt**: copy-pasted tab-delimited table of sample IDs and experimental group (ex. "Blood\_WT\_Case")
3. **41467\_2019\_10601\_MOESM4\_ESM\_minimal.csv**: an edited version of Supplementary Data 2, "Genes in the lung modules with average normalised read counts for each group", used by customGO\_v3\_Lmodules to map ENSEMBL IDs to custom gene sets
4. **41467\_2019\_10601\_MOESM5\_ESM\_minimal.csv**: an edited version of Supplementary Data 3, "Genes in the blood modules with average normalised read counts for each group", used by customGO\_v3 to map ENSEMBL IDs to custom gene sets
5. **41467\_2019\_10601\_MOESM6\_ESM\_minimal.csv**: an edited version of Supplementary Data 4, "Annotation of the lung modules", used by customGO\_v3\_Lmodules to map module names to human-readable functional annotations
6. **41467\_2019\_10601\_MOESM7\_ESM\_minimal.csv**: an edited version of Supplementary Data 5, "Annotation of the blood modules", used by customGO\_v3 to map module names to human-readable functional annotations

### Scripts:

1. **deseq\_v4.R**: DESeq2 2-factor design (WT/DKO, Case/Control): type + condition + type:condition. Saves DESeqResultsObjects for all contrasts into /downstream\_inputs that can be used for downstream scripts
2. **GO\_v4.R**: GO enrichment analysis and dotplots for DE\_DKO\_adj. Creates three plots: one for all genes, one within subset of downregulated genes, one within subset of upregulated genes
3. **customGO\_v3.R:** enrichment analysis of custom gene sets, representing the custom blood modules defined by Singhania et al. 2019.
      - Runs a similar enrichment analysis to GO\_v3.R on DE\_DKO\_adj, DE\_healthy, and DE\_diseased. 
      - It also exports customGO\_typeI\_ranks.txt and customGO\_typeII\_ranks.txt, which ranks each module by its percentage of type I IFN and type II IFN genes (as annotated in Supp. data 3).
      - It also exports IFNtype\_genelist\_inputs.RData
4. **customGO\_v3\_Lmodules.R:** enrichment analysis of custom gene sets, representing the custom lung modules defined by Singhania et al. 2019.
5. **volcanoplot\_v2.R:** creates volcano plots for DE\_diseased. Also creates an alternate version which overlays the significant type I IFN genes on the plot.

### /downstream\_inputs

1. **DE\_healthy\_inputs.RData**: "effect of DKO vs WT, for healthy case". type\_DKO\_vs\_WT
2. **DE\_diseased\_inputs.RData**: "effect of DKO vs WT, for diseased case". c("type\_DKO\_vs\_WT", "typeDKO.conditionCase")
3. **DE\_WT\_inputs.RData:** "effect of case vs control, for WT". condition\_Case\_vs\_Control
4. **DE\_DKO\_inputs.RData**: "effect of case vs control, for DKO". this is the BULK difference and is not adjusted by the observations in WT. c("condition\_Case\_vs\_Control", "typeDKO.conditionCase")
5. **DE\_DKO\_adj\_inputs.RData:** "effect of case vs healthy control, for DKO, ADJUSTED by effect in WT". this is the ADDITIONAL difference, beyond what changes in WT. typeDKO.conditionCase
6. **IFNtype\_genelist\_inputs.RData**: list created by customGO\_v3.R which allows for labeling of type I or type II genes in volcanoplot\_v2.R



### /plots

deseq.pdf:
1. DispEsts plot to assess good fit after variance stabilizing transform (VST)
2. PCA

GO.pdf:
1. GO enrichment analysis dotplot of all DE genes in DE\_diseased
2. GO enrichment analysis dotplot of only upregulated genes
3. GO enrichment analysis dotplot of only downregulated genes

customGO.pdf:
1. Custom enrichment analysis dotplot of all DE genes in DE\_diseased, annotated with blood modules
2. Custom enrichment analysis dotplot of only upregulated genes
3. Custom enrichment analysis dotplot of only downregulated genes

customGO_Lmodules.pdf:
1. Custom enrichment analysis dotplot of all DE genes in DE\_diseased, annotated with lung modules
2. Custom enrichment analysis dotplot of only upregulated genes
3. Custom enrichment analysis dotplot of only downregulated genes

volcanoplot.pdf:
1. volcano plot of DE\_DKO\_adj
2. volcano plot of DE\_DKO\_adj with type I IFN genes labeled
3. volcano plot of DE\_diseased
4. volcano plot of DE\_diseased with type I IFN genes labeled
