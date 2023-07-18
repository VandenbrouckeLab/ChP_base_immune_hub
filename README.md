# ChP_base_immune_hub
“Identification of the choroid plexus base as immune hub at the CSF interface”

## About

The choroid plexus (ChP) remains one of the most understudied tissues of the CNS, despite its ability to influence the function of the entire CNS. Positioned within the brain ventricles, the choroid plexus is primarily comprised of  specialized epithelial cells that produce the cerebrospinal fluid (CSF). We recently discovered the ChP base barrier cells (BBCs) at the ChP base and attachments sites with the brain. Both epithelial- and BBCs express tight junctions and together form the blood‑CSF barrier. They encapsulate the ChP stroma which contains fenestrated capillaries, fibroblasts and a broad range of immune cells. Previous reports have demonstrated, during disease, immune cell infiltration and accumulation in the ChP stroma, with suggestions of subsequent translocation across the epithelial layer to access the CSF, although robust evidence is limited. Here, we report that, in a mouse model for systemic inflammation, the ChP is infiltrated by innate immune cells and resident DCs adopt a migratory phenotype. DCs preferentially reside near the ChP base. In CCR7 KO mice, DCs accumulate in both naïve and systemically inflamed conditions, potentially unable to sense CCR7 ligands in the CSF and therefore suggesting a continuous egress of DCs near/across the ChP BBCs. In the EAE mouse model for multiple sclerosis, we observed T-cell accumulation adjacent to the ChP BBCs with frequent examples suggesting migration into the CSF. Altogether, we have identified the ChP base as an immune hub which might serve as gateway to reach the CNS.


## Overview scripts

Here's an overview of the various R scripts used in processing the scRNA-Seq data in the manuscript Verhaege et al.:
- 1_script_LPS_kinetics_Bulk_RNAseq.R: edgeR-Limma workflow for processing LPS kinetics bulk RNA-seq data from adult ChP comparing all three timepoints (1h, 6h and 24h post-injection) individually with the untreated condition.
- 2_script_scRNAseq.R: Seurat workflow for processing scRNA-Seq data of adult mouse lateral and fourth ventricle ChP in naive and systemic inflammation conditions.

## Citation

...
