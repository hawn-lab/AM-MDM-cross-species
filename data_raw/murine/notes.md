# AM EXPERIMENT

- ex vivo murine AM
- H37Rv stim (6 hr) vs untreated, 3 AM conditions (Unvax_Mtb, scBCG_Mtb, CoTB_Mtb)
	- AM isolated 8 weeks after either scBCG or coMtb or ctrls

## AM FILES

`GSEA_AM_H37Rv_RNAseq.csv` : AM GSEA Hallmark pathways for AM Ex Vivo Stim RNA-seq data

`AM_H37Rv_RNAseq.xlsx` : counts and DEGs
	- DEG defined at abs log2 fold change > 1 & FDR < 0.05
	- Pairwise contrasts within AM conditions
 
# BMDM EXPERIMENT

- ex vivo WT BMDM
- H37Rv stim (4 hr) vs untreated

## BMDM FILES

`BMDM_H37Rv_RNA-seq.xlsx` : counts and DEGs
	- DEG defined at abs log2 fold change > 1 & FDR < 0.01
`BMDM_H37Rv_allLogCPM_GSO30.xlsx`
  - log CPM counts
`BMDM_H37Rv_FoldChange_FDR_GSO30.xlsx`
  - edgeR derived fold change and FDR (GSO)
  - GSEA from `gso30.Bl6.N.dRv.4` contrast of Mtb infected vs uninfected WT BMDMs, filtered with ave logCPM > 1.0. I also included the equivalent data for IFNAR-/- BMDMs that were part of that study, in case they are useful to us down the road, although they are not part of the analysis so far. 
GSEA output folders:
  - Mtb ex vivo stim AM: unvax
  - Mtb ex vivo stim AM: scBCG
  - Mtb ex vivo stim AM: coMtb
  - Mtb stim BMDM (Greg Olson-GSO)
