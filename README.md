# SCLC_rapid_research_autopsy

Written by Bea Szeitz.

Description: scripts to reproduce the analyses reported in the manuscript entitled "Comprehensive molecular profiling of inter-lesion heterogeneity in small cell lung cancer" by Megyesfalvi et al.

To run the code, clone the GitHub repository on your computer and download Data S1-3 accompanying the manuscript. Open the project in R (e.g., in RStudio), and run renv::restore() to automatically install the specified package versions into the project's private library.

## **Main scripts:**

The scripts are listed in the suggested order of execution.

<details>
  <summary><b>clean_input_data.Rmd</b></summary>
  
  - Description: Cleaning the input data for downstream scripts through QC and generation of RData objects.
  - Dependencies:
    - helpers/load_packages.R
    - helpers/utility_functions.R
    - helpers/color_list.R
    - files/Gene_counts_tximport.rds (note: this file will likely be available on GEO)
    - files/Data S1.xlsx
    - files/Data S2.xlsx
    - files/Data S3.xlsx
  - Outputs:
    - rdata/sample_annotations.RData
    - rdata/rna_seq.RData
    - rdata/proteomics.RData
    - rdata/ihc.RData
    - rdata/gene_protein_annotations.RData
    - figs/qc_prot_cohort1_outlier.png/.svg/.pdf
    - figs/qc_prot_cohort2_outlier.png/.svg/.pdf
    - figs/qc_prot_batchcorr.png/svg/pdf
</details>


<details>
  <summary><b>single_sample_scores.Rmd</b></summary>
  
  - Description: Generating single-sample scores from the RNA-Seq data.
  - Dependencies:
    - helpers/load_packages.R
    - helpers/utility_functions.R
    - helpers/color_list.R
    - rdata/gene_protein_annotations.RData
    - rdata/proteomics.RData
    - rdata/rna_seq.RData
    - rdata/gene_protein_annotations.RData
    - rdata/sample_annotations.RData
    - rdata/msigdb_v7.5.1.RData (note: see manuscript for instructions on how to get access to this file)
    - files/NE_and_non-NE_genes_Zhang2018.xlsx (note: see manuscript for instructions on how to get access to this file)
    - files/EMT_genes_Kohn2014.xlsx (note: see manuscript for instructions on how to get access to this file)
  - Outputs:
    - rdata/single_sample_scores_rnaseq.RData
    - supplementary_tables/Table S3 - Single-sample scores.xlsx
    - figs/NE_nonNE_scores_prot_correlation.png/.svg/.pdf
    - figs/Inflamed_scores_prot_correlation.png/.svg/.pdf
</details>


<details>
  <summary><b>genewise_corr.Rmd</b></summary>
  
  - Description: Performing the gene-wise correlation analysis between transcripts and proteins.
  - Dependencies:
    - helpers/load_packages.R
    - helpers/utility_functions.R
    - helpers/color_list.R
    - rdata/gene_protein_annotations.RData
    - rdata/proteomics.RData
    - rdata/rna_seq.RData
    - files/proteinatlas.tsv (note: see manuscript for instructions on how to get access to this file)
  - Outputs:
    - rdata/genewise_corr_ora.RData
    - supplementary_tables/Table S2 - Gene-wise correlations.xlsx
    - figs/genewise_corr_ctrl_vs_tumors_tissue_specific.png/.svg/.pdf
    - figs/genewise_corr_top_pathways.png/.svg/.pdf
    - figs/genewise_corr_tumors.png/.svg/.pdf
    - figs/genewise_corr_ctrls.png/.svg/.pdf
</details>



<details>
  <summary><b>omics_exploratory.Rmd</b></summary>
  
  - Description: Exploratory analysis of the proteomics and RNA-Seq data, including the detection of sample-composition associated proteins.
  - Dependencies:
    - helpers/load_packages.R
    - helpers/utility_functions.R
    - helpers/color_list.R
    - rdata/gene_protein_annotations.RData
    - rdata/proteomics.RData
    - rdata/rna_seq.RData
    - rdata/gene_protein_annotations.RData
    - rdata/single_sample_scores_rnaseq.RData
    - rdata/sample_overview_table.RData
    - files/Desai24_DataS3.xlsx (note: see manuscript for instructions on how to get access to this file)
  - Outputs:
    - rdata/histology_ora.RData
    - rdata/histology_linregr.RData
    - supplementary_tables/Table S4 - Histology associations.xlsx
    - figs/sample_overview.png/.svg/.pdf
    - figs/nr_proteins_transcripts_boxplot.png/.svg/.pdf
    - figs/unsupervised_clustering_prot_with_quality_comments.png/.svg/.pdf
    - figs/unsupervised_clustering_prot.png/.svg/.pdf
    - figs/overlap_with_desai24.png/.svg/.pdf
    - figs/tumor_associated_proteins_heatmap.png/.svg/.pdf
    - figs/stroma_associated_proteins_heatmap.png/.svg/.pdf
    - figs/necrosis_associated_proteins_heatmap.png/.svg/.pdf
    - figs/histology_associated_proteins_ORA_verbose.png/.svg/.pdf
    - figs/histology_associated_proteins_ORA.png/.svg/.pdf
</details>



<details>
  <summary><b>tumor_vs_control.Rmd</b></summary>
  
  - Description: Differential expression analysis between tumors and normals.
  - Dependencies:
    - helpers/load_packages.R
    - helpers/utility_functions.R
    - helpers/color_list.R
    - rdata/gene_protein_annotations.RData
    - rdata/proteomics.RData
    - rdata/rna_seq.RData
    - rdata/gene_protein_annotations.RData
    - rdata/single_sample_scores_rnaseq.RData
    - rdata/reactome_hierarchy.RData (note: instructions on how to create this RData object will be available later)
    - files/Gene_counts_tximport.rds (note: this file will likely be available on GEO)
    - files/Liu2024_TableS3.xlsx (note: see manuscript for instructions on how to get access to this file)
    - files/SCLC_RRA protein list_KB_MP_BS.xlsx (note: this file isn't available at the moment)
    - files/proteinatlas.tsv (note: see manuscript for instructions on how to get access to this file)
  - Outputs:
    - rdata/tumor_vs_ctrl_results.RData
    - rdata/tumorspec_pGSEA.RData
    - rdata/tumor_DEA_proteins.RData
    - supplementary_tables/Table S5 - Tumor vs ctrl multiomic comparison.xlsx
    - figs/tumor_vs_ctrl_upsetplot_onlyUP.png/.svg/.pdf
    - figs/tumor_vs_ctrl_upsetplot.png/.svg/.pdf
    - figs/ours_vs_Liu_log2fcs_tumor_vs_ctrl.png/.svg/.pdf
    - figs/tumor_vs_ctrl_gsea.png/.svg/.pdf
    - figs/tumor_vs_ctrl_volcano.png/.svg/.pdf
    - figs/wishlist_proteins_tumor_vs_ctrl.png/.svg/.pdf
</details>


<details>
  <summary><b>tumor_wise_similarities.Rmd</b></summary>
  
  - Description: Tumor-wise similarity calculations to highlight within-patient homogeneity.
  - Dependencies:
    - helpers/load_packages.R
    - helpers/utility_functions.R
    - helpers/color_list.R
    - rdata/gene_protein_annotations.RData
    - rdata/proteomics.RData
    - rdata/rna_seq.RData
    - rdata/gene_protein_annotations.RData
    - rdata/single_sample_scores_rnaseq.RData
    - rdata/histology_linregr.RData
    - rdata/tumor_vs_ctrl_results.RData
    - rdata/diagnosis_time_info.RData (note: this information has not been made publicly available about the patients, should be skipped when rerunning the script)
  - Outputs:
    - supplementary_tables/rna_scores_per_patient.xlsx
    - figs/patient_summary_distances.png/.svg/.pdf
    - figs/rna_scores_per_patient.png/.svg/.pdf
    - figs/corr_heatmap_patient_summary_avrtumor.png/.svg/.pdf
</details>

<details>
  <summary><b>clinically_relevant_proteins.Rmd</b></summary>
  
  - Description: Overview and analysis of clinically relevant proteins.
  - Dependencies:
    - helpers/load_packages.R
    - helpers/utility_functions.R
    - helpers/color_list.R
    - rdata/gene_protein_annotations.RData
    - rdata/proteomics.RData
    - rdata/rna_seq.RData
    - rdata/gene_protein_annotations.RData
    - rdata/single_sample_scores_rnaseq.RData
    - rdata/histology_linregr.RData
    - rdata/tumor_vs_ctrl_results.RData
    - files/SCLC_RRA protein list_KB_MP_BS.xlsx (note: this file isn't available at the moment)
  - Outputs:
    - supplementary_tables/Table S6 - Clinically relevant proteins.xlsx
    - figs/clin_rel_proteins_MAD_density.png/.svg/.pdf
    - figs/wishlist_proteins_variability.png/.svg/.pdf
</details>


<details>
  <summary><b>stable_expression.Rmd</b></summary>
  
  - Description: Detection of proteins with stable expression across the tumor sites.
  - Dependencies:
    - helpers/load_packages.R
    - helpers/utility_functions.R
    - helpers/color_list.R
    - rdata/gene_protein_annotations.RData
    - rdata/proteomics.RData
    - rdata/rna_seq.RData
    - rdata/gene_protein_annotations.RData
    - rdata/single_sample_scores_rnaseq.RData
    - rdata/histology_linregr.RData
    - rdata/tumor_vs_ctrl_results.RData
    - rdata/tumor_DEA_proteins.RData
  - Outputs:
    - rdata/median_FCs.RData
    - rdata/stable_expr_ora.RData
    - supplementary_tables/Table S7 - Stably expressed proteins.xlsx
    - figs/stable_proteins_targetable_boxplot.png/.svg/.pdf
    - figs/stable_proteins_heatmap_targetable.png/.svg/.pdf
    - figs/stable_proteins_heatmap_all.png/.svg/.pdf
    - figs/stably_expressed_proteins_ORA.png/.svg/.pdf
</details>


<details>
  <summary><b>tumor_site_specific_proteins.Rmd</b></summary>
  
  - Description: Detection of tumor-site specific proteins.
  - Dependencies:
    - helpers/load_packages.R
    - helpers/utility_functions.R
    - helpers/color_list.R
    - rdata/gene_protein_annotations.RData
    - rdata/proteomics.RData
    - rdata/rna_seq.RData
    - rdata/gene_protein_annotations.RData
    - rdata/single_sample_scores_rnaseq.RData
    - rdata/histology_linregr.RData
    - rdata/tumor_vs_ctrl_results.RData
  - Outputs:
    - supplementary_tables/site_specific_proteins.xlsx
    - figs/site_specific_boxplot.png/.svg/.pdf
    - figs/site_specific_heatmap.png/.svg/.pdf
</details>


<details>
  <summary><b>ihc_exploratory.Rmd</b></summary>
  
  - Description: Analysis of the IHC data.
  - Dependencies:
    - helpers/load_packages.R
    - helpers/utility_functions.R
    - helpers/color_list.R
    - rdata/gene_protein_annotations.RData
    - rdata/proteomics.RData
    - rdata/rna_seq.RData
    - rdata/gene_protein_annotations.RData
    - rdata/ihc.RData
  - Outputs:
    - rdata/IHC_per_tumor.RData
    - supplementary_tables/IHC_summary_statistics.xlsx
    - figs/cd3_high_vs_low_heatmap.png/.svg/.pdf
    - figs/h_score_overview.png/.svg/.pdf
    - figs/ihc_subtype_summary_per_tumor.png/.svg/.pdf
    - figs/tf_vs_cd3_per_tumor.png/.svg/.pdf
    - figs/ihc_per_tumor_summary.png/.svg/.pdf
    - figs/ihc_vs_rnaseq_subtype.png/.svg/.pdf
    - figs/ihc_subtype_vs_other_marker_ihc.png/.svg/.pdf
    - figs/ihc_subtype_vs_protein_expr.png/.svg/.pdf
    - figs/ihc_h_scores_per_tumor.png/.svg/.pdf
    - figs/circular_ihc_subtype_summary_per_tumor.png/.svg/.pdf
</details>


<details>
  <summary><b>immune_hot_vs_cold_markers.Rmd</b></summary>
  
  - Description: Differential protein expression analysis between immune hot vs cold tumors.
  - Dependencies:
    - helpers/load_packages.R
    - helpers/utility_functions.R
    - helpers/color_list.R
    - rdata/gene_protein_annotations.RData
    - rdata/proteomics.RData
    - rdata/rna_seq.RData
    - rdata/gene_protein_annotations.RData
    - rdata/single_sample_scores_rnaseq.RData
    - rdata/IHC_per_tumor.RData
    - rdata/reactome_hierarchy.RData (note: instructions on how to create this RData object will be available later)
  - Outputs:
    - rdata/ihc_markers_pGSEA.RData
    - supplementary_tables/Table S8 - Immune-hot vs -cold comparison.xlsx
    - figs/cd3_high_vs_low_heatmap.png/.svg/.pdf
    - figs/cd3_high_vs_low_gsea.png/.svg/.pdf
</details>



## **Helper scripts:**

These scripts are used to load all necessary packages, define custom functions and color lists.

<details>
  <summary><b>helpers/load_packages.R</b></summary>
  
  - Description: Used to load all required packages at the beginning of the script.
  - Dependencies: -
  - Outputs: -
</details>

<details>
  <summary><b>helpers/color_list.R</b></summary>
  
  - Description: Colors specified for figures.
  - Dependencies: -
  - Outputs: -
</details>

<details>
  <summary><b>helpers/utility_functions.R</b></summary>
  
  - Description: Custom functions used in the scripts.
  - Dependencies: -
  - Outputs: -
</details>