#!/usr/bin/env dsc

%include module/simulate.dsc
%include module/brv.dsc

DSC:
  define:
    simulate_phenotype: bin_phenotype
    brv: brv_exome, brv_hrc, brv_topmed, brv_hrc_topmed, brv_hrc_topmed_exome
  run: parse_input * simulate_phenotype * brv
  exec_path: code
  global:
    exome_folder:               "~/project/imputation-rvtest/workflows/imputation_aggregated_analysis/gene/exome_maf001"
    hrc_topmed_exome_folder:    "~/project/imputation-rvtest/workflows/imputation_aggregated_analysis/gene/hrc_topmed_v3_exome_rsq03_maf001"
    hrc_topmed_folder:          "~/project/imputation-rvtest/workflows/imputation_aggregated_analysis/gene/hrc_topmed_v3_rsq03_maf001"
    hrc_folder:                 "~/project/imputation-rvtest/workflows/imputation_aggregated_analysis/gene/hrc_rsq03_maf001"
    topmed_folder:              "~/project/imputation-rvtest/workflows/imputation_aggregated_analysis/gene/topmed_v3_rsq03_maf001"
    out_folder:                 "/mnt/mfs/statgen/tl3031/git/imputation_brv/workflow/dsc_pipeline_rsq03/bin_pheno_168206id_brv_rsq03/result_csv"
