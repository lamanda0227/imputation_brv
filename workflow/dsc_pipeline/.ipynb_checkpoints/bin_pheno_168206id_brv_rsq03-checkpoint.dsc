#!/usr/bin/env dsc

%include module/simulate.dsc
%include module/brv.dsc

DSC:
  define:
    simulate_effect: fixed_effect_prop
    simulate_phenotype: bin_phenotype
    method: brv_pre_exome, brv_pre_exome_cadd, brv_pre_hrc, brv_pre_hrc_cadd, brv_pre_topmed, brv_pre_topmed_cadd, brv_pre_hrc_topmed, brv_pre_hrc_topmed_cadd, brv_pre_hrc_topmed_exome, brv_pre_hrc_topmed_exome_cadd
  run: parse_input * simulate_effect * simulate_phenotype * method * brv
  exec_path: code
  global:
    exome_folder:               "~/project/imputation-rvtest/workflows/imputation_aggregated_analysis/gene/exome"
    hrc_topmed_exome_folder:    "~/project/imputation-rvtest/workflows/imputation_aggregated_analysis/maf001/hrc_topmed_exome_rsq03"
    hrc_topmed_folder:          "~/project/imputation-rvtest/workflows/imputation_aggregated_analysis/maf001/hrc_topmed_rsq03"
    hrc_folder:                 "~/project/imputation-rvtest/workflows/imputation_aggregated_analysis/gene/hrc_rsq03"
    topmed_folder:              "~/project/imputation-rvtest/workflows/imputation_aggregated_analysis/gene/topmed_rsq03"
    out_folder:                 "~/project/git/imputation_brv/workflow/dsc_pipeline/bin_pheno_168206id_brv_rsq03/result_csv"
