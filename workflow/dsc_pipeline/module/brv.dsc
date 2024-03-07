# BRV module

## This module performs the Burden of Rare Variant method

brv_pre_exome: brv_preprocess.R + \
    R(out = brv_preprocess(fname, case_control_ref, sample))
    fname:             $exome
    case_control_ref:  $case_control_ref
    sample:            "exome"
    $X:                out$X
    $sample_name:      out$sample
    
brv_pre_exome_cadd (brv_pre_exome):
    fname:             $exome_cadd
    sample:       "exome_cadd"
    
brv_pre_hrc: brv_preprocess.R + \
    R(out = brv_preprocess(fname, case_control_ref, sample))
    fname:             $hrc_fname
    case_control_ref:  $case_control_ref
    sample:            "hrc"
    $X:                out$X
    $sample_name:      out$sample
    
brv_pre_topmed (brv_pre_hrc):
    fname:             $topmed_fname
    sample:            "topmed"
    
brv_pre_hrc_topmed (brv_pre_hrc):
    fname:             $hrc_topmed_fname
    sample:            "hrc_topmed"    
    
brv_pre_hrc_topmed_exome (brv_pre_hrc):
    fname:             $X_hrc_topmed_exome
    sample:            "hrc_topmed_exome"    
    
brv_pre_hrc_cadd (brv_pre_hrc):
    fname:             $hrc_cadd_fname
    sample:            "hrc_cadd"  
    
brv_pre_topmed_cadd (brv_pre_hrc):
    fname:             $topmed_cadd_fname
    sample:            "topmed_cadd"  
     
brv_pre_hrc_topmed_cadd (brv_pre_hrc):
    fname:             $hrc_topmed_cadd_fname
    sample:            "hrc_topmed_cadd" 
    
brv_pre_hrc_topmed_exome_cadd (brv_pre_hrc):
    fname:             $hrc_topmed_exome_cadd_fname
    sample:            "hrc_topmed_exome_cadd" 
    
brv: brv.R + \
    R(out = BRV(X, y, sample, gene, out_folder, es, causal_prop, prev))
    X:                 $X
    y:                 $y
    sample:            $sample_name
    gene:              $gene
    es:                $es
    causal_prop:       $causal_prop
    prev:              $prev
    out_folder:        ${out_folder}
    $result_df:        out$result_df
    $sample_name:      out$sample_name
    
