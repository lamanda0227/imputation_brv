# BRV module

## This module performs the Burden of Rare Variant method

brv_exome: brv.R + \
    R(out = BRV(fname, case_control_ref, y, sample, gene, es, causal_prop, prev, out_folder))
    fname:             $exome_fname
    case_control_ref:  $case_control_ref
    y:                 $y
    sample:            "exome"
    gene:              $gene
    es:                $es
    causal_prop:       $causal_prop
    prev:              $prev
    out_folder:        ${out_folder}
    $result_df:        out$result_df
    $sample_name:      out$sample_name
    
    
brv_hrc (brv_exome):
    fname:             $hrc_fname
    sample:            "hrc"
    
brv_topmed (brv_exome):
    fname:             $topmed_fname
    sample:            "topmed"
    
brv_hrc_topmed (brv_exome):
    fname:             $hrc_topmed_fname
    sample:            "hrc_topmed"    
    
brv_hrc_topmed_exome (brv_exome):
    fname:             $hrc_topmed_exome_fname
    sample:            "hrc_topmed_exome"    
    