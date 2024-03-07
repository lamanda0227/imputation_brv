# Simulation module

## Generate Separate Input File
## This module passes a separate input file for each gene separately
parse_input: input.R + \
    R(out = parse_input(exome_folder, 
                        hrc_topmed_exome_folder, hrc_topmed_folder, 
                        hrc_folder, topmed_folder, 
                        idx))
    exome_folder:                    ${exome_folder}
    hrc_topmed_exome_folder:         ${hrc_topmed_exome_folder}
    hrc_topmed_folder:               ${hrc_topmed_folder}
    hrc_folder:                      ${hrc_folder}
    topmed_folder:                   ${topmed_folder}
    idx:                             R{1:2}
    $gene:                           out$gene
    $bim_fname:                      out$bim_fname
    $exome:                          out$exome
    $exome_cadd:                     out$exome_cadd
    $hrc_topmed_exome_fname:         out$hrc_topmed_exome_fname
    $hrc_topmed_exome_cadd_fname:    out$hrc_topmed_exome_cadd_fname
    $hrc_fname:                      out$hrc_fname
    $hrc_cadd_fname:                 out$hrc_cadd_fname
    $topmed_fname:                   out$topmed_fname
    $topmed_cadd_fname:              out$topmed_cadd_fname
    $hrc_topmed_fname:               out$hrc_topmed_fname
    $hrc_topmed_cadd_fname:          out$hrc_topmed_cadd_fname
    
## Generate Effect Sizes
### This module simulate fixed effect by proportion of variants
fixed_effect_prop: simulate_effect.R + \
    R(out = generate_fixed_effect_prop(bim_fname, hrc_topmed_exome_fname, option, deletrious_effect, protective_effect, causal_del[1], causal_del[2]))
    bim_fname:              $bim_fname
    hrc_topmed_exome_fname: $hrc_topmed_exome_fname
    option:                 "hrc_topmed_exome"
    deletrious_effect:      1.2, 1.5, 1.8
    protective_effect:      1
    causal_del:             (0,0), (1,1), (0.75, 1), (0.5, 1)
    $effect:                out$effect
    $es:                    out$es
    $causal_prop:           out$causal_prop

## Generate Phenotype
### This module simulate binary phenotype for case-control study
bin_phenotype: simulate_phenotype.R + \
    R(out = simulate_bin_phenotype(effect, prevalence, case_control[1], case_control[2], hrc_topmed_exome_fname))

    effect:                 $effect
    prevalence:             0.1, 0.2
    case_control:           (40000, 60000)
    hrc_topmed_exome_fname: $hrc_topmed_exome_fname
    $prev:                  out$prev
    $y:                     out$y
    $X_hrc_topmed_exome:    out$X
    $case_control_ref:      out$case_control_ref
