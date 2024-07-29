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
    idx:                             R{1:3203}
    $gene:                           out$gene
    $exome_fname:                    out$exome_fname
    $hrc_fname:                      out$hrc_fname
    $topmed_fname:                   out$topmed_fname
    $hrc_topmed_fname:               out$hrc_topmed_fname
    $hrc_topmed_exome_fname:         out$hrc_topmed_exome_fname

## Generate Phenotype
### This module simulate binary phenotype for case-control study
bin_phenotype: simulate_phenotype.R + \
    R(out = simulate_bin_phenotype(hrc_topmed_exome_fname, deletrious_effect, protective_effect, causal_del[1], causal_del[2], prevalence, case_control[1], case_control[2]))
    hrc_topmed_exome_fname: $hrc_topmed_exome_fname
    deletrious_effect:      1.2, 1.5, 1.8
    protective_effect:      1
    causal_del:             (1,1), (0.75, 1), (0.5, 1)
    prevalence:             0.1, 0.2
    case_control:           (40000, 60000)
    $prev:                  out$prev
    $es:                    out$es
    $y:                     out$y
    $causal_prop:           out$causal_prop
    $case_control_ref:      out$case_control_ref
