# Rare Variant Aggregate Association Analysis of Imputed Variants

### Aims

------

This project aims to explore whether imputed data is a viable alternative to detect rare-variant aggregate associations when exome sequence data is unavailable

### Folders

-----

- `analysis`: analysis of exome sequence data and imputed genotypes

    - `exome_sequencde_processing.ipynb`: annotation and filtering steps for exome sequencing data

    - `imputed_data_processing.ipynb`: annotation and filtering steps for imputed data (HRC, TOPMed)

    - `data_summary_analysis.ipynb`: notebook that generates data summary tables and plots

    - `pcsk9_apoc3_brv_analysis.ipynb`: analysis of *PCSK9* and *APOC3* versus LDL and TG respectively, using BRV

- `workflow`: simulations and applications of rare variant association studies
    - `simulation_pipeline.ipynb`: notebook documenting simulation steps and result analysis, accompanied by `dsc_pipeline` folder, where code for simulations can be found
    - `additional_queries.ipynb`: notebook investigating $R^2$ (imputation quality) and $r^2$ (correlation) between imputed dataset and simulation data

### Collaborators

---

This project is a joint effort between Suzanne (CU), Gao (CU), Andy (Yale) and Paul (MCW), with contributions from 

- Amanda/Tianyi (CU)
- Elnaz (MSK)
