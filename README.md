Mechanistic Pathology (MechPath)

As the prevalence of non-alcoholic liver disease (NAFLD) and non-alcoholic steatohepatitis (NASH) are rapidly increasing, there is an urgent need for precise patient stratification to understand molecular profiles that can aid the drug development process. Currently, available patient stratification is solely based on fibrosis staging given after biopsy examination. However that may not consider biological complexity and heterogeneity in the disease manifestations. Furthermore, this lack of patient stratification may explain the lack of successful pharmacotherapies. Our study investigates the biological patterns of the NAFLD-NASH continuum and identifies patient subgroups independently of pathology scores (NAS score, Fibrosis score) in the discovery dataset (Govaere) and is validated in two independent datasets.

This work encompasses several comprehensive analyses. Firstly, the identification of 14 gene modules using Weighted Gene Correlation Network Analysis (WGCNA) in the discovery dataset that reflect well-known pathology related mechanisms such as cholesterol biosynthesis, immune pathways and hepatic fibrosis signaling. From 14 gene modules, the corresponding hub gene signature which embedded the biological complexity was used to stratify the patient population in the discovery dataset into 6 patient subgroups. Interestingly, DESeq2 analysis allowed us to identify fibrosis stage-related and subgroup-specific genes (one vs rest approach), which translated into pathways. Using subgroup-specific genes, pathway directionality (z-score) showed relevant biological differences across subgroups. For instance, subgroup 1 and subgroup 6 show opposite effects in both extracellular matrix organization and integrin signaling. Taken together, this shows that patient subgroups specific molecular profiles may aid in the assessment of potential pharmacotherapies as well as in target discovery for early virtual screening in drug development.

Patient stratification workflow
1.0_preprocessing_RNAseq.Rmd 
2.0_WGCNA_govaere.Rmd                                                         
2.1 DESeq2_Disease_Stage.Rmd
2.2_DESeq2_DEG_analysis.Rmd                                                    
3.0_Visualization Heatmaps_Patient Stratification.Rmd 
3.1 Visualization Heatmaps_Gene modules.Rmd          
3.2_Circular heatmap patient subgroups.Rmd
4.0_Multinomial regression 6 patient subgroups.R              