<h3>The materials in this repository are associated with the following manuscript:</h3> 
Hannah L. Noel, Hannah M. Maddox, Owen J. Singleton, Stephanie K. Drumheller, Jennifer M. DeBruyn. “Persistence of tegu lizard skin during decomposition selects for a unique fungal decomposer community” (in prep)

<h3>Below is a list of folders and explanations of the files within each:</h3>

**MOTHUR_Pipelines:** <br>
Includes detailed pipelines for processing raw amplicon sequences of both the 16S rRNA and ITS gene regions. The file “MOTHUR 16S Pipeline v 1.48.3 (SILVA v138.2)” corresponds to the bacterial 16S rRNA pipeline and the file “MOTHUR v1.48.3 ITS Pipeline (UNITE 10.0)” corresponds to the fungal ITS pipeline. Raw reads can be found on the NCBI Short Read Archive under PRJNA1537017. <br>

**MOTHUR_Output:** <br>
This folder contains the shared and taxonomy files from both the 16S rRNA and ITS pipelines. <br>
                                                                                        
- 16S: “tegu16s2.shared” and “tegu16s2.taxonomy”
- ITS: “TEGUITS2.shared” and “TEGUITS2.taxonomy”
  
**R_Analysis:** <br>
This folder contains the metadata as a .csv that is needed for all chemical and sequencing data analysis, along with R code. <br>

- "Tegu_FinalMetadata.csv”:Metadata file that is needed for all chemical and sequencing data analysis. 
- “Tegu_Chemical_Data_Analysis”: R file with code for statistical analysis and data visualization of soil chemical parameters. Metadata file needed. 
- “miseqR.R”: Contains functions to normalize and scale sequence data and allow for tests like adonis to be performed on phyloseq objects. Originally created by the Denef Lab (http://deneflab.github.io/MicrobeMiseq) 
- “16S_ClusteredAnalysis”: Quarto documents that has code for analysis of 16S rRNA data. Metadata, shared, taxonomy, and miseqR files needed for analysis. 
- “ITS_ClusterAnalysis”: Quarto document that has code for the analysis of ITS data. Metadata, shared, taxonomy, and miseqR files needed for analysis. 
