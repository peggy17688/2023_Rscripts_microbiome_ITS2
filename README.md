# 2023_Rscripts_microbiome_ITS2
R scripts for manuscript submission to Ocean Microbiology. Will update paper title later.

# File description:
1. 20230301.R: Contains R codes for main statistical analyses including LME (Linear Mixed Models) for alpha diversity matrices (Chao1, Shannon Diversity, Inverse Simpson), PERMANOVA/PERMDISP for coral microbiome samples.
   <p> *Disclaimer* The initial steps (DADA2/phyloseq object creation) and DESeq2 analyses were run in the RConsole instead of the script window, hence the incompleteness. Details of the parameters are found in the Materials and Methods section of the manuscript.</p>
2. 20230401_ITS2.R: Contains R codes for main statistical analyses as mentioned (LME, PERMANOVA/PERMDISP) for coral Symbiodiniaceae ITS2 samples.
3. for_initial_samples.R: Contains R codes for PERMANOVA/PERMDISP analyses and NMDS (Bray-Curtis) plots for <i>Pocillopora acuta</i> and <i>Merulina ampliata</i> initial samples (before 2 weeks acclimation).
4. alleopathy_symbiont.R: This script was developed by my collaborator Dr. Lee Li Keat for the coral Symbiodiniaceae ITS2 figures in both manuscript and supplementary sections. 
<p>Raw sequences are found on the Sequence Read Archive (SRA) under submission number SUB14601872, BioProject PRJNA1135796 : AN INTEGRATED APPROACH TO STUDY THE EFFECTS OF MACROALGAL INTERACTIONS ON THE CORAL MICROBIOME IN AN URBANIZED REEF SYSTEM. </p>
