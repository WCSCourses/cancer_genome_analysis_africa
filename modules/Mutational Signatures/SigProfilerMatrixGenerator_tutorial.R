########### Tutorial script for Module 4 Mutational Signatures #############
# Running SigProfilerMatrixGenerator to generate a mutational matrix and 
# mutational profiles for each sample

# Two examples are included here - one for the 21 breast cancers and one for 
# the TCGA BRCA HER2 breast cancers from CBIO portal

# load the library needed for calling python code from R
# NOTE: the install.packages line will be needed if running this script
# outside of the virtual machine
#install.packages("devtools")
#install.packages("reticulate")
library("reticulate")

# now we need to point to the version of Python that we want to use
use_python("/home/manager/miniconda3/bin/") 
py_config()

# load the library needed for installing packages from Github
# and then install SigProfilerMatrixGenerator
# and the different versions of the human genome
# NOTE: this is done for you on the virtual machine but will be needed 
# if running this script off of the virtual machine
#library("devtools")
#install_github("AlexandrovLab/SigProfilerMatrixGeneratorR")
#library("SigProfilerMatrixGeneratorR")
#install('GRCh37', rsync=FALSE, bash=TRUE)
#install('GRCh38', rsync=FALSE, bash=TRUE)

# load the SigProfilerMatrixGeneratorR
# This is necessary so that we can use the tool from our script
library("SigProfilerMatrixGeneratorR")

# load tidyverse
# this is a library that helps with data handling 
#install.packages("tidyverse")
library(tidyverse)


########################  FIRST EXAMPLE ########################  
#### Using the 21 breast cancers that come with the tool #### 
################################################################


# First, lets load the data
# We will be starting with VCF files for the 21 breast cancers for this example
path_to_vcf_files = "./datasets/21BRCA_vcf/" 

# Next, lets set the name of our project
project_name_example1 = 'MatrixGenerator_21_breast_cancers'

# Finally, it's time to run SigProfilerMatrixGenerator
# This will generate the mutational matrices 
# and the mutational profiles for each sample
matrices <- SigProfilerMatrixGeneratorR(project_name_example1, "GRCh37", path_to_vcf_files, plot=T, exome=F)

# NOTE: output will be stored in the path that you specify for path_to_vcf_files
# This make take a few minutes to run
# If everything runs successfully, you will see the following message:

#Starting matrix generation for SNVs and DINUCs...Completed! Elapsed time: 110.53 seconds.
#Matrices generated for 21 samples with 0 errors. Total of 183916 SNVs, 969 DINUCs, and 0 INDELs were successfully analyzed.

# 3 new folders (input, output, and logs) will be generated
# You can go to output/plots/SBS_96_plots_MatrixGenerator_21_breast_cancers to
# see the same 21 mutational profile plots that we saw in the slides




########################  SECOND EXAMPLE ########################  
#### Using the HER2 data from CBIO Portal #### 
################################################################

# First, we need to convert the .csv file to a maf file
# and format it correctly

location_of_HER2_cbioportal_csv = "./datasets/tcga_brca/original_data/tcga_brca_her2.csv" 

# formatting
dataframe_for_maf=read.csv(location_of_HER2_cbioportal_csv) %>%
  select(Hugo=Hugo_Symbol,	Entrez=Entrez_Gene_Id,	Center,	Genome=NCBI_Build,
         Chrom=Chromosome,	Start=Start_Position,	End=End_Position,	Strand,	Classification=Variant_Classification,
         Type=Variant_Type,	Ref=Reference_Allele,	Alt1=Tumor_Seq_Allele1,	Alt2=Tumor_Seq_Allele2,	
         dbSNP=dbSNP_RS,	SNP_Val_status=dbSNP_Val_Status,	Tumor_sample=Tumor_Sample_Barcode)

# Next, we save the dataframe as a MAF file to the folder where we want the output
# first, let's make the directory
dir.create("./outputs/SigProfilerMatrixGenerator_TCGA_BRCA_HER2/")
write_tsv(dataframe_for_maf,"./outputs/SigProfilerMatrixGenerator_TCGA_BRCA_HER2/tcga_brca_her2.maf") 

# Now, lets save the path to our new MAF file and set our project name
path_to_her2_maf_file = "./outputs/SigProfilerMatrixGenerator_TCGA_BRCA_HER2/" 
project_name_example2 = "MatrixGenerator_TCGA_BRCA_HER2"

# Finally, let's run SigProfilerMatrixGeneratorR again
# We need to set exome=T because these are exome samples.
# This make take a few minutes to run 
matrices2 <- SigProfilerMatrixGeneratorR(project_name_example2, "GRCh37", path_to_her2_maf_file, plot=T, exome=T)

# If everything runs successfully, you will see the following message:

#Starting matrix generation for SNVs and DINUCs...Completed! Elapsed time: 395.08 seconds.
#Starting matrix generation for INDELs...Completed! Elapsed time: 82.7 seconds.
#Matrices generated for 78 samples with 18 errors. Total of 19123 SNVs, 47 DINUCs, and 407 INDELs were successfully analyzed.


