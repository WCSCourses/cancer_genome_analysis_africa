########### Tutorial script for Module 4 Mutational Signatures #############
# Running SigProfilerExtractor to extract de novo signatures from 
# the HER2 data from CBIO portal

# see SigProfilerMatrixGenerator_tutorial.R for more 
# information on running the scripts if NOT running on the virtual machine

# load the library needed for calling python code from R
library("reticulate")

# now we need to point to the version of Python that we want to use
use_python("/home/manager/miniconda3/bin/") 
py_config()

# load the library needed for installing packages from Github
# and then install SigProfilerExtractorR
# NOTE: this is done for you on the virtual machine but will be needed 
# if running this script off of the virtual machine
#library(devtools)
#install_github("AlexandrovLab/SigProfilerExtractorR")

# load the SigProfilerExtractorR
# this is necessary so that we can use the tool from our script
library(SigProfilerExtractorR)

# we would normally install the different versions of the genome,
# but this should already be done from our tutorial using 
# SigProfilerMatrixGenerator
#install("GRCh37", rsync=FALSE, bash=TRUE)
#install("GRCh38", rsync=FALSE, bash=TRUE)


# There are many different ways that we can run SigProfilerExtractor.
# Here, we will go through 3 different ways of running the tool
  # EXAMPLE1: starting from the 21 VCF breast cancer files included with the tool
  # EXAMPLE2: starting from a MAF file for the TCGA BRCA HER2 samples
  # EXAMPLE3: starting from a mutational matrix for the TCGA BRCA HER2 samples
# NOTE: SigProfilerExtractor can only do extraction if SigProfilerMatrixGenerator
# has been run. In EXAMPLE1 and EXAMPLE2, when we use SigProfilerExtractor the tool
# will actually run SigProfilerMatrixGenerator as the first step. In EXAMPLE3, we 
# will give SigProfilerExtractor a mutational matrix that we generated as part of
# the SigProfilerMatrixGenerator_tutorial, in which case SigProfilerExtractor will
# not run SigProfilerMatrixGenerator and will instead go directly into the extraction.



########################  FIRST EXAMPLE ########################  
#### Using the 21 breast cancers that come with the tool #### 
################################################################

# First, we need to specify the location of the input data.
# Because the 21 breast cancers come with the tool, our
# path will just be the data itself
path_to_21_breast_cancers = importdata("vcf")

# Also, let's specify where the output will go
output_21_breast_cancers = "./outputs/SigProfilerExtractor_output_21_breast_cancers" 

# Next, we need to specify what type of input we are using
# The input can either be VCF or matrix format
# VCF files refer to an individual file per sample (like what we used for the 21 breast cancers in the SigProfilerMatrixGenerator tutorial)
# The matrix refers to the output from SigProfilerMatrixGenerator
# We will use the matrix output from the SigProfilerMatrixGenerator tutorial
input_type = "vcf"

# Finally, to run SigProfilerExtractorR we use the below command
# NOTE: We are setting nmf_replicates = 5 to decrease how long the command takes to run
# For accurate results, you will want to have nmf_replicates = 100 (or greater)
sigprofilerextractor(input_type, output_21_breast_cancers, path_to_21_breast_cancers, minimum_signatures=1, maximum_signatures=5,
                     nmf_replicates = 5)


# This command will print updates on the processes it's running and the time that these take
# Also - expect that this command will take a substantial amount of time to run

# This process will generate a lot of output in the specified output folder
# We'll go through how to interpret the output during the tutorial, and additional
# information is provided in the pdf accompanying the tutorials




########################  SECOND EXAMPLE ########################  
#### Using the MAF file for the TCGA BRCA HER2 samples #### 
################################################################

# First, we need to set the path to the MAF file
# This MAF file was generated in the SigProfilerMatrixGenerator_tutorial
path_to_her2_maf_file = "./outputs/SigProfilerMatrixGenerator_TCGA_BRCA_HER2/" 

# Also we also need to name our output folder
output_HER2_MAF = "./outputs/SigProfilerExtractor_output_TCGA_BRCA_HER2_MAF" 


# Next, we need to specify what type of input we are using
# Although the input is MAF, we need to set it VCF 
input_type = "vcf"

# Finally, to run SigProfilerExtractorR we use the below command
# NOTE: We are setting nmf_replicates = 5 to decrease how long the command takes to run
# For accurate results, you will want to have nmf_replicates = 100 (or greater)
# Since these are exomes, we want to do exome=T
sigprofilerextractor(input_type, output_HER2_MAF, path_to_her2_maf_file, minimum_signatures=1, maximum_signatures=5,
                     nmf_replicates = 5,exome=T)




########################  THIRD EXAMPLE ########################  
#### Using a mutational matrix for TCGA BRCA HER2 samples that #### 
#### we generated as part of SigProfilerMatrixGenerator_tutorial.R #### 
################################################################

# First, we need to set the path for where our mutational matrix is located
# This is what we generated in the SigProfilerMatrixGenerator_tutorial.R script
# Each type of mutational context format gets a different matrix 
# Let's start with SBS96
path_to_sbs96_HER2_matrix = './outputs/SigProfilerMatrixGenerator_TCGA_BRCA_HER2/output/SBS/MatrixGenerator_TCGA_BRCA_HER2.SBS96.exome' 

# Next, we need to specify what type of input we are using
# We use "matrix" because we're starting with a mutational matrix
# that we got as output from SigProfilerMatrixGenerator
input_type = "matrix"

# And we also need to decide what we want to name the output folder
# This will be created in our current directory unless a path is specified before the folder name
output_HER2_matrix = './outputs/SigProfilerExtractor_output_TCGA_BRCA_HER2_matrix' 

# Finally, to run SigProfilerExtractorR we use the below command
# NOTE: We are setting nmf_replicates = 5 to decrease how long the command takes to run
# For accurate results, you will want to have nmf_replicates = 100 (or greater)
# Since these are exomes, we want to do exome=T
sigprofilerextractor(input_type, output_HER2_matrix, path_to_sbs96_HER2_matrix, minimum_signatures=1, maximum_signatures=5,
                     nmf_replicates = 5,exome=T)


# NOTE: You'll notice that the extraction for the sbs96 context looks similar
# between EXAMPLE2 and EXAMPLE3 because the data that we're using is the same
# (the TCGA BRCA HER2 samples) even though the data type that we give SigProfilerExtractor 
# (MAF in EXAMPLE2 and a mutational matrix in EXAMPLE3) is different. 


