########### Tutorial script for Module 4 Mutational Signatures #############
# Running SigProfilerAssignment to assign reference mutational signatures
# to individual samples in the BRCA HER2 dataset

# see SigProfilerMatrixGenerator_tutorial.R for more 
# information on running the scripts if NOT running on the virtual machine

# load the library needed for calling python code from R
library("reticulate")

# now we need to point to the version of Python that we want to use
use_python("/home/manager/miniconda3/bin/") 
py_config()

# import SigProfilerAssignment
# this is necessary so that we can use the tool from our script
# NOTE: previously we just used library(NAME) to load the tools we need
# Here, we're using import because SigProfilerAssignment is ONLY 
# available as a Python package
SigProfilerAssignment <- import("SigProfilerAssignment")


# We will be running the assignment tool with our two datasets
# EXAMPLE1 will be using the 21 breast cancer data that comes with the SigProfiler tools
# EXAMPLE2 will be using the TCGA BRCA HER2 samples that we've been working with

########################  FIRST EXAMPLE ########################  
#### Using the 21 breast cancers that come with the tool #### 
################################################################

# First, we need to specify where the mutational matrix is for the 21 breast cancers.
# Remember, we generated this as part of the SigProfilerMatrixGenerator_tutorial. This
# was EXAMPLE1

path_to_sbs96_21_breast_cancers_matrix = './datasets/21BRCA_vcf/output/SBS/MatrixGenerator_21_breast_cancers.SBS96.all'

# Next, lets specify where we want the output to be
output_folder_21_breast_cancers = "./outputs/SBS96_assignment_output_21_breast_cancers" 

# Finally, let's run the assignment command
# cosmic_fit will assign the reference mutational signatures 
# from COSMIC to our 21 breast cancers
SigProfilerAssignment$Analyzer$cosmic_fit(path_to_sbs96_21_breast_cancers_matrix,output_folder_21_breast_cancers,
                                          genome_build="GRCh37", 
                                          make_plots=T,
                                          exome=F)

# If this runs correctly, you will see something like:

#Assigning COSMIC sigs or Signature Database ...... 
#|████████████████████████████████████████| 21/21 [100%] in 5.8s (3.61/s) 



########################  FIRST EXAMPLE ########################  
          #### Using the TCGA BRCA HER2 samples #### 
################################################################


# First, we need to specify where our mutational matrix is stored
# Remember: this is what we generated as part of SigProflerMatrixGenerator_tutorial.R,
# and this is the same matrix that we used in SigProfilerExtractor_tutorial.R
# Let's start with SBS96
path_to_sbs96_HER2_matrix = './outputs/SigProfilerMatrixGenerator_TCGA_BRCA_HER2/output/SBS/MatrixGenerator_TCGA_BRCA_HER2.SBS96.exome'

# Next, lets specify where we want the output to be
output_folder = "./outputs/SBS96_assignment_output_TCGA_BRCA_HER2" 

# Finally, let's run the assignment command
# cosmic_fit will assign the reference mutational signatures 
# from COSMIC to our input BRCA HER2 samples
# Since these are exomes, we want to do exome=T
SigProfilerAssignment$Analyzer$cosmic_fit(path_to_sbs96_HER2_matrix,output_folder,
                                          genome_build="GRCh37", 
                                          make_plots=T,
                                          exome=T)

# If this runs successfully, you will see:

# Assigning COSMIC sigs or Signature Database ...... 
# |████████████████████████████████████████| 78/78 [100%] in 14.9s (5.25/s) 


