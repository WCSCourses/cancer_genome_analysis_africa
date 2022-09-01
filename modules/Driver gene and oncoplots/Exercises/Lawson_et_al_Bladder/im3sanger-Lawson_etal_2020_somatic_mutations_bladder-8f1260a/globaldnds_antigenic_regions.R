# Inigo Martincorena - 01.05.2020
# Script to calculate global dN/dS in antigenic regions (based on PMID:31768072)
# applied to the bladder whole-exome data

globaldnds_gr = function(mutations, gr, gene_list = NULL, L = NULL, sm = "192r_3w") {
    
    # Restricting the mutation table to the gr regions
    gr_muts = GenomicRanges::GRanges(mutations$chr, IRanges::IRanges(mutations$pos,mutations$pos))
    ol = as.data.frame(GenomicRanges::findOverlaps(gr_muts, gr, type="any", select="all"))
    mutations = mutations[unique(ol[,1]), ]
    
    # If the L matrix is not provided, we calculate it
    if (is.null(L)) {
        
        # Generating the RefCDSnew object with codon annotation if it does not exist
        if (!exists("RefCDSnew")) {
            library(GenomicRanges)
            data("refcds_hg19", package="dndscv") # Load RefCDS
            RefCDSnew = buildcodon(RefCDS) # Adding codon information to RefCDS
        }
        
        ## 1. Find antigenic regions overlapping RefCDS sequences 
        
        ind = setNames(1:length(RefCDS), sapply(RefCDS,function(x) x$gene_name))
        gr_genes_ind = ind[gr_genes$names]
        ol = as.data.frame(GenomicRanges::findOverlaps(gr, gr_genes, type="any", select="all"))
        gr = gr[ol[,1]]
        geneind = gr_genes_ind[ol[,2]]
        
        ## 2. Calculating the L matrix: using codondnds code
        
        for (j in 1:length(RefCDS)) {
            if (RefCDSnew[[j]]$strand==1) {
                RefCDSnew[[j]]$cdspos = unlist(apply(RefCDSnew[[j]]$intervals_cds, 1, function(x) x[1]:x[2]))
            } else { # -1 strand
                RefCDSnew[[j]]$cdspos = unlist(rev(apply(RefCDSnew[[j]]$intervals_cds, 1, function(x) x[2]:x[1])))
            }
        }
        
        L = array(0, dim = c(192,4))
        for (j in 1:length(geneind)) {
            ind = rep(RefCDSnew[[geneind[j]]]$cdspos, each=3) %in% (start(gr[j]):end(gr[j]))
            cod_impt = RefCDSnew[[geneind[j]]]$codon_impact[ind]
            cod_rate = RefCDSnew[[geneind[j]]]$codon_rate[ind]
            for (h in 1:length(cod_impt)) {
                L[cod_rate[h],cod_impt[h]] = L[cod_rate[h],cod_impt[h]] + 1
            }
            if (round(j/5000)==(j/5000)) { message(sprintf('    %0.3g%% ...', round(j/length(geneind),2)*100)) }
        }
    }
    
    dndsout = dndscv(mutations, gene_list=gene_list, outmats=T, outp=1, max_muts_per_gene_per_sample=Inf, max_coding_muts_per_sample=Inf)
    N = apply(dndsout$N, c(1,2), sum) # N matrix
    
    ## 3. Global dN/dS ratios using Poisson regression
    
    # Substitution model (The user can also input a custom substitution model as a matrix)
    if (length(sm)==1) {
        data(list=sprintf("submod_%s",sm), package="dndscv")
    } else {
        substmodel = sm
    }

    # Subfunction: fitting substitution model
    
    fit_substmodel = function(N, L, substmodel) {
        
        l = c(L); n = c(N); r = c(substmodel)
        n = n[l!=0]; r = r[l!=0]; l = l[l!=0]
        
        params = unique(base::strsplit(x=paste(r,collapse="*"), split="\\*")[[1]])
        indmat = as.data.frame(array(0, dim=c(length(r),length(params))))
        colnames(indmat) = params
        for (j in 1:length(r)) {
            indmat[j, base::strsplit(r[j], split="\\*")[[1]]] = 1
        }
        
        model = glm(formula = n ~ offset(log(l)) + . -1, data=indmat, family=poisson(link=log))
        mle = exp(coefficients(model)) # Maximum-likelihood estimates for the rate params
        ci = exp(confint.default(model)) # Wald confidence intervals
        par = data.frame(name=gsub("\`","",rownames(ci)), mle=mle[rownames(ci)], cilow=ci[,1], cihigh=ci[,2])
        return(list(par=par, model=model))
    }
    
    # Fitting all mutation rates and the 3 global selection parameters
    
    par = fit_substmodel(N, L, substmodel)$par # Original substitution model
    globaldnds = par[c("wmis","wnon"),]
    
    return(list(N=N, L=L, globaldnds=globaldnds))
}


## 0. Preparing the environment for the bladder analysis

# Annotating RefCDS with codon information (optional, if absent it will be generated inside the function)

if (!exists("RefCDSnew")) {
    library(dndscv)
    library(GenomicRanges)
    data("refcds_hg19", package="dndscv") # Load RefCDS
    RefCDSnew = buildcodon(RefCDS) # Adding codon information to RefCDS
}

# Loading the bladder mutations

exome_mutations = "bld_exome_caveman_and_pindel_calls.tsv"
lcm_file = "2019-10-01_LCM_database.rds"
lcmdata = readRDS(lcm_file)
patient_file = "bladder_patient_info_2019-10-30.csv"
patientdata = read.table(patient_file, header=1, sep=",", stringsAsFactors=F)

muts.exome = read.table(exome_mutations, sep="\t", header=1, stringsAsFactors=F)
urot_ids = as.vector(lcmdata$SupplierSampleName[lcmdata$Feature=="Urothelium" & as.vector(lcmdata$Donor) %in% patientdata$external_id[patientdata$patient_type=="transplant"]])
mutations = muts.exome[which(muts.exome$sampleID %in% urot_ids),]
mutations$sampleID = substr(mutations$sampleID,1,7)
mutations = unique(mutations)

# Gene list (bladder passengers)

bladder_cancer_genes_file = "Bladder_cancer_genes_PMID_28988769_29056346.tsv"
bladder_cancer_genes = read.table(bladder_cancer_genes_file, header=0, sep="\t", stringsAsFactors=F)[,1]
non_driver_genes = setdiff(sapply(RefCDSnew, function(x) x$gene_name), bladder_cancer_genes)


## 1. Global dN/dS in the entire exome

dndsout = dndscv(mutations, gene_list=non_driver_genes, outmats=T, outp=1, max_muts_per_gene_per_sample=Inf, max_coding_muts_per_sample=Inf)
dnds1 = dndsout$globaldnds[c("wmis","wnon"), -1]


## 2. Antigenic regions based on the 6 most common HLA alleles (PMID:31768072)

load("antigenic_regions_PMID31768072.rda")
chrs = c(1:22,"X","Y")
gr_hla6 = renameSeqlevels(annotation_gr_HLA_Mean, setNames(chrs, paste("chr",chrs,sep="")))

out2 = globaldnds_gr(mutations = mutations, gr = gr_hla6, gene_list = non_driver_genes) # If L is not provided it will generate it (slow)
L = out2$L; save(L, file="antigenic_regions_PMID31768072_Lmatrix.rda")
#load("antigenic_regions_PMID31768072_Lmatrix.rda")
#out2 = globaldnds_gr(mutations = mutations, gr = gr_hla6, gene_list = non_driver_genes, L = L) # Using pre-saved L
dnds2 = out2$globaldnds[c("wmis","wnon"), -1]

## 3. IEDB epitopes (PMID:31768072)

gr_iedb = readRDS("IEDB_epitopes_gr_PMID31768072.rds")

out3 = globaldnds_gr(mutations = mutations, gr = gr_iedb, gene_list = non_driver_genes) # If L is not provided it will generate it (slow)
L = out3$L; save(L, file="IEDB_epitopes_gr_PMID31768072_Lmatrix.rda")
#load("IEDB_epitopes_gr_PMID31768072_Lmatrix.rda")
#out3 = globaldnds_gr(mutations = mutations, gr = gr_iedb, gene_list = non_driver_genes, L = L) # Using pre-saved L
dnds3 = out3$globaldnds[c("wmis","wnon"), -1]

## 4. IEDB epitopes (PMID:31768072) intersected with antigenic regions using 6 common HLA (PMID:31768072)

gr_hla6 = GenomicRanges::GRanges(gsub("chr","",as.vector(seqnames(annotation_gr_HLA_Mean))), IRanges(start(annotation_gr_HLA_Mean),end(annotation_gr_HLA_Mean)))
gr_iedb = readRDS("IEDB_epitopes_gr_PMID31768072.rds")
gr_iedb = GenomicRanges::GRanges(as.vector(seqnames(gr_iedb)), IRanges(start(gr_iedb),end(gr_iedb)))
gr_antigens = intersect(gr_iedb, gr_hla6)
out4 = globaldnds_gr(mutations = mutations, gr = gr_antigens, gene_list = non_driver_genes) # If L is not provided it will generate it (slow)
dnds4 = out4$globaldnds[c("wmis","wnon"), -1]

## 5. Plotting the results

dev.new(width=6, height=4)
b = barplot(cbind(dnds1[,1],dnds2[,1],dnds3[,1],dnds4[,1]), beside=T, ylim=c(0,max(cbind(dnds1[,3],dnds2[,3],dnds3[,3],dnds4[,3]))),
            ylab="dN/dS", col=c("grey30","grey60"), border=NA, las=1, names.arg = c("All sites","6HLA","IEDB","6HLA+IEDB"))
segments(x0=b, y0=cbind(dnds1[,2],dnds2[,2],dnds3[,2],dnds4[,2]), y1=cbind(dnds1[,3],dnds2[,3],dnds3[,3],dnds4[,3]))
abline(h=1, col="cadetblue")
legend(x=min(b), y=2.5, pch=15, col=c("grey30","grey60"), legend=c("Missense","Nonsense"), box.col="white")
dev.copy(pdf, "fig_dNdS_putative_antigenic_regions.pdf", width=6, height=4)
dev.off(); dev.off()

# Testing the code
if (0) { # This should be highly enriched in missense mutations
    mis = dndsout$annotmuts[dndsout$annotmuts$impact=="Missense", ]
    gr_mis = GenomicRanges::GRanges(mis$chr, IRanges(mis$pos, mis$pos))
    outtest = globaldnds_gr(mutations = mutations, gr = gr_mis)
}

# Zapata's substitution model
if (0) {
    data("submod_192r_3w", package = "dndscv")
    s = substmodel
    s[,1] = "t"; s[,2] = "t*wmis"; s[,3] = "t*wnon"; s[,4] = "t*wspl"
    ntcomp = c(A="T",C="G",G="C",T="A")
    cpg = (substr(rownames(s),2,3)=="CG" | substr(rownames(s),1,2)=="CG")
    for (h in 1:nrow(s)) {
        basefrom = substr(rownames(s)[h],2,2)
        baseto = substr(rownames(s)[h],6,6)
        if (basefrom %in% c("A","G")) {
            basefrom = ntcomp[basefrom]
            baseto = ntcomp[baseto]
        }
        if (cpg[h]) { # CpG site
            s[h,] = paste(s[h,],"CpG>X",sep="*")
        } else { # No CpG
            s[h,] = paste(s[h,],sprintf("%s>%s",basefrom,baseto),sep="*")
        }
    }
    s[nrow(s),] = c("t","t*wmis","t*wnon","t*spl")
    # Running dN/dS
    gr_antigens = intersect(gr_iedb, gr_hla6)
    outz = globaldnds_gr(mutations = mutations, gr = gr_antigens, gene_list = non_driver_genes, sm = s) # If L is not provided it will generate it (slow)
    dndsz = outz$globaldnds[c("wmis","wnon"), -1]
}