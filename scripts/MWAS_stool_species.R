setwd("~/Documents/BGI Work/GC/ms/CRM/Revision2/speTest/scripts/")
library(data.table)
library(Maaslin2)
library(ANCOMBC)
library(mia)
library(ggvenn)

# take first column as the rownames
dt.to.matrix <- function(x) {
  x <- as.data.frame(x)
  x <- x[!is.na(x[, 1]), ]
  rownames(x) <- x[,1]
  x <- as.matrix(x[,-1])
  x
}

stool <- fread("../data/Cohort1_stool_dt.tsv")
qc_summary <- fread("../data/Cohort1_reads_qc_summary.txt")

# summary statistics
spe_summ <- stool[, .(prevalence=round(sum(abundance>0)/.N*100, digits=2)), by=.(Taxon)]
spe_summ_group <- stool[abundance>0, .(meanAbd=mean(abundance), N=.N), by=.(Taxon, group)]
spe_summ_group_N <- dcast(spe_summ_group, Taxon~group, value.var = "N")
spe_summ_group_meanAbd <- dcast(spe_summ_group, Taxon~group, value.var = "meanAbd")
spe_summ_group_dt <- merge(spe_summ_group_N[, .(Taxon, N_GC=cancer, N_ChG=gastritis)], 
                           spe_summ_group_meanAbd[, .(Taxon, RA_GC=cancer, RA_ChG=gastritis)])

# filter out species present in <10% individuals of cohort 1
spe_kept <- spe_summ[!(prevalence<10)] # 485 species left
stool1 <- stool[Taxon %in% spe_kept$Taxon]

# wilcoxon test based on CLR transform relative abundance
wilcoxon_res <- stool1[, Filter(Negate(is.null), kruskal.test(as.numeric(clr_abundance)~group)), by="Taxon"]
wilcoxon_res[, p.value.fdr:=p.adjust(p.value, method = "fdr")]

# MaAslin2 analysis, no adjustment for any covariate
clr_mat <- dt.to.matrix(dcast(stool, ID~Taxon, value.var = "clr_abundance"))
phe <- unique(stool[, .(ID, group, Age, Sex, BMI)])
phe <- dt.to.matrix(phe)
phe <- as.data.frame(phe)
phe$group <- factor(phe$group, levels = c("gastritis", "cancer"))
maaslin <- Maaslin2(input_data     = as.data.frame(clr_mat), 
                    input_metadata = phe, 
                    min_prevalence = 0.1,
                    normalization  = "NONE",
                    transform      = "NONE",
                    output         = "masslin", 
                    fixed_effects  = c("group"),
                    correction     = "BH",
                    plot_heatmap   = FALSE,
                    plot_scatter   = FALSE,
                    max_pngs       = 100
)
maaslin_res <- as.data.table(maaslin$results)

# ANCOM-BC analysis
# taxon names in row, samples in column
abd_mat <- t(dt.to.matrix(dcast(stool, ID~Taxon, value.var = "abundance")))
qc_summary <- qc_summary[match(colnames(abd_mat), Library_ID)]
# put back the sequencing depth (effective library size) information for each sample
# in my application, the sequencing depth did not change the primary test results, 
# but changed the sensitivity analysis for pseudo-count addition.
abd_mat_count <- sweep(abd_mat, 2, qc_summary$Raw_reads_count/100, `*`)

# taxonomic class information is for place-holder only
tax_tab <- data.table(taxon=rownames(abd_mat), Kindom="Bacteria", 
                      Phylum=rownames(abd_mat), Class=rownames(abd_mat), 
                      Order=rownames(abd_mat), Family=rownames(abd_mat),
                      Genus=rownames(abd_mat), Species=rownames(abd_mat))
tax_tab <- dt.to.matrix(tax_tab)

# sample metadata
smd <- unique(stool1[, .(ID, group, Age, Sex, BMI)])
smd <- smd[, .(sampleID=ID, subjectID=ID, group, Age, Sex, BMI)]
# follow the same order as abd_mat
smd <- smd[match(colnames(abd_mat), sampleID)]
smd$group <- factor(smd$group, levels = c("gastritis", "cancer"))

tse <- TreeSummarizedExperiment(assays=SimpleList(counts = abd_mat_count), 
                                colData=smd, 
                                rowData=DataFrame(tax_tab))

out <- ancombc2(data = tse, assay_name = "counts", tax_level = "Species",
                fix_formula = "group", p_adj_method = "BH", pseudo_sens = TRUE,
                prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                group = NULL, struc_zero = FALSE, neg_lb = FALSE,
                alpha = 0.05, n_cl = 1, verbose = TRUE)
ancombc_res <- as.data.table(out$res)

# combine wilcoxon_res, masslin_res
res <- merge(wilcoxon_res[, .(Taxon, P_wilcoxon=p.value, FDR_wilcoxon=p.value.fdr)],
             maaslin_res[, .(Taxon=feature, Beta=coef, SE=stderr, P_maaslin=pval, FDR_maaslin=qval, 
                             FC=2**(coef), FC_lower=2**(coef-stderr), FC_upper=2**(coef+stderr))])
res <- merge(spe_summ_group_dt, res, all.y = TRUE)
res <- merge(res, ancombc_res[, .(Taxon=taxon, Beta_ancombc=lfc_groupcancer, SE_ancombc=se_groupcancer, 
                                  P_ancombc=p_groupcancer, FDR_ancombc=q_groupcancer, 
                                  FC_ancombc=exp(lfc_groupcancer),
                                  FC_lower_ancombc=exp(lfc_groupcancer-se_groupcancer), 
                                  FC_upper_ancombc=exp(lfc_groupcancer+se_groupcancer),
                                  passed_ss_groupcancer)])

# Venn plot
venn_data <- data.frame(
  Wilcoxon = res$FDR_wilcoxon<0.05,
  MaAslin2 = res$FDR_maaslin<0.05,
  "ANCOM-BC" = res$FDR_ancombc<0.05 & res$passed_ss_groupcancer
)
p_venn <- ggvenn(venn_data, 
                 fill_color = c("grey", "lightblue", "darkblue"), # alpha(brewer.pal(11,"BrBG")[c(3, 2, 10)], 1),
                 stroke_color = "black",  text_color = "black",
                 stroke_size = 0.5, text_size = 8, show_percentage = F
) #+ theme(plot.margin = margin(10, 20, 10, 10, "pt"))
ggsave("../results/MWAS_stool_species_res_pie.pdf", p_venn, width = 2.6*1.2, height = 2.8*1.6, units = c("in"), limitsize = FALSE)


### validation in Cohort2
stool_v <- fread("../data/Cohort2_stool_dt.tsv")
o2g <- fread("../data/Schmidt_2019_oral_gut_transmitter.tsv")

# only worked on differential species identified in step1 (Cohort1)
stool_v <- stool_v[Taxon %in% res[FDR_wilcoxon<0.05 & FDR_maaslin<0.05 & FDR_ancombc<0.05 & passed_ss_groupcancer]$Taxon]

# summary statistics
#spe_summ <- stool1[, .(prevalence=round(sum(abundance>0)/.N*100, digits=2)), by=.(Taxon)]
spe_summ_group_v <- stool_v[abundance>0, .(meanAbd=mean(abundance), N=.N), by=.(Taxon, group)]
spe_summ_group_N_v <- dcast(spe_summ_group_v, Taxon~group, value.var = "N")
spe_summ_group_meanAbd_v <- dcast(spe_summ_group_v, Taxon~group, value.var = "meanAbd")
spe_summ_group_dt_v <- merge(spe_summ_group_N_v[, .(Taxon, N_GC_V=cancer, N_ChG_V=gastritis)], 
                           spe_summ_group_meanAbd_v[, .(Taxon, RA_GC_V=cancer, RA_ChG_V=gastritis)])

# wilcoxon test based on CLR transform relative abundance
wilcoxon_res_v <- stool_v[, Filter(Negate(is.null), kruskal.test(as.numeric(clr_abundance)~group)), by="Taxon"]
wilcoxon_res_v[, p.value.fdr:=p.adjust(p.value, method = "fdr")]

# MaAslin2 analysis, no adjustment for any covariate
clr_mat_v <- dt.to.matrix(dcast(stool_v, ID~Taxon, value.var = "clr_abundance"))
phe_v <- unique(stool_v[, .(ID, group, Age, Sex, BMI)])
phe_v <- dt.to.matrix(phe_v)
phe_v <- as.data.frame(phe_v)
phe_v$group <- factor(phe_v$group, levels = c("gastritis", "cancer"))
maaslin_v <- Maaslin2(input_data     = as.data.frame(clr_mat_v), 
                    input_metadata = phe_v, 
                    min_prevalence = 0,
                    normalization  = "NONE",
                    transform      = "NONE",
                    output         = "masslin", 
                    fixed_effects  = c("group"),
                    correction     = "BH",
                    plot_heatmap   = FALSE,
                    plot_scatter   = FALSE,
                    max_pngs       = 100
)
maaslin_res_v <- as.data.table(maaslin_v$results)

# combine wilcoxon_res, masslin_res
res <- merge(res, spe_summ_group_dt_v[, .(Taxon, N_GC_V, N_ChG_V, RA_GC_V, RA_ChG_V)], all.x = TRUE)

res <- merge(res, wilcoxon_res_v[, .(Taxon, P_wilcoxon_V=p.value)], all.x = TRUE)

res <- merge(res, 
             maaslin_res_v[, .(Taxon=feature, Beta_V=coef, SE_V=stderr, P_maaslin_V=pval, 
                               FC_V=2**(coef), FC_lower_V=2**(coef-stderr), FC_upper_V=2**(coef+stderr))],
             all.x = TRUE)

res[, Concordance:=FALSE]
res[Beta*Beta_V>0, Concordance:=TRUE]
res[is.na(Beta_V), Concordance:=NA]
res <- merge(o2g, res, all.y = TRUE)
fwrite(res[order(P_wilcoxon)], "../results/MWAS_stool_species_res.tsv", quote = F, sep = "\t", na = "-")


