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

pathway <- fread("../data/Cohort2_saliva_pathway_dt.tsv")
qc_summary <- fread("../data/Cohort2_saliva_reads_qc_summary.txt")

# summary statistics
pathway_summ <- pathway[, .(prevalence=round(sum(abundance>0)/.N*100, digits=2)), by=.(pathway)]
pathway_summ_group <- pathway[abundance>0, .(meanAbd=mean(abundance), N=.N), by=.(pathway, group)]
pathway_summ_group_N <- dcast(pathway_summ_group, pathway~group, value.var = "N")
pathway_summ_group_meanAbd <- dcast(pathway_summ_group, pathway~group, value.var = "meanAbd")
pathway_summ_group_dt <- merge(pathway_summ_group_N[, .(pathway, N_GC=GC, N_ChG=ChG)], 
                               pathway_summ_group_meanAbd[, .(pathway, RA_GC=GC, RA_ChG=ChG)])

# remove pathway that only detectable in one group
pathway <- pathway[!(pathway %in% pathway_summ_group_dt[is.na(N_GC) | is.na(N_ChG)]$pathway)]

# filter out pathway present in <10% individuals of cohort 1
pathway_kept <- pathway_summ[!(prevalence<10)] 
pathway1 <- pathway[pathway %in% pathway_kept$pathway] # 388 pathways left

# wilcoxon test based on CLR transform relative abundance
wilcoxon_res <- pathway1[, Filter(Negate(is.null), kruskal.test(as.numeric(clr_abundance)~group)), by="pathway"]
wilcoxon_res[, p.value.fdr:=p.adjust(p.value, method = "fdr")]

# MaAslin2 analysis, no adjustment for any covariate
# rename long pathway name to a simple alias and put it back in the output
pathway_summ[, og_name:=pathway]
pathway_summ[, alias_name:=paste0("P", 1:nrow(pathway_summ))]
pathway <- merge(pathway, pathway_summ[, .(pathway, alias_name)])
pathway[, pathway:=alias_name]

clr_mat <- dt.to.matrix(dcast(pathway, sampleID~pathway, value.var = "clr_abundance"))
phe <- unique(pathway[, .(sampleID, group, Age, Sex, BMI)])
phe <- dt.to.matrix(phe)
phe <- as.data.frame(phe)
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
maaslin_res <- merge(pathway_summ[, .(pathway=alias_name, og_name)], maaslin_res, 
                     by.x="pathway", by.y="feature", all.y=TRUE)
maaslin_res[, pathway:=og_name]

# ANCOM-BC analysis
# pathway names in row, samples in column
abd_mat <- t(dt.to.matrix(dcast(pathway, sampleID~pathway, value.var = "abundance")))
qc_summary <- qc_summary[match(colnames(abd_mat), Library_ID)]
# put back the sequencing depth (effective library size) information for each sample
# in my application, the sequencing depth did not change the primary test results, 
# but changed the sensitivity analysis for pseudo-count addition.
abd_mat_count <- sweep(abd_mat, 2, qc_summary$Raw_reads_count/100, `*`)

# level information is for place-holder only
tax_tab <- data.table(pathway=rownames(abd_mat), Kindom="Bacteria", 
                      Phylum=rownames(abd_mat), Class=rownames(abd_mat), 
                      Order=rownames(abd_mat), Family=rownames(abd_mat),
                      Genus=rownames(abd_mat), Species=rownames(abd_mat))
tax_tab <- dt.to.matrix(tax_tab)

# sample metadata
smd <- unique(pathway1[, .(sampleID, group, Age, Sex, BMI)])
smd <- smd[, .(sampleID, subjectID=sampleID, group, Age, Sex, BMI)]
smd <- smd[match(colnames(abd_mat), sampleID)] # follow the same order as abd_mat

tse <- TreeSummarizedExperiment(assays=SimpleList(counts = abd_mat_count), 
                                colData=smd, 
                                rowData=DataFrame(tax_tab))

out <- ancombc2(data = tse, assay_name = "counts", tax_level = "Species",
                fix_formula = "group", p_adj_method = "BH", pseudo_sens = TRUE,
                prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                group = NULL, struc_zero = FALSE, neg_lb = FALSE,
                alpha = 0.05, n_cl = 1, verbose = TRUE)
ancombc_res <- as.data.table(out$res)
ancombc_res <- merge(pathway_summ[, .(pathway=alias_name, og_name)], ancombc_res, 
                     by.x="pathway", by.y="taxon", all.y=TRUE)
ancombc_res[, pathway:=og_name]

# combine wilcoxon_res, masslin_res
res <- merge(wilcoxon_res[, .(pathway, P_wilcoxon=p.value, FDR_wilcoxon=p.value.fdr)],
             maaslin_res[, .(pathway, Beta=coef, SE=stderr, P_maaslin=pval, FDR_maaslin=qval, 
                             FC=2**(coef), FC_lower=2**(coef-stderr), FC_upper=2**(coef+stderr))])
res <- merge(pathway_summ_group_dt, res, all.y = TRUE)
res <- merge(res, ancombc_res[, .(pathway, Beta_ancombc=lfc_groupGC, SE_ancombc=se_groupGC, 
                                  P_ancombc=p_groupGC, FDR_ancombc=q_groupGC, 
                                  FC_ancombc=exp(lfc_groupGC),
                                  FC_lower_ancombc=exp(lfc_groupGC-se_groupGC), 
                                  FC_upper_ancombc=exp(lfc_groupGC+se_groupGC),
                                  passed_ss_groupGC)])

# Venn plot
venn_data <- data.frame(
  Wilcoxon = res$FDR_wilcoxon<0.05,
  MaAslin2 = res$FDR_maaslin<0.05,
  "ANCOM-BC" = res$FDR_ancombc<0.05 & res$passed_ss_groupGC
)
p_venn <- ggvenn(venn_data, 
                 fill_color = c("grey", "lightblue", "darkblue"), # alpha(brewer.pal(11,"BrBG")[c(3, 2, 10)], 1),
                 stroke_color = "black",  text_color = "black",
                 stroke_size = 0.5, text_size = 8, show_percentage = F
) #+ theme(plot.margin = margin(10, 20, 10, 10, "pt"))
ggsave("../results/MWAS_oral_pathway_res_pie.pdf", p_venn, width = 2.6*1.2, height = 2.8*1.6, units = c("in"), limitsize = FALSE)

### validation in Harbin tongue cohort
pathway_v <- fread("../data/Harbin_tongue_pathway_dt.tsv")
# only worked on 76 differential pathways identified in Cohort2 saliva metagenomes
pathway_v <- pathway_v[pathway %in% res[FDR_wilcoxon<0.05 & FDR_maaslin<0.05 & FDR_ancombc<0.05 & passed_ss_groupGC]$pathway]
pathwayV_summ_group <- pathway_v[abundance>0, .(meanAbd=mean(abundance), N=.N), by=.(pathway, group)]
pathwayV_summ_group_N <- dcast(pathwayV_summ_group, pathway~group, value.var = "N")
pathwayV_summ_group_meanAbd <- dcast(pathwayV_summ_group, pathway~group, value.var = "meanAbd")
pathwayV_summ_group_dt <- merge(pathwayV_summ_group_N[, .(pathway, N_GC_V=GC, N_ChG_V=ChG)], 
                                pathwayV_summ_group_meanAbd[, .(pathway, RA_GC_V=GC, RA_ChG_V=ChG)])

# wilcoxon test based on CLR transform relative abundance
pathwayV_wilcoxon_res <- pathway_v[, Filter(Negate(is.null), kruskal.test(as.numeric(clr_abundance)~group)), by="pathway"]
pathwayV_wilcoxon_res[, p.value.fdr:=p.adjust(p.value, method = "fdr")]

# MaAslin2 analysis, no adjustment for any covariate
# rename long pathway name to a simple alias and put it back in the output
pathwayV_summ_group_dt[, og_name:=pathway]
pathwayV_summ_group_dt[, alias_name:=paste0("P", 1:nrow(pathwayV_summ_group_dt))]
pathway_v <- merge(pathway_v, pathwayV_summ_group_dt[, .(pathway, alias_name)])
pathway_v[, pathway:=alias_name]

pathwayV_clr_mat <- dt.to.matrix(dcast(pathway_v, sampleID~pathway, value.var = "clr_abundance"))
pathwayV_phe <- unique(pathway_v[, .(sampleID, group, age, sex)])
pathwayV_phe <- dt.to.matrix(pathwayV_phe)
pathwayV_phe <- as.data.frame(pathwayV_phe)
pathwayV_maaslin <- Maaslin2(input_data = as.data.frame(pathwayV_clr_mat), 
                             input_metadata = pathwayV_phe, 
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
pathwayV_maaslin_res <- as.data.table(pathwayV_maaslin$results)
pathwayV_maaslin_res <- merge(pathwayV_summ_group_dt[, .(pathway=alias_name, og_name)], pathwayV_maaslin_res, 
                              by.x="pathway", by.y="feature", all.y=TRUE)
pathwayV_maaslin_res[, pathway:=og_name]

# combine discovery and validation results 
res <- merge(res, pathwayV_summ_group_dt[, .(pathway, N_GC_V, N_ChG_V, RA_GC_V, RA_ChG_V)], all.x = TRUE)

res <- merge(res, pathwayV_wilcoxon_res[, .(pathway, P_wilcoxon_V=p.value)], all.x = TRUE)

res <- merge(res, 
             pathwayV_maaslin_res[, .(pathway, Beta_V=coef, SE_V=stderr, P_maaslin_V=pval,
                                      FC_V=2**(coef), FC_lower_V=2**(coef-stderr), FC_upper_V=2**(coef+stderr))],
             all.x = TRUE)

res[, Concordance:=FALSE]
res[Beta*Beta_V>0, Concordance:=TRUE]
res[is.na(Beta_V), Concordance:=NA]
fwrite(res[order(P_wilcoxon)], "../results/MWAS_oral_pathway_res.tsv", quote = F, sep = "\t", na = "-")

