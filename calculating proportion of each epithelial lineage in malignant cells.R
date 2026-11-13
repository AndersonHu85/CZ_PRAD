
library(dplyr)
library(ggplot2)
library(Seurat)
library(ggthemes)

marker_list2 <- list(pan_dLum=c("MMP7","RARRES1","PIGR","SLPI","WFDC2","CRABP2","TNFAIP2","CP","KRT7","PPP1R1B","RHOV","SCUBE2","CLDN7","TSPAN3",
              "RARRES3","SELENBP1","PNPLA2","LY6E","LXN","C19orf33","CRACR2B","VTCN1","QSOX1","TMEM45B"),
  pan_basal=c("KRT5","KRT15","DST","TIMP3","TRIM29","ITGA6","GJA1","TP63","ISYNA1","IGFBP4","DSC3","CSTA","ITGA2","SPARCL1"),
  pan_Lum=c("ACPP","KLK3","MSMB","KLK2","KLK4","SLC45A3","NKX3-1","SORD","AZGP1","DHRS7","CPE","DBI","SEC11C","PMEPA1","SLC30A4","RAB27A","LCP1","FOLH1","PAK1IP1","DPP4","RCAN3","ANPEP","NIPAL3","SPOCK1",
             "ALOX15B","RAB3B","SCD","PTN","SPDEF","HIST1H2AC","ZNF761","SYT7","MBOAT2","NAAA","GREB1","AFF3","ANO7","CLGN","HGD","PLPP1","CPLX3","ACLY","SOAT1","LIFR","CNTNAP2","FASN","BMPR1B","EPHX2","CPNE4",
             "CREB3L4","CWH43","PCAT4"),
  club=c("NT5E","SCGB3A1","SRD5A2","SCGB1A1"),
  KRT4_Hillock=c("OLFM4","KRT13","SNCG","KRT4","MUC4","UPK1A","ALOX5","ERP27","CEACAM6","UPK2","MUC15","LPCAT4","VGLL1","MUC1"),
  LTF_dLum=c("LTF","TGM2","C1R","MGP","C1S","CP","MUC3A","PDE4B","C3","MMP7","NPC2","KLK10","PAPSS1",'CHI3L2','CLDN10','DNAJC12','SERPINA1','TF'),
  basal_KRT13=c("KRT13","SERPINB5",'LYPD3','OLFM4','PLAT','TMPRSS4','MUC4','CLCA4','ZNF750','RAPGEFL1','DUOX2','CLEC2B',"FGFR3"),
  basal_mature=c("SFRP1", "EFEMP1","COL17A1", "TNC",  "LTBP2", "SPARC", "PDPN","FLRT3", "JAG2", "PDGFC", "HTRA1", "MXRA5","TMEM237" ),
  basal_ion=c( "FOXI1", "ATP6V1C2", "RCAN2",  "KIT", "ATP6V0D2",  "ATP6V0A4", "CELF2", "CA12", "ATP6V1B1", "ADAMTSL3", "SFRP4",  "KIAA1210", "RARRES2", "SIX2"),
  LTF_Lum=c("CPB1","TMEM178A","OLFML3","ASRGL1","GUSB","PTP4A3","AMOT","SLC9A9","ACSM3","COL22A1","CRAT"),
  type1_lum=c("CD177", "LSAMP", "NEFH", "CCK", "CD38","SLC2A5", "ZNF827", "FADS2","PRDM8", "FADS1","AOC1","SV2B"),
  type2_lum=c("RLN1","TMEFF2","PCGEM1","OR51E2","SPON2","PRIM2","TFF3","TRGC1","TRPM8","PRIM2","RLN2","OR51E1"),
  type3_lum=c( "SFTPA2","MIPEP","TRPM4","HSD17B6","SCIN","REPS2","SLC38A11","GCSAML","CST1","CST4","MIPEP","CST2","TGM3"),
  type4_lum=c("TGM4","PGC","CDK19","GATA2","CLSTN2","CYP1B1","PRUNE2","PLA2G4D","NXPE4","SLC39A2","SYT16","GALNTL6","DNAH8"),
  Ciliated=c("C1orf194","CCDC170","FOXJ1","ZMYND10","TPPP3","ARMC3","CFAP126","ZBBX","C9orf24","CAPSL","CDHR3","FAM81B","RSPH4A","CFAP73","CFAP45","DYDC2","LRRC46","ROPN1L","CCDC114","WDR63",
             "ENKUR","C6orf118","PRR29","HYDIN","MAP3K19","CDHR4","SPAG17","DRC1","CFAP54","CFAP77","FAM216B","ERICH3","MORN5","WDR38","CCDC39","C20orf85","CFAP57","ANKRD66","SNTN",
             "TEKT1","CFAP52","LRRC10B","DNAI1","CCDC65","TUBA4B","CFAP157","VWA3A","TTC29","C7orf57"),
  NE=c("SCG2","SCGN","CHGA","CHGB","SCG3","ASCL1","INSM1","DDC", "VGF", "CELF3","APLP1","TPH1","SLC18A1","SYT4","CACNA1A","NOL4","SCG5","KCNH2","MAOB","NFASC","GRP",
       "RIMBP2","CPLX2","CALCA","RBP4","ANK2","SNAP25","PPFIA2","TOX", "ELAVL3","PHGR1","TMEM59L","TAGLN3","RUNDC3A","PCSK2","SYT5","PEX5L","GCG", "CDK5R2","NKX2-2","ELAVL4",
       "VWA5B2","EFR3B","PTPRN","PCSK1")
)
merge_metacell <- AddModuleScore(merge_metacell,marker_list2)
colnames(merge_metacell@meta.data)[12:27] <- names(marker_list2)
metadata <- merge_metacell@meta.data
major_cols <- c("pan_Lum", "pan_basal", "pan_dLum")
metadata$major_class <- major_cols[
  max.col(metadata[, major_cols], ties.method = "first")
]
subclass_cols <- list(
  pan_basal = c("basal_KRT13", "basal_mature", "basal_ion"),
  pan_dLum = c("club", "KRT4_Hillock", "LTF_dLum"),
  pan_Lum = c("type1_lum", "type2_lum", "type3_lum", "type4_lum", "LTF_Lum")
)
metadata$identity[metadata$NE>0.25] <- "NE"
metadata$identity[metadata$Ciliated>0.15] <- "Ciliated"
metadata <- metadata %>%
  group_by(major_class) %>%
  group_modify(function(df, key) {
    current_sub_cols <- subclass_cols[[key$major_class]]
    df$identity <- current_sub_cols[
      max.col(df[, current_sub_cols], ties.method = "first")
    ]
    return(df)
  }) %>%
  ungroup()
merge_metacell@meta.data <- metadata
identity_sample_counts <- metadata %>%
  group_by(sample, group4, identity) %>%
  summarise(total_cells = sum(size), .groups = "drop")
sample_total <- metadata %>%
  group_by(sample) %>%
  summarise(total_sample_cells = sum(size), .groups = "drop")


identity_proportion <- identity_sample_counts %>%
  left_join(sample_total, by = "sample") %>%
  mutate(proportion = total_cells / total_sample_cells) %>%
  select(sample, group4, identity, total_cells, total_sample_cells, proportion)
identity_proportion <- identity_proportion %>%
  arrange(sample, desc(proportion))
type2_lum_ratio <- identity_proportion %>%
  group_by(sample) %>%
  summarise(
    type2_lum_prop = ifelse(
      "type2_lum" %in% identity,
      proportion[identity == "type2_lum"],
      0
    ),
    .groups = "drop"
  )
identity_proportion_sorted <- identity_proportion %>%
  left_join(type2_lum_ratio, by = "sample") %>%
  arrange(desc(type2_lum_prop), sample, group4, desc(proportion)) %>%
  select(-type2_lum_prop) 

p1 <- ggplot(identity_proportion_sorted,aes(factor(sample,levels = unique(identity_proportion_sorted$sample)),weight=total_cells,fill=factor(identity,levels = names(marker_list2))))+geom_bar(position="fill")+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  scale_fill_manual(values = celltype)
p1 + ggthemes::theme_base()+facet_grid(~factor(group4,levels = c("boundary","localized", "Locally advanced PCa",
                                                                 "mHSPC", "post_ADT", "CRPC","NEPC","PBCC","LN_M","BM_M")),scales = "free", space = "free")+theme(
                                                                    panel.spacing.x = unit(1, "mm")
                                                                 )
ggsave("barplot_lineage_v2.pdf", width = 20, height = 6)


shannon_data <- identity_proportion %>%
  filter(proportion > 0, !is.na(proportion)) %>%
  select(sample, group4, identity, proportion)

sample_shannon <- shannon_data %>%
  group_by(sample, group4) %>%
  summarise(
    shannon_index = -sum(proportion * log(proportion)),
    celltype_count = n_distinct(identity),
    evenness = shannon_index / log(n_distinct(identity)),
    .groups = "drop"
  )

sample_shannon_sorted <- sample_shannon %>%
  arrange(group4, desc(shannon_index))

ggplot(sample_shannon_sorted, aes(
  x = group4,
  y = shannon_index,
  fill = group4
)) +
  geom_boxplot()+geom_jitter()
