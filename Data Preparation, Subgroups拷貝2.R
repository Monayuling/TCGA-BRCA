# 參考教學影片：https://www.youtube.com/watch?v=wPzeea1Do18
# Raw Subgroups: Tumor, Normal, Meta
# Young: <= 50 years old, Old: >50 years old
# Menopause: <=45 years old, non-menopause: >= 55 years old
# N0, N1; R0, Rx

library(DESeq2)
library(apeglm)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
#install.packages("openxlsx")
library(openxlsx)
library(dplyr)
library(survival)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("pcaExplorer")
#browseVignettes("pcaExplorer")
library("pcaExplorer")
#install.packages("markdown")
library(markdown)
library(clusterProfiler)
library(AnnotationDbi)
library(enrichplot)
library(org.Hs.eg.db)
library(topGO)
library(Rgraphviz)


# ================ 1.前置作業: 準備expression matrix 及 clinical matrix
# ==== 讀入expression matrix，並新增一欄sample種類
raw_exp <- read.csv('/Users/mona/Desktop/Subgroup/Rawcount/luminal_subgroup_exp_matrix_nogroup.csv', 
                    header = T, row.names = 1)
samples_expression <- t(raw_exp)
# 把sample type 讀出來
type <- lapply(rownames(samples_expression), function(rowname) {
  substr(rowname, 14, 16)
})
type_matrix <- matrix(unlist(type), ncol = 1)
# 看到底有幾種
type_vector <- unlist(type)
unique_types <- unique(type_vector)
cat("Unique types:", toString(unique_types), "\n")
# 01: Primary Solid Tumor; 11: Solid Tissue Normal; 06: Metastatic Site
# 把Sample資訊加回去
samples_expression <- cbind(samples_expression, type_matrix)
# ==== 建立新的expression matrix # 這些沒有加上metadata所以不能直接用
# Create a matrix for "tumor" samples (values with "01A" or "01B")
tumor_matrix <- samples_expression[grep("01A|01B", type), ]
tumor_barcode <- rownames(tumor_matrix)
# Create a matrix for "normal" samples (values with "11A" or "11B")
normal_matrix <- samples_expression[grep("11A|11B", type), ]
normal_barcode <- rownames(normal_matrix)
# Create a matrix for "Metastatic" samples (values with "06A")
metastatic_matrix <- samples_expression[grep("06A", type), ]
meta_barcode <- rownames(metastatic_matrix)
meta_patients <- substr(meta_barcode, 1, 12)
# ==== 加入臨床資料
# ====== 讀入之前已經分類好的clinical information, 並加入新資訊：
complete_info <- read.csv('/Users/mona/Desktop/Subgroup/Rawcount/true_age_info.csv'
                          , header = T, row.names = 1)
# ======== 新增Race, Stage, N stage, M, New Tumor event
additional <- read.csv('/Users/mona/Desktop/Subgroup/Annotated_metadata.csv'
                       , header = T)
rownames(additional) <- additional$bcr_patient_barcode
additional$bcr_patient_barcode <- NULL
# ========= 加到complete_info
# Extract the first 12 characters of row names in complete_info
partial_match <- substr(rownames(complete_info), 1, 12)
# Find common partial row names
common_partial_rows <- partial_match %in% substr(rownames(additional), 1, 12)
# Update 'Race' column in 'complete_info' for common partial row names
complete_info$Race <- NA
complete_info$Race[common_partial_rows] <- additional$race[match(partial_match[common_partial_rows], substr(rownames(additional), 1, 12))]
# 1.White 2.Black 3.Asian 4.Hispatic
# Update other columns
complete_info$Stage <- NA
complete_info$Stage[common_partial_rows] <- additional$Stage[match(partial_match[common_partial_rows], substr(rownames(additional), 1, 12))]
# Stage 1,2,3,4
complete_info$Node <- NA
complete_info$Node[common_partial_rows] <- additional$N.Stage[match(partial_match[common_partial_rows], substr(rownames(additional), 1, 12))]
# 0: N0 1:Nx
complete_info$Metastasis <- NA
complete_info$Metastasis[common_partial_rows] <- additional$pathologic_M[match(partial_match[common_partial_rows], substr(rownames(additional), 1, 12))]
# No: No metastasis Yes: Metastasis 
complete_info$Event <- NA
complete_info$Event[common_partial_rows] <- additional$new_tumor_event_after_initial_treatment[match(partial_match[common_partial_rows], substr(rownames(additional), 1, 12))]
# New tumor event
#==目前raw_exp 的colname 有 '.', info 的 rowname則是有'-'
#Extract row names from the info matrix
complete_row_names <- rownames(complete_info)
# Compare and replace row names
for (i in 1:length(complete_row_names)) {
  # Replace "-" with "."
  c_new_row_name <- gsub("\\-", ".", complete_row_names[i])
  #print(new_row_name)
  if (c_new_row_name == colnames(raw_exp)[i]) {
    # Substitute row name in the info matrix
    rownames(complete_info)[i] <- c_new_row_name
  }
}
# ========== 新增sample information (現在rowname跟barcode一樣)
complete_info$Sample <- NA
# Create new column 'Sample' in complete_info
complete_info$Sample <- ifelse(rownames(complete_info) %in% tumor_barcode, "Tumor",
                               ifelse(rownames(complete_info) %in% normal_barcode, "Normal",
                                      ifelse(rownames(complete_info) %in% meta_barcode, "Meta", NA)))
# 先把raw_exp 和 complete_info merge 在一起
t_raw <- as.data.frame(t(raw_exp))
complete_info$ID <- rownames(complete_info)
t_raw$ID <- rownames(t_raw)
large <- merge(complete_info, t_raw, by = 'ID')
# 按照年齡排序
large$age_true <- as.integer(large$age_true)
large <- large[order(large$age_true),]
row.names(large) <- large[,1]
large <- large[,-1]
write.xlsx(large, "Metadata_and_expression.xlsx", rowNames = TRUE)

# Stage Subgroups
get_stage_barcode <- function(data, barcodes, stage) {
  row_names <- rownames(data)
  selected_barcodes <- row_names[row_names %in% barcodes & !is.na(data$Stage) & data$Stage == stage]
  return(selected_barcodes)
}
# Usage for tumor
S1_tumor_barcode <- get_stage_barcode(large, tumor_barcode, 1)
S2_tumor_barcode <- get_stage_barcode(large, tumor_barcode, 2)
S3_tumor_barcode <- get_stage_barcode(large, tumor_barcode, 3)
S4_tumor_barcode <- get_stage_barcode(large, tumor_barcode, 4)
# Usage for normal
S1_normal_barcode <- get_stage_barcode(large, normal_barcode, 1)
S2_normal_barcode <- get_stage_barcode(large, normal_barcode, 2)
S3_normal_barcode <- get_stage_barcode(large, normal_barcode, 3)
S4_normal_barcode <- get_stage_barcode(large, normal_barcode, 4)
# N0, Nx subgroups
get_node_barcode <- function(data, barcodes, node) {
  row_names <- rownames(data)
  selected_barcodes <- row_names[row_names %in% barcodes & !is.na(data$Node) & data$Node == node]
  return(selected_barcodes)
}
# Tumor
N0_tumor_barcode <- get_node_barcode(large, tumor_barcode, 0)
N1_tumor_barcode <- get_node_barcode(large, tumor_barcode, 1)
N2_tumor_barcode <- get_node_barcode(large, tumor_barcode, 2)
N3_tumor_barcode <- get_node_barcode(large, tumor_barcode, 3)
Nx_tumor_barcode <-c(N1_tumor_barcode, N2_tumor_barcode, N3_tumor_barcode)
# N0 patients with metastasis:
N0_patients <- substr(N0_tumor_barcode, 1, 12)
N0_Mx_patients <- intersect(N0_patients, meta_patients) # None
# Normal
N0_normal_barcode <- get_node_barcode(large, normal_barcode, 0)
N1_normal_barcode <- get_node_barcode(large, normal_barcode, 1)
N2_normal_barcode <- get_node_barcode(large, normal_barcode, 2)
N3_normal_barcode <- get_node_barcode(large, normal_barcode, 3)
Nx_normal_barcode <- c(N1_normal_barcode, N2_normal_barcode, N3_normal_barcode)
# Tumor events (R0, Rx)
get_event_barcode <- function(data, barcodes, event) {
  row_names <- rownames(data)
  selected_barcodes <- row_names[row_names %in% barcodes & !is.na(data$Event) & data$Event == event]
  return(selected_barcodes)
}
R0_tumor_barcode <- get_event_barcode(large, tumor_barcode, 'NO')
Rx_tumor_barcode <- get_event_barcode(large, tumor_barcode, 'YES')
R0_normal_barcode <- get_event_barcode(large, normal_barcode, 'NO')
Rx_normal_barcode <- get_event_barcode(large, normal_barcode, 'YES')

# Find Menopause barcodes: <= 45 years: No, >= 55 years: Yes,
nonmeno_barcodes <- rownames(large)[which(large$age_true <= 45)]
menopause_barcodes <- rownames(large)[which(large$age_true >= 55)]

large <- cbind(large[, 1], Menopause_True = NA, large[, -1])
large$Menopause_True[which(large$age_true <= 45)] <- "No"  # Assign "Yes" to rows <= 45
large$Menopause_True[which(large$age_true >= 55)] <- "Yes" # Assign "Yes" to rows >= 55

# Drop the first column, 已確定資料無誤，刪除原本抓下來的menopause
large <- large[, -1]

large_tumor_matrix <- large[row.names(large) %in% tumor_barcode, ]
large_normal_matrix <- large[row.names(large) %in% normal_barcode, ]

# ======= TN Match ============
# Find Tumor-Normal Pair DEGs
# Extract the first 12 characters of row names in normal_matrix
TN_match <- substr(rownames(large[row.names(large) %in% normal_barcode, ]), 1, 12)
#TN_match <- substr(rownames(normal_matrix_large), 1, 12)
# Find common partial row names
TN_match_matrix <- large[substr(rownames(large), 1, 12) %in% TN_match, ]
# Remove meta samples
TN_match_matrix <- TN_match_matrix[!grepl("Meta", TN_match_matrix$Sample), ] # 剩下都是T_N match, 多一塊tumor 檢體
TN_match_barcode <- rownames(TN_match_matrix)
TN_tumor_barcode <- rownames(TN_match_matrix[grepl("Tumor", TN_match_matrix$Sample), ])
TN_normal_barcode <- rownames(TN_match_matrix[grepl("Normal", TN_match_matrix$Sample), ])

TNM <- substr(rownames(large[row.names(large) %in% meta_barcode, ]), 1, 12)
TNM_matrix <- large[substr(rownames(large), 1, 12) %in% TNM, ]
TNM_matrix <- TNM_matrix[, -c(3, 14)] # tumor & meta pair => 6 pairs, should be LN, age: 31,35,40,45,49,58

# Subgroups
TN_match_young_matrix <- TN_match_matrix[grepl("Young", TN_match_matrix$age), ]
TN_match_old_matrix <- TN_match_matrix[grepl("Old", TN_match_matrix$age), ]

# Survival 
# Adding Survival Data
# ======== 新增days to new tumor event, days to death
# Load the survival package if not already loaded

additional_2 <- read.csv('/Users/mona/Desktop/Subgroup/Fresh Re Run/survival metadata_2.csv'
                         , header = T)
additional <- read.csv('/Users/mona/Desktop/Subgroup/Annotated_metadata.csv'
                       , header = T)
rownames(additional_2) <- additional$bcr_patient_barcode
additional_2$row_names <- rownames(additional_2)
additional_2$row_names <- gsub("-", ".", additional$bcr_patient_barcode)
rownames(additional_2) <- additional_2$row_names

# ========= 新增survival_info
# Create survival_info dataframe
survival_info <- complete_info
# Extract the first 12 characters of row names in complete_info
partial_match <- substr(rownames(survival_info), 1, 12)
# Find common partial row names
common_partial_rows <- partial_match %in% substr(rownames(additional_2), 1, 12)
# Update 'Days to Death' column in 'complete_info' for common partial row names
survival_info$Days_to_Death <- NA
survival_info$Days_to_Death[common_partial_rows] <- additional_2$days_to_death[match(partial_match[common_partial_rows], substr(rownames(additional_2), 1, 12))]
# Update "Days to Recurrence"
survival_info$Days_to_Recurrence <- NA
survival_info$Days_to_Recurrence[common_partial_rows] <- additional_2$days_to_new_tumor_event_after_initial_treatment[match(partial_match[common_partial_rows], substr(rownames(additional_2), 1, 12))]
# Update "Recurrence Type"
survival_info$Recurrence_Type <- NA
survival_info$Recurrence_Type[common_partial_rows] <- additional_2$new_neoplasm_event_type[match(partial_match[common_partial_rows], substr(rownames(additional_2), 1, 12))]
# Update "Recurrence Site"
survival_info$Recurrence_Site <- NA
survival_info$Recurrence_Site[common_partial_rows] <- additional_2$new_neoplasm_event_occurrence_anatomic_site[match(partial_match[common_partial_rows], substr(rownames(additional_2), 1, 12))]
# Update "Days to last follow up "
survival_info$Last_follow <- NA
survival_info$Last_follow[common_partial_rows] <- additional_2$days_to_last_followup[match(partial_match[common_partial_rows], substr(rownames(additional_2), 1, 12))]
survival_info$Therapy <- NA
survival_info$Therapy[common_partial_rows] <- additional_2$therapy_type[match(partial_match[common_partial_rows], substr(rownames(additional_2), 1, 12))]
survival_info$Drug <- NA
survival_info$Drug[common_partial_rows] <- additional_2$drug_name[match(partial_match[common_partial_rows], substr(rownames(additional_2), 1, 12))]

library(dplyr)
# Update "Others" in Recurrence Site
survival_info$Recurrence_Site_sp <- NA
survival_info$Recurrence_Site_sp[common_partial_rows] <- additional_2$new_neoplasm_occurrence_anatomic_site_text[match(partial_match[common_partial_rows], substr(rownames(additional_2), 1, 12))]
# Replace the Other, specify with actual text, then delete this column
survival_info$Recurrence_Site <- ifelse(survival_info$Recurrence_Site == "Other, specify", survival_info$Recurrence_Site_sp, survival_info$Recurrence_Site)
# Remove the Recurrence_site_sp column
survival_info <- survival_info %>%
  dplyr::select(-Recurrence_Site_sp)

  
  # =======
# 把surivival info 和 large merge 在一起
# 按照年齡排序
survival_info <- survival_info[order(survival_info$age_true),]
survival_info$Recurrence_code <- NA
# Create new column Recurrence_code based on Recurrence_Type
survival_info$Recurrence_code <- ifelse(
  survival_info$Recurrence_Type == "", 0,
  ifelse(
    survival_info$Recurrence_Type %in% c("Locoregional Disease", "Locoregional Recurrence"), 1,
    ifelse(
      survival_info$Recurrence_Type == "Distant Metastasis", 2,
      ifelse(
        survival_info$Recurrence_Type == "New Primary Tumor", 9,
        NA  # Default value for unmatched cases
      )
    )
  )
)

survival_info$Death_code <- ifelse(
  is.na(survival_info$Days_to_Death), 0, 1
)

large_survival <- cbind(survival_info, large[, -(1:10)])

# Update ER,PR,HER2
IHC_info <- complete_info
# ER
IHC_info$ER <- NA
IHC_info$ER[common_partial_rows] <- additional_2$breast_carcinoma_estrogen_receptor_status[match(partial_match[common_partial_rows], substr(rownames(additional_2), 1, 12))]
# PR
IHC_info$PR <- NA
IHC_info$PR[common_partial_rows] <- additional_2$breast_carcinoma_progesterone_receptor_status[match(partial_match[common_partial_rows], substr(rownames(additional_2), 1, 12))]
# HER2
IHC_info$HER2 <- NA
IHC_info$HER2[common_partial_rows] <- additional_2$her2_immunohistochemistry_level_result[match(partial_match[common_partial_rows], substr(rownames(additional_2), 1, 12))]
# ER %
IHC_info$ER_per <- NA
IHC_info$ER_per[common_partial_rows] <- additional_2$er_level_cell_percentage_category[match(partial_match[common_partial_rows], substr(rownames(additional_2), 1, 12))]
# PR %
IHC_info$PR_per <- NA
IHC_info$PR_per[common_partial_rows] <- additional_2$progesterone_receptor_level_cell_percent_category[match(partial_match[common_partial_rows], substr(rownames(additional_2), 1, 12))]

# 把IHC info 和 large merge 在一起
# 按照年齡排序
IHC_info <- IHC_info[order(IHC_info$age_true),]
IHC_info$menopause <- large$Menopause_True
large_IHC <- cbind(IHC_info, large_survival[, -(1:10)])

#Transform large_IHC dataframe => N0Nx, M0Mx, R0Rx
# Replace all 0 values with "N0" in the 'Node' column
large_IHC$Node[large_IHC$Node == 0] <- "N0"
# Replace all other non-NA values with "Nx" in the 'Node' column
large_IHC$Node[!is.na(large_IHC$Node) & large_IHC$Node != "N0"] <- "Nx"
# Replace 0 values with "No", 1 values with "Yes", and leave NAs unchanged in the 'Metastasis' column
large_IHC$Metastasis[large_IHC$Metastasis == 0] <- "No"
large_IHC$Metastasis[large_IHC$Metastasis == 1] <- "Yes"
large_IHC$Event[!(large_IHC$Event %in% c("NO", "YES"))] <- NA

IHC_info$age_true
large_IHC$age_true
large$age_true