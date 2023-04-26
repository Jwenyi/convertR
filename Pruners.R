#-- set work environments
ps <- c("data.table","tidyverse","foreach","doFuture","progressr","purrr","survival")
lapply(ps, function(x){
  suppressMessages(library(x, character.only = T))}) 
handlers(handler_progress(
  format   = ":spin :current/:total (:message) [:bar] :percent in :elapsed ETA: :eta",
  width    = 80,
  complete = "="
))
options(future.globals.maxSize= 5*1024*1024^2)

#--set function

Filter_SL <- function(surv_panca = "./data/GDC-PANCAN.survival.tsv",
                      surv_os ="./data/TARGET-OS.survival.tsv",
                      exprs_dir = "./expression/",
                      cas9_effect = "./data/CRISPR_gene_effect.csv",
                      rnai_effect = "./data/Achilles_gene_effect.csv",
                      ccle_exprs = "./data/CCLE_expression.csv",
                      ccle_scna = "./data/CCLE_gene_cn.csv",
                      ccle_mutation = "./data/CCLE_mutations.csv",
                      ccl_info = "./data/sample_info2.csv",
                      bp_file = "./data/BP_sim.csv",
                      cc_file = "./data/CC_sim.csv",
                      mf_file = "./data/MF_sim.csv",
                      cor_method = "spearman",
                      cor_method_cas9_rnai = "pearson",
                      type = c("BLCA", "SKCM", "SARC", "PRAD", "PAAD",  "OV", "NBL", "DLBC", "LUSC", "LUAD", "LIHC", "AML", "ALL", "KIRC", "HNSC", "STAD", "UVM", "ESCA", "UCEC", "COAD", "CESC", "BRCA", "GBM", "OS", "CHOL"),
                      cut_inactive = 0.2,
                      co_dependencies = 0.3,
                      neg_cor_cut_abs = 0.10,
                      surv_pos_p = 0.05,
                      surv_neg_p = 0.1,
                      amp_scna = 0.3, 
                      del_scna = -0.3,
                      cas9_RNAi_effect_p_pos = 0.05,
                      cas9_RNAi_effect_p_neg = 0.05,
                      Fitler_neg = TRUE,
                      pos_SL = "./data/Human_SL.csv",
                      Fitler_pos = FALSE,
                      Filter_both = FALSE
){
  
  message("Loading data...")
  
  meta_data <- list()
  for(i in c(ccle_mutation,ccl_info)){meta_data[[i]] <- fread(i,header = T)}
  
  for(i in c(cas9_effect,
             rnai_effect,
             ccle_exprs,
             ccle_scna,
             bp_file,
             cc_file,
             mf_file
  )){
    meta_data[[i]] <- fread(i,header = T)
    ids_row <- meta_data[[i]][[1]]
    meta_data[[i]] <- meta_data[[i]][,-1] %>% as.data.frame()
    row.names(meta_data[[i]]) <- ids_row
  }
  names(meta_data) <- c("ccle_mutation",
                        "ccl_info",
                        "cas9_effect",
                        "rnai_effect",
                        "ccle_exprs",
                        "ccle_scna",
                        "cc",
                        "mf",
                        "bp")
  surv <- read.table(surv_panca,header = T)
  surv$sample <- substr(surv$sample,1,nchar(surv$sample)-1)
  surv <- rbind(surv,
                read.table(surv_os,header = T))
  surv <- surv %>% dplyr::distinct(sample, .keep_all = TRUE)
  surv <- surv %>% dplyr::select("sample","OS","OS.time")
  surv$OS.time <- surv$OS.time/365
  
  genes_tcga <- fread(paste0(exprs_dir,"OS",".csv"),header = T) %>% colnames()
  message("Done.")
  message("Data pre-processing...")
  
  # removeRowsAllNa  <- function(x){x[apply(x, 1, function(y) any(!is.na(y))),]}
  # removeColsAllNa  <- function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}
  # for(i in c("cc","mf","bp")){
  #   meta_data[[i]] <- removeRowsAllNa(meta_data[[i]])
  #   meta_data[[i]] <- removeColsAllNa(meta_data[[i]])
  # }
  
  for(i in c("cc","mf","bp")){
    meta_data[[i]] <- janitor::remove_empty(meta_data[[i]],which = c("rows", "cols"))
  }
  
  colnames(meta_data[["cas9_effect"]]) <- gsub('\\s+\\(+[a-zA-Z0-9]+\\)','',colnames(meta_data[["cas9_effect"]]))
  colnames(meta_data[["rnai_effect"]]) <- gsub('\\s+\\(+[a-zA-Z0-9]+\\)','',colnames(meta_data[["rnai_effect"]]))
  colnames(meta_data[["ccle_exprs"]]) <- gsub('\\s+\\(+[a-zA-Z0-9]+\\)','',colnames(meta_data[["ccle_exprs"]]))
  colnames(meta_data[["ccle_scna"]]) <- gsub('\\s+\\(+[a-zA-Z0-9]+\\)','',colnames(meta_data[["ccle_scna"]]))
  
  go_cas9_nodes <- Reduce(intersect,list(colnames(meta_data[["bp"]]),
                                         colnames(meta_data[["cc"]]),
                                         colnames(meta_data[["mf"]]),
                                         genes_tcga,
                                         colnames(meta_data[["ccle_exprs"]]),
                                         colnames(meta_data[["ccle_scna"]]),
                                         colnames(meta_data[["cas9_effect"]]),
                                         colnames(meta_data[["rnai_effect"]]),
                                         unique(meta_data[["ccle_mutation"]][["Hugo_Symbol"]])
  ))#[1:200]##should remove [1:200]
  nodes_to_predict <- go_cas9_nodes
  #---------------------------------------
  if(type == "ALL" | Fitler_neg == TRUE){
    message("Extracting Pos SLs...(For Filtering Negative Samples only)")
    pos_SLs <- fread(pos_SL,header = T)
    # pos_SLs <- pos_SLs %>% subset(pos_SLs$SL.source %in% c("GenomeRNAi", "High Throughput", "RNAi Screen", 
    #                                                        "Synlethality;GenomeRNAi", "Synlethality", 
    #                                                        "Drug Screen", "CRISPR/CRISPRi","Decipher",
    #                                                        "Low Throughput","Synlethality;Decipher","High Throughput|Low Throughput",
    #                                                        "GenomeRNAi;Decipher"))
    pos_SL_nodes <- c(pos_SLs[["gene_a.name"]],
                      pos_SLs[["gene_b.name"]]) %>% unique()
    
    go_cas9_nodes <- intersect(go_cas9_nodes,pos_SL_nodes)
    pos_SLs <- pos_SLs %>% subset(pos_SLs[["gene_a.name"]] %in% go_cas9_nodes &
                                    pos_SLs[["gene_b.name"]] %in% go_cas9_nodes) %>% 
      dplyr::select("gene_a.name","gene_b.name")
    pos_SLs$label <- 1
    colnames(pos_SLs) <- c("node_1","node_2","label")
  }
  #---------------------------------------
  for(i in c("cc","mf","bp")){
    meta_data[[i]] <- meta_data[[i]] %>% dplyr::select(all_of(go_cas9_nodes))
    meta_data[[i]] <- meta_data[[i]][go_cas9_nodes,]
  }
  
  for(i in c("ccle_exprs","ccle_scna")){
    meta_data[[i]] <- meta_data[[i]] %>% dplyr::select(all_of(go_cas9_nodes))
  }
  
  message("All matrix intersected")
  #-- pre-processing ccle data
  #cnv
  meta_data[["ccle_scna"]] <- meta_data[["ccle_scna"]] %>% as.matrix()
  meta_data[["ccle_scna"]] <- 2^meta_data[["ccle_scna"]]-1
  meta_data[["ccle_scna"]] <- log2(meta_data[["ccle_scna"]])
  meta_data[["ccle_scna"]] <- ifelse(meta_data[["ccle_scna"]] > amp_scna,1,meta_data[["ccle_scna"]])
  meta_data[["ccle_scna"]] <- ifelse(meta_data[["ccle_scna"]] < del_scna,-1,meta_data[["ccle_scna"]])
  meta_data[["ccle_scna"]] <- ifelse(meta_data[["ccle_scna"]] >= del_scna & meta_data[["ccle_scna"]] <= amp_scna,0,meta_data[["ccle_scna"]])
  #mutation
  
  #meta_data[["ccle_mutation"]] = fread(ccle_mutation,header = T)
  meta_data[["ccle_mutation"]] <- meta_data[["ccle_mutation"]] %>% dplyr::select("Hugo_Symbol","Variant_Classification","DepMap_ID")
  meta_data[["ccle_mutation"]]$Variant_Classification <- ifelse(
    meta_data[["ccle_mutation"]]$Variant_Classification %in% c(
      "Nonsense_Mutation","Frame_Shift_Ins", "Frame_Shift_Del"),-1,0)
  meta_data[["ccle_mutation"]] <- meta_data[["ccle_mutation"]] %>% distinct(Hugo_Symbol,DepMap_ID,.keep_all = T)
  meta_data[["ccle_mutation"]] <- meta_data[["ccle_mutation"]] %>% pivot_wider(names_from = Hugo_Symbol, values_from = Variant_Classification, values_fill = 0)
  mut_meta_ids <- meta_data[["ccle_mutation"]][[1]]
  meta_data[["ccle_mutation"]] <-  meta_data[["ccle_mutation"]][,-1] %>% as.data.frame()
  rownames(meta_data[["ccle_mutation"]]) <- mut_meta_ids
  
  #Removing ccl containing neither exprs scna and mut data in cas9 and rnai
  NAs_ccl <- unique(
    c(
      Reduce(setdiff,list(rownames(meta_data[["cas9_effect"]]),
                          rownames(meta_data[["ccle_exprs"]]),
                          rownames(meta_data[["ccle_mutation"]]),
                          rownames(meta_data[["ccle_scna"]])
      )),
      Reduce(setdiff,list(rownames(meta_data[["rnai_effect"]]),
                          rownames(meta_data[["ccle_exprs"]]),
                          rownames(meta_data[["ccle_mutation"]]),
                          rownames(meta_data[["ccle_scna"]])
      ))
    )
  )
  
  meta_data[["cas9_effect"]] <- meta_data[["cas9_effect"]] %>%
    subset(rownames(meta_data[["cas9_effect"]]) %in% NAs_ccl == FALSE)
  meta_data[["rnai_effect"]] <- meta_data[["rnai_effect"]] %>%
    subset(rownames(meta_data[["rnai_effect"]]) %in% NAs_ccl == FALSE)
  
  
  #--------------------Filter only POS in each cancer-----
  if(Fitler_neg == FALSE & Fitler_pos == TRUE & Filter_both == FALSE){
    res_data <- list()
    for(i in type){
      # dark_SL_clone <- unite(dark_SL, "Pair", 
      #                        node_1, node_2, 
      #                        sep = " ", 
      #                        remove = T) %>% as.data.table() %>% setkey("Pair")
      
      message("Detecting SLs in ",i,"...")
      message("Examining ",i," Cas9/RNAi-included ccl status (inactive and active of given genes)...")
      
      subtype_ccl_id <- subset(meta_data[["ccl_info"]], meta_data[["ccl_info"]][["Subtype"]] == i)[["DepMap_ID"]]
      subtype_cas9 <- meta_data[["cas9_effect"]] %>% subset(rownames(meta_data[["cas9_effect"]]) %in% subtype_ccl_id)
      subtype_rnai <- meta_data[["rnai_effect"]] %>% subset(rownames(meta_data[["rnai_effect"]]) %in% subtype_ccl_id)
      subtype_exprs <- meta_data[["ccle_exprs"]] %>% subset(rownames(meta_data[["ccle_exprs"]]) %in% subtype_ccl_id)
      subtype_mut <- meta_data[["ccle_mutation"]] %>% subset(rownames(meta_data[["ccle_mutation"]]) %in% subtype_ccl_id)
      subtype_scna <- meta_data[["ccle_scna"]] %>% subset(rownames(meta_data[["ccle_scna"]]) %in% subtype_ccl_id) %>% as.data.frame()
      subtype_cas9_ids <- rownames(subtype_cas9)
      subtype_rnai_ids <- rownames(subtype_rnai)
      
      ccl_cut_low <- subtype_exprs %>% 
        apply(2,function(x){return(quantile(x, prob = seq(0,1,cut_inactive))[[2]])})
      
      registerDoFuture()
      plan(multisession,workers = 76)
      with_progress({
        p <- progressor(along = go_cas9_nodes) 
        ccl_status_cas9 <- foreach(ccl_node = go_cas9_nodes,
                                   .packages = c("data.table","tidyverse","purrr"),
                                   .combine = cbind) %dopar% {
                                     p()
                                     # Testing mutation
                                     mut_ccl_node <- subtype_mut[[ccl_node]]
                                     names(mut_ccl_node) <- rownames(subtype_mut)
                                     
                                     # Testing scna
                                     scna_ccl_node <- subtype_scna[[ccl_node]]
                                     names(scna_ccl_node) <- rownames(subtype_scna)
                                     
                                     # Testing exprs
                                     exprs_ccl_node <- subtype_exprs[[ccl_node]]
                                     exprs_ccl_node <- ifelse(exprs_ccl_node <= ccl_cut_low[[ccl_node]],-1,1)
                                     names(exprs_ccl_node) <- rownames(subtype_exprs)
                                     
                                     # grouping
                                     group_ccl_node <- subtype_cas9_ids %>% purrr::map(function(cas_id){
                                       out = 0
                                       names(out) = cas_id
                                       # same ccl has censored data
                                       judge_group <- c(ifelse(cas_id %in% names(mut_ccl_node),mut_ccl_node[[cas_id]],99),
                                                        ifelse(cas_id %in% names(scna_ccl_node),scna_ccl_node[[cas_id]],99),
                                                        ifelse(cas_id %in% names(exprs_ccl_node),exprs_ccl_node[[cas_id]],99)
                                       )
                                       
                                       if(-1 %in% judge_group){
                                         return(out)
                                       }else{
                                         out[[cas_id]] = 1
                                         return(out)
                                       } 
                                     })%>% unlist()
                                     df_out <- group_ccl_node %>% as.data.frame()
                                     colnames(df_out) <- all_of(ccl_node)
                                     return(df_out)
                                   }
      })
      
      with_progress({
        p <- progressor(along = go_cas9_nodes) 
        ccl_status_rnai <- foreach(ccl_node = go_cas9_nodes,
                                   .packages = c("data.table","tidyverse","purrr"),
                                   .combine = cbind) %dopar% {
                                     p()
                                     # Testing mutation
                                     mut_ccl_node <- subtype_mut[[ccl_node]]
                                     names(mut_ccl_node) <- rownames(subtype_mut)
                                     
                                     # Testing scna
                                     scna_ccl_node <- subtype_scna[[ccl_node]]
                                     names(scna_ccl_node) <- rownames(subtype_scna)
                                     
                                     # Testing exprs
                                     exprs_ccl_node <- subtype_exprs[[ccl_node]]
                                     exprs_ccl_node <- ifelse(exprs_ccl_node <= ccl_cut_low[[ccl_node]],-1,1)
                                     names(exprs_ccl_node) <- rownames(subtype_exprs)
                                     
                                     # grouping
                                     group_ccl_node <- subtype_rnai_ids %>% purrr::map(function(cas_id){
                                       out = 0
                                       names(out) = cas_id
                                       # same ccl has censored data
                                       judge_group <- c(ifelse(cas_id %in% names(mut_ccl_node),mut_ccl_node[[cas_id]],99),
                                                        ifelse(cas_id %in% names(scna_ccl_node),scna_ccl_node[[cas_id]],99),
                                                        ifelse(cas_id %in% names(exprs_ccl_node),exprs_ccl_node[[cas_id]],99)
                                       )
                                       
                                       if(-1 %in% judge_group){
                                         return(out)
                                       }else{
                                         out[[cas_id]] = 1
                                         return(out)
                                       } 
                                     })%>% unlist()
                                     df_out <- group_ccl_node %>% as.data.frame()
                                     colnames(df_out) <- all_of(ccl_node)
                                     return(df_out)
                                   }
      })
      
      
      message(paste0("Curating ",i," clinical data..."))
      
      cli_exprs <- fread(paste0(exprs_dir,i,".csv"),header = T)
      cli_exprs_matrix <- cli_exprs %>% dplyr::select(-1)
      
      
      #--Curating Survival Data
      cli_exprs_matrix <- cli_exprs_matrix %>% as.data.frame()
      row.names(cli_exprs_matrix) <- cli_exprs[[1]]
      genes_cut_surv <- cli_exprs_matrix %>% 
        apply(2,function(x){return(quantile(x, prob = seq(0,1,cut_inactive))[[2]])})
      genes_cut_high <- cli_exprs_matrix %>% 
        apply(2,function(x){return(quantile(x, prob = seq(0,1,cut_inactive))[[5]])})
      patients <- intersect(row.names(cli_exprs_matrix),surv$sample)
      surv_used <- surv %>% 
        subset(surv$sample %in% patients)
      patients <- intersect(row.names(cli_exprs_matrix),surv$sample)
      surv_used <- surv %>% 
        subset(surv$sample %in% patients)
      
      message("Filtering positive samples within ",i,"...")
      
      with_progress({
        p <- progressor(along = go_cas9_nodes)
        res_data[[i]][["All_filtered_res"]] <- foreach(k = go_cas9_nodes,
                                .packages = c("data.table","tidyverse","purrr","survival"),
                                .combine = rbind) %dopar% {
                                  p()
                                  
                                  partners <- purrr::map(go_cas9_nodes,function(x){
                                    
                                    # Set Tree-structured filter pipeline
                                    if(k == x | 
                                       length(table(ccl_status_cas9[[k]])) == 1 |
                                       length(table(ccl_status_rnai[[k]])) == 1){
                                      #Discard self-paired SLs
                                      
                                      #Seting all params as NA!!!!!!
                                      temp <- c(NA, #Pair
                                                NA, #Type
                                                NA, #node_1
                                                NA, #node_2
                                                NA, #cor
                                                NA, #abs_cor
                                                NA, #p
                                                NA, #statistic
                                                NA, #wilcox_cas9_p
                                                NA, #wilcox_cas9_stats
                                                NA, #wilcox_cas9_meanscore_inactive
                                                NA, #wilcox_cas9_meanscore_active
                                                NA, #wilcox_cas9_diff
                                                NA, #wilcox_rnai_p
                                                NA, #wilcox_rnai_stats
                                                NA, #wilcox_rnai_meanscore_inactive
                                                NA, #wilcox_rnai_meanscore_active
                                                NA, #wilcox_rnai_diff
                                                NA, #cas_9_pass
                                                NA, #rnai_pass
                                                NA #Surv p 
                                      )
                                      
                                      
                                    }else{
                                      #Let OTHER SLs passed to pipeline
                                      
                                      #Step.1 Do clinical correlation detection
                                      cor_res <- cor.test(cli_exprs_matrix[[k]],
                                                          cli_exprs_matrix[[x]],
                                                          method = cor_method)
                                      
                                      # Allocating SLs into branches
                                      
                                      if(cor_res[["p.value"]] < 0.05 & cor_res[["estimate"]][[1]] > co_dependencies){
                                        #Let POS SLs meeting our clinical cor requirement pass
                                        
                                        #Step.2_POS Do Cas9/RNAi-based functional examination procedure
                                        
                                        #Cas9
                                        group_cas9_k <- ccl_status_cas9[[k]]
                                        #removing group sample size < 2
                                        
                                        if(table(group_cas9_k)[[1]] >= 2 &
                                           table(group_cas9_k)[[2]] >= 2){
                                          exprs_cas9_x <- subtype_cas9[[x]]
                                          
                                          wilcox_cas9_x <- exprs_cas9_x[group_cas9_k == 0]
                                          wilcox_cas9_y <- exprs_cas9_x[group_cas9_k == 1]
                                          
                                          wilcox_cas9 <- wilcox.test(
                                            wilcox_cas9_x,
                                            wilcox_cas9_y
                                          )
                                          wilcox_cas9_p = wilcox_cas9$p.value
                                          wilcox_cas9_stats = wilcox_cas9$statistic
                                          wilcox_cas9_meanscore_inactive = mean(wilcox_cas9_x)
                                          wilcox_cas9_meanscore_active = mean(wilcox_cas9_y)
                                          wilcox_cas9_diff = ifelse(
                                            wilcox_cas9_meanscore_inactive < wilcox_cas9_meanscore_active,
                                            "Important_in_inactive_ccl",
                                            "Not_important_in_inactive_ccl")
                                          
                                        }else{
                                          wilcox_cas9_p = NA
                                          wilcox_cas9_stats = NA
                                          wilcox_cas9_meanscore_inactive = NA
                                          wilcox_cas9_meanscore_active = NA
                                          wilcox_cas9_diff = NA
                                        }
                                        #RNAi
                                        group_rnai_k <- ccl_status_rnai[[k]]
                                        #removing group sample size < 2
                                        if(table(group_rnai_k)[[1]] >= 2 &
                                           table(group_rnai_k)[[2]] >= 2){
                                          exprs_rnai_x <- subtype_rnai[[x]]
                                          
                                          wilcox_rnai_x <- exprs_rnai_x[group_rnai_k == 0]
                                          wilcox_rnai_y <- exprs_rnai_x[group_rnai_k == 1]
                                          
                                          wilcox_rnai <- wilcox.test(
                                            wilcox_rnai_x,
                                            wilcox_rnai_y
                                          )
                                          wilcox_rnai_p = wilcox_rnai$p.value
                                          wilcox_rnai_stats = wilcox_rnai$statistic
                                          wilcox_rnai_meanscore_inactive = mean(wilcox_rnai_x)
                                          wilcox_rnai_meanscore_active = mean(wilcox_rnai_y)
                                          wilcox_rnai_diff = ifelse(
                                            wilcox_rnai_meanscore_inactive < wilcox_rnai_meanscore_active,
                                            "Important_in_inactive_ccl",
                                            "Not_important_in_inactive_ccl")
                                          
                                        }else{
                                          wilcox_rnai_p = NA
                                          wilcox_rnai_stats = NA
                                          wilcox_rnai_meanscore_inactive = NA
                                          wilcox_rnai_meanscore_active = NA
                                          wilcox_rnai_diff = NA
                                        }
                                        
                                        #judgement_cas9
                                        # 0 represents NO pass, 1 represents PASS
                                        if(anyNA(wilcox_cas9_p)){
                                          judgement_cas9 = 0
                                        }else{
                                          if(wilcox_cas9_p < cas9_RNAi_effect_p_pos & 
                                             wilcox_cas9_diff == "Important_in_inactive_ccl"){
                                            judgement_cas9 = 1
                                          }else{
                                            judgement_cas9 = 0
                                          }
                                        }
                                        
                                        #judgement_rnai
                                        # 0 represents NO pass, 1 represents PASS
                                        if(anyNA(wilcox_rnai_p)){
                                          judgement_rnai = 0
                                        }else{
                                          if(wilcox_rnai_p < cas9_RNAi_effect_p_pos & 
                                             wilcox_rnai_diff == "Important_in_inactive_ccl"){
                                            judgement_rnai = 1
                                          }else{
                                            judgement_rnai = 0
                                          }
                                        }
                                        
                                        if(judgement_cas9 == 1 | judgement_rnai == 1){
                                          #Let POS SLs meeting RNAi/cas9 requirement pass
                                          
                                          #Step.3_POS Do surv analyses (It return log_rank_p)
                                          #Grouping patients and removing those without inactive/active SL
                                          groups <- purrr::map(patients, function(w){
                                            node1_cli_exprs <- cli_exprs_matrix[w,k]
                                            node2_cli_exprs <- cli_exprs_matrix[w,x]
                                            if(node1_cli_exprs < genes_cut_surv[k][[1]] & node2_cli_exprs < genes_cut_surv[x][[1]]){
                                              return(0)
                                            }else if(node1_cli_exprs > genes_cut_high[k][[1]] & node2_cli_exprs > genes_cut_high[x][[1]]){
                                              return(1)
                                            }else{return(NA)}
                                          }) %>% unlist()
                                          
                                          #survival 
                                          if(length(unique(na.omit(groups))) <= 1){
                                            log_rank_p = NA
                                          }else if(table(na.omit(groups))[[1]] < 3 | table(na.omit(groups))[[2]] < 3){
                                            log_rank_p = NA
                                          }else{
                                            surv_df <- cbind(sample = patients,
                                                             group = groups)
                                            surv_df <- merge(surv_used,surv_df,by = "sample")
                                            if(length(unique(na.omit(surv_df[["OS"]]))) <= 1){
                                              log_rank_p = NA
                                            }else{
                                              surv_df <- surv_df[complete.cases(surv_df),]
                                              fit <- survfit(Surv(OS.time, OS) ~group, data = surv_df)
                                              diff <- survdiff(Surv(OS.time, OS) ~group,data = surv_df)
                                              log_rank_p = 1-pchisq(diff$chisq,df=1)
                                            }
                                          }
                                          
                                        }else{
                                          #Drop POS SLs NOT meeting RNAi/cas9 requirement
                                          log_rank_p = NA
                                        }
                                        
                                        #outputing
                                        if(anyNA(log_rank_p)){
                                          temp <- c(paste0(k," ",x),
                                                    "NS",
                                                    k,
                                                    x,
                                                    cor_res[["estimate"]],
                                                    abs(cor_res[["estimate"]]),
                                                    cor_res[["p.value"]],
                                                    cor_res[["statistic"]],
                                                    wilcox_cas9_p,
                                                    wilcox_cas9_stats,
                                                    wilcox_cas9_meanscore_inactive,
                                                    wilcox_cas9_meanscore_active,
                                                    wilcox_cas9_diff,
                                                    wilcox_rnai_p,
                                                    wilcox_rnai_stats,
                                                    wilcox_rnai_meanscore_inactive,
                                                    wilcox_rnai_meanscore_active,
                                                    wilcox_rnai_diff,
                                                    judgement_cas9,
                                                    judgement_rnai,
                                                    log_rank_p)
                                        }else if(log_rank_p < surv_pos_p){
                                          temp <- c(paste0(k," ",x),
                                                    "POS",
                                                    k,
                                                    x,
                                                    cor_res[["estimate"]],
                                                    abs(cor_res[["estimate"]]),
                                                    cor_res[["p.value"]],
                                                    cor_res[["statistic"]],
                                                    wilcox_cas9_p,
                                                    wilcox_cas9_stats,
                                                    wilcox_cas9_meanscore_inactive,
                                                    wilcox_cas9_meanscore_active,
                                                    wilcox_cas9_diff,
                                                    wilcox_rnai_p,
                                                    wilcox_rnai_stats,
                                                    wilcox_rnai_meanscore_inactive,
                                                    wilcox_rnai_meanscore_active,
                                                    wilcox_rnai_diff,
                                                    judgement_cas9,
                                                    judgement_rnai,
                                                    log_rank_p)
                                        }else{
                                          temp <- c(paste0(k," ",x),
                                                    "NS",
                                                    k,
                                                    x,
                                                    cor_res[["estimate"]],
                                                    abs(cor_res[["estimate"]]),
                                                    cor_res[["p.value"]],
                                                    cor_res[["statistic"]],
                                                    wilcox_cas9_p,
                                                    wilcox_cas9_stats,
                                                    wilcox_cas9_meanscore_inactive,
                                                    wilcox_cas9_meanscore_active,
                                                    wilcox_cas9_diff,
                                                    wilcox_rnai_p,
                                                    wilcox_rnai_stats,
                                                    wilcox_rnai_meanscore_inactive,
                                                    wilcox_rnai_meanscore_active,
                                                    wilcox_rnai_diff,
                                                    judgement_cas9,
                                                    judgement_rnai,
                                                    log_rank_p)
                                        }
                                        
                                      }else{
                                        #Step.1 Discard other SLs
                                        
                                        temp <- c(paste0(k," ",x),
                                                  "NS",
                                                  k,
                                                  x,
                                                  cor_res[["estimate"]],
                                                  abs(cor_res[["estimate"]]),
                                                  cor_res[["p.value"]],
                                                  cor_res[["statistic"]],
                                                  NA, #wilcox_cas9_p
                                                  NA, #wilcox_cas9_stats
                                                  NA, #wilcox_cas9_meanscore_inactive
                                                  NA, #wilcox_cas9_meanscore_active
                                                  NA, #wilcox_cas9_diff
                                                  NA, #wilcox_rnai_p
                                                  NA, #wilcox_rnai_stats
                                                  NA, #wilcox_rnai_meanscore_inactive
                                                  NA, #wilcox_rnai_meanscore_active
                                                  NA, #wilcox_rnai_diff
                                                  NA, #cas9_pass
                                                  NA, #rnai_pass
                                                  NA #Surv p 
                                        )
                                        
                                      }
                                    }
                                    names(temp) <- c("Pair",
                                                     "Type",
                                                     "node_1",
                                                     "node_2",
                                                     "cor",
                                                     "abs_cor",
                                                     "p",
                                                     "statistic",
                                                     "wilcox_cas9_p",
                                                     "wilcox_cas9_stats",
                                                     "wilcox_cas9_meanscore_inactive",
                                                     "wilcox_cas9_meanscore_active",
                                                     "wilcox_cas9_diff",
                                                     "wilcox_rnai_p",
                                                     "wilcox_rnai_stats",
                                                     "wilcox_rnai_meanscore_inactive",
                                                     "wilcox_rnai_meanscore_active",
                                                     "wilcox_rnai_diff",
                                                     "cas9_pass",
                                                     "rnai_pass",
                                                     "log_rank_p")
                                    
                                    return(temp)
                                  }) %>% bind_rows() 
                                  
                                  return(partners)
                                } %>% as.data.table() %>% setkey("Pair")
        
      })
      rm(p)
      #pos_SLs <- Filtered_res %>% subset(Filtered_res$Type == "POS")
      res_data[[i]][["pos_SLs"]] <- res_data[[i]][["All_filtered_res"]] %>% 
        subset(res_data[[i]][["All_filtered_res"]][["Type"]] == "POS")
      #res_data[[i]][["All_filtered_res"]] <- Filtered_res

      #hold grouped data
      res_data[[i]][["group_info"]] <- list(ccl_cut_low,
                                            ccl_status_cas9,
                                            ccl_status_rnai,
                                            subtype_cas9,
                                            subtype_rnai,
                                            subtype_mut,
                                            subtype_scna,
                                            subtype_exprs,
                                            cli_exprs_matrix,
                                            genes_cut_surv,
                                            genes_cut_high,
                                            surv_used)
      names(res_data[[i]][["group_info"]]) <- c("ccl_cut_low",
                                                "ccl_status_cas9",
                                                "ccl_status_rnai",
                                                "subtype_cas9",
                                                "subtype_rnai",
                                                "subtype_mut",
                                                "subtype_scna",
                                                "subtype_exprs",
                                                "cli_exprs_matrix",
                                                "genes_cut_surv",
                                                "genes_cut_high",
                                                "surv_used")
    }
  }
  
  #--------------------Filter only POS in each cancer-----  
  return(res_data)
}
