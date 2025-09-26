#Collected by: Bea Szeitz
#Contains custom functions used across multiple R Notebooks

export_heatmap_png <- function(filename, width, height, ht_object, annot_list = NULL, not_heatmap=FALSE){
  png(filename = paste0(filename,".png"),
      width = width, height = height, units = "in",res = 600)
  
  if (not_heatmap){
    print(ht_object)
  } else {
    draw(ht_object)
    if (!is.null(annot_list)){
      for (i in 1:length(annot_list)){
        draw(annot_list[[i]][[1]], x = unit(annot_list[[i]][[2]], "npc"), y = unit(annot_list[[i]][[3]], "npc"))
      }
    }
  }
  
  dev.off()
  
  
  svg(filename = paste0(filename,".svg"),
      width = width, height = height)
  
  if (not_heatmap){
    print(ht_object)
  } else {
    draw(ht_object)
    if (!is.null(annot_list)){
      for (i in 1:length(annot_list)){
        draw(annot_list[[i]][[1]], x = unit(annot_list[[i]][[2]], "npc"), y = unit(annot_list[[i]][[3]], "npc"))
      }
    }
  }
  
  dev.off()
  
  pdf(file = paste0(filename,".pdf"),
      width = width, height = height)
  
  if (not_heatmap){
    print(ht_object)
  } else {
    draw(ht_object)
    if (!is.null(annot_list)){
      for (i in 1:length(annot_list)){
        draw(annot_list[[i]][[1]], x = unit(annot_list[[i]][[2]], "npc"), y = unit(annot_list[[i]][[3]], "npc"))
      }
    }
  }
  
  dev.off()
}

filter_missingvalues <- function (df, percentage) {
  df <- as.data.frame(df)
  perc <- ceiling(ncol(df)*(percentage/100))
  vv <- rowSums(!is.na(df))
  to_be_kept <- unlist(lapply(vv, function(x){ 
    if (x < perc) { y <- F}
    else y <- T}))
  return(as.data.frame(df[to_be_kept,]))
}


select_gene_with_highest_expr <- function(expr_mat, rowannot){
  # rowannot 1st column: rowname, 2nd column: genename
  rowannot <- rowannot[rowannot[,1] %in% row.names(expr_mat),]
  expr_mat <- expr_mat[rowannot[,1],]
  
  original_rownames <- row.names(expr_mat)
  
  row.names(expr_mat) <- paste(rowannot[,2], seq(1,nrow(rowannot)), sep="_")
  
  Gene_summary <- data.frame(gene = rowannot[,2],
                             rowname = row.names(expr_mat),
                             original_rowname = original_rownames,
                             sum_expr = rowSums(expr_mat, na.rm = TRUE),
                             row.names = row.names(expr_mat))
  Gene_summary <- Gene_summary[order(Gene_summary$sum_expr, decreasing = T),]
  
  dupl.genes <- unique(Gene_summary[(duplicated(Gene_summary$gene)),"gene"])
  print(length(dupl.genes))
  row_to_remove <- vector()
  
  i=1
  for (i in 1:length(dupl.genes)){
    Gene_summary.sub <- Gene_summary[Gene_summary$gene == dupl.genes[i],]
    row_to_remove <- c(row_to_remove, Gene_summary.sub[-1,"rowname"])
  }
  
  reduced_tab <- expr_mat[! row.names(expr_mat) %in% row_to_remove,]
  row.names(reduced_tab) <- sapply(strsplit(row.names(reduced_tab), split="_"),"[[",1)
  Gene_summary$removed <- ifelse(Gene_summary$rowname %in% row_to_remove, "*","")
  
  return(list(newtab=as.data.frame(reduced_tab),summary=Gene_summary))
  
}



# Function to create unique expression matrices from replicate data
average_expr_data <- function(expr_dat, annot_tab) {
  #annot_tab should contain 2 columns: original sample names and the larger IDs to which we need to collapse
  
  unique_specimens <- unique(annot_tab[,2])
  
  # Initialize output matrices
  new_expr_dat <- matrix(nrow = nrow(expr_dat), ncol = length(unique_specimens))
  colnames(new_expr_dat) <- unique_specimens
  rownames(new_expr_dat) <- rownames(expr_dat)
  
  # Fill matrices
  i = 1
  for (i in seq_along(unique_specimens)) {
    
    specimen <- unique_specimens[i]
    
    # RNA part
    corresponding_samples <- annot_tab[annot_tab[,2] == specimen, 1]
    
    if (length(corresponding_samples) == 1) {
      new_expr_dat[, specimen] <- expr_dat[, corresponding_samples]
    } else if (length(corresponding_samples) > 1) {
      new_expr_dat[, specimen] <- apply(expr_dat[, corresponding_samples, drop = FALSE], 1, mean, na.rm = TRUE)
    }
    
  }
  
  new_expr_dat[is.nan(new_expr_dat)] <- NA
  return(new_expr_dat)
}


# Retrieved files from the Reactome website
#reactome_ids <- read.delim("files/reactome_ids_2025-04-27.txt", header = FALSE,sep = "\t")
#reactome_ids <- reactome_ids[grepl("HSA",reactome_ids[,1]),c(1,2)]
#colnames(reactome_ids) = c("ID", "Description")
#row.names(reactome_ids) <- reactome_ids$ID
#reactome_hierarchy <- read.delim("files/reactome_hierarchy_2025-04-27.txt", header = FALSE)
#reactome_hierarchy <- reactome_hierarchy[grepl("HSA",reactome_hierarchy[,1]),]
#colnames(reactome_hierarchy) <- c("parent","child")
#reactome_hierarchy$parent_description <- reactome_ids[reactome_hierarchy$parent,"Description"]
#reactome_hierarchy$child_description <- reactome_ids[reactome_hierarchy$child,"Description"]
#save(reactome_hierarchy,
#     file = "rdata/reactome_hierarchy.RData")
load("rdata/reactome_hierarchy.RData")

#pw_id <- "R-HSA-111447"
map_reactome_hierarchy <- function(pw_id) {
  
  explore_paths <- function(child_id, path_so_far = c()) {
    # If child_id is not found among children, we've reached the top
    if (!(child_id %in% reactome_hierarchy$child)) {
      return(list(path_so_far))
    }
    
    # Find all parent entries
    subdf <- reactome_hierarchy[reactome_hierarchy$child == child_id, ]
    paths <- list()
    
    for (i in seq_len(nrow(subdf))) {
      parent_id <- subdf$parent[i]
      parent_desc <- subdf$parent_description[i]
      
      # Recurse for each parent
      new_paths <- explore_paths(parent_id, c(path_so_far, parent_desc))
      paths <- c(paths, new_paths)
    }
    
    return(paths)
  }
  
  all_paths <- explore_paths(pw_id)
  
  # Process the results into categories and subcategories
  categories <- sapply(all_paths, function(path) rev(path)[1])  # top-level parent
  subcategories <- sapply(all_paths, function(path) paste(rev(path[-length(path)]), collapse = " -> "))
  
  # Collapse multiple categories and subcategories into combined strings
  category_result <- paste(unique(categories), collapse = " // ")
  subcategory_result <- paste(unique(subcategories), collapse = " // ")
  
  if (length(category_result) ==1 & category_result =="NULL"){
    #print(pw_id)
    category_result = unique(reactome_hierarchy[reactome_hierarchy$parent == pw_id,"parent_description"])
    if (length(category_result)==0){
      category_result = "---"
    }
    
    subcategory_result = "---"
    
  }
  
  subcategory_result[subcategory_result==""] <- "---" 
  
  return(list(category = category_result,
              subcategory = subcategory_result))
}



#genes_to_test = unique(pos_correlating_tumors)
#universe_genes = unique(tumors_tpm_spearman_hpa$gene)

kegg_and_reactome_ora <- function(genes_to_test, universe_genes){
  kegg_res <- enrichKEGG(gene = genes_to_test,
                         organism = 'hsa',
                         universe = universe_genes,
                         minGSSize = 10,
                         maxGSSize = 500,
                         pvalueCutoff = 1)
  reac_res <- enrichPathway(gene=genes_to_test, 
                            organism = 'human',
                            universe = universe_genes,
                            minGSSize = 10,
                            maxGSSize = 500,
                            pvalueCutoff = 1, 
                            readable=FALSE)
  
  colnames(kegg_res@result)
  reac_res_df <- reac_res@result
  #grep("R-HSA-4839726",reac_res_df$ID)
  rec_mappings <- lapply(reac_res@result$ID, map_reactome_hierarchy)
  reac_res_df$category <- sapply(rec_mappings,"[[",1)
  reac_res_df$subcategory <- sapply(rec_mappings,"[[",2)
  reac_res_df$Database <- "Reactome"
  
  both_res <- kegg_res
  both_res@result$Database <- "KEGG"
  both_res@result <- rbind(both_res@result,
                           reac_res_df[,colnames(both_res@result)])
  
  #gd = "1845/3654/4790/51135/5594/5595/5598/5608/6195/6416/7186/9261/929"
  entrez_map_sub = suppressMessages(bitr(genes_to_test, 
                                         fromType="ENTREZID", 
                                         toType = "SYMBOL", 
                                         OrgDb=org.Hs.eg.db,
                                         drop = FALSE))
  entrez_map_sub$SYMBOL <- ifelse(is.na(entrez_map_sub$SYMBOL),"-",entrez_map_sub$SYMBOL)
  entrez_map_sub = entrez_map_sub %>% 
    group_by(ENTREZID) %>% 
    summarise(SYMBOL = paste(SYMBOL,collapse=","))%>% as.data.frame()
  row.names(entrez_map_sub) <- entrez_map_sub$ENTREZID
  
  #entrez_map_sub[duplicated(entrez_map_sub$SYMBOL),]
  
  
  both_res@result$geneID_symbols <- sapply(both_res@result$geneID, function(gd){
    symbols = entrez_map_sub[strsplit(gd, split = "/",fixed = T)[[1]],"SYMBOL"]
    return(paste(symbols,collapse = "/"))
  })
  
  both_res@result <- both_res@result[order(both_res@result$pvalue),c("Database",
                                                                     "Description",
                                                                     "pvalue",
                                                                     "p.adjust",
                                                                     "category",
                                                                     "subcategory",
                                                                     "geneID_symbols",
                                                                     "Count",
                                                                     "ID",
                                                                     "GeneRatio",
                                                                     "BgRatio",
                                                                     "RichFactor",
                                                                     "FoldEnrichment",
                                                                     "zScore",
                                                                     "qvalue",
                                                                     "geneID")]
  
  both_res@geneSets <- append(both_res@geneSets,
                              reac_res@geneSets)
  
  return(both_res)
  
}



# taken from https://stackoverflow.com/questions/44785229/change-chart-correlation-defaults-to-produce-best-fit-line-rather-than-smoothed
chart.Correlation.linear <- function (R, histogram = TRUE, method=c("pearson", "kendall", "spearman"), ...){ # @author R Development Core Team
    # @author modified by Peter Carl & Marek Lahoda
    # Visualization of a Correlation Matrix. On top the (absolute) value of the correlation plus the result 
    # of the cor.test as stars. On botttom, the bivariate scatterplots, with a linear regression fit. 
    # On diagonal, the histograms with probability, density and normal density (gaussian) distribution.
    
    x = PerformanceAnalytics::checkData(R, method="matrix")
    
    if(missing(method)) method=method[1] #only use one
    cormeth <- method
    
    # Published at http://addictedtor.free.fr/graphiques/sources/source_137.R
    panel.cor <- function(x, y, digits=2, prefix="", use="pairwise.complete.obs", method=cormeth, cex.cor, ...)
    {
      usr <- par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      r <- cor(x, y, use=use, method=method) # MG: remove abs here
      txt <- format(c(r, 0.123456789), digits=digits)[1]
      txt <- paste(prefix, txt, sep="")
      if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
      
      test <- cor.test(as.numeric(x),as.numeric(y), method=method)
      # borrowed from printCoefmat
      Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
                       cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                       symbols = c("***", "**", "*", ".", " "))
      # MG: add abs here and also include a 30% buffer for small numbers
      text(0.5, 0.5, txt, cex = cex * (abs(r) + .3) / 1.3)
      text(.8, .8, Signif, cex=cex, col=2)
    }
    
    #remove method from dotargs
    dotargs <- list(...)
    dotargs$method <- NULL
    rm(method)
    
    hist.panel = function (x, ...=NULL ) {
      par(new = TRUE)
      hist(x,
           col = "snow3",
           probability = TRUE,
           axes = FALSE,
           main = "",
           breaks = "FD")
      #lines(density(x, na.rm=TRUE),
      #      col = "red",
      #      lwd = 1)
      # adding line representing density of normal distribution with parameters correponding to estimates of mean and standard deviation from the data 
      #ax.x = seq(min(x), max(x), 0.1)                                                  # ax.x containts points corresponding to data range on x axis
      #density.est = dnorm(ax.x, mean = mean(x), sd = sd(x))   # density corresponding to points stored in vector ax.x 
      #lines(ax.x, density.est, col = "blue", lwd = 1, lty = 1)                                # adding line representing density into histogram
      #rug(x)
    }
    
    # Linear regression line fit over points
    reg <- function(x, y, ...) {
      points(x,y, ...)
      abline(lm(y~x), col = "red") 
    }
    
    # Draw the chart
    if(histogram)
      pairs(x, gap=0, lower.panel=reg, upper.panel=panel.cor, diag.panel=hist.panel)
    else
      pairs(x, gap=0, lower.panel=reg, upper.panel=panel.cor) 
  }



create_rank <- function(dea_res, coef_col, p_col, gene_col, use_p=TRUE){
  
  rank_df <- data.frame(gene = dea_res[,gene_col],
                        coef = as.numeric(dea_res[,coef_col]),
                        p = as.numeric(dea_res[,p_col]))
  
  if (use_p){
    rank_df$rank <- rank_df$coef #* -log10(rank_df$p)
  } else {
    rank_df$rank <- rank_df$coef * -log10(rank_df$p)
  }
  
  
  
  rank_df <- rank_df %>%
    group_by(gene) %>%
    summarise(avr_rank = mean(rank))
  
  #print(hist(rank_df$avr_rank))
  
  rank_df <- rank_df[order(rank_df$avr_rank, decreasing = TRUE),]
  
  return(setNames(rank_df$avr_rank,rank_df$gene))
  
  
}



extract_kegg_reactome_genesets <- function(universe_genes){
  kegg_res <- enrichKEGG(gene = universe_genes,
                         organism = 'hsa',
                         universe = universe_genes,
                         minGSSize = 0,
                         maxGSSize = Inf,
                         qvalueCutoff = 1,
                         pvalueCutoff = 1)
  reac_res <- enrichPathway(gene=universe_genes, 
                            organism = 'human',
                            universe = universe_genes,
                            minGSSize = 0,
                            maxGSSize = Inf,
                            pvalueCutoff = 1,
                            qvalueCutoff = 1,
                            readable=FALSE)
  
  kegg_info = unique(kegg_res@result[,c("category", "subcategory", "ID",
                                        "Description")])
  kegg_info$Database = "KEGG"
  reac_info = unique(reac_res@result[,c("ID",
                                        "Description")])
  rec_mappings <- lapply(reac_info$ID, map_reactome_hierarchy)
  reac_info$category <- sapply(rec_mappings,"[[",1)
  reac_info$subcategory <- sapply(rec_mappings,"[[",2)
  reac_info$Database <- "Reactome"
  
  
  kegg_gs_df <- enframe(kegg_res@geneSets, name = "term", value = "gene") %>%
    unnest(gene)
  reac_gs_df <- enframe(reac_res@geneSets, name = "term", value = "gene") %>%
    unnest(gene)
  
  
  #gd = "1845/3654/4790/51135/5594/5595/5598/5608/6195/6416/7186/9261/929"
  entrez_map_sub = suppressMessages(bitr(universe_genes[universe_genes!="" &
                                                          !is.na(universe_genes)], 
                                         fromType="ENTREZID", 
                                         toType = "SYMBOL", 
                                         OrgDb=org.Hs.eg.db,
                                         drop = FALSE))
  entrez_map_sub$SYMBOL <- ifelse(is.na(entrez_map_sub$SYMBOL),"-",entrez_map_sub$SYMBOL)
  entrez_map_sub = entrez_map_sub %>% 
    group_by(ENTREZID) %>% 
    summarise(SYMBOL = paste(SYMBOL,collapse=","))%>% as.data.frame()
  row.names(entrez_map_sub) <- entrez_map_sub$ENTREZID
  
  #entrez_map_sub[duplicated(entrez_map_sub$SYMBOL),]
  
  
  return(list(geneset_info = rbind(kegg_info[,c("Database","Description","category","subcategory","ID")],
                                   reac_info[,c("Database","Description","category","subcategory","ID")]),
              genesets_for_gsea = rbind(kegg_gs_df,
                                        reac_gs_df),
              entrez_to_symbol = entrez_map_sub))
  
}



run_pgsea <- function(rank_list, geneset_info_for_gsea){
  
  GSEA_list <- list()
  
  GSEA_results_matrix <- as.data.frame(matrix(nrow=0, ncol=16))
  colnames(GSEA_results_matrix) <- c("Data", "Description" ,
                                     "NES" ,"pvalue" , "p.adjust", 
                                     "Database", "category","subcategory",
                                     "ID" ,"core_enrichment_symbols",
                                     "setSize" ,"enrichmentScore" ,
                                     "qvalue" ,
                                     "rank","leading_edge","core_enrichment")
  
  i=1
  for (i in 1:length(rank_list)){
    
    rank <- na.omit(rank_list[[i]])
    rank <- rank[names(rank) !=""]
    rank <- sort(rank, decreasing = TRUE)
    #print(hist(rank))
    set.seed(12345)
    Res_gsea <- GSEA(rank, TERM2GENE=geneset_info_for_gsea$genesets_for_gsea, 
                     verbose=FALSE, pvalueCutoff = 1.0, eps = 0, seed = TRUE)
    Res_gsea@result$Data <- names(rank_list)[i]
    
    Res_gsea@result = merge(Res_gsea@result[,-grep("Description",
                                                   colnames(Res_gsea@result))],
                            geneset_info_for_gsea$geneset_info,
                            by="ID", all.x=T)
    Res_gsea@result = Res_gsea@result[order(Res_gsea@result$pvalue, decreasing = FALSE),]
    
    Res_gsea@result$core_enrichment_symbols <- sapply(Res_gsea@result$core_enrichment, function(gd){
      symbols = geneset_info_for_gsea$entrez_to_symbol[strsplit(gd, split = "/",fixed = T)[[1]],"SYMBOL"]
      return(paste(symbols,collapse = "/"))
    })
    
    
    Res_gsea@result = Res_gsea@result[,colnames(GSEA_results_matrix)]
    
    GSEA_list <- append(GSEA_list, Res_gsea)
    
    GSEA_results_matrix <- rbind(GSEA_results_matrix, 
                                 Res_gsea@result[,colnames(GSEA_results_matrix)])
    
  }
  
  return(list(GSEA_list = GSEA_list,
              GSEA_matrix = GSEA_results_matrix))
  
}


histology_correction <- function(expression_df, named_histology_values) {
  
  named_histology_values <- named_histology_values[colnames(expression_df)]
  
  corrected_residuals <- t(apply(expression_df, 1, function(gene_expr) {
    complete_idx <- which(!is.na(gene_expr) & !is.na(named_histology_values))
    
    res <- rep(NA, length(gene_expr))
    
    if (length(complete_idx) >= 2) {  
      model <- lm(gene_expr[complete_idx] ~ named_histology_values[complete_idx])
      res[complete_idx] <- resid(model)
    }
    
    return(res)
  }))
  
  corrected_residuals <- as.data.frame(corrected_residuals)
  rownames(corrected_residuals) <- rownames(expression_df)
  colnames(corrected_residuals) <- colnames(expression_df)
  
  return(corrected_residuals)
}


calculate_pooled_mad = function(expr_tab, patient_vector, sample_patient_match){
  
  mad_per_patient_list <- lapply(patient_vector,function(x){
    apply(expr_tab[,colnames(expr_tab) %in% sample_patient_match[sample_patient_match$patient ==x,1]],1,mad,na.rm=T)
  })
  names(mad_per_patient_list) <- patient_vector
  
  mad_per_patient_df <- as.data.frame(do.call(cbind, mad_per_patient_list))
  
  mad_per_patient_df$pooled_mad <- sapply(seq(1,nrow(mad_per_patient_df)), function(x){
    
    gene <- row.names(mad_per_patient_df)[x]
    
    nr_samples <- sapply(patient_vector,function(z){
      length(na.omit(expr_tab[gene,colnames(expr_tab) %in% sample_patient_match[sample_patient_match$patient ==z,1]]))
    })
    
    mads <- nr_samples * mad_per_patient_df[x,]
    
    sum(mads, na.rm=T) / sum(nr_samples)
    
  })
  
  mad_per_patient_df$Row.names = row.names(mad_per_patient_df)
  return(mad_per_patient_df)
}
