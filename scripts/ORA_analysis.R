library(clusterProfiler)
#BiocManager::install("enrichplot")
library(enrichplot)
library(ReactomePA)
organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
library(DOSE)
library(ggplot2)
##https://yulab-smu.top/biomedical-knowledge-mining-book/025-do-enrichment.html

gene_list <- "/home/yolanda/tblab/yolanda/GLOWgenes/berta/networks_exploration/networks_exploration/louvain_clustering/phenotype_clusteing.csv"
gene_list <- "/home/yolanda/tblab/yolanda/GLOWgenes/berta/networks_exploration/networks_exploration/louvain_clustering/coessenciality_clusteing.csv"
gene_list <- "/home/yolanda/tblab/yolanda/GLOWgenes/berta/networks_exploration/networks_exploration/louvain_clustering/genetics_clusteing.csv"
gene_list <- "/home/yolanda/tblab/yolanda/GLOWgenes/berta/networks_exploration/networks_exploration/louvain_clustering/phenotype_clusteing_weights.csv"
# gene_list <- "/home/yolanda/tblab/yolanda/GLOWgenes/berta/networks_exploration/networks_exploration/louvain_clustering/phenotype_subclusteing_weights_cluster7.csv"
# gene_list <- "/home/yolanda/tblab/yolanda/GLOWgenes/berta/networks_exploration/networks_exploration/louvain_clustering/phenotype_subclusteing_weights_subcluster_11.csv"
# gene_list <- "/home/yolanda/tblab/yolanda/GLOWgenes/berta/networks_exploration/networks_exploration/louvain_clustering/phenotype_subclusteing_weights_subcluster_8.csv"
gene_list <- "/home/yolanda/tblab/yolanda/GLOWgenes/berta/networks_exploration/networks_exploration/louvain_clustering/leiden_clustering.csv"
gene_list <- "/home/yolanda/tblab/yolanda/GLOWgenes/berta/networks_exploration/networks_exploration/louvain_clustering/leiden_clustering_no_weight.csv"
# gene_list <- "/home/yolanda/tblab/yolanda/GLOWgenes/berta/networks_exploration/networks_exploration/louvain_clustering/mcl_clustering_2_0.csv"
# gene_list <- "/home/yolanda/tblab/yolanda/GLOWgenes/berta/networks_exploration/networks_exploration/louvain_clustering/phenotype_clusteing_NO_weights.csv"


#background <- "/home/yolanda/tblab/yolanda/GLOWgenes/berta/background_genes.csv"

# reading in input
df_genes = read.csv(gene_list, header=TRUE)

# si quiero todos los genes de la red de fenotipo para el background
#df_bg_phenotype = read.csv("/home/yolanda/tblab/yolanda/GLOWgenes/berta/networks_exploration/networks_exploration/louvain_clustering/phenotype_clusteing.csv", header=FALSE)
#gb_gene_list = df_bg_phenotype[["V1"]]

gb_gene_list = df_genes$gen

cluster_table <- as.data.frame(table(df_genes[2]))

for (cluster_id in cluster_table[cluster_table["Freq"] > 10, 1]) {
  
  gene_list <- df_genes[df_genes[2] == cluster_id,][["gen"]] # by cluster_* column
  cat("############### \n")
  cat("Cluster", cluster_id, "contains", length(gene_list), "genes \n")
  
  go_enrich <- enrichGO(gene = gene_list,
                        universe = gb_gene_list,
                        OrgDb = organism, 
                        keyType = 'SYMBOL',
                        readable = T,
                        ont = "BP",
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.1)
  
  
  ## términos significativos
  res_go <- go_enrich@result
  sig_go <- res_go[res_go$p.adjust < 0.05, ]
  sig_go
  
  cat(c("number of GO-BP significant terms:", as.character(nrow(sig_go)), "\n"))
  if (nrow(sig_go) > 0){
    b <- barplot(go_enrich, 
                 drop = TRUE, 
                 showCategory = 20, 
                 title = paste("GO Biological Pathways barplot cluster:", cluster_id),
                 font.size = 8)
    #print(b)
    
    d <- dotplot(go_enrich, showCategory = 20,
                 title = paste("GO Biological Pathways dotplot cluster:", cluster_id)
    )
    #print(d)
    
    u <- upsetplot(go_enrich,
                   title = paste("GO Biological Pathways upsetplot cluster:", cluster_id)
    )
    print(u)
    
    print(b+d)
    
    
    go_enrich <- pairwise_termsim(go_enrich) # pairwise_termsim() calculates semantic similarity between enriched GO terms, producing a term–term similarity matrix, which emapplot uses to draw networks.
    enrichplot::emapplot(go_enrich) # A network graph of enriched GO terms based on semantic similarity.
    emaplot_50 <- enrichplot::emapplot(go_enrich, 
                                       showCategory = 50) # A network graph of enriched GO terms based on semantic similarity.
    print(emaplot_50)
    # si aumento el número (por defoult es 30), se me juntan los clusters
    
    go_simple <- clusterProfiler::simplify(
      go_enrich,
      cutoff = 0.7,       # similarity threshold for merging redundant terms
      by = "p.adjust",
      select_fun = min
    ) # This merges highly redundant GO terms
    
    cat(c("number of significant simplified terms:", as.character(nrow(go_simple@result)), "\n"))
    
    emapplot_50_simpl <- emapplot(go_simple, showCategory = 50, 
                                  layout = "fr")
    print(emapplot_50_simpl)
    
    emapplot(go_simple, layout = "fr")
    
  }
  
  ## REACTOME: 
  
  entrez_gene_list <- as.data.frame(mapIds(org.Hs.eg.db, keys = gene_list,
                                           column = "ENTREZID", keytype = "SYMBOL"))
  
  entrez_bg_gene_list <- as.data.frame(mapIds(org.Hs.eg.db, keys = gb_gene_list,
                                              column = "ENTREZID", keytype = "SYMBOL"))
  
  reactome_enrich <- enrichPathway(gene= entrez_gene_list[[1]], 
                                   universe =  entrez_bg_gene_list[[1]] ,
                                   pvalueCutoff = 0.05, 
                                   readable=TRUE)
  head(reactome_enrich)
  
  res_react <- reactome_enrich@result
  sig_react <- res_react[res_react$p.adjust < 0.05, ]
  sig_react
  
  cat(c("number of significant reactome terms:", as.character(nrow(sig_react)), "\n"))
  
  b <- barplot(reactome_enrich, 
               drop = TRUE, 
               showCategory = 20, 
               title = paste("REACTOME Pathways barplot cluster:", cluster_id),
               font.size = 8)
  
  d <- dotplot(reactome_enrich, showCategory = 20,
               title = paste("REACTOME Pathways dotplot cluster:", cluster_id)
  )
   
  print(b+d)
  
  
  
  ### DISEASE
  
  hdo_enrich <- enrichDO(gene          = entrez_gene_list[[1]],
                ont           = "HDO",
                pvalueCutoff  = 0.05,
                pAdjustMethod = "BH",
                universe      = entrez_bg_gene_list[[1]],
                minGSSize     = 5,
                maxGSSize     = 500,
                qvalueCutoff  = 0.05,
                readable      = FALSE)
  head(hdo_enrich)
  
  res_hdo <- hdo_enrich@result
  sig_hdo <- res_hdo[res_hdo$p.adjust < 0.05, ]
  sig_hdo
  
  cat(c("number of significant hdo terms:", as.character(nrow(sig_hdo)), "\n"))
  
  b <- barplot(hdo_enrich, 
               drop = TRUE, 
               showCategory = 20, 
               title = paste("HDO Pathways barplot cluster:", cluster_id),
               font.size = 8)
  
  d <- dotplot(hdo_enrich, showCategory = 20,
               title = paste("hdo Pathways dotplot cluster:", cluster_id)
  )
  
  print(b+d)
  
  hpo_enrich <- enrichDO(gene          = entrez_gene_list[[1]],
                         ont           = "HPO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         universe      = entrez_bg_gene_list[[1]],
                         minGSSize     = 5,
                         maxGSSize     = 500,
                         qvalueCutoff  = 0.05,
                         readable      = FALSE)
  head(hpo_enrich)
  
  res_hpo <- hpo_enrich@result
  sig_hpo <- res_hpo[res_hpo$p.adjust < 0.05, ]
  sig_hpo
  
  cat(c("number of significant hpo terms:", as.character(nrow(sig_hpo)), "\n"))
  
  b <- barplot(hpo_enrich, 
               drop = TRUE, 
               showCategory = 20, 
               title = paste("HPO Pathways barplot cluster:", cluster_id),
               font.size = 8)
  
  d <- dotplot(hpo_enrich, showCategory = 20,
               title = paste("hpo Pathways dotplot cluster:", cluster_id)
  )
  
  print(b+d)
  
  dgn_enrich <- enrichDGN(gene          = entrez_gene_list[[1]],
                         pvalueCutoff  = 0.05,
                         universe      = entrez_bg_gene_list[[1]],
                         qvalueCutoff  = 0.05,
                         readable      = FALSE)
  head(dgn_enrich)
  
  res_dng <- dgn_enrich@result
  sig_dng <- res_dng[res_dng$p.adjust < 0.05, ]
  sig_dng
  
  cat(c("number of significant dng terms:", as.character(nrow(sig_dng)), "\n"))
  
  b <- barplot(dgn_enrich, 
               drop = TRUE, 
               showCategory = 20, 
               title = paste("dng Pathways barplot cluster:", cluster_id),
               font.size = 8)
  
  d <- dotplot(dgn_enrich, showCategory = 20,
               title = paste("dng Pathways dotplot cluster:", cluster_id)
  )
  
  print(b+d)
  
  ##### resumen de n significativos de cada enrichment
  
  resume_enrichr <- data.frame(enrich_type=c("go", "reactome", "hdo", "hpo", "dgn"), 
                               n_sig=c(nrow(sig_go), nrow(sig_react), nrow(sig_hdo), nrow(sig_hpo), nrow(sig_dng)))
  
  resume_enrichr$enrich_type <- factor(resume_enrichr$enrich_type, level = resume_enrichr[order(resume_enrichr$n_sig, decreasing = TRUE), ][["enrich_type"]])
  
  p <- ggplot(data=resume_enrichr, aes(x=enrich_type, y=n_sig)) +
    geom_bar(stat="identity", fill="steelblue")+
    theme_minimal()
  
  print(p)

 }



