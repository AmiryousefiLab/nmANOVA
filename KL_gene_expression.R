library( readxl )
library( pheatmap )
library( maditr )
library( dplyr )
library( bladderbatch ) # BiocManager::install("bladderbatch")
library( ggplot2 )
library( ggrepel ) 

data( bladderdata )

pheno <- pData(bladderEset)
edata <- exprs(bladderEset)

#now load Kullback Leibler divergence matrix
#the code to obtain KL_data.RData has been executed using CSC clusters - full code for CSC environment upon request
load("KL_data.RData")
bl <- bladderKL
colnames( bl )[18:29] <- paste( 'mTCC', seq(1:12), sep = '_' )
colnames( bl )[30:41] <- paste( 'sTCC+CIS', seq(1:12), sep = '_' )
colnames( bl )[42:57] <- paste( 'sTCC-CIS', seq(1:16), sep = '_' )

#Figure S2
pheatmap( log( bl + 0.01), labels_row = sub('_[^_]*$', '', colnames( bl )), labels_col =  sub('_[^_]*$', '', colnames( bl )), fontsize_row = 12, fontsize_col = 13 )


#for PCA
prc <- prcomp( t( edata ), center = T, scale = T )
prc_res <- as.data.frame( prc$x )
prc_res <- merge( prc_res, pheno, by = 'row.names' )

#Figure 2C 
ggplot( prc_res , aes( x =  PC1, y = PC2, col = as.factor( outcome ) ) ) + geom_point( size = 10, alpha = 0.7 ) +
  theme(panel.grid.minor = element_line(colour = "grey", size = 0.1), 
        panel.grid.major = element_line(colour = "grey", size = 0.5)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_bw()+ xlab( "PC1" ) + ylab( "PC2" ) + guides(col=guide_legend(title="Sample"), shape = guide_legend(title='concentrations'))+ theme(plot.title = element_text(size = 25),axis.text.x = element_text(angle = 0, hjust = 1, size = 25), axis.title=element_text(size=25), axis.text.y = element_text(angle = 0, hjust = 1, size = 25), legend.title=element_text(size=25), legend.text=element_text(size=25) )+ annotate("text", x=-15, y=95, label= expression( paste( italic( "p-value = 0.041" )  ) ), size = 10) 


#collect results for all the partitions in nm_df 
#change the column/row names
better_names <- c( paste0( 'Normal', '_', c( 1:8 ) ), paste0( 'Biopsy', '_', c( 1:9 ) ), paste0( 'mTCC', '_', c( 1:12 ) ), paste0( 'sTCC+CIS', '_', c( 1:12 ) ),  paste0( 'sTCC-CIS', '_', c( 1:16 ) ) ) 
colnames( bladderKL ) <- better_names
rownames( bladderKL ) <- better_names
diss_df <- bladderKL

#select partitions 
partition <- list()
partition[[1]] <-c( "Normal"   ,         "Biopsy" )
partition[[2]] <-c( "Normal"   ,         "mTCC" )
partition[[3]] <-c( "Normal"   ,         "sTCC+CIS" )
partition[[4]] <-c( "Normal"   ,         "sTCC-CIS" )
partition[[5]] <-c( "Biopsy"   ,         "mTCC" )
partition[[6]] <-c( "Biopsy"   ,         "sTCC+CIS" )
partition[[7]] <-c( "Biopsy"   ,         "sTCC-CIS" )
partition[[8]] <-c( "mTCC"   ,         "sTCC+CIS" )
partition[[9]] <-c( "mTCC"   ,         "sTCC-CIS" )
partition[[10]] <-c( "sTCC+CIS"   ,         "sTCC-CIS" )

partition[[11]] <- c( "mTCC"  ,     "sTCC+CIS" ,  "sTCC-CIS" )

partition[[12]] <- c( "Normal"   , "mTCC"  ,     "sTCC+CIS" ,  "sTCC-CIS"   )
partition[[13]] <- c( "Biopsy"   , "mTCC"  ,     "sTCC+CIS" ,  "sTCC-CIS"   )
partition[[14]] <- c( "Normal"   ,      "mTCC+sTCC+CIS+sTCC-CIS" ) 
partition[[15]] <- c( "Biopsy"   ,       "mTCC+sTCC+CIS+sTCC-CIS" )

partition[[16]] <- c( "Normal"   ,         "Biopsy" ,   "sTCC+CIS" )
partition[[17]] <- c( "Normal+Biopsy"   , "sTCC+CIS" )

partition[[18]] <- c( "Normal",         "Biopsy" ,   "sTCC-CIS" )
partition[[19]] <- c( "Normal+Biopsy", "sTCC-CIS" )

partition[[20]] <- c( "Normal",         "Biopsy" ,   "mTCC" )
partition[[21]] <- c( "Normal+Biopsy", "mTCC" )


partition[[22]] <- c( "Normal+Biopsy"  , "mTCC"  ,     "sTCC+CIS" ,  "sTCC-CIS" )
partition[[23]] <- c( "Normal+Biopsy"  ,"mTCC+sTCC+CIS+sTCC-CIS" )

partition[[24]] <- c( "Biopsy"   ,"Normal", "mTCC"  ,     "sTCC+CIS" ,  "sTCC-CIS"   )



d <- as.matrix( diss_df )


joints <- c( "Normal+Biopsy", "mTCC+sTCC+CIS+sTCC-CIS" )
ind_Norm_Biopsy <-  which( sub('_[^_]*$', '', colnames( d ) ) %in% c( 'Normal', 'Biopsy' ) )
ind_Cancers <-  which( sub('_[^_]*$', '', colnames( d ) ) %in% c( "mTCC"  ,     "sTCC+CIS" ,  "sTCC-CIS" ) )


#run npANOVA 

res_df <- data.frame( N_partitions = rep( NA_character_, 24 ),
                      partition = rep( NA_character_, 24 ),
                      n_partition = rep( NA_character_, 24 ),
                      F_stat = rep( NA_real_, 24 ),
                      p_value = rep( NA_real_, 24 ) )

res_m <- list()

for( i in 1:length( partition ) ){
 
  
  if( length( intersect( partition[[i]], joints ) ) > 0 ){
    
    sel_partition <- c( rep('Normal', 8), rep('Biopsy', 9), rep("mTCC" , 12 ), rep("sTCC+CIS" , 12 ), rep("sTCC-CIS" , 16 ) )
    l <- unique( partition[[i]] )
    
    sub1 <- vector()
    sub2 <- vector()
    
    if( "Normal+Biopsy" %in% l )  { 
      sel_partition[ ind_Norm_Biopsy ] <- 'Normal+Biopsy' 
      sub1 <- union( l, c( 'Normal', 'Biopsy') )
    }
    if( "mTCC+sTCC+CIS+sTCC-CIS" %in% l )  { 
      sel_partition[ ind_Cancers ] <- "mTCC+sTCC+CIS+sTCC-CIS" 
      sub2 <- union( l, c( 'mTCC', 'sTCC+CIS', 'sTCC-CIS') )
    }
    
    sub <- union( sub1 , sub2 )
    subpartition <- partition[[i]]
    res_m[[i]] <- nmANOVA(d, sel_partition,  10,  subpartition )
    res_m[[i]] <- res_m[[i]]$summary
    
    subpartition2 <- paste0( sel_partition, '_', seq(1, length( sel_partition ) ) )
    
    
    dd <- as.data.frame( d )
    colnames( dd ) <- subpartition2
    rownames( dd ) <- subpartition2
    
    pq_vec <- vector()
    cur_partition <- vector() 
    n_partition <- vector()
    
    for( j in 1:length( partition[[i]] ) ){
      
      pq <- which( grepl( partition[[i]][j] , colnames( dd ), fixed = TRUE ) == TRUE )
      pq_vec <- c( pq_vec, pq )
      cur_partition <- c( cur_partition, rep( partition[[i]][j], length( pq ) ) )
      n_partition <- c( n_partition, length( pq ) )
      
    }
    
    dd <- dd[pq_vec, pq_vec]
    
    ad_res <- adonis2( as.dist( dd ) ~ cur_partition )
    res_df$p_value[i] <- as.numeric( ad_res$`Pr(>F)` )[1]
    res_df$F_stat[i] <- as.numeric( ad_res$F )[1]
    cur_partition <- as.factor( cur_partition )
    res_df$N_partitions[i] = length( levels( cur_partition ) )
    res_df$partition[i] = paste( levels( cur_partition ) , collapse = ',')
    res_df$n_partition[i] = paste( n_partition , collapse= ','  )
  
  
    
    
  } 
  else{
    
    dd <- as.data.frame( d )
    
    pq_vec <- vector()
    cur_partition <- vector() 
    n_partition <- vector()
    
    for( j in 1:length( partition[[i]] ) ){
      
      pq <- which( grepl( partition[[i]][j] , colnames( dd ), fixed = TRUE  ) == TRUE )
      pq_vec <- c( pq_vec, pq )
      cur_partition <- c( cur_partition, rep( partition[[i]][j], length( pq ) ) )
      n_partition <- c( n_partition, length( pq ) )
      
    }
    
    dd <- dd[pq_vec, pq_vec]
  
    ad_res <- adonis2( as.dist( dd ) ~ cur_partition )
    res_df$p_value[i] <- as.numeric( ad_res$`Pr(>F)` )[1]
    res_df$F_stat[i] <- as.numeric( ad_res$F )[1]
    cur_partition <- as.factor( cur_partition )
    res_df$N_partitions[i] = length( levels( cur_partition ) )
    res_df$partition[i] = paste(  levels( cur_partition ) , collapse = ',')
    res_df$n_partition[i] = paste( n_partition , collapse= ','  )
    
    
  }
 
  

  
}  

#final data frame with p-values, data for Table S3 has been taken from nm_df
np_df <-  res_df


load("/Users/amalyuti/Desktop/nmANOVA/revision_codes/KL_divergence/JSD_np_results.RData")
load("/Users/amalyuti/Desktop/nmANOVA/revision_codes/KL_divergence/KL_results.RData")

#Figure S3
options(ggrepel.max.overlaps = Inf)
com_df <- 
  np_df %>%
  left_join( nm_df, by = c( 'partition' ), suffix = c( '.np', '.nm' ) )

ggplot( com_df, aes( x =  log10( p_value.np ), y = log10( p_value.nm ) ) ) +  
  geom_point( alpha = 0.7,  size = 10, col = 'darkorchid')  +
  geom_hline( yintercept = log10( 0.05) ) + geom_vline( xintercept = log10( 0.05 ) ) +
  
  theme(panel.grid.minor = element_line(colour = "grey", size = 0.1), 
        panel.grid.major = element_line(colour = "grey", size = 0.5)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))  +
  theme_bw()+ 
  theme(plot.title = element_text(size = 25),axis.text.x = element_text(angle = 0, hjust = 1, size = 25), axis.title=element_text(size=28), axis.text.y = element_text(angle = 0, hjust = 1, size =25), legend.title=element_text(size=28), legend.text=element_text(size=28) ) +
  guides( col=guide_legend(title="Species"), shape = guide_legend(title='Species')) + geom_label_repel(  aes(  log10( p_value.np ),   log10( p_value.nm ) , label = partition  ), size = 3, data = com_df) 






