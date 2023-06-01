#nmANOVA function
nmANOVA <- function( diss_data, partition, N_sampling = 1, subpartition = NULL ) {
  
  #diss_data - a numeric matrix or data.frame with dissimilarity measures,columns and rows should be in the same order
  #partition - a vector specifying the groups, it should be same length and in the same order as diss_data columns/rows 
  #N_sampling - number of proportional samplings to perform 
  #subpartition - a vector of groups present in partition for which the method has to be applied to 
  #if subpartition is NULL, the method will be applied to all the groups present in partition ( unique( partition ) )
  
  if ( nrow( diss_data ) != ncol( diss_data ) )
    stop( "matrix has to be square" )
  
  if ( nrow( diss_data ) != length( partition ) )
    stop( "partition has to be same length as the number of row/columns in the matrix" )
  
  if ( N_sampling <= 0 )
    stop( "N_sampling has to be a positive integer" )
  
  if ( !is.null( subpartition ) & length( unique( subpartition ) ) < 2 )
    stop( "there should be at least two unique subpartitions" )
  
  
  
  partition <- as.factor( as.vector( partition ) )
  partition_sorted <- partition[ order( partition ) ]
  
  inner_partition_names <- vector()
  for ( i in 1:length( levels( partition_sorted ) ) ){
    inner_partition_names <- c( inner_partition_names, paste0( levels( partition_sorted )[i],'_', seq( 1:table( partition_sorted )[i] ) ) )
  }
  
  #set the names for each column/row based on a partition they belong to
  inner_partition_names <- inner_partition_names[ order( order ( partition ) ) ]
  partition_df <- as.data.frame( diss_data )
  colnames( partition_df ) <- inner_partition_names
  rownames( partition_df ) <- inner_partition_names
  
  #check if there is a subpartition provided  
  if ( !is.null( subpartition ) ){
    if ( length( unique( intersect( partition, subpartition ) ) ) !=  length( unique( subpartition )  )  )
      stop( "some elements from the subpartition do not match with partition" )
    subpartition <- as.factor( as.vector( subpartition ) )
    partition <- partition[ which( partition %in% intersect( partition, subpartition ) ) ]
  }
  
  selected_partitions <- which( sub('_[^_]*$', '', colnames( partition_df ) ) %in% partition )
  selected_partitions <- which( sub('_[^_]*$', '', rownames( partition_df ) ) %in% partition )
  partition_df <- partition_df[ selected_partitions, selected_partitions ]
  partition_df <- partition_df[ order( row.names( partition_df ) ), order( names( partition_df ) ) ]
  
  #the actual partition
  partition <- as.factor( sub('_[^_]*$', '', colnames( partition_df ) ) )
  #number of elements in each partition
  n_partition <- as.vector( table( partition ) )
  
  #calculate delta_w and delta_jj
  vec_within <- vector()
  delta_jj <- vector()
  SNM_w <- 0
  for( i in 1:length( levels( partition ) ) ){
    
    ind_partition <- which(  sub('_[^_]*$', '', colnames( partition_df ) ) == levels( partition )[i] )
    cur_partition <-  as.numeric( as.matrix( partition_df[ ind_partition, ind_partition ] ) )
    vec_within <- c( vec_within, cur_partition)
    delta_jj[i] <- mean( cur_partition, na.rm = TRUE )
    
  }
  
  delta_w <- mean( vec_within, na.rm = TRUE )
  
  #calculate SNM_w
  SNM_w <- 0
  for( i in 1:length( levels( partition ) ) ){
    
    ind_partition <- which(  sub('_[^_]*$', '', colnames( partition_df ) ) == levels( partition )[i] )
    cur_partition <-  as.numeric( as.matrix( partition_df[ ind_partition, ind_partition ] ) )
    
    n_cur_partition <- sum( is.na( cur_partition ) == FALSE )
    
    
    SNM_w <- SNM_w + n_cur_partition * (  ( delta_jj[i] - delta_w ) )^2
  }
  
  
  #calculate delta_jj', delta_prop and SNM_b
  
  off_diag_new <- 0
  for( i in 1:length( levels( partition ) ) ){
    for( j in 1:length( levels( partition ) ) ){
      if( i != j ){ 
        
        ind_row_partition <- which(  sub('_[^_]*$', '', colnames( partition_df ) ) == levels( partition )[i] )
        ind_col_partition <- which(  sub('_[^_]*$', '', colnames( partition_df ) ) == levels( partition )[j] )
        cur_partition <-  as.numeric( as.matrix( partition_df[ ind_row_partition, ind_col_partition ] ) )
        
        off_diag_new <- off_diag_new +  sum( is.na( cur_partition ) == FALSE )
        
        
        
      }
    }
  }
  
  
  num_off_diag_groups <- length( levels( partition ) )^2 - length( levels( partition ) )
  num_off_diag <- off_diag_new
  
  F_stat_array <- vector()
  prop_sampl_df_list <- list()
  
  for( ii in 1:N_sampling ){
    
    diagonal_sampling <- vec_within[!is.na(vec_within)]
    
    
    prop_sampl_df <- data.frame( N_sampling_num = rep( ii, num_off_diag_groups ), #which N_sampling is performed
                                 sampling_num = numeric( num_off_diag_groups ),#how many elements are sampled
                                 partitions = character( num_off_diag_groups ),# for what between group partition the sampling is performed
                                 part_mean = numeric( num_off_diag_groups ),# mean of the current between group partition
                                 prop_mean = numeric( num_off_diag_groups ),# mean of the sampled elements  
                                 sampling_ind = character( num_off_diag_groups), #indexes for the sampled elements from the vector of the whole
                                 sampled_from = character(num_off_diag_groups), 
                                 stringsAsFactors=FALSE )  
    
    non_na_vec_within <- sum( is.na( vec_within ) == FALSE )
    
    
    k <- 0
    delta_jj_prime <- vector()
    SNM_b <- 0
    SNM_b_array <- vector()
    for( i in 1:length( levels( partition ) ) ){
      for( j in 1:length( levels( partition ) ) ){
        if( i != j ){ 
          
          ind_row_partition <- which(  sub('_[^_]*$', '', colnames( partition_df ) ) == levels( partition )[i] )
          ind_col_partition <- which(  sub('_[^_]*$', '', colnames( partition_df ) ) == levels( partition )[j] )
          cur_partition <-  as.numeric( as.matrix( partition_df[ ind_row_partition, ind_col_partition ] ) )
          non_na_cur_partition <-  sum( is.na( cur_partition ) == FALSE )
          delta_jj_prime <- mean( cur_partition,  na.rm = TRUE )
          
          prop_num <- floor( non_na_cur_partition  *  ( non_na_vec_within  / as.numeric( num_off_diag ) )  ) 
          
          if( prop_num == 0 ){
            print( 'prop_num = 0 -> changed to 1' )
            prop_num <- 1 }
          
          sample_ind <- sample( c( 1:length( diagonal_sampling ) ), prop_num, replace = FALSE )
          k <- k + 1
          prop_sampl_df$sampling_num[k] <- prop_num
          prop_sampl_df$partitions[k] <- paste0( levels( partition )[i], ',', levels( partition )[j] )
          prop_sampl_df$sampling_ind[k] <- paste0( sample_ind, collapse = ',' )
          prop_sampl_df$sampled_from[k] <- paste0( round( diagonal_sampling, 3) , collapse = ',' )
          delta_prop <- mean( diagonal_sampling[ sample_ind ],  na.rm = TRUE )
          
          
          if( floor( non_na_cur_partition  * ( non_na_vec_within  / as.numeric( num_off_diag ) )  )  != 0 ){ diagonal_sampling <- diagonal_sampling[ -sample_ind ] }
          prop_sampl_df$prop_mean[k] <- delta_prop
          
          
          delta_jj_prime <- mean( cur_partition,  na.rm = TRUE )
          
          prop_sampl_df$part_mean[k] <- delta_jj_prime
          SNM_b <- SNM_b + ( ( delta_jj_prime - delta_prop )/( sqrt( ( 1/( non_na_cur_partition ) + 1/prop_num ) ) ) )^2
          
        }
      }
    }
    
    #calculate statistics
    prop_sampl_df_list[[ii]] <- prop_sampl_df 
    SNM_b_array[ii] <- SNM_b
    F_stat_array[ii] <- ( SNM_b/length( levels( partition ) ) ) / SNM_w
    
  }
  
  
  
  F_stat <- mean( F_stat_array, na.rm = TRUE )
  
  df1 <- length( levels( partition ) )^2 - length( levels( partition ) )
  df2 <-  length( levels( partition ) ) - 1
  p_value <- pf( F_stat, df1, df2, lower.tail = F ) 
  partition_df1 <- partition_df
  
  
  res <- list()
  
  res[[1]] <- data.frame( N_partitions = length( levels( partition ) ), partition = paste(  levels( partition ) , collapse = ','), n_partition = paste( n_partition , collapse= ','  ), delta_w = delta_w,  delta_jj = paste( delta_jj, collapse= ',' ),  SNM_w = SNM_w,  SNM_b = mean( SNM_b_array, na.rm = TRUE ), F_stat =  F_stat, p_value = p_value  )
  
  res[[2]] <- do.call( rbind.data.frame, prop_sampl_df_list )
  
  names( res ) <- c( 'summary', 'proportional_sampling' ) 
  
  res
  
}

library( readxl )
library( pheatmap )
library( pvclust )
library( maditr )
library( dplyr )
library( plyr )

oneil <- read_excel("drug_combination_sensitivity_data.xlsx" )


all_cl <- unique( oneil$cell_line_name )
partition_df_list <- list()
all_res <- list()
k <- 0

#YOU DO NOT NEED TO RUN THIS CODE TO REPRODUCE THE FIGURE S1 (FIGURE 4B has 2 heatmaps with smallest and highest p-values which are visualized at figure S1)
###################

for (ii in 1:length( all_cl ) ){
  
  
  k <- k+1
  cur_cl <- all_cl[ii]
  oneil_data <- oneil[oneil$cell_line_name==cur_cl,]
  
  oneil_data_single = ddply(oneil_data, .(drug_row,drug_col), summarise, 
                            css_row = mean(css_row), css_col = mean(css_col),
                            ri_row = mean(ri_row), ri_col = mean(ri_col))
  # combo
  oneil_data_1 = oneil_data_single[,c('drug_row','drug_col','css_row')]
  oneil_data_2 = oneil_data_single[,c('drug_col','drug_row','css_col')]
  
  # single drug
  oneil_data_3 = oneil_data_single[,c('drug_col','drug_row','ri_row')]
  oneil_data_4 = oneil_data_single[,c('drug_col','drug_row','ri_col')]
  oneil_data_3$drug_col=oneil_data_3$drug_row
  oneil_data_4$drug_row=oneil_data_4$drug_col
  colnames(oneil_data_2) = colnames(oneil_data_1)
  colnames(oneil_data_3) = colnames(oneil_data_1)
  colnames(oneil_data_4) = colnames(oneil_data_1)
  oneil_data_all = rbind(oneil_data_1, oneil_data_2, oneil_data_3, oneil_data_4)
  # css_mat: asymetric
  
  #ri values on the diagonal
  css_mat = dcast(oneil_data_all, drug_row ~ drug_col, value.var = "css_row", fun.aggregate = mean)
  rownames(css_mat) = css_mat$drug_row
  css_mat = css_mat[,-1]
  
  
  #impute NA with the average of RI of single drugs
  na_index = which(is.na(css_mat)==T, arr.ind = T)
  for (i in 1:nrow(na_index)){
    css_mat[na_index[i,1], na_index[i,2]] = (css_mat[na_index[i,1],na_index[i,1]] + css_mat[na_index[i,2],na_index[i,2]])/2 
  }
  
  
  # to be used by hclust, force it to be symmetric by taking lower triangle
  # css_mat2 = 1 - abs(cor(css_mat)) # use correlation as distance
  css_mat2 =  1 - abs(cor(css_mat))
  css_mat2 = 1 - abs(cor(css_mat)) # use correlation as distance, column wise
  css_mat3 = 1 - abs(cor(t(css_mat))) # row wise, you have two versions of correlation matrix
  # then merge these two correlation matrix, lower diagnoal as css_mat2, upper diagonal as css_mat3
  css_mat4 = css_mat2
  css_mat4[upper.tri(css_mat4)] <- t(css_mat3)[upper.tri(css_mat3)] # take lower triangle
  clusters <- hclust(as.dist(css_mat2))
  clusterCut <- cutree(clusters, 2)
  
  
  res.pv <- pvclust(css_mat, nboot=1000,  parallel=FALSE)
  
  cl_res <- pvpick( res.pv, alpha=0.99 )
  
  partition_df <- 
    data.frame( drug = toupper( colnames( css_mat4 ) ) )%>%
    mutate( partition = 0 , cell_line = all_cl[ii] )
  
 
  for ( i in 1:length( cl_res$clusters  ) ){
    
 
  ind <- which( partition_df$drug %in%  toupper(  cl_res$clusters[[i]] ) )
  partition_df$partition[ind] <- i
  
  
  }
  
  
partition_df_list[[ii]] <- partition_df


}



#select significant partitions with cluster size > 3 
check_partitions <- function( df ){
  
  
  df_stat <- 
    df%>%
    filter( partition != 0 )%>%
    group_by( partition )%>%
    dplyr::summarise( n = n() )%>%
    filter( n > 3 )

  
  if ( nrow( df_stat ) == 1 ){
    
    df$nmA_analysis <- 'no'
    
  }else {
    
    
    df <- 
      df %>%
      mutate( nmA_analysis = ifelse(  partition %in% df_stat$partition , 'yes', nmA_analysis) )
    
    
  }
  
  df
  
}

partition_res <- 
  ldply( partition_df_list, data.frame)%>%
  mutate( nmA_analysis = 'no' )%>%
  group_by( cell_line )%>%
  do( { check_partitions(.) } )

#select significant partitions
partition_res_sel <- 
  partition_res %>%
  filter( nmA_analysis == 'yes' )

#save( partition_res_sel, file = 'pvclust_significant_partitions.RData')

load( "pvclust_significant_partitions.RData" )

all_cl <- unique( partition_res_sel$cell_line )
all_res <- list()


#check the other way around - select the clusters significant by the function
for (ii in 1:length( all_cl ) ){

  cur_cl <- all_cl[ii]
  oneil_data <- oneil[oneil$cell_line_name==cur_cl,]
  
  oneil_data_single = ddply(oneil_data, .(drug_row,drug_col), summarise, 
                            css_row = mean(css_row), css_col = mean(css_col),
                            ri_row = mean(ri_row), ri_col = mean(ri_col))
  # combo
  oneil_data_1 = oneil_data_single[,c('drug_row','drug_col','css_row')]
  oneil_data_2 = oneil_data_single[,c('drug_col','drug_row','css_col')]
  
  # single drug
  oneil_data_3 = oneil_data_single[,c('drug_col','drug_row','ri_row')]
  oneil_data_4 = oneil_data_single[,c('drug_col','drug_row','ri_col')]
  oneil_data_3$drug_col=oneil_data_3$drug_row
  oneil_data_4$drug_row=oneil_data_4$drug_col
  colnames(oneil_data_2) = colnames(oneil_data_1)
  colnames(oneil_data_3) = colnames(oneil_data_1)
  colnames(oneil_data_4) = colnames(oneil_data_1)
  oneil_data_all = rbind(oneil_data_1, oneil_data_2, oneil_data_3, oneil_data_4)
  # css_mat: asymetric
  
  #ri values on the diagonal
  css_mat = dcast(oneil_data_all, drug_row ~ drug_col, value.var = "css_row", fun.aggregate = mean)
  rownames(css_mat) = css_mat$drug_row
  css_mat = css_mat[,-1]
  
  
  #impute NA with the average of RI of single drugs
  na_index = which(is.na(css_mat)==T, arr.ind = T)
  for (i in 1:nrow(na_index)){
    css_mat[na_index[i,1], na_index[i,2]] = (css_mat[na_index[i,1],na_index[i,1]] + css_mat[na_index[i,2],na_index[i,2]])/2 
  }
  
  css_mat2 = 1 - abs(cor(css_mat)) # use correlation as distance, column wise
  css_mat3 = 1 - abs(cor(t(css_mat))) # row wise, you have two versions of correlation matrix
  # then merge these two correlation matrix, lower diagnoal as css_mat2, upper diagonal as css_mat3
  css_mat4 = css_mat2
  css_mat4[upper.tri(css_mat4)] <- t(css_mat3)[upper.tri(css_mat3)] # take lower triangle
  
  
  
  # run nmanova
  
  partition <- vector()
  dr <- toupper( colnames( css_mat4 ) )
  
  cur_part <- 
    partition_res_sel %>%
    filter( cell_line == cur_cl )
  
  for( j in 1:length( dr ) ){
    
    ind <- which( cur_part$drug == dr[j] )
    
    if( length( ind ) > 0 ){
      
      partition[j] <- paste0( 'partition', '_', cur_part$partition[ind] )
    }else{
      
      partition[j] <- 'no_partition'
      
    }
    
  }
    
    ind <- which( partition ==  'no_partition' )
    if (length( ind ) > 0 ){
      
      css_mat4 <- css_mat4[-ind, -ind]
      partition <- partition[-ind]
    }
    
    colnames( css_mat4 ) <- toupper( colnames( css_mat4 ) )
    ind <- which( colnames( css_mat4 ) == toupper( '7-Ethyl-10-hydroxycamptothecin' ) )
    colnames( css_mat4 )[ind ] <- 'SN-38'
    
    
    rownames( css_mat4 ) <- toupper( rownames( css_mat4 ) )
    ind <- which( rownames( css_mat4 ) == toupper( '7-Ethyl-10-hydroxycamptothecin' ) )
    rownames( css_mat4 )[ind ] <- 'SN-38'
    
  res <- nmANOVA( ( css_mat4 ), partition, 10 )
  
  all_res[[ii]] <- res$summary
  all_res[[ii]]$cell_line <- cur_cl
  
  ind <- which( partition == partition[1] ) 
  
  css_mat_pheatmap <- css_mat4[c(ind, setdiff( 1:nrow( css_mat4 ), ind ) ), c(ind, setdiff( 1:nrow( css_mat4 ), ind ) )]
 
  pheatmap( css_mat_pheatmap, cluster_rows = FALSE, cluster_cols =FALSE, main = paste0( cur_cl, ', p-value = ',  round( res$summary$p_value, 4 ) ), fontsize_row = 7, fontsize_col = 7 )
  
    
  
}
#get all the results in one data frame; data for Table 2 has been taken from  res
final_res <- ldply( all_res, data.frame ) 




