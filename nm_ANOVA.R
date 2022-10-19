nmANOVA <- function( diss_data, partition, N_sampling = 1, subpartition = NULL ) {
  
  #diss_data - a numeric matrix or data.frame with dissimilarity measures,columns and rows should be in the same order
  #partition - a vector specifying the groups, it should be same length and in the same order as diss_data columns/rows 
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
          
          prop_num <- floor( non_na_cur_partition  *  non_na_vec_within  / as.numeric( num_off_diag )  ) 
          
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
          
          
          if( floor( non_na_cur_partition  *  non_na_vec_within  / as.numeric( num_off_diag )  )  != 0 ){ diagonal_sampling <- diagonal_sampling[ -sample_ind ] }
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


