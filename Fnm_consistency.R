#Codes to replicate figure 1 from the manuscript
#Note that because of random sampling in the nmANOVA algorithm and system-specific numerical instabilities, executing the 
#provided replication materials within a different environment will likely
#produce minor differences from what has been shown in the paper

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
          
          prop_num <- floor( non_na_cur_partition  *  ( non_na_vec_within  / as.numeric( num_off_diag ) ) ) 
          
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
          
          
          if( floor( non_na_cur_partition  * ( non_na_vec_within  / as.numeric( num_off_diag ) ) )  != 0 ){ diagonal_sampling <- diagonal_sampling[ -sample_ind ] }
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


#Figure 1 A

#YOU DO NOT NEED TO RUN THIS CODE TO REPRODUCE THE FIGURE 1 A, AS WE SAVED THE RESULTS IN figure1A_p_values.RData, figure1A_p_values_points.RData ( run from load("figure1A_p_values_points.RData") )
###################

#functions to generate partitions based on selected distributions
generate_diss <- function( distribution, group_N, df_dist_par, show_plots ){
  
  if ( distribution == 'uniform' ) {
    
    dist_function <- function( N, par1, par2 ){
      r <- runif( N, par1, par2 )
      r
    } 
  }
  
  if ( distribution == 'normal' ) {
    
    dist_function <- function( N, mean, sd ){
      r <- rnorm( N, mean, sd )
      r
    } 
  }
  
  if ( distribution == 'gamma' ) {
    
    dist_function <- function( N, shape, rate ){
      r <- rgamma( N, shape, rate )
      r
    } 
  }
  
  if ( distribution == 'chisq' ) {
    
    dist_function <- function( N, df, ignored ){
      r <- rchisq( N, df )
      r
    } 
  }
  
  block_list <- list()
  block_list_row <- list()
  k <- 0
  l <- 0
  for ( i in 1:length( group_N ) ){
    l <- l + 1
    for ( j in 1:length( group_N ) ){
      k <- k + 1
      ind <- which( df_dist_par$groups == paste0( 'group', i, ',', 'group', j ) )
      par <- as.numeric( unlist( strsplit( df_dist_par$dist_par[ind], split=',', fixed=TRUE ) ) )
      block_list[[k]] <- matrix( dist_function( group_N[i]*group_N[j], par[1], par[2] ), ncol = group_N[j] )
      
      
    }
    
    r <- block_list[[1]]
    
    for( m in 2:length( group_N ) ){
      r <- cbind( r, block_list[[m]])
      
    }
    
    block_list_row[[l]] <- r
    k <- 0
    
  }
  
  r <- block_list_row[[1]]
  for( m in 2:length( group_N ) ){
    r <- rbind( r, block_list_row[[m]])
    
  }
  
  if( show_plots == TRUE)  pheatmap(r, cluster_rows = F, cluster_cols = F )
  
  r
  
}
nmANOVA_sim <- function( distribution, group_N, par_vec, H1_group = NULL, show_plots = FALSE, N_sampling = 1  ){
  
  
  
  par1 <- par_vec[1]
  par2 <- par_vec[2]
  
  df_dist_par <- data.frame( groups = character( length( group_N )^2 ), 
                             dist_par = character( length( group_N )^2 ), 
                             stringsAsFactors=FALSE )
  k <- 0
  for ( i in 1:length( group_N ) ){
    for ( j in 1:length( group_N ) ){
      
      k <- k + 1 
      df_dist_par$groups[k] <- paste0( 'group', i, ',', 'group', j )
      
      if( ( !is.null( H1_group ) ) & ( length( which( H1_group$groups ==  df_dist_par$groups[k] ) ) != 0 )  ){
        
        par1_new <- H1_group$dist_par1[ which( H1_group$groups ==  df_dist_par$groups[k] )]
        par2_new <- H1_group$dist_par2[ which( H1_group$groups ==  df_dist_par$groups[k] )]
        df_dist_par$dist_par[k] <- paste0( par1_new,',', par2_new )
        
      } else{
        df_dist_par$dist_par[k] <- paste0( par1,',', par2 )
      }  
      
      
      
    }
    
  }
  
  
  #generate datasets based on selection
  diss_data <- generate_diss( distribution, group_N, df_dist_par, show_plots )
  
  partition <- vector()
  for( i in 1:length( group_N ) ){
    partition <- c( partition, rep( paste0( 'group', i ), group_N[i] ) )
    
  }
  N_sampling <- 1
  res <-list( nmANOVA( diss_data, partition, N_sampling ), diss_data )
  
  
  res
  
  
}
distribution <- 'uniform' #distribution for the simulation matrix
group_N <- c( 10, 10, 10 )#group sizes 
par_vec <- c( 0, 1 )# distribution parameters

#in case you want to introduce the difference - select some groups and change parameters in H1_group
#dist_par1 is a new parameter to use instead of par1 for the specified groups
#dist_par2 is a new parameter to use instead of par2 for the specified groups


#perform simulations and check distribution of p-values
Nsim <- 1
p_value <- vector()
delta <- seq( 0, 10, length.out = 1000 )



for( i in 1:( 3*length( delta ) ) ){
  print(i)
  if ( i <= length( delta ) ){
    
    H1_group <- data.frame( groups = c( 'group2,group1', 'group1,group2' ), dist_par1 = c( 0 + delta[i], 0 + delta[i] ), dist_par2 = c( 1, 1 ) )
    sim_results <- nmANOVA_sim( distribution, group_N, par_vec, H1_group, show_plots = FALSE )
    p_value[i] <- sim_results[[1]]$summary$p_value 
  }
  
  if ( ( i > length( delta ) )  & ( i <= 2*length( delta ) ) ){
    
    H1_group <- data.frame( groups = c( 'group2,group1', 'group1,group2', 'group1,group3', 'group3,group1'  ), dist_par1 = c( 0 + delta[i-length( delta )], 0 + delta[i-length( delta )],  0 + delta[i-length( delta )], 0 + delta[i-length( delta )] ), dist_par2 = c( 1 , 1 , 1 , 1  ) )
    sim_results <- nmANOVA_sim( distribution, group_N, par_vec, H1_group, show_plots = FALSE )
    p_value[i] <- sim_results[[1]]$summary$p_value 
    
    
  }
  
  if ( ( i > 2*length( delta ) )   ){
    
    H1_group <- data.frame( groups = c( 'group2,group1', 'group1,group2', 'group1,group3', 'group3,group1' ,  'group2,group3', 'group3,group2'  ), dist_par1 = c( 0 + delta[i-2*length( delta )], 0 + delta[i-2*length( delta )], 0 + delta[i-2*length( delta )], 0 + delta[i-2*length( delta )],  0 + delta[i-2*length( delta )], 0 + delta[i-2*length( delta )] ), dist_par2 = c( 1, 1, 1, 1, 1, 1  ) )
    sim_results <- nmANOVA_sim( distribution, group_N, par_vec, H1_group, show_plots = FALSE )
    p_value[i] <- sim_results[[1]]$summary$p_value 
    
    
  }
  
  
}

p_value_points <- p_value





for( i in 1:( 3*length( delta ) ) ){
  print(i)
  if ( i <= length( delta ) ){
    set.seed(411)
    H1_group <- data.frame( groups = c( 'group2,group1', 'group1,group2' ), dist_par1 = c( 0 + delta[i], 0 + delta[i] ), dist_par2 = c( 1, 1 ) )
    sim_results <- nmANOVA_sim( distribution, group_N, par_vec, H1_group, show_plots = FALSE )
    p_value[i] <- sim_results[[1]]$summary$p_value 
  }
  
  if ( ( i > length( delta ) )  & ( i <= 2*length( delta ) ) ){
    set.seed(411)
    H1_group <- data.frame( groups = c( 'group2,group1', 'group1,group2', 'group1,group3', 'group3,group1'  ), dist_par1 = c( 0 + delta[i-length( delta )], 0 + delta[i-length( delta )],  0 + delta[i-length( delta )], 0 + delta[i-length( delta )] ), dist_par2 = c( 1 , 1 , 1 , 1  ) )
    sim_results <- nmANOVA_sim( distribution, group_N, par_vec, H1_group, show_plots = FALSE )
    p_value[i] <- sim_results[[1]]$summary$p_value 
    
    
  }
  
  if ( ( i > 2*length( delta ) )   ){
    set.seed(411)
    H1_group <- data.frame( groups = c( 'group2,group1', 'group1,group2', 'group1,group3', 'group3,group1' ,  'group2,group3', 'group3,group2'  ), dist_par1 = c( 0 + delta[i-2*length( delta )], 0 + delta[i-2*length( delta )], 0 + delta[i-2*length( delta )], 0 + delta[i-2*length( delta )],  0 + delta[i-2*length( delta )], 0 + delta[i-2*length( delta )] ), dist_par2 = c( 1, 1, 1, 1, 1, 1  ) )
    sim_results <- nmANOVA_sim( distribution, group_N, par_vec, H1_group, show_plots = FALSE )
    p_value[i] <- sim_results[[1]]$summary$p_value 
    
    
  }
  
  
}


###################


load("figure1A_p_values_points.RData")
p_value_points <- p_value

load("figure1A_p_values.RData")
delta <- seq( 0, 10, length.out = 1000 )

library(dplyr)
library( ggplot2 )



df <- data.frame( delta = rep( delta, 3 ), num = c( seq(1,1000, length.out = 1000), seq( 1100,2100, length.out = 1000),seq( 2200,3200, length.out = 1000) ) ,  group_num = c( rep( 1, 1000 ), rep( 2, 1000 ), rep( 3, 1000 ) ),  x = p_value )
df2 <- data.frame( delta = rep( delta, 3 ), num = c( seq(1,1000, length.out = 1000), seq( 1100,2100, length.out = 1000),seq( 2200,3200, length.out = 1000) ) ,  group_num = c( rep( 1, 1000 ), rep( 2, 1000 ), rep( 3, 1000 ) ),  x = p_value_points )


#data frame with mean and median values used for the figure
sum_stat <- 
  df %>%
  group_by( group_num ) %>%
  dplyr::summarise( mean = mean(x), median = median( x) )%>%
  mutate( log_mean = log( mean ), log_median = log( median ) )


ggplot( df, aes( x =  num, y = ( log(x) ) ) ) +  
  
  geom_point( data = df2, aes( x =  num, y = ( log(x) ) ),  alpha = 0.05, size = 4, col = 'red')  +
  geom_line( alpha = 0.7,  size = 2, col = 'blue')  +
  
  geom_segment( aes(x = 1, xend=200,  y = log( sum_stat$mean[1] ) ), yend=log( sum_stat$mean[1]), size = 1) +
  geom_segment( aes(x = 1100, xend=1300,  y = log( sum_stat$mean[2] ) ), yend=log( sum_stat$mean[2]), size = 1) +
  geom_segment( aes(x = 2200, xend=2400,  y = log( sum_stat$mean[3]) ), yend=log( sum_stat$mean[3]), size = 1) +
  
  geom_segment( aes(x = 400, xend=600,  y = log( sum_stat$median[1] ) ), yend=log( sum_stat$median[1]), size = 1) +
  geom_segment( aes(x = 1500, xend=1700,  y = log( sum_stat$median[2] ) ), yend=log( sum_stat$median[2]), size = 1) +
  geom_segment( aes(x = 2600, xend=2800,  y = log( sum_stat$median[3] ) ), yend=log( sum_stat$median[3]), size = 1) +
  
  scale_x_continuous(breaks=c(1, 1000, 1100, 2100, 2200,3200) ,
                     labels=c("1",'10',  "1", "10", '1', '10')) +
  
  theme(panel.grid.minor = element_line(colour = "grey", size = 0.1), 
        panel.grid.major = element_line(colour = "grey", size = 0.5)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme_bw()+ xlab( expression(delta) ) + ylab( "log(p-value)" ) +  
  theme(plot.title = element_text(size = 25),axis.text.x = element_text(angle = 0, hjust = 1, size = 25), axis.title=element_text(size=25), axis.text.y = element_text(angle = 0, hjust = 1, size = 25), legend.title=element_text(size=25), legend.text=element_text(size=25) )  + theme(legend.position = "none")

###################


#Figure 1 B

#YOU DO NOT NEED TO RUN THIS CODE TO REPRODUCE THE FIGURE 1 B, AS WE SAVED THE RESULTS IN figure1B_p_values.RData, figure1B_p_values_points.RData ( run from load("figure1B_p_values_points.RData") )
#########################

distribution <- 'uniform' #distribution for the simulation matrix
group_N <- c( 10, 10, 10 )#group sizes 
par_vec <- c( 0, 1 )# distribution parameters

#in case you want to introduce the difference - select some groups and change parameters in H1_group
#dist_par1 is a new parameter to use instead of par1 for the specified groups
#dist_par2 is a new parameter to use instead of par2 for the specified groups


#perform simulations and check distribution of p-values
Nsim <- 1
p_value <- vector()
delta <- seq( 0, 10, length.out = 1000 )



for( i in 1:( 3*length( delta ) ) ){
  print(i)
  if ( i <= length( delta ) ){
    
    H1_group <- data.frame( groups = c( 'group2,group1', 'group1,group2' ), dist_par1 = c( 0 + delta[i], 0 + delta[i] ), dist_par2 = c( 1 + delta[i], 1 + delta[i] ) )
    sim_results <- nmANOVA_sim( distribution, group_N, par_vec, H1_group, show_plots = FALSE )
    p_value[i] <- sim_results[[1]]$summary$p_value 
  }
  
  if ( ( i > length( delta ) )  & ( i <= 2*length( delta ) ) ){
    
    H1_group <- data.frame( groups = c( 'group2,group1', 'group1,group2', 'group1,group3', 'group3,group1'  ), dist_par1 = c( 0 + delta[i-length( delta )], 0 + delta[i-length( delta )],  0 + delta[i-length( delta )], 0 + delta[i-length( delta )] ), dist_par2 = c( 1 + delta[i-length( delta )], 1 + delta[i-length( delta )], 1 + delta[i-length( delta ) ], 1 + delta[i-length( delta )]  ) )
    sim_results <- nmANOVA_sim( distribution, group_N, par_vec, H1_group, show_plots = FALSE )
    p_value[i] <- sim_results[[1]]$summary$p_value 
    
    
  }
  
  if ( ( i > 2*length( delta ) )   ){
    
    H1_group <- data.frame( groups = c( 'group2,group1', 'group1,group2', 'group1,group3', 'group3,group1' ,  'group2,group3', 'group3,group2'  ), dist_par1 = c( 0 + delta[i-2*length( delta )], 0 + delta[i-2*length( delta )], 0 + delta[i-2*length( delta )], 0 + delta[i-2*length( delta )],  0 + delta[i-2*length( delta )], 0 + delta[i-2*length( delta )] ), dist_par2 = c( 1 + delta[i-2*length( delta )], 1 + delta[i-2*length( delta )], 1 + delta[i-2*length( delta ) ], 1 + delta[i-2*length( delta )], 1 + delta[i-2*length( delta )], 1 + delta[i-2*length( delta )]  ) )
    sim_results <- nmANOVA_sim( distribution, group_N, par_vec, H1_group, show_plots = FALSE )
    p_value[i] <- sim_results[[1]]$summary$p_value 
    
    
  }
  
  
}

p_value_points <- p_value



for( i in 1:( 3*length( delta ) ) ){
  print(i)
  if ( i <= length( delta ) ){
    set.seed(412)
    H1_group <- data.frame( groups = c( 'group2,group1', 'group1,group2' ), dist_par1 = c( 0 + delta[i], 0 + delta[i] ), dist_par2 = c( 1 + delta[i], 1 + delta[i] ) )
    sim_results <- nmANOVA_sim( distribution, group_N, par_vec, H1_group, show_plots = FALSE )
    p_value[i] <- sim_results[[1]]$summary$p_value 
  }
  
  if ( ( i > length( delta ) )  & ( i <= 2*length( delta ) ) ){
    set.seed(412)
    H1_group <- data.frame( groups = c( 'group2,group1', 'group1,group2', 'group1,group3', 'group3,group1'  ), dist_par1 = c( 0 + delta[i-length( delta )], 0 + delta[i-length( delta )],  0 + delta[i-length( delta )], 0 + delta[i-length( delta )] ), dist_par2 = c( 1 + delta[i-length( delta )], 1 + delta[i-length( delta )], 1 + delta[i-length( delta ) ], 1 + delta[i-length( delta )]  ) )
    sim_results <- nmANOVA_sim( distribution, group_N, par_vec, H1_group, show_plots = FALSE )
    p_value[i] <- sim_results[[1]]$summary$p_value 
    
    
  }
  
  if ( ( i > 2*length( delta ) )   ){
    set.seed(412)
    H1_group <- data.frame( groups = c( 'group2,group1', 'group1,group2', 'group1,group3', 'group3,group1' ,  'group2,group3', 'group3,group2'  ), dist_par1 = c( 0 + delta[i-2*length( delta )], 0 + delta[i-2*length( delta )], 0 + delta[i-2*length( delta )], 0 + delta[i-2*length( delta )],  0 + delta[i-2*length( delta )], 0 + delta[i-2*length( delta )] ), dist_par2 = c( 1 + delta[i-2*length( delta )], 1 + delta[i-2*length( delta )], 1 + delta[i-2*length( delta ) ], 1 + delta[i-2*length( delta )], 1 + delta[i-2*length( delta )], 1 + delta[i-2*length( delta )]  ) )
    sim_results <- nmANOVA_sim( distribution, group_N, par_vec, H1_group, show_plots = FALSE )
    p_value[i] <- sim_results[[1]]$summary$p_value 
    
    
  }
  
  
}

#########################


load("figure1B_p_values_points.RData")
p_value_points <- p_value

load("figure1B_p_values.RData")
delta <- seq( 0, 10, length.out = 1000 )

library(dplyr)
library( ggplot2 )


df <- data.frame( delta = rep( delta, 3 ), num = c( seq(1,1000, length.out = 1000), seq( 1100,2100, length.out = 1000),seq( 2200,3200, length.out = 1000) ) ,  group_num = c( rep( 1, 1000 ), rep( 2, 1000 ), rep( 3, 1000 ) ),  x = p_value )

df2 <- data.frame( delta = rep( delta, 3 ), num = c( seq(1,1000, length.out = 1000), seq( 1100,2100, length.out = 1000),seq( 2200,3200, length.out = 1000) ) ,  group_num = c( rep( 1, 1000 ), rep( 2, 1000 ), rep( 3, 1000 ) ),  x = p_value_points )

#data frame with mean and median values used for the figure
sum_stat <- 
  df %>%
  group_by( group_num ) %>%
  dplyr::summarise( mean = mean(x), median = median( x) )%>%
  mutate( log_mean = log( mean ), log_median = log( median ) )


ggplot( df, aes( x =  num, y = ( log(x) ) ) ) +  
  
  geom_point( data = df2, aes( x =  num, y = log( (x) ) ),  alpha = 0.05, size = 4, col = 'red')  +
  geom_line( alpha = 0.7,  size = 2, col = 'blue')  +
  
  geom_segment( aes(x = 200, xend=400,  y = log( sum_stat$mean[1])), yend=log( sum_stat$mean[1]), size = 1) +
  geom_segment( aes(x = 1300, xend=1500,  y = log( sum_stat$mean[2])), yend=log( sum_stat$mean[2]), size = 1) +
  geom_segment( aes(x = 2400, xend=2600,  y = log( sum_stat$mean[3])), yend=log( sum_stat$mean[3]), size = 1) +
  
  geom_segment( aes(x = 400, xend=600,  y = log( sum_stat$median[1])), yend=log( sum_stat$median[1]), size = 1) +
  geom_segment( aes(x = 1500, xend=1700,  y = log( sum_stat$median[2])), yend=log( sum_stat$median[2]), size = 1) +
  geom_segment( aes(x = 2600, xend=2800,  y = log( sum_stat$median[3])), yend=log( sum_stat$median[3]), size = 1) +
  
  scale_x_continuous(breaks=c(1, 1000, 1100, 2100, 2200,3200) ,
                     labels=c("1",'10',  "1", "10", '1', '10')) +
  
  theme(panel.grid.minor = element_line(colour = "grey", size = 0.1), 
        panel.grid.major = element_line(colour = "grey", size = 0.5)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme_bw()+ xlab( expression(delta) ) + ylab( "log( p-value )" ) +  theme(plot.title = element_text(size = 25),axis.text.x = element_text(angle = 0, hjust = 1, size = 25), axis.title=element_text(size=25), axis.text.y = element_text(angle = 0, hjust = 1, size = 25), legend.title=element_text(size=25), legend.text=element_text(size=25) )  + theme(legend.position = "none")

