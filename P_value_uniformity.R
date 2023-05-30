#Codes to replicate figure 1 from the manuscript
#Note that because of random sampling in the nmANOVA algorithm and system-specific numerical instabilities, executing the 
#provided replication materials within a different environment will likely
#produce minor differences from what has been shown in the paper

#nmANOVA function
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
          
          
          if( floor( non_na_cur_partition  *  ( non_na_vec_within  / as.numeric( num_off_diag ) )  )  != 0 ){ diagonal_sampling <- diagonal_sampling[ -sample_ind ] }
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

library( ggplot2 )
library( pheatmap )

#functions to generate partitions based on selected distributions
generate_diss <- function( distribution, group_N, df_dist_par, show_plots ){
  
  if ( distribution == 'Fisher' ) {
    
    dist_function <- function( N, f1, f2 ){
      r <- rf( N, f1, f2 )
      r
    } 
  }
  if ( distribution == 'Exp' ) {
    
    dist_function <- function( N, e1 ){
      r <- rexp( N, e1 )
      r
    } 
  }
  if ( distribution == 'Poisson' ) {
    
    dist_function <- function( N, p1 ){
      r <- rpois( N, p1 )
      r
    } 
  }
  if ( distribution == 'Binomial' ) {
    
    dist_function <- function( N, bi1, bi2 ){
      r <- rbinom( N, bi1, bi2 )
      r
    } 
  }
  if ( distribution == 'PD' ) {
    
    dist_function <- function( N, pd1 ){
      r <- rPD( N, pd )
      r
    } 
  }
  if ( distribution == 'Pareto' ) {
    
    dist_function <- function( N, pa1, pa2 ){
      r <- rpareto( N, pa1, pa2 )
      r
    } 
  }
  if ( distribution == 'Rayleigh' ) {
    
    dist_function <- function( N, re1 ){
      r <- rrayleigh( N, re1 )
      r
    } 
  }
  if ( distribution == 'Vonmises' ) {
    
    dist_function <- function( N, vm1, vm2 ){
      r <- rvonmises( N, vm1, vm2 )
      r
    } 
  }
  if ( distribution == 'Bort' ) {
    
    dist_function <- function( N, bo1, bo2 ){
      r <- rbort( N, bo1, bo2 )
      r
    } 
  }
  
  if ( distribution == 'Hermite' ) {
    
    dist_function <- function( N, h1, h2 ){
      r <- rhermite( N, h1, h2 )
      r
    } 
  }
  if ( distribution == 'Rice' ) {
    
    dist_function <- function( N, r1, r2 ){
      r <- rrice( N, r1, r2 )
      r
    } 
  }
  if ( distribution == 'Zipf' ) {
    
    dist_function <- function( N, z1, z2 ){
      r <- rzipf( N, z1, z2 )
      r
    } 
  }
  if ( distribution == 'Uniform' ) {
    
    dist_function <- function( N, u1, u2 ){
      r <- runif( N, u1, u2 )
      r
    } 
  }
  
  if ( distribution == 'Normal' ) {
    
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
  
  if ( distribution == 'Cauchy' ) {
    
    dist_function <- function( N, c1, c2 ){
      r <- rcauchy( N, c1, c2 )
      r
    } 
  }
  
  if ( distribution == 'Beta' ) {
    
    dist_function <- function( N, b1, b2 ){
      r <- rbeta( N, b1, b2 )
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

mr<- c("Normal", "Beta", "Cauchy", "Fisher", "Uniform", "Exp", "Poisson", "Binomial", "PD","Pareto", "Reiyghlee", "vonMies", "Borel", "Hermite", "Rice", "Zipf")

#####parameter assignment 
set.seed(123456)
s<-sample(1:400, 28)
#normal
n1<-s[1] ; n2<-s[2] 
#poisson-dirichlete
pd<-s[3]
#beta
b1<-s[4]; b2<-s[5]
#cuachy
c1<-s[6];c2<-s[7]
#borel
bo1<-s[8] ; bo2<-s[9]/400
#zipf
z1<-s[10];z2<-s[11]
#hermite
h1<-s[12]; h2<-s[13]
#rice
r1<-s[14]; r2<-s[15]
#unif
u1<-min(s[16], s[17]); u2<-max(s[17], s[16])+1; 
#vonMies
vm1<-s[18];vm2<-s[19]
#reighlee
re1<-s[20]
#poisson
p1<-s[21]
#exp
e1<-s[22]
#pareto
pa1<-s[23]; pa2<-s[24]
#bionomial
bi1<-s[25]; bi2<-s[26]/400
#Fisher
f1<-s[27]; f2<-s[28]


#YOU DO NOT NEED TO RUN THIS CODE TO REPRODUCE THE FIGURE 1 A, AS WE SAVED THE RESULTS IN figure1A_p_values.RData ( run from load("figure1A_p_values.RData") )
###################

distribution <- 'normal' #distribution for the simulation matrix
group_N <- c( 10, 10, 10 )#group sizes = 10
par_vec <- c( 0, 1 )# distribution parameters

#in case you want to introduce the difference - select some groups and change parameters in H1_group
#dist_par1 is a new parameter to use instead of par1 for the specified groups
#dist_par2 is a new parameter to use instead of par2 for the specified groups
H1_group <- data.frame( groups = c( 'group2,group1', 'group1,group2', 'group3,group1', 'group1,group3', 'group2,group3', 'group3,group2' ), dist_par1 = rep( 0, 6 ), dist_par2 = rep( 1, 6 ) )
# here we say that for all 6 between group matrices we have distribution parameters 0 and 1 (in this case mean and sd)
# you can change their distribution eg dist_par1 = rep( 1, 6 ), dist_par2 = rep( 2, 6 ) means all between group matrices will be from U(1,2)
# dist_par1 = c(1, 2, 3, 4, 5, 5), dist_par2 = rep( 1, 6 ) here sd is 1 but means are different for the groups
N_rep <- 1
k <- 0 
p_value <- vector()
norm_p<- numeric(16)
hist(runif(10000, 0,1), probability = T, col = "honeydew", breaks = 20, border = "burlywood")
for (j in 1:16){
  distribution <- mr[j]
for( i in 1:1000 ){
  
  sim_results <- nmANOVA_sim( distribution, group_N, par_vec, H1_group, show_plots = FALSE, N_sampling = 1  )
  diss_mat <- sim_results[[2]]
  partition <- c( rep( 'group_1', 10 ), rep( 'group_2', 10 ), rep( 'group_3', 10 ) )
  prop_sampl_df <-  nmANOVA( diss_mat, partition, N_sam  )
  p_value[i] <- prop_sampl_df$summary$p_value 
  
}

norm_p[j]<- ks.test( p_value, "punif", 0, 1 )$p.value
#hist(p_value,breaks = 30, probability = T)
lines(density(na.omit(p_value)), col = rainbow(16)[j], lwd=1.1)

}

mean(norm_p)
sd(norm_p)

load("figure1A_p_values.RData")

df <- data.frame( x = p_value_h0 )
ks.test( p_value_h0, "punif", 0, 1 )

ggplot(df, aes( x=x ) )+
  geom_histogram(color="black", fill="lightblue", alpha = 0.6, bins = 50) + theme_bw()  + xlab( 'p-value' ) + ylab( "Count" ) + theme(plot.title = element_text(size = 25),axis.text.x = element_text(angle = 0, hjust = 1, size = 25), axis.title=element_text(size=25), axis.text.y = element_text(angle = 0, hjust = 1, size = 25), legend.title=element_text(size=25), legend.text=element_text(size=25) )  + theme(legend.position = "none") +  annotate("text", x=0.52, y=31, label= expression( paste( italic( "Kolmogorov-Smirnov test p-value = 0.5074" )  ) ), size = 7) 


