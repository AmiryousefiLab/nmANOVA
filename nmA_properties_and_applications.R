#Codes to replicate figure 2 A-D from the manuscript
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

#Figure 2 A
###########################

library( stringr )
library( dplyr )
library( plyr )
library( ggplot2 )

#the code to obtain figure2A_data.RData has been executed using CSC clusters - full code for CSC environment upon request

load( "figure2A_data.RData" )
df_list  <- list()
gg_list  <- list()

for( i in 3:5 ){
  
  df2 <- 
    df %>%
    filter( clust_num == i )
  
  
  if( i == 3 )
    partition <- c( rep( 'group_1' ,34 ), rep( 'group_2' ,34 ),  rep( 'group_3' ,55 ) )
  
  if( i == 4 )
    partition <- c( rep( 'group_1' ,13 ), rep( 'group_2' ,21),  rep( 'group_3' ,34 ), rep( 'group_4' ,55 ) )
 
  
  if( i == 5 )
    partition <- c( rep( 'group_1', 13 ), rep( 'group_2' ,21 ), rep( 'group_3' ,34),  rep( 'group_4' ,21 ), rep( 'group_5' ,34 ) )
  
  
  indi <- which( ( df2$N_sampl == 0 ) &  ( df2$num_change == 0 ) )
  
  df2_null <-
    df2 %>%
    filter( N_sampl == 0, num_change == 0   )
  p_value <-  median( df2_null$p_value )
  

  d <- data.frame( num_change = 0, mean = ( p_value), min = (p_value ), max = (p_value), med = (p_value))
  
  ind <- which(df2$part == paste0( partition, collapse = ',' ) )
  df2 <- df2[-ind, ]
  
  
  df2 <- 
    df2 %>%
    filter( N_sampl != 0, num_change != 0   )%>%
    group_by( num_change, part)%>%
    dplyr::summarise(   p_value = mean( p_value) )
  
  
  
  
  d <- data.frame( num_change = 0, mean = ( p_value), min = (p_value ), max = (p_value), med = (p_value))
  
  gg <- 
    df2 %>%
    group_by(num_change )%>%
    dplyr::summarise( mean = mean( p_value), min = min(p_value ), max = max(p_value), med = median(p_value) )
  
  gg <- rbind( d, gg)
  
  df_list[[i-1]]  <-
    df2 %>%
    mutate( clust_num =  i )
  
  gg_list[[i-1]] <- 
    gg %>%
    mutate( clust_num =  i )  
  
  
}



df_all <- ldply( df_list, as.data.frame )
gg_all <- ldply( gg_list, as.data.frame )

df_all <- 
  df_all %>%
  mutate( hist_factor = paste( num_change, clust_num, sep = '_' ) )

ggplot( df_all, aes(   x =  ( num_change ), y = (p_value ), col = as.factor( clust_num ) ) ) +  

  geom_ribbon(data=gg_all,aes(ymin=( min ),ymax= max, fill =as.factor( clust_num )) ,alpha=0.1, colour = NA ) +
  geom_line(data=gg_all ,aes( ( num_change ), ( med ), col = as.factor( clust_num ), group = as.factor( clust_num ) ), size = 1, alpha = 0.7) +
  geom_point(data=gg_all,aes( ( num_change ), (med), col = as.factor( clust_num ), group = as.factor( clust_num ) ), size =1.4) +
  
  scale_x_continuous(breaks=c(0,  10,  30, 60 ) ,
                     labels=c('0',  '10',  '30',  '60' ) ) + ylim( 0,1 ) +
  theme(panel.grid.minor = element_line(colour = "grey", size = 0.1), 
        panel.grid.major = element_line(colour = "grey", size = 0.5)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + guides( col=guide_legend( title="Number of clusters" ), fil= guide_legend( title="Number of clusters" ) )
  + ylab( 'p-value' ) + xlab( "Number of pairs switched" ) +  theme(plot.title = element_text(size = 25),axis.text.x = element_text(angle = 0, hjust = 1, size = 20), axis.title=element_text(size=20), axis.text.y = element_text(angle = 0, hjust = 1, size = 20), legend.title=element_text(size=20), legend.text=element_text(size=20) ) # + theme(legend.position = "none") #+theme(axis.title=element_blank() )# +  annotate("text", x = 3, y = 11.5, label= expression( paste( italic( "p value for original matrix" )  ) ), size = 7



###########################



#Figure 2 B
###########################
library( ggplot2 )

#the code to obtain 2 million p-values has been executed using CSC clusters - full code for CSC environment upon request
#load("figure2B_diss_mat.RData") - dissimilarity matrix has been used to obtain p-values

load( "figure2B_p_values.RData" )
m <- mean( ( df$p_value ) )
df2 <- df[c(50:2000000),]# only for visualization purposes - to highlight convergence

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

#takes long time to execute
ggplot(df2, aes(x=log10( prop_sampl_num ), y=( -log10( p_value ) )))+ 
  geom_line()+
  scale_x_continuous(breaks=c( log10( 50 ), log10( 10000),  log10( 100000), log10( 500000 ),    log10( 2000000 ) ),
                     labels= scientific_10( c( 50, 10000, 100000,  500000, 2000000 ) ) ) +
  geom_hline( yintercept = -log10( mean( ( df$p_value ) ) ), col = 'red' ) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1) ) +
  xlab( "Number of proportional samplings performed" ) + ylab( "-log10(p-value)" ) + 
  theme(plot.title = element_text(size = 25),axis.text.x = element_text(angle = 0, hjust = 1, size = 20), axis.title=element_text(size=25), axis.text.y = element_text(angle = 0, hjust = 1, size = 25), legend.title=element_text(size=20), legend.text=element_text(size=20) ) 

###########################



#Figure 2 C
###########################
#the code to obtain figure2C_data.RData has been executed using CSC clusters - full code for CSC environment upon request
load( "figure2C_data.Rdata" )

library( stringr )
library( dplyr )
library( plyr )
library( fmsb )
library( RColorBrewer )

df_thr <- 
  df %>%
  mutate( nm_stat_0.05 = ifelse( nm_pval <= 0.05, 'yes', 'no' ), nm_stat_0.1 = ifelse( nm_pval <= 0.1, 'yes', 'no' ),
          np_stat_0.05 = ifelse( np_pval <= 0.05, 'yes', 'no' ), np_stat_0.1 = ifelse( np_pval <= 0.1, 'yes', 'no' ),
          cl_stat_0.05 = ifelse( cl_pval <= 0.05, 'yes', 'no' ), cl_stat_0.1 = ifelse( cl_pval <= 0.1, 'yes', 'no' ) )%>%
  mutate( group = paste( dist_method, vector_num, sep = '_' ) )

all_groups <- unique( df_thr$group )

df_5 <- as.data.frame( matrix(0, ncol = length( all_groups ), nrow = 5 ) )
colnames( df_5 ) <- all_groups
rownames( df_5 ) <- c( 'max', 'min', 'cl', 'np', 'nm' )
df_5[1, ] <- 0
df_5[2, ] <- -2

df_1 <- df_5

for( j in 1:length( all_groups ) ){
  
  cur_df <- 
    df_thr %>%
    filter( group == all_groups[j] ) %>%
    arrange( delta )
  
  
  
  df_1[3,j] <- -cur_df$delta[ which( cur_df$cl_stat_0.1 == 'yes')[1] ]
  print( 'cl' ) 
  print( length( which( cur_df$cl_stat_0.1[ ( which( cur_df$cl_stat_0.1 == 'yes' ) +1 ):nrow( cur_df ) ] == 'no' ) ) ) 
  
  
  df_1[4,j] <- -cur_df$delta[ which( cur_df$np_stat_0.1 == 'yes')[1] ]
  print( 'np' ) 
  print( length( which( cur_df$np_stat_0.1[ ( which( cur_df$np_stat_0.1 == 'yes' ) +1 ):nrow( cur_df ) ] == 'no' ) ) ) 
  
  
  
  if(  length( which( cur_df$nm_stat_0.1 == 'yes') ) >0 ){
    df_1[5,j] <- -cur_df$delta[ which( cur_df$nm_stat_0.1 == 'yes')[1] ]
    print( 'nm' ) 
    print( length( which( cur_df$nm_stat_0.1[ ( which( cur_df$nm_stat_0.1 == 'yes' ) +1 ):nrow( cur_df ) ] == 'no' ) ) ) 
    
  }
  
  
}

# Set graphic colors

colors_border=c(  rgb(0,1,0.6,0.9), rgb(1,0.2,0.4,0.6) ,  rgb(0,0.6,1,0.7) )
colors_in=c( rgb(0,1,0.6,0.4), rgb(1,0.2,0.4,0.4) , rgb(0,0.6,1,0.2) )


# plot with default options:
radarchart( df_1  , axistype=1 , 
            #custom polygon
            pcol=colors_border , pfcol=colors_in , plwd=2 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="black",  cglwd=0.8, caxislabels=c( 2, 1.5, 1, 0.5, 0),
            #custom labels
            vlcex=1 
)

# Add a legend
legend(x=1.5, y=1, legend = rownames(df_1[-c(1,2),]), bty = "n", pch=20 , col=colors_in , text.col = "grey", cex=1.2, pt.cex=3)



###########################



#Figure 2 D
###########################
#the code to obtain figure2D_data.RData has been executed using CSC clusters - full code for CSC environment upon request
load( "figure2D_data.Rdata" )

library( stringr )
library( dplyr )
library( plyr )
library( fmsb )
library( RColorBrewer )


df_thr <- 
  df %>%
  mutate( nm_stat_0.05 = ifelse( nm_pval <= 0.05, 'yes', 'no' ), nm_stat_0.1 = ifelse( nm_pval <= 0.1, 'yes', 'no' ),
          np_stat_0.05 = ifelse( np_pval <= 0.05, 'yes', 'no' ), np_stat_0.1 = ifelse( np_pval <= 0.1, 'yes', 'no' )
  )%>%
  mutate( group = paste( dist_method, vector_num, sep = '_' ) )

all_groups <- unique( df_thr$group )

df_5 <- as.data.frame( matrix(0, ncol = length( all_groups ), nrow = 4) )
colnames( df_5 ) <- all_groups
rownames( df_5 ) <- c( 'max', 'min',  'np', 'nm' )
df_5[1, ] <- 0
df_5[2, ] <- -2

df_1 <- df_5

for( j in 1:length( all_groups ) ){
  
  cur_df <- 
    df_thr %>%
    filter( group == all_groups[j] ) %>%
    arrange( delta )
  
  
  
  df_1[3,j] <- -cur_df$delta[ which( cur_df$np_stat_0.1 == 'yes')[1] ]
  print( 'np' ) 
  print( length( which( cur_df$np_stat_0.1[ ( which( cur_df$np_stat_0.1 == 'yes' ) +1 ):nrow( cur_df ) ] == 'no' ) ) ) 
  
  df_1[4,j] <- -cur_df$delta[ which( cur_df$nm_stat_0.1 == 'yes')[1] ]
  print( 'nm' ) 
  print( length( which( cur_df$nm_stat_0.1[ ( which( cur_df$nm_stat_0.1 == 'yes' ) +1 ):nrow( cur_df ) ] == 'no' ) ) ) 
  
  
  
}

df_1 <- df_1[c(1,2, 4, 3),]


# Color vector
colors_border=c( rgb(0,0.6,1,0.7), rgb(1,0.2,0.4,0.6) , rgb(0.7,0.5,0.1,0.9) )
colors_in=c( rgb(0,0.6,1,0.2), rgb(1,0.2,0.4,0.4) , rgb(0.7,0.5,0.1,0.4) )

# plot with default options:
radarchart( df_1  , axistype=1 , 
            #custom polygon
            pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=c(20, 15, 10, 5, 1 ), cglwd=0.8,
            #custom labels
            vlcex=0.8 
)

# Add a legend
legend(x=1.5, y=1, legend = rownames(df_1[-c(1,2),]), bty = "n", pch=20 , col=colors_in , text.col = "grey", cex=1.2, pt.cex=3)

###########################


