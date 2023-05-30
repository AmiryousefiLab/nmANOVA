##nmA asymptotics

#libraries with specific functions

library(PEkit)
library(bbl)
library(Rfast)
library(VGAM)
library(EnvStats)
library(hermite)
library(sads)
library(vioplot)
library(moments)

hist(rvonmises(1000, 2, 25, rads = TRUE))

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
            #print( 'prop_num = 0 -> changed to 1' )
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

#sampling from each dist
n<- 100000000
mr<- matrix(0, 100000000, 16)
mr<- as.data.frame(mr)
colnames(mr)<- c("Normal", "Beta", "Cauchy", "Fisher", "Uniform", "Exp", "Poisson", "Binomial", "PD","Pareto", "Reiyghlee", "vonMies", "Borel", "Hermite", "Rice", "Zipf")
mr$Normal<- rnorm(n, n1,n2) 
mr$Beta<- rbeta(n, b1, b2)
mr$Cauchy<- rcauchy(n, c1, c2)
mr$Fisher <- rf(n, f1, f2)
mr$Unif<- runif(n, u1, u2)
mr$Exp<- rexp(n, e1)
mr$Poisson<- rpois(n, p1)
mr$Bionomial<- rbinom(n, bi1, bi2)
mr$PD <- sample(rPD(n/10000, pd), n, T)
mr$Pareto <- rpareto(n, pa1, pa2)
mr$Reiyghlee<- rrayleigh(n, re1)
mr$vonMies<- rvonmises(n, vm1, vm2)
mr$Borel<- sample(rbort(n/10000, bo1, bo2), n, T)
mr$Hermite<- rhermite(n, h1, h2)
mr$Rice <- rrice(n, r1,r2)
mr$Zipf <- rzipf(n, z1, 1)

#forming the blocks of each distribution to be used for nmA input

#B names
Bname<- as.character()
for (i in 1:16){Bname[i]<-paste("B", i, sep="")}

# Ball<- list()
# for (aName in Bname){
#   Ball[[aName]] <- list(aName)
# }
# B<- list()
# for (aName in Bname){
#   B[[aName]] <- get(Bname[aName])
#   for (i in 1:16){
#     Ball[[aName]][i]<- matrix(mr[1:2^(which(Bname==aName)+10),i], sqrt(2^(which(Bname==aName)+10)), sqrt(2^(which(Bname==aName)+10)))
#   }
#   names(Ball[[aName]])<-colnames(mr)
# }



B1<- list()
for (i in 1:16){
  B1[[i]]<- matrix(mr[1:2^11,i], sqrt(2^11), sqrt(2^11))
}
names(B1)<- colnames(mr)

B2<- list()
for (i in 1:16){
  B2[[i]]<- matrix(mr[1:2^12,i], sqrt(2^12), sqrt(2^12))
}
names(B2)<- colnames(mr)

B3<- list()
for (i in 1:16){
  B3[[i]]<- matrix(mr[1:2^13,i], sqrt(2^13), sqrt(2^13))
}
names(B3)<- colnames(mr)

B4<- list()
for (i in 1:16){
  B4[[i]]<- matrix(mr[1:2^14,i], sqrt(2^14), sqrt(2^14))
}
names(B4)<- colnames(mr)

B5<- list()
for (i in 1:16){
  B5[[i]]<- matrix(mr[1:2^15,i], sqrt(2^15), sqrt(2^15))
}
names(B5)<- colnames(mr)

B6<- list()
for (i in 1:16){
  B6[[i]]<- matrix(mr[1:2^16,i], sqrt(2^16), sqrt(2^16))
}
names(B6)<- colnames(mr)

B7<- list()
for (i in 1:16){
  B7[[i]]<- matrix(mr[1:2^17,i], sqrt(2^17), sqrt(2^17))
}
names(B7)<- colnames(mr)

B8<- list()
for (i in 1:16){
  B8[[i]]<- matrix(mr[1:2^18,i], sqrt(2^18), sqrt(2^18))
}
names(B8)<- colnames(mr)


B9<- list()
for (i in 1:16){
  B9[[i]]<- matrix(mr[1:2^19,i], sqrt(2^19), sqrt(2^19))
}
names(B9)<- colnames(mr)

B10<- list()
for (i in 1:16){
  B10[[i]]<- matrix(mr[1:2^20,i], sqrt(2^20), sqrt(2^20))
}
names(B10)<- colnames(mr)

B11<- list()
for (i in 1:16){
  B11[[i]]<- matrix(mr[1:2^21,i], sqrt(2^21), sqrt(2^21))
}
names(B11)<- colnames(mr)

B12<- list()
for (i in 1:16){
  B12[[i]]<- matrix(mr[1:2^22,i], sqrt(2^22), sqrt(2^22))
}
names(B12)<- colnames(mr)


B13<- list()
for (i in 1:16){
  B13[[i]]<- matrix(mr[1:2^23,i], sqrt(2^23), sqrt(2^23))
}
names(B13)<- colnames(mr)

B14<- list()
for (i in 1:16){
  B14[[i]]<- matrix(mr[1:2^24,i], sqrt(2^24), sqrt(2^24))
}
names(B14)<- colnames(mr)

B15<- list()
for (i in 1:16){
  B15[[i]]<- matrix(mr[1:2^25,i], sqrt(2^25), sqrt(2^25))
}
names(B15)<- colnames(mr)

B16<- list()
for (i in 1:16){
  B16[[i]]<- matrix(mr[1:2^26,i], sqrt(2^26), sqrt(2^26))
}
names(B16)<- colnames(mr)


#for ming the ranodm sectioning for the inner matrices from 2:1024
#these would be the lower tiangel indices for each column

I1<- matrix(0, 16, 16)
for (i in 1:16){
  I1[,i]<- c(sample(2:(sqrt(2^11)-1), i, F), sample(0, 16-i, T))
}
I1<-apply(I1, 2, sort)

I2<- matrix(0, 16, 16)
for (i in 1:16){
  I2[,i]<- c(sample(2:(sqrt(2^12)-1), i, F), sample(0, 16-i, T))
}
I2<-apply(I2, 2, sort)

I3<- matrix(0, 16, 16)
for (i in 1:16){
  I3[,i]<- c(sample(2:(sqrt(2^13)-1), i, F), sample(0, 16-i, T))
}
I3<-apply(I3, 2, sort)

I4<- matrix(0, 16, 16)
for (i in 1:16){
  I4[,i]<- c(sample(2:(sqrt(2^14)-1), i, F), sample(0, 16-i, T))
}
I4<-apply(I4, 2, sort)

I5<- matrix(0, 16, 16)
for (i in 1:16){
  I5[,i]<- c(sample(2:(sqrt(2^15)-1), i, F), sample(0, 16-i, T))
}
I5<-apply(I5, 2, sort)

I6<- matrix(0, 16, 16)
for (i in 1:16){
  I6[,i]<- c(sample(2:(sqrt(2^16)-1), i, F), sample(0, 16-i, T))
}
I6<-apply(I6, 2, sort)

I7<- matrix(0, 16, 16)
for (i in 1:16){
  I7[,i]<- c(sample(2:(sqrt(2^17)-1), i, F), sample(0, 16-i, T))
}
I7<-apply(I7, 2, sort)

I8<- matrix(0, 16, 16)
for (i in 1:16){
  I8[,i]<- c(sample(2:(sqrt(2^18)-1), i, F), sample(0, 16-i, T))
}
I8<-apply(I8, 2, sort)

I9<- matrix(0, 16, 16)
for (i in 1:16){
  I9[,i]<- c(sample(2:(sqrt(2^19)-1), i, F), sample(0, 16-i, T))
}
I9<-apply(I9, 2, sort)

I10<- matrix(0, 16, 16)
for (i in 1:16){
  I10[,i]<- c(sample(2:(sqrt(2^20)-1), i, F), sample(0, 16-i, T))
}
I10<-apply(I10, 2, sort)

I11<- matrix(0, 16, 16)
for (i in 1:16){
  I11[,i]<- c(sample(2:(sqrt(2^21)-1), i, F), sample(0, 16-i, T))
}
I11<-apply(I11, 2, sort)

I12<- matrix(0, 16, 16)
for (i in 1:16){
  I12[,i]<- c(sample(2:(sqrt(2^22)-1), i, F), sample(0, 16-i, T))
}
I12<-apply(I12, 2, sort)

I13<- matrix(0, 16, 16)
for (i in 1:16){
  I13[,i]<- c(sample(2:(sqrt(2^23)-1), i, F), sample(0, 16-i, T))
}
I13<-apply(I13, 2, sort)

I14<- matrix(0, 16, 16)
for (i in 1:16){
  I14[,i]<- c(sample(2:(sqrt(2^24)-1), i, F), sample(0, 16-i, T))
}
I14<-apply(I14, 2, sort)

I15<- matrix(0, 16, 16)
for (i in 1:16){
  I15[,i]<- c(sample(2:(sqrt(2^25)-1), i, F), sample(0, 16-i, T))
}
I15<-apply(I15, 2, sort)

I16<- matrix(0, 16, 16)
for (i in 1:16){
  I16[,i]<- c(sample(2:(sqrt(2^26)-1), i, F), sample(0, 16-i, T))
}
I16<-apply(I16, 2, sort)



p_maker<- function(inner, outer, B, I ,N){
  #inner//outer==  c("Normal", "Beta", "Cauchy", "Fisher", "Unif", "Exp", "Poisson", "Bionomial",
  #"PD","Pareto", "Reiyghlee", "vonMies", "Borel", "Hermite", "Rice", "Zipf")
  #b--> matrix index=1...16
  #N--> N_partition== 1...16
  base<- B[[which(names(B)==outer)]]
  part<-c(I[(16-N+1):16,N])
  leg<- dim(base)[1]
  partition<- list()
  partition[[1]]<- rep(1,part[1])
  if (length(part) > 1) {
    for(i in 2:length(part)){
      partition[[i]]<- rep(i, part[i]-part[(i-1)])
    }
  }
  partition[[(length(part)+1)]]<- rep((length(part)+1),(leg-part[length(part)]))
  
  partition<- unlist(partition)
  
  dist<-base
  inn<-B[[which(names(B)==inner)]]
  dist[1:part[1], 1:part[1]]<-inn[1:part[1], 1:part[1]]
  if(length(part) > 1 ){
    for(i in 1:(length(part)-1)){
      dist[part[i]:part[(i+1)], part[i]:part[(i+1)]] <- inn[part[i]:part[(i+1)], part[i]:part[(i+1)]]
    }
  }
  dist[part[length(part)]:(leg-part[length(part)]), part[length(part)]:(leg-part[length(part)])]<- inn[part[length(part)]:(leg-part[length(part)]), part[length(part)]:(leg-part[length(part)])]
  return(nmANOVA(dist, partition)$summary$p_value)
}

##forming the input for the final plot

MAP<- matrix(0, 16, 16)

#making the four vecotrs of the names as anti-clockwise strating from the bottom
ax1<- dens
ax2<- dens[c(1:4,9:12, 5:8, 13:16)]
ax3<- dens[c(1,9,2,10,3,11,4,12,5,13,6,14,7,15,8,16)] 
ax4<- dens[c(1,2,9,10,3,4,11,12,5,6,13,14,7,8,15,16)]

B<- list(B1, B2,B3,B4,B5,B6, B7,B8,B9,B10,B11,B12,B13,B14,B15,B16)
I<- list(I1, I2,I3,I4,I5,I6, I7,I8,I9,I10,I11,I12,I13,I14,I15,I16)
#lower right-tirangle 
for(j in 1:16){
  for (i in (17-j):16){
    MAP[i,j]<- p_maker(ax1[j], ax2[17-i], B[[j]], I[[j]], (17-i))
  }
}

#upper left-tiangle
for(j in 1:16){
  for (i in 1:(16-j)){
    MAP[i,j]<- p_maker(ax4[i], ax3[j], B[[j]], I[[j]], (17-i))
  }
}

for (i in 1:16){
  MAP[i, 17-i]<- 0.5
}

MAP<- as.data.frame(MAP)
row.names(MAP)<- ax2
colnames(MAP)<- ax1

#Plotting
library(ComplexHeatmap)

# Define some graphics to display the distribution of columns
.hist = anno_histogram(MAP, gp = gpar(fill = "chocolate1"))
.density = anno_density(MAP, type = "line", gp = gpar(col = "blue"))
ha_mix_top = HeatmapAnnotation(
  Histogram = .hist, Density = .density,
  height = unit(3.8, "cm")
)
# Define some graphics to display the distribution of rows
.violin = anno_density(MAP, type = "violin", 
                       gp = gpar(fill = "lightblue"), which = "row")
.boxplot = anno_boxplot(MAP, which = "row")
ha_mix_right = HeatmapAnnotation(Violin = .violin, Boxplot = .boxplot,
                                 which = "row", width = unit(4, "cm"))

pdf("heatmapPub2.pdf", width = 8)

# Combine annotation with heatmap
Heatmap(MAP, name = "p values", show_row_dend = TRUE, cluster_rows = FALSE,
        cluster_row_slices = FALSE, column_dend_reorder = FALSE, cluster_columns = FALSE,
        cluster_column_slices = FALSE, row_dend_reorder = FALSE, show_column_dend = TRUE, 
        
        row_labels = ax4,
        row_names_side = "left",
        show_row_names = TRUE,
        row_names_max_width = unit(6, "cm"),
        row_names_gp = gpar(fontsize = 10),
        row_names_rot = 45,
        row_names_centered = FALSE,
        column_labels = ax1,
        column_names_side = "bottom",
        show_column_names = TRUE,
        column_names_max_height = unit(6, "cm"),
        column_names_rot = 45,
        column_names_centered = FALSE,
        
        column_names_gp = gpar(fontsize = 10),
        top_annotation = ha_mix_top) + ha_mix_right 
dev.off()


MAP1<-MAP





MAP



#sanity check 
pM<- matrix(1, 16, 16)
dens<- c("Normal", "Beta", "Cauchy", "Fisher", "Uniform", "Exp", "Poisson", "Binomial", "PD","Pareto", "Reiyghlee", "vonMies", "Borel", "Hermite", "Rice", "Zipf")
for (i in dens){
  for (j in dens){
    pM[which(dens==i),which(dens==j)]<- p_maker(i, j, B1, I1, 8)
  }
}
##survey the mean and sds of the samples
MEAN<- numeric(16)
STD<- numeric(16)
for (i in 1:16){
  MEAN[i]<-mean(mr[,i])
  STD[i]<- sd(mr[,i])
}


saveRDS(dens, file = "Main axis of variables names.Rds")
saveRDS(s, file = "Parameter template.Rds")
saveRDS(MAP, file = "P-value table.Rds")
saveRDS(I, file = "Blocks of Indices.Rds")
saveRDS(B, file = "Blocks of random variables.Rds")
saveRDS(mr, file = "Master board of ranndom variables.Rds")



mr<- readRDS("Master board of ranndom variables.Rds")
s<- readRDS("Parameter template.Rds")
dens<- readRDS("Main axis of variables names.Rds")
titname<- dens
titname[1]<-paste("Normal(", n1, ",", n2, ")", sep="")
titname[2]<-paste("Beta(", b1, ",", b2, ")", sep="")
titname[3]<-paste("Cauchy(", c1, ",", c2, ")", sep="")
titname[4]<-paste("Fisher(", f1, ",", f2, ")", sep="")
titname[5]<-paste("Uniform(", u1, ",", u2, ")", sep="")
titname[6]<-paste("Exponential(", e1, ")", sep="")
titname[7]<-paste("Poisson(", p1, ")", sep="")
titname[8]<-paste("Binomial(", bi1, ",", bi2, ")", sep="")
titname[9]<-paste("Poisson-Dirichlete(", pd, ")", sep="")
titname[10]<-paste("Pareto(", pa1, ",", pa2, ")", sep="")
titname[11]<-paste("Rayleigh(", re1, ")", sep="")
titname[12]<-paste("von Mises(", vm1, ",", vm2, ")", sep="")
titname[13]<-paste("Borel(", bo1, ",", bo2, ")", sep="")
titname[14]<-paste("Hermite(", h1, ",", h2, ")", sep="")
titname[15]<-paste("Rice(", r1, ",", r2, ")", sep="")
titname[16]<-paste("Zipf(", z1, ",", z2, ")", sep="")


round(skewness(data), 2)
round(mean(data), 2)
round(sd(data), 2)
round(kurtosis(data), 2)

mmr<- as.matrix(mr)

pdf("Distributions.pdf", height = 8, width = 8)
par(mfrow=c(4,4), mar = c(5.1, 3, 4.1, 2))
for (i in 1:16){
  data<- mmr[1:1000000,i]
  hist(data, probability = TRUE, col = "lightgrey", axes = FALSE,
       main = "", xlab = "",  ylab = "")
  
  # X-axis
  axis(1)
  axis(2)
  
  # Density
  lines(density(data), lwd = 2, col = rainbow(16)[i])
  
  # Add violin plot
  par(new = TRUE)
  vioplot(data, horizontal = TRUE, yaxt = "n", axes = FALSE,
          col = rgb(0, 1, 1, alpha = 0.15))
  title(paste(titname[i]), xlab = paste("M=", round(mean(data), 2),", ",  "SD=", round(sd(data), 2), sep=""), paste("S=", round(skewness(data), 2),", ", "K=", round(kurtosis(data), 2), sep=""))
}
dev.off()
dev.off()
