#### Visualize Data ####

library("oro.nifti")

library(ggplot2)

library(rgl)

####################################################
### Glioma segmentation
###################################################
#line 572 "dicom_nifti3.Rnw"
fname = "C:/Users/matth/OneDrive/Master Statistics/Masterarbeit/2. Datenübergabe 2023-02-20/Gliomsimulation/Allocortex/1/wSEG^FLAIR_20040513_533^533^_533.nii"
(data_pid_533_t1 <- readNIfTI(fname))
aux.file(data_pid_533_t1)
descrip(data_pid_533_t1)

file_pid_533_t2 = "C:/Users/matth/OneDrive/Master Statistics/Masterarbeit/2. Datenübergabe 2023-02-20/Gliomsimulation/Allocortex/2/wSEG^FLAIR_20061204_533^533^_533"
(data_pid_533_t2 <- readNIfTI(file_pid_533_t2))
aux.file(data_pid_533_t2)
descrip(data_pid_533_t2)

dim(data_pid_533_t1@.Data)



dim_1_from = 50
dim_1_to = 170
dim_2_from = 30
dim_2_to = 160
dim_3_from = 100
dim_3_to = 100


image(data_pid_533_t1@.Data[dim_1_from:dim_1_to,
                            dim_2_from:dim_2_to,
                            dim_3_from:dim_3_to])



longData<-reshape2::melt(data_pid_533_t1@.Data[,,100])
longData<-longData[longData$value!=0,]

ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="Var2", y="Var1", title="Visualization of Glioma") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11)) +
  xlim(0,200) +
  ylim(0,180)





#### FUNCTION ####

    plot_glioma = function(file_NIfTI,
                           fixed_axis = 3, 
                           fixed_axis_val = 100,
                           omit_zeros = TRUE){
      
      if(fixed_axis == 1) longData<-reshape2::melt(file_NIfTI@.Data[fixed_axis_val,,])
      if(fixed_axis == 2) longData<-reshape2::melt(file_NIfTI@.Data[,fixed_axis_val,])
      if(fixed_axis == 3) longData<-reshape2::melt(file_NIfTI@.Data[,,fixed_axis_val])
      
      if(omit_zeros) longData<-longData[longData$value!=0,]
      
      ggplot(longData, aes(x = Var2, y = Var1)) + 
        geom_raster(aes(fill=value)) + 
        scale_fill_gradient(low="grey90", high="red") +
        labs(x="Var2", y="Var1", title="Visualization of Glioma") +
        theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                           axis.text.y=element_text(size=9),
                           plot.title=element_text(size=11)) +
        xlim(0,200) +
        ylim(0,180)
      
      
    }
    
    plot_glioma(data_pid_533_t1)
    plot_glioma(data_pid_533_t2)
    
    plot_glioma(data_pid_533_t1,fixed_axis_val = 106)
    plot_glioma(data_pid_533_t2,fixed_axis_val = 106)


####################################################
### MNI
###################################################
#line 572 "dicom_nifti3.Rnw"
fname = "C:/Users/matth/OneDrive/Master Statistics/Masterarbeit/2. Datenübergabe 2023-02-20/Gliomsimulation/MNI/mni_icbm152_nlin_sym_09c/mni_icbm152_csf_tal_nlin_sym_09c.nii"
(mni <- readNIfTI(fname))
aux.file(mni)
descrip(mni)


dim(mni@.Data)



dim_1_from = 50
dim_1_to = 170
dim_2_from = 30
dim_2_to = 160
dim_3_from = 100
dim_3_to = 100


image(mni@.Data[dim_1_from:dim_1_to,
                            dim_2_from:dim_2_to,
                            dim_3_from:dim_3_to])



longData<-reshape2::melt(mni@.Data[,,100])
longData<-longData[longData$value!=0,]

ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="Var2", y="Var1", title="Visualization of Glioma") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11)) +
  xlim(0,200) +
  ylim(0,180)
#  xlim(0,229) +
#  ylim(0,193)










####################################################
### ChatGPT
###################################################



library(rgl)

matrix_3d <- array(1:27, dim=c(3, 3, 3))

rnomr

# Create a 3D scatter plot with colored points
plot3d(matrix_3d, col=rainbow(length(matrix_3d)), size=10, type="s")

# Add labels for each point
text3d(x=rep(1:3, each=9), y=rep(1:3, times=3), z=rep(1:3, each=9),
       labels=matrix_3d, adj=c(0.5, 0.5), cex=2)




library(rgl)

# Generate random multivariate normal data
set.seed(123)
n <- 40
mu <- c(0, 0, 0)
sigma <- matrix(c(1, .5, 0.2,
                  0.5, 1, 0.4,
                  0.2, .5, 1), nrow=3, ncol=3)
data <- MASS::mvrnorm(n, mu, sigma)
data

# Plot the data using 3D scatter plot
plot3d(data, col="blue", size=2, type="s")





dim_1_from = 50
dim_1_to = 170
dim_2_from = 30
dim_2_to = 160
dim_3_from = 100
dim_3_to = 100

longData2<-reshape2::melt(data_pid_533_t1@.Data[50:170,30:160,50:150])
longData2
longData2<-longData2[longData2$value!=0,]
longData2
plot3d(longData2, col="blue", size=.1, type="s")




    
# reduce the resolution
    
    matrix3d_downsample <- function(matrix, factor) {
      # Compute the new dimensions of the downsampled matrix
      nx <- ceiling(dim(matrix)[1]/factor)
      ny <- ceiling(dim(matrix)[2]/factor)
      nz <- ceiling(dim(matrix)[3]/factor)
      
      # Downsample the matrix
      downsampled_matrix <- array(0, dim=c(nx, ny, nz))
      for (i in 1:nx) {
        for (j in 1:ny) {
          for (k in 1:nz) {
            x <- (i-1)*factor + 1
            y <- (j-1)*factor + 1
            z <- (k-1)*factor + 1
            downsampled_matrix[i,j,k] <- mean(matrix[x:min(x+factor-1, dim(matrix)[1]),
                                                     y:min(y+factor-1, dim(matrix)[2]),
                                                     z:min(z+factor-1, dim(matrix)[3])])
          }
        }
      }
      
      # Return the downsampled matrix
      return(downsampled_matrix)
    }    



# plot in 3d
        
    
    plot_glioma_3d = function(data, resolution_factor=3){
      red_matrix = matrix3d_downsample(data,resolution_factor)
      
      longData2<-reshape2::melt(red_matrix)#data_pid_533_t1@.Data[50:170,30:160,50:150])
      longData2
      longData2<-longData2[longData2$value!=0,]
      longData2
      plot3d(longData2, col="red", size=.5, type="s")
    }
    
    plot_glioma_3d(data_pid_533_t1@.Data)
    plot_glioma_3d(data_pid_533_t2@.Data)
    plot_glioma_3d((data_pid_533_t1@.Data-data_pid_533_t2@.Data))
    


