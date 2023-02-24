#### Computations ####

library("oro.nifti")


fname = "C:/Users/matth/OneDrive/Master Statistics/Masterarbeit/2. Datenübergabe 2023-02-20/Gliomsimulation/Allocortex/1/wSEG^FLAIR_20040513_533^533^_533.nii"
(data_pid_533_t1 <- readNIfTI(fname))
aux.file(data_pid_533_t1)
descrip(data_pid_533_t1)

file_pid_533_t2 = "C:/Users/matth/OneDrive/Master Statistics/Masterarbeit/2. Datenübergabe 2023-02-20/Gliomsimulation/Allocortex/2/wSEG^FLAIR_20061204_533^533^_533"
(data_pid_533_t2 <- readNIfTI(file_pid_533_t2))
aux.file(data_pid_533_t2)
descrip(data_pid_533_t2)


voxel_matrix_data_pid_533_t1 = data_pid_533_t1@.Data

# loop over 
sum_voxel = 0
for(dim_1 in 1:dim(voxel_matrix_data_pid_533_t1)[1]){
  for(dim_2 in 1:dim(voxel_matrix_data_pid_533_t1)[2]){
    for(dim_3 in 1:dim(voxel_matrix_data_pid_533_t1)[3]){
      sum_voxel = sum_voxel+1
      #print(voxel_matrix_data_pid_533_t1[dim_1,1:3,1:3])
    }
  }
}

sum_voxel
181 *217 *181

loop_through_matrix <- function(matrix_3d) {
  for (i in 1:nrow(matrix_3d)) {
    for (j in 1:ncol(matrix_3d)) {
      for (k in 1:dim(matrix_3d)[3]) {
        # Do something with matrix_3d[i, j, k]
        print(matrix_3d[i, j, k])

      }
    }
  }
}


loop_through_matrix_surrounding <- function(matrix_3d) {
  for (i in 1:nrow(matrix_3d)) {
    for (j in 1:ncol(matrix_3d)) {
      for (k in 1:dim(matrix_3d)[3]) {
        # Get the indices of the surrounding entries
        i_surrounding <- max(1, i-1):min(nrow(matrix_3d), i+1)
        j_surrounding <- max(1, j-1):min(ncol(matrix_3d), j+1)
        k_surrounding <- max(1, k-1):min(dim(matrix_3d)[3], k+1)
        
        print(matrix_3d[i, j, k])
        print(":")
        
        # Loop through the surrounding entries
        for (i2 in i_surrounding) {
          for (j2 in j_surrounding) {
            for (k2 in k_surrounding) {
              # Do something with matrix_3d[i2, j2, k2]
              print(matrix_3d[i2, j2, k2])
              
              
            }
          }
        }
        cat("\n\n")
        
      }
    }
  }
}

loop_through_matrix_surrounding <- function(matrix_3d) {
  for (i in 1:nrow(matrix_3d)) {
    for (j in 1:ncol(matrix_3d)) {
      for (k in 1:dim(matrix_3d)[3]) {
        # Get the indices of the surrounding entries
        i_surrounding <- max(1, i-1):min(nrow(matrix_3d), i+1)
        j_surrounding <- max(1, j-1):min(ncol(matrix_3d), j+1)
        k_surrounding <- max(1, k-1):min(dim(matrix_3d)[3], k+1)
        
        # Print out the current entry and its surrounding entries
        cat("Current entry:", matrix_3d[i, j, k], "\n")
        cat("Surrounding entries:\n")
        for (i2 in i_surrounding) {
          for (j2 in j_surrounding) {
            for (k2 in k_surrounding) {
              cat(matrix_3d[i2, j2, k2], " ")
           #   cat("i2: ",i2,", j2: ", j2, ", k2:", k2, " ")
            }
            cat("\n")
          }
          cat("\n")
        }
        cat("--------------------------\n")
      }
    }
  }
}

matrix_3d <- array(1:27, dim=c(3, 3, 3))
loop_through_matrix_surrounding(matrix_3d)


subset = voxel_matrix_data_pid_533_t1[105:110,105:110,55:60]

any(subset==1)
loop_through_matrix_surrounding(subset)
















######## Function for Reducing the Resolution by Downsampling ########


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


red_matrix = matrix3d_downsample(data_pid_533_t1@.Data,3)
red_matrix
