
source('kernelridge.R')


x <- seq(0, 1, by = 0.005)
y <- seq(0, 1, by = 0.005)
z = matrix(NA, length(x), length(y))
for(i in 1:nrow(z)){
  for(j in 1:nrow(z)){
    z[i,j] = laplacian_kernel(x[i], x[j], 0.1)
  }
}
filled.contour(x, y, z, main = "Laplacian(0.1)", cex.main = 2)



x <- seq(0, 1, by = 0.005)
y <- seq(0, 1, by = 0.005)
z = matrix(NA, length(x), length(y))
for(i in 1:nrow(z)){
  for(j in 1:nrow(z)){
    z[i,j] = laplacian_kernel(x[i], x[j], 1)
  }
}
filled.contour(x, y, z, main = "Laplacian(1)", cex.main = 2)



x <- seq(0, 1, by = 0.005)
y <- seq(0, 1, by = 0.005)
z = matrix(NA, length(x), length(y))
for(i in 1:nrow(z)){
  for(j in 1:nrow(z)){
    z[i,j] = gaussian_kernel(x[i], x[j], 0.1)
  }
}
filled.contour(x, y, z, main = "Gaussian(0.1)", cex.main = 2)



x <- seq(0, 1, by = 0.005)
y <- seq(0, 1, by = 0.005)
z = matrix(NA, length(x), length(y))
for(i in 1:nrow(z)){
  for(j in 1:nrow(z)){
    z[i,j] = gaussian_kernel(x[i], x[j], 0.05)
  }
}
filled.contour(x, y, z, main = "Gaussian(0.05)", cex.main = 2)





# x <- seq(0, 1, by = 0.005)
# y <- seq(0, 1, by = 0.005)
# z = matrix(NA, length(x), length(y))
# for(i in 1:nrow(z)){
#   for(j in 1:nrow(z)){
#     z[i,j] = polynomial_kernel(x[i], x[j], 2)
#   }
# }
# filled.contour(x, y, z, main = "Polynomial(degree = 2)")
# 
# 
# 
# x <- seq(0, 1, by = 0.005)
# y <- seq(0, 1, by = 0.005)
# z = matrix(NA, length(x), length(y))
# for(i in 1:nrow(z)){
#   for(j in 1:nrow(z)){
#     z[i,j] = polynomial_kernel(x[i], x[j], 5)
#   }
# }
# filled.contour(x, y, z, main = "Polynomial(degree = 5)")





x <- seq(0, 1, by = 0.005)
y <- seq(0, 1, by = 0.005)
z = matrix(NA, length(x), length(y))
for(i in 1:nrow(z)){
  for(j in 1:nrow(z)){
    z[i,j] = periodic_kernel(x[i], x[j], 1)
  }
}
filled.contour(x, y, z, main = "Periodic(1)", cex.main = 2)



x <- seq(0, 1, by = 0.005)
y <- seq(0, 1, by = 0.005)
z = matrix(NA, length(x), length(y))
for(i in 1:nrow(z)){
  for(j in 1:nrow(z)){
    z[i,j] = periodic_kernel(x[i], x[j], 0.5)
  }
}
filled.contour(x, y, z, main = "Periodic(0.5)", cex.main = 2)
