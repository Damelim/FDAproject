source('kernelridge.R')

load('Simulationdata.rdata')
n = 800  
set.seed(1) ; x = sort(runif(n, 0, 1)) ; x_pred = seq(0.001, 0.999, by = 0.002) ## function evaluation points





#### gamma = 0.1 ####

K = generate_kernel_matrix(x, kernel_type = "laplacian", params = 1e-1)
K_pred = generate_prediction_kernel_matrix(trainx = x, testx = x_pred, kernel_type = "laplacian", params = 1e-1)

lam_vec = c(1e-10, 0.0001, 0.1, 1, 100)

Smoothermat1 = eigenMapMatMult(eigenMapMatMult(K, armaInv(eigenMapMatMult(t(K), K) + lam_vec[1] * diag(1, n))) , t(K) )
Smoothermat2 = eigenMapMatMult(eigenMapMatMult(K, armaInv(eigenMapMatMult(t(K), K) + lam_vec[2] * diag(1, n))) , t(K) )
Smoothermat3 = eigenMapMatMult(eigenMapMatMult(K, armaInv(eigenMapMatMult(t(K), K) + lam_vec[3] * diag(1, n))) , t(K) )
Smoothermat4 = eigenMapMatMult(eigenMapMatMult(K, armaInv(eigenMapMatMult(t(K), K) + lam_vec[4] * diag(1, n))) , t(K) )
Smoothermat5 = eigenMapMatMult(eigenMapMatMult(K, armaInv(eigenMapMatMult(t(K), K) + lam_vec[5] * diag(1, n))) , t(K) )

eigS1 = eigen(Smoothermat1);eigS2 = eigen(Smoothermat2);eigS3 = eigen(Smoothermat3) ;eigS4 = eigen(Smoothermat4) ;eigS5 = eigen(Smoothermat5)

evalS1 = as.numeric(eigS1$values);evalS2 = as.numeric(eigS2$values);evalS3 = as.numeric(eigS3$values) ;evalS4 = as.numeric(eigS4$values) ;evalS5= as.numeric(eigS5$values)

evecS1 = eigS1$vectors ; evecS2 = eigS2$vectors ; evecS3 = eigS3$vectors ; evecS4 = eigS4$vectors ; evecS5 = eigS5$vectors 
for(i in 1:n){
  for(j in 1:n){
    evecS1[i,j] = as.numeric(evecS1[i,j]) # in case the eigenvectors contain imaginary parts (by numerical error)
    evecS2[i,j] = as.numeric(evecS2[i,j]) # in case the eigenvectors contain imaginary parts (by numerical error)
    evecS3[i,j] = as.numeric(evecS3[i,j]) # in case the eigenvectors contain imaginary parts (by numerical error)
    evecS4[i,j] = as.numeric(evecS4[i,j]) # in case the eigenvectors contain imaginary parts (by numerical error)
    evecS5[i,j] = as.numeric(evecS5[i,j]) # in case the eigenvectors contain imaginary parts (by numerical error)
  }
}



png('G:/My Drive/FunctionalDataAnalysis/Project/kernelridge/Laplacian(0.1)_prediction.png', width = 15, height = 8, units = "in", res = 150)

par(mfrow = c(1,3))
lwd_vec = c(4, 0.3, 2, 2, 2, 2)

train_sigma = train_mediumsigma
test_Y_signal = f8(x_pred)
y = f8(x) + rnorm(n, mean = 0, sd = train_sigma)
mny = mean(y)
ycentered = y - mny

ypred_1 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[1] * diag(1, n))) , t(K) )  , ycentered)
ypred_2 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[2] * diag(1, n))) , t(K) )  , ycentered)
ypred_3 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[3] * diag(1, n))) , t(K) )  , ycentered)
ypred_4 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[4] * diag(1, n))) , t(K) )  , ycentered)
ypred_5 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[5] * diag(1, n))) , t(K) )  , ycentered)
# ypred_1 = eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[1] * diag(1, n))) , t(K) )  , y)

plot(x, y, xlab = "", ylab ="", cex = 0.4, main = paste("Laplacian(0.1)- prediction, rank(K) =",Matrix::rankMatrix(K)[1],', n=',n,', sigma=',train_sigma), cex.main = 1.25)
lines(x_pred, test_Y_signal, col = "black", lwd = 4)
lines(x_pred, ypred_1 , col = "rosybrown1", lwd = 0.3)
lines(x_pred, ypred_2, col = "orange", lwd = 2)
lines(x_pred, ypred_3, col = "green", lwd = 2)
lines(x_pred, ypred_4, col = "blue", lwd = 2)
lines(x_pred, ypred_5, col = "purple", lwd = 2)
legend('topright',c('true',paste('lam = ', lam_vec)), col = c('black',"rosybrown1",'orange','green','blue','purple'), lwd = lwd_vec)



train_sigma = train_mediumsigma
test_Y_signal = f17(x_pred)
y = f17(x) + rnorm(n, mean = 0, sd = train_sigma)
mny = mean(y)
ycentered = y - mny

ypred_1 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[1] * diag(1, n))) , t(K) )  , ycentered)
ypred_2 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[2] * diag(1, n))) , t(K) )  , ycentered)
ypred_3 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[3] * diag(1, n))) , t(K) )  , ycentered)
ypred_4 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[4] * diag(1, n))) , t(K) )  , ycentered)
ypred_5 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[5] * diag(1, n))) , t(K) )  , ycentered)
# ypred_1 = eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[1] * diag(1, n))) , t(K) )  , y)

plot(x, y, xlab = "", ylab ="", cex = 0.4, main = paste("Laplacian(0.1)- prediction, rank(K) =",Matrix::rankMatrix(K)[1],', n=',n,', sigma=',train_sigma), cex.main = 1.25)
lines(x_pred, test_Y_signal, col = "black", lwd = 4)
lines(x_pred, ypred_1 , col = "rosybrown1", lwd = 0.3)
lines(x_pred, ypred_2, col = "orange", lwd = 2)
lines(x_pred, ypred_3, col = "green", lwd = 2)
lines(x_pred, ypred_4, col = "blue", lwd = 2)
lines(x_pred, ypred_5, col = "purple", lwd = 2)
legend('topright',c('true',paste('lam = ', lam_vec)), col = c('black',"rosybrown1",'orange','green','blue','purple'), lwd = lwd_vec)




train_sigma = train_mediumsigma
test_Y_signal = f12(x_pred)
y = f12(x) + rnorm(n, mean = 0, sd = 0.1)
mny = mean(y)
ycentered = y - mny

ypred_1 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[1] * diag(1, n))) , t(K) )  , ycentered)
ypred_2 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[2] * diag(1, n))) , t(K) )  , ycentered)
ypred_3 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[3] * diag(1, n))) , t(K) )  , ycentered)
ypred_4 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[4] * diag(1, n))) , t(K) )  , ycentered)
ypred_5 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[5] * diag(1, n))) , t(K) )  , ycentered)
# ypred_1 = eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[1] * diag(1, n))) , t(K) )  , y)

plot(x, y, xlab = "", ylab ="", cex = 0.4, main = paste("Laplacian(0.1)- prediction, rank(K) =",Matrix::rankMatrix(K)[1],', n=',n,', sigma=', 0.1), cex.main = 1.25)
lines(x_pred, test_Y_signal, col = "black", lwd = 4)
lines(x_pred, ypred_1 , col = "rosybrown1", lwd = 0.3)
lines(x_pred, ypred_2, col = "orange", lwd = 2)
lines(x_pred, ypred_3, col = "green", lwd = 2)
lines(x_pred, ypred_4, col = "blue", lwd = 2)
lines(x_pred, ypred_5, col = "purple", lwd = 2)
legend('topright',c('true',paste('lam = ', lam_vec)), col = c('black',"rosybrown1",'orange','green','blue','purple'), lwd = lwd_vec)

dev.off()














png('G:/My Drive/FunctionalDataAnalysis/Project/kernelridge/Laplacian(0.1)_eigenanalysis.png', width = 15, height = 8, units = "in", res = 150)

m1 = matrix(c(
  1,1,2,3,4,
  1,1,5,6,7,
  1,1,8,9,10), nrow = 3, ncol = 5, byrow = T
)
layout(m1)



plot(1:25, seq(-0.2,1.2, length = 25), type = "n", xaxt = "n", xlab = "order", cex.axis = 1.5, cex.lab = 1.5, ylab = "eigenvalues", main = paste("Laplacian(0.1)-First 25 eigenvalues of smoother matrix, rank(K) =", Matrix::rankMatrix(K)[1]), cex.main = 1.25)
abline(h = 1, lty = 2)
abline(h = 0, lty = 2)
points(1:25,evalS1[1:25], type = "o", col = "rosybrown1", pch = 16)
points(1:25,evalS2[1:25], type = "o", col = "orange", pch = 16)
points(1:25,evalS3[1:25], type = "o", col = "green", pch = 16)
points(1:25,evalS4[1:25], type = "o", col = "blue", pch = 16)
points(1:25,evalS5[1:25], type = "o", col = "purple", pch = 16)
# legend('topright',paste('lam = ', lam_vec), col = c("rosybrown1",'orange','green','blue','purple'), lwd = 2)
axis(1, xaxp = c(1,25,24), las = 2)





orders = 1:9

plot(x, evalS1[orders[1]] * evecS1[,orders[1]] + rep(-0.001, n), col = "rosybrown1", lwd = 0.1, type = "l", xlab = "", ylab = "", main = "1st eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, -evalS2[orders[1]] * evecS2[,orders[1]] + rep(-0.001, n), lwd = 1, col = "orange")
lines(x, evalS3[orders[1]] * evecS3[,orders[1]] + rep(-0.002, n), lwd = 1, col = "green")
lines(x, evalS4[orders[1]] * evecS4[,orders[1]] + rep(0.002, n), lwd = 1,col = "blue")
lines(x, evalS5[orders[1]] * evecS5[,orders[1]], col = "purple")
# legend('topright',paste('lam = ', lam_vec), col = c("rosybrown1",'orange','green','blue','purple'), lwd = 2)

plot(x, evalS1[orders[2]] * evecS1[,orders[2]] + rep(0.002, n), col = "rosybrown1", type = "l", xlab = "", ylab = "", main = "2nd eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, -evalS2[orders[2]] * evecS2[,orders[2]] , col = "orange")
lines(x, -evalS3[orders[2]] * evecS3[,orders[2]] , col = "green")
lines(x, evalS4[orders[2]] * evecS4[,orders[2]] , col = "blue")
lines(x, evalS5[orders[2]] * evecS5[,orders[2]], col = "purple")
# legend('topright',paste('lam = ', lam_vec), col = c("rosybrown1",'orange','green','blue','purple'), lwd = 2)

plot(x, -evalS1[orders[3]] * evecS1[,orders[3]] + rep(0.002, n), col = "rosybrown1", type = "l", xlab = "", ylab = "", main = "3rd eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, -evalS2[orders[3]] * evecS2[,orders[3]] + rep(-0.001, n), col = "orange")
lines(x, -evalS3[orders[3]] * evecS3[,orders[3]] , col = "green")
lines(x, -evalS4[orders[3]] * evecS4[,orders[3]] , col = "blue")
lines(x, -evalS5[orders[3]] * evecS5[,orders[3]], col = "purple")
# legend('topright',paste('lam = ', lam_vec), col = c("rosybrown1",'orange','green','blue','purple'), lwd = 2)

plot(x, -evalS1[orders[4]] * evecS1[,orders[4]] + rep(0.002, n), col = "rosybrown1", type = "l", xlab = "",ylab = "", main = "4th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, evalS2[orders[4]] * evecS2[,orders[4]] + rep(-0.001, n), col = "orange")
lines(x, -evalS3[orders[4]] * evecS3[,orders[4]], col = "green")
lines(x, -evalS4[orders[4]] * evecS4[,orders[4]], col = "blue")
lines(x, -evalS5[orders[4]] * evecS5[,orders[4]], col = "purple")
# legend('topright',paste('lam = ', lam_vec), col = c("rosybrown1",'orange','green','blue','purple'), lwd = 2)

plot(x, evalS1[orders[5]] * evecS1[,orders[5]] + rep(0.002, n), col = "rosybrown1", type = "l", xlab = "",ylab = "", main = "5th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, evalS2[orders[5]] * evecS2[,orders[5]] + rep(-0.001, n), col = "orange")
lines(x, evalS3[orders[5]] * evecS3[,orders[5]], col = "green")
lines(x, -evalS4[orders[5]] * evecS4[,orders[5]], col = "blue")
lines(x,-evalS5[orders[5]] * evecS5[,orders[5]], col = "purple")
# legend('topright',paste('lam = ', lam_vec), col = c("rosybrown1",'orange','green','blue','purple'), lwd = 2)

plot(x, -evalS1[orders[6]] * evecS1[,orders[6]] + rep(0.002, n), col = "rosybrown1", type = "l", xlab = "",ylab = "", main = "6th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, evalS2[orders[6]] * evecS2[,orders[6]] + rep(-0.001, n), col = "orange")
lines(x,evalS3[orders[6]] * -evecS3[,orders[6]], col = "green")
lines(x, -evalS4[orders[6]] * evecS4[,orders[6]], col = "blue")
lines(x, -evalS5[orders[6]] * evecS5[,orders[6]], col = "purple")
# legend('topright',paste('lam = ', lam_vec), col = c("rosybrown1",'orange','green','blue','purple'), lwd = 2)

plot(x, evalS1[orders[7]] * evecS1[,orders[7]] + rep(0.002, n), col = "rosybrown1", type = "l", xlab = "",ylab = "", main = "7th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, evalS2[orders[7]] * evecS2[,orders[7]] + rep(-0.001, n), col = "orange")
lines(x, -evalS3[orders[7]] * evecS3[,orders[7]], col = "green")
lines(x, evalS4[orders[7]] * evecS4[,orders[7]], col = "blue")
lines(x, evalS5[orders[7]] * evecS5[,orders[7]], col = "purple")
# legend('topright',paste('lam = ', lam_vec), col = c("rosybrown1",'orange','green','blue','purple'), lwd = 2)

plot(x, -evalS1[orders[8]] * evecS1[,orders[8]] + rep(0.002, n), col = "rosybrown1", type = "l", xlab = "",ylab = "", main = "8th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, evalS2[orders[8]] * evecS2[,orders[8]] + rep(-0.001, n), col = "orange")
lines(x, evalS3[orders[8]] * evecS3[,orders[8]], col = "green")
lines(x, -evalS4[orders[8]] * evecS4[,orders[8]], col = "blue")
lines(x, -evalS5[orders[8]] * evecS5[,orders[8]], col = "purple")
# legend('topright',paste('lam = ', lam_vec), col = c("rosybrown1",'orange','green','blue','purple'), lwd = 2)

plot(x, -evalS1[orders[9]] * evecS1[,orders[9]] + rep(0.002, n), col = "rosybrown1", type = "l", xlab = "",ylab = "", main = "9th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, evalS2[orders[9]] * evecS2[,orders[9]] + rep(-0.001, n), col = "orange")
lines(x, evalS3[orders[9]] * evecS3[,orders[9]], col = "green")
lines(x, -evalS4[orders[9]] * evecS4[,orders[9]], col = "blue")
lines(x, evalS5[orders[9]] * evecS5[,orders[9]], col = "purple")
# legend('topright',paste('lam = ', lam_vec), col = c("rosybrown1",'orange','green','blue','purple'), lwd = 2)

dev.off()





























































K = generate_kernel_matrix(x, kernel_type = "laplacian", params = 1)
K_pred = generate_prediction_kernel_matrix(trainx = x, testx = x_pred, kernel_type = "laplacian", params = 1)

lam_vec = c(1e-10, 0.0001, 0.1, 1, 100)

Smoothermat1 = eigenMapMatMult(eigenMapMatMult(K, armaInv(eigenMapMatMult(t(K), K) + lam_vec[1] * diag(1, n))) , t(K) )
Smoothermat2 = eigenMapMatMult(eigenMapMatMult(K, armaInv(eigenMapMatMult(t(K), K) + lam_vec[2] * diag(1, n))) , t(K) )
Smoothermat3 = eigenMapMatMult(eigenMapMatMult(K, armaInv(eigenMapMatMult(t(K), K) + lam_vec[3] * diag(1, n))) , t(K) )
Smoothermat4 = eigenMapMatMult(eigenMapMatMult(K, armaInv(eigenMapMatMult(t(K), K) + lam_vec[4] * diag(1, n))) , t(K) )
Smoothermat5 = eigenMapMatMult(eigenMapMatMult(K, armaInv(eigenMapMatMult(t(K), K) + lam_vec[5] * diag(1, n))) , t(K) )

eigS1 = eigen(Smoothermat1);eigS2 = eigen(Smoothermat2);eigS3 = eigen(Smoothermat3) ;eigS4 = eigen(Smoothermat4) ;eigS5 = eigen(Smoothermat5)

evalS1 = as.numeric(eigS1$values);evalS2 = as.numeric(eigS2$values);evalS3 = as.numeric(eigS3$values) ;evalS4 = as.numeric(eigS4$values) ;evalS5= as.numeric(eigS5$values)

evecS1 = eigS1$vectors ; evecS2 = eigS2$vectors ; evecS3 = eigS3$vectors ; evecS4 = eigS4$vectors ; evecS5 = eigS5$vectors 
for(i in 1:n){
  for(j in 1:n){
    evecS1[i,j] = as.numeric(evecS1[i,j]) # in case the eigenvectors contain imaginary parts (by numerical error)
    evecS2[i,j] = as.numeric(evecS2[i,j]) # in case the eigenvectors contain imaginary parts (by numerical error)
    evecS3[i,j] = as.numeric(evecS3[i,j]) # in case the eigenvectors contain imaginary parts (by numerical error)
    evecS4[i,j] = as.numeric(evecS4[i,j]) # in case the eigenvectors contain imaginary parts (by numerical error)
    evecS5[i,j] = as.numeric(evecS5[i,j]) # in case the eigenvectors contain imaginary parts (by numerical error)
  }
}



png('G:/My Drive/FunctionalDataAnalysis/Project/kernelridge/Laplacian(1)_prediction.png', width = 15, height = 8, units = "in", res = 150)

par(mfrow = c(1,3))
lwd_vec = c(4, 0.3, 2, 2, 2, 2)

train_sigma = train_mediumsigma
test_Y_signal = f8(x_pred)
y = f8(x) + rnorm(n, mean = 0, sd = train_sigma)
mny = mean(y)
ycentered = y - mny

ypred_1 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[1] * diag(1, n))) , t(K) )  , ycentered)
ypred_2 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[2] * diag(1, n))) , t(K) )  , ycentered)
ypred_3 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[3] * diag(1, n))) , t(K) )  , ycentered)
ypred_4 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[4] * diag(1, n))) , t(K) )  , ycentered)
ypred_5 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[5] * diag(1, n))) , t(K) )  , ycentered)
# ypred_1 = eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[1] * diag(1, n))) , t(K) )  , y)

plot(x, y, xlab = "", ylab ="", cex = 0.4, main = paste("Laplacian(1)- prediction, rank(K) =",Matrix::rankMatrix(K)[1],', n=',n,', sigma=',train_sigma), cex.main = 1.25)
lines(x_pred, test_Y_signal, col = "black", lwd = 4)
lines(x_pred, ypred_1 , col = "rosybrown1", lwd = 0.3)
lines(x_pred, ypred_2, col = "orange", lwd = 2)
lines(x_pred, ypred_3, col = "green", lwd = 2)
lines(x_pred, ypred_4, col = "blue", lwd = 2)
lines(x_pred, ypred_5, col = "purple", lwd = 2)
legend('topright',c('true',paste('lam = ', lam_vec)), col = c('black',"rosybrown1",'orange','green','blue','purple'), lwd = lwd_vec)



train_sigma = train_mediumsigma
test_Y_signal = f17(x_pred)
y = f17(x) + rnorm(n, mean = 0, sd = train_sigma)
mny = mean(y)
ycentered = y - mny

ypred_1 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[1] * diag(1, n))) , t(K) )  , ycentered)
ypred_2 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[2] * diag(1, n))) , t(K) )  , ycentered)
ypred_3 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[3] * diag(1, n))) , t(K) )  , ycentered)
ypred_4 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[4] * diag(1, n))) , t(K) )  , ycentered)
ypred_5 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[5] * diag(1, n))) , t(K) )  , ycentered)
# ypred_1 = eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[1] * diag(1, n))) , t(K) )  , y)

plot(x, y, xlab = "", ylab ="", cex = 0.4, main = paste("Laplacian(1)- prediction, rank(K) =",Matrix::rankMatrix(K)[1],', n=',n,', sigma=',train_sigma), cex.main = 1.25)
lines(x_pred, test_Y_signal, col = "black", lwd = 4)
lines(x_pred, ypred_1 , col = "rosybrown1", lwd = 0.3)
lines(x_pred, ypred_2, col = "orange", lwd = 2)
lines(x_pred, ypred_3, col = "green", lwd = 2)
lines(x_pred, ypred_4, col = "blue", lwd = 2)
lines(x_pred, ypred_5, col = "purple", lwd = 2)
legend('topright',c('true',paste('lam = ', lam_vec)), col = c('black',"rosybrown1",'orange','green','blue','purple'), lwd = lwd_vec)




train_sigma = train_mediumsigma
test_Y_signal = f12(x_pred)
y = f12(x) + rnorm(n, mean = 0, sd = 0.1)
mny = mean(y)
ycentered = y - mny

ypred_1 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[1] * diag(1, n))) , t(K) )  , ycentered)
ypred_2 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[2] * diag(1, n))) , t(K) )  , ycentered)
ypred_3 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[3] * diag(1, n))) , t(K) )  , ycentered)
ypred_4 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[4] * diag(1, n))) , t(K) )  , ycentered)
ypred_5 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[5] * diag(1, n))) , t(K) )  , ycentered)
# ypred_1 = eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[1] * diag(1, n))) , t(K) )  , y)

plot(x, y, xlab = "", ylab ="", cex = 0.4, main = paste("Laplacian(1)- prediction, rank(K) =",Matrix::rankMatrix(K)[1],', n=',n,', sigma=', 0.1), cex.main = 1.25)
lines(x_pred, test_Y_signal, col = "black", lwd = 4)
lines(x_pred, ypred_1 , col = "rosybrown1", lwd = 0.3)
lines(x_pred, ypred_2, col = "orange", lwd = 2)
lines(x_pred, ypred_3, col = "green", lwd = 2)
lines(x_pred, ypred_4, col = "blue", lwd = 2)
lines(x_pred, ypred_5, col = "purple", lwd = 2)
legend('topright',c('true',paste('lam = ', lam_vec)), col = c('black',"rosybrown1",'orange','green','blue','purple'), lwd = lwd_vec)

dev.off()














png('G:/My Drive/FunctionalDataAnalysis/Project/kernelridge/Laplacian(1)_eigenanalysis.png', width = 15, height = 8, units = "in", res = 150)

m1 = matrix(c(
  1,1,2,3,4,
  1,1,5,6,7,
  1,1,8,9,10), nrow = 3, ncol = 5, byrow = T
)
layout(m1)



plot(1:25, seq(-0.2,1.2, length = 25), type = "n", xaxt = "n", xlab = "order", cex.axis = 1.5, cex.lab = 1.5, ylab = "eigenvalues", main = paste("Laplacian(0.1)-First 25 eigenvalues of smoother matrix, rank(K) =", Matrix::rankMatrix(K)[1]), cex.main = 1.25)
abline(h = 1, lty = 2)
abline(h = 0, lty = 2)
points(1:25,evalS1[1:25], type = "o", col = "rosybrown1", pch = 16)
points(1:25,evalS2[1:25], type = "o", col = "orange", pch = 16)
points(1:25,evalS3[1:25], type = "o", col = "green", pch = 16)
points(1:25,evalS4[1:25], type = "o", col = "blue", pch = 16)
points(1:25,evalS5[1:25], type = "o", col = "purple", pch = 16)
# legend('topright',paste('lam = ', lam_vec), col = c("rosybrown1",'orange','green','blue','purple'), lwd = 2)
axis(1, xaxp = c(1,25,24), las = 2)





orders = 1:9

plot(x, -evalS1[orders[1]] * evecS1[,orders[1]] + rep(-0.001, n), col = "rosybrown1", lwd = 0.1, type = "l", xlab = "", ylab = "", main = "1st eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, -evalS2[orders[1]] * evecS2[,orders[1]] + rep(-0.001, n), lwd = 1, col = "orange")
lines(x, evalS3[orders[1]] * evecS3[,orders[1]] + rep(-0.002, n), lwd = 1, col = "green")
lines(x, evalS4[orders[1]] * evecS4[,orders[1]] + rep(0.002, n), lwd = 1,col = "blue")
lines(x, evalS5[orders[1]] * evecS5[,orders[1]], col = "purple")
# # legend('topright',paste('lam = ', lam_vec), col = c("rosybrown1",'orange','green','blue','purple'), lwd = 2)

plot(x, evalS1[orders[2]] * evecS1[,orders[2]] + rep(0.002, n), col = "rosybrown1", type = "l", xlab = "", ylab = "", main = "2nd eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, -evalS2[orders[2]] * evecS2[,orders[2]] , col = "orange")
lines(x, -evalS3[orders[2]] * evecS3[,orders[2]] , col = "green")
lines(x, -evalS4[orders[2]] * evecS4[,orders[2]] , col = "blue")
lines(x, evalS5[orders[2]] * evecS5[,orders[2]], col = "purple")
# # legend('topright',paste('lam = ', lam_vec), col = c("rosybrown1",'orange','green','blue','purple'), lwd = 2)

plot(x, -evalS1[orders[3]] * evecS1[,orders[3]] + rep(0.002, n), col = "rosybrown1", type = "l", xlab = "", ylab = "", main = "3rd eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, -evalS2[orders[3]] * evecS2[,orders[3]] + rep(-0.001, n), col = "orange")
lines(x, evalS3[orders[3]] * evecS3[,orders[3]] , col = "green")
lines(x, -evalS4[orders[3]] * evecS4[,orders[3]] , col = "blue")
lines(x, -evalS5[orders[3]] * evecS5[,orders[3]], col = "purple")
# # legend('topright',paste('lam = ', lam_vec), col = c("rosybrown1",'orange','green','blue','purple'), lwd = 2)

plot(x, -evalS1[orders[4]] * evecS1[,orders[4]] + rep(0.002, n), col = "rosybrown1", type = "l", xlab = "",ylab = "", main = "4th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, evalS2[orders[4]] * evecS2[,orders[4]] + rep(-0.001, n), col = "orange")
lines(x, evalS3[orders[4]] * evecS3[,orders[4]], col = "green")
lines(x, evalS4[orders[4]] * evecS4[,orders[4]], col = "blue")
lines(x, -evalS5[orders[4]] * evecS5[,orders[4]], col = "purple")
# # legend('topright',paste('lam = ', lam_vec), col = c("rosybrown1",'orange','green','blue','purple'), lwd = 2)

plot(x, evalS1[orders[5]] * evecS1[,orders[5]] + rep(0.002, n), col = "rosybrown1", type = "l", xlab = "",ylab = "", main = "5th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, -evalS2[orders[5]] * evecS2[,orders[5]] + rep(-0.001, n), col = "orange")
lines(x, evalS3[orders[5]] * evecS3[,orders[5]], col = "green")
lines(x, evalS4[orders[5]] * evecS4[,orders[5]], col = "blue")
lines(x,-evalS5[orders[5]] * evecS5[,orders[5]], col = "purple")
# # legend('topright',paste('lam = ', lam_vec), col = c("rosybrown1",'orange','green','blue','purple'), lwd = 2)

plot(x, -evalS1[orders[6]] * evecS1[,orders[6]] + rep(0.002, n), col = "rosybrown1", type = "l", xlab = "",ylab = "", main = "6th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, -evalS2[orders[6]] * evecS2[,orders[6]] + rep(-0.001, n), col = "orange")
lines(x, -evalS3[orders[6]] * -evecS3[,orders[6]], col = "green")
lines(x, -evalS4[orders[6]] * evecS4[,orders[6]], col = "blue")
lines(x, -evalS5[orders[6]] * evecS5[,orders[6]], col = "purple")
# legend('topright',paste('lam = ', lam_vec), col = c("rosybrown1",'orange','green','blue','purple'), lwd = 2)

plot(x, evalS1[orders[7]] * evecS1[,orders[7]] + rep(0.002, n), col = "rosybrown1", type = "l", xlab = "",ylab = "", main = "7th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, -evalS2[orders[7]] * evecS2[,orders[7]] + rep(-0.001, n), col = "orange")
lines(x, -evalS3[orders[7]] * evecS3[,orders[7]], col = "green")
lines(x, evalS4[orders[7]] * evecS4[,orders[7]], col = "blue")
lines(x, -evalS5[orders[7]] * evecS5[,orders[7]], col = "purple")
# legend('topright',paste('lam = ', lam_vec), col = c("rosybrown1",'orange','green','blue','purple'), lwd = 2)

plot(x, -evalS1[orders[8]] * evecS1[,orders[8]] + rep(0.002, n), col = "rosybrown1", type = "l", xlab = "",ylab = "", main = "8th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, evalS2[orders[8]] * evecS2[,orders[8]] + rep(-0.001, n), col = "orange")
lines(x, -evalS3[orders[8]] * evecS3[,orders[8]], col = "green")
lines(x, -evalS4[orders[8]] * evecS4[,orders[8]], col = "blue")
lines(x, -evalS5[orders[8]] * evecS5[,orders[8]], col = "purple")
# legend('topright',paste('lam = ', lam_vec), col = c("rosybrown1",'orange','green','blue','purple'), lwd = 2)

plot(x, -evalS1[orders[9]] * evecS1[,orders[9]] + rep(0.002, n), col = "rosybrown1", type = "l", xlab = "",ylab = "", main = "9th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, evalS2[orders[9]] * evecS2[,orders[9]] + rep(-0.001, n), col = "orange")
lines(x, evalS3[orders[9]] * evecS3[,orders[9]], col = "green")
lines(x, -evalS4[orders[9]] * evecS4[,orders[9]], col = "blue")
lines(x, evalS5[orders[9]] * evecS5[,orders[9]], col = "purple")
# legend('topright',paste('lam = ', lam_vec), col = c("rosybrown1",'orange','green','blue','purple'), lwd = 2)

dev.off()




































# K = generate_kernel_matrix(x, kernel_type = "laplacian", params = 0.01)
# K_pred = generate_prediction_kernel_matrix(trainx = x, testx = x_pred, kernel_type = "laplacian", params = 0.01)
# 
# lam_vec = c(1e-10, 0.0001, 0.1, 1, 100)
# 
# Smoothermat1 = eigenMapMatMult(eigenMapMatMult(K, armaInv(eigenMapMatMult(t(K), K) + lam_vec[1] * diag(1, n))) , t(K) )
# Smoothermat2 = eigenMapMatMult(eigenMapMatMult(K, armaInv(eigenMapMatMult(t(K), K) + lam_vec[2] * diag(1, n))) , t(K) )
# Smoothermat3 = eigenMapMatMult(eigenMapMatMult(K, armaInv(eigenMapMatMult(t(K), K) + lam_vec[3] * diag(1, n))) , t(K) )
# Smoothermat4 = eigenMapMatMult(eigenMapMatMult(K, armaInv(eigenMapMatMult(t(K), K) + lam_vec[4] * diag(1, n))) , t(K) )
# Smoothermat5 = eigenMapMatMult(eigenMapMatMult(K, armaInv(eigenMapMatMult(t(K), K) + lam_vec[5] * diag(1, n))) , t(K) )
# 
# eigS1 = eigen(Smoothermat1);eigS2 = eigen(Smoothermat2);eigS3 = eigen(Smoothermat3) ;eigS4 = eigen(Smoothermat4) ;eigS5 = eigen(Smoothermat5)
# 
# evalS1 = as.numeric(eigS1$values);evalS2 = as.numeric(eigS2$values);evalS3 = as.numeric(eigS3$values) ;evalS4 = as.numeric(eigS4$values) ;evalS5= as.numeric(eigS5$values)
# 
# evecS1 = eigS1$vectors ; evecS2 = eigS2$vectors ; evecS3 = eigS3$vectors ; evecS4 = eigS4$vectors ; evecS5 = eigS5$vectors 
# for(i in 1:n){
#   for(j in 1:n){
#     evecS1[i,j] = as.numeric(evecS1[i,j]) # in case the eigenvectors contain imaginary parts (by numerical error)
#     evecS2[i,j] = as.numeric(evecS2[i,j]) # in case the eigenvectors contain imaginary parts (by numerical error)
#     evecS3[i,j] = as.numeric(evecS3[i,j]) # in case the eigenvectors contain imaginary parts (by numerical error)
#     evecS4[i,j] = as.numeric(evecS4[i,j]) # in case the eigenvectors contain imaginary parts (by numerical error)
#     evecS5[i,j] = as.numeric(evecS5[i,j]) # in case the eigenvectors contain imaginary parts (by numerical error)
#   }
# }
# 
# 
# 
# png('G:/My Drive/FunctionalDataAnalysis/Project/kernelridge/Laplacian(0.01)_prediction.png', width = 15, height = 8, units = "in", res = 150)
# 
# par(mfrow = c(1,3))
# lwd_vec = c(4, 0.3, 2, 2, 2, 2)
# 
# train_sigma = train_mediumsigma
# test_Y_signal = f8(x_pred)
# y = f8(x) + rnorm(n, mean = 0, sd = train_sigma)
# mny = mean(y)
# ycentered = y - mny
# 
# ypred_1 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[1] * diag(1, n))) , t(K) )  , ycentered)
# ypred_2 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[2] * diag(1, n))) , t(K) )  , ycentered)
# ypred_3 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[3] * diag(1, n))) , t(K) )  , ycentered)
# ypred_4 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[4] * diag(1, n))) , t(K) )  , ycentered)
# ypred_5 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[5] * diag(1, n))) , t(K) )  , ycentered)
# # ypred_1 = eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[1] * diag(1, n))) , t(K) )  , y)
# 
# plot(x, y, xlab = "", ylab ="", cex = 0.4, main = paste("Laplacian(0.01)- prediction, rank(K) =",Matrix::rankMatrix(K)[1],', n=',n,', sigma=',train_sigma), cex.main = 1.25)
# lines(x_pred, test_Y_signal, col = "black", lwd = 4)
# lines(x_pred, ypred_1 , col = "rosybrown1", lwd = 0.3)
# lines(x_pred, ypred_2, col = "orange", lwd = 2)
# lines(x_pred, ypred_3, col = "green", lwd = 2)
# lines(x_pred, ypred_4, col = "blue", lwd = 2)
# lines(x_pred, ypred_5, col = "purple", lwd = 2)
# legend('topright',c('true',paste('lam = ', lam_vec)), col = c('black',"rosybrown1",'orange','green','blue','purple'), lwd = lwd_vec)
# 
# 
# 
# train_sigma = train_mediumsigma
# test_Y_signal = f17(x_pred)
# y = f17(x) + rnorm(n, mean = 0, sd = train_sigma)
# mny = mean(y)
# ycentered = y - mny
# 
# ypred_1 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[1] * diag(1, n))) , t(K) )  , ycentered)
# ypred_2 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[2] * diag(1, n))) , t(K) )  , ycentered)
# ypred_3 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[3] * diag(1, n))) , t(K) )  , ycentered)
# ypred_4 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[4] * diag(1, n))) , t(K) )  , ycentered)
# ypred_5 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[5] * diag(1, n))) , t(K) )  , ycentered)
# # ypred_1 = eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[1] * diag(1, n))) , t(K) )  , y)
# 
# plot(x, y, xlab = "", ylab ="", cex = 0.4, main = paste("Laplacian(0.01)- prediction, rank(K) =",Matrix::rankMatrix(K)[1],', n=',n,', sigma=',train_sigma), cex.main = 1.25)
# lines(x_pred, test_Y_signal, col = "black", lwd = 4)
# lines(x_pred, ypred_1 , col = "rosybrown1", lwd = 0.3)
# lines(x_pred, ypred_2, col = "orange", lwd = 2)
# lines(x_pred, ypred_3, col = "green", lwd = 2)
# lines(x_pred, ypred_4, col = "blue", lwd = 2)
# lines(x_pred, ypred_5, col = "purple", lwd = 2)
# legend('topright',c('true',paste('lam = ', lam_vec)), col = c('black',"rosybrown1",'orange','green','blue','purple'), lwd = lwd_vec)
# 
# 
# 
# 
# train_sigma = train_mediumsigma
# test_Y_signal = f12(x_pred)
# y = f12(x) + rnorm(n, mean = 0, sd = 0.1)
# mny = mean(y)
# ycentered = y - mny
# 
# ypred_1 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[1] * diag(1, n))) , t(K) )  , ycentered)
# ypred_2 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[2] * diag(1, n))) , t(K) )  , ycentered)
# ypred_3 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[3] * diag(1, n))) , t(K) )  , ycentered)
# ypred_4 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[4] * diag(1, n))) , t(K) )  , ycentered)
# ypred_5 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[5] * diag(1, n))) , t(K) )  , ycentered)
# # ypred_1 = eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[1] * diag(1, n))) , t(K) )  , y)
# 
# plot(x, y, xlab = "", ylab ="", cex = 0.4, main = paste("Laplacian(0.01)- prediction, rank(K) =",Matrix::rankMatrix(K)[1],', n=',n,', sigma=', 0.1), cex.main = 1.25)
# lines(x_pred, test_Y_signal, col = "black", lwd = 4)
# lines(x_pred, ypred_1 , col = "rosybrown1", lwd = 0.3)
# lines(x_pred, ypred_2, col = "orange", lwd = 2)
# lines(x_pred, ypred_3, col = "green", lwd = 2)
# lines(x_pred, ypred_4, col = "blue", lwd = 2)
# lines(x_pred, ypred_5, col = "purple", lwd = 2)
# legend('topright',c('true',paste('lam = ', lam_vec)), col = c('black',"rosybrown1",'orange','green','blue','purple'), lwd = lwd_vec)
# 
# dev.off()
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# png('G:/My Drive/FunctionalDataAnalysis/Project/kernelridge/Laplacian(0.01)_eigenanalysis.png', width = 15, height = 8, units = "in", res = 150)
# 
# m1 = matrix(c(
#   1,1,2,3,4,
#   1,1,5,6,7,
#   1,1,8,9,10), nrow = 3, ncol = 5, byrow = T
# )
# layout(m1)
# 
# 
# 
# plot(1:25, seq(-0.2,1.2, length = 25), type = "n", xaxt = "n", xlab = "order", cex.axis = 1.5, cex.lab = 1.5, ylab = "eigenvalues", main = paste("Laplacian(0.1)-First 25 eigenvalues of smoother matrix, rank(K) =", Matrix::rankMatrix(K)[1]), cex.main = 1.25)
# abline(h = 1, lty = 2)
# abline(h = 0, lty = 2)
# points(1:25,evalS1[1:25], type = "o", col = "rosybrown1", pch = 16)
# points(1:25,evalS2[1:25], type = "o", col = "orange", pch = 16)
# points(1:25,evalS3[1:25], type = "o", col = "green", pch = 16)
# points(1:25,evalS4[1:25], type = "o", col = "blue", pch = 16)
# points(1:25,evalS5[1:25], type = "o", col = "purple", pch = 16)
# # legend('topright',paste('lam = ', lam_vec), col = c("rosybrown1",'orange','green','blue','purple'), lwd = 2)
# axis(1, xaxp = c(1,25,24), las = 2)
# 
# 
# 
# 
# 
# orders = 1:9
# 
# plot(x, evalS1[orders[1]] * evecS1[,orders[1]] + rep(-0.001, n), col = "rosybrown1", lwd = 0.1, type = "l", xlab = "", ylab = "", main = "1st eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
# abline(h = 0, lty = 2)
# lines(x, -evalS2[orders[1]] * evecS2[,orders[1]] + rep(-0.001, n), lwd = 1, col = "orange")
# lines(x, evalS3[orders[1]] * evecS3[,orders[1]] + rep(-0.002, n), lwd = 1, col = "green")
# lines(x, evalS4[orders[1]] * evecS4[,orders[1]] + rep(0.002, n), lwd = 1,col = "blue")
# lines(x, evalS5[orders[1]] * evecS5[,orders[1]], col = "purple")
# # legend('topright',paste('lam = ', lam_vec), col = c("rosybrown1",'orange','green','blue','purple'), lwd = 2)
# 
# plot(x, evalS1[orders[2]] * evecS1[,orders[2]] + rep(0.002, n), col = "rosybrown1", type = "l", xlab = "", ylab = "", main = "2nd eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
# abline(h = 0, lty = 2)
# lines(x, -evalS2[orders[2]] * evecS2[,orders[2]] , col = "orange")
# lines(x, -evalS3[orders[2]] * evecS3[,orders[2]] , col = "green")
# lines(x, evalS4[orders[2]] * evecS4[,orders[2]] , col = "blue")
# lines(x, evalS5[orders[2]] * evecS5[,orders[2]], col = "purple")
# # legend('topright',paste('lam = ', lam_vec), col = c("rosybrown1",'orange','green','blue','purple'), lwd = 2)
# 
# plot(x, -evalS1[orders[3]] * evecS1[,orders[3]] + rep(0.002, n), col = "rosybrown1", type = "l", xlab = "", ylab = "", main = "3rd eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
# abline(h = 0, lty = 2)
# lines(x, -evalS2[orders[3]] * evecS2[,orders[3]] + rep(-0.001, n), col = "orange")
# lines(x, -evalS3[orders[3]] * evecS3[,orders[3]] , col = "green")
# lines(x, -evalS4[orders[3]] * evecS4[,orders[3]] , col = "blue")
# lines(x, -evalS5[orders[3]] * evecS5[,orders[3]], col = "purple")
# # legend('topright',paste('lam = ', lam_vec), col = c("rosybrown1",'orange','green','blue','purple'), lwd = 2)
# 
# plot(x, -evalS1[orders[4]] * evecS1[,orders[4]] + rep(0.002, n), col = "rosybrown1", type = "l", xlab = "",ylab = "", main = "4th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
# abline(h = 0, lty = 2)
# lines(x, evalS2[orders[4]] * evecS2[,orders[4]] + rep(-0.001, n), col = "orange")
# lines(x, -evalS3[orders[4]] * evecS3[,orders[4]], col = "green")
# lines(x, -evalS4[orders[4]] * evecS4[,orders[4]], col = "blue")
# lines(x, -evalS5[orders[4]] * evecS5[,orders[4]], col = "purple")
# # legend('topright',paste('lam = ', lam_vec), col = c("rosybrown1",'orange','green','blue','purple'), lwd = 2)
# 
# plot(x, evalS1[orders[5]] * evecS1[,orders[5]] + rep(0.002, n), col = "rosybrown1", type = "l", xlab = "",ylab = "", main = "5th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
# abline(h = 0, lty = 2)
# lines(x, evalS2[orders[5]] * evecS2[,orders[5]] + rep(-0.001, n), col = "orange")
# lines(x, evalS3[orders[5]] * evecS3[,orders[5]], col = "green")
# lines(x, -evalS4[orders[5]] * evecS4[,orders[5]], col = "blue")
# lines(x,-evalS5[orders[5]] * evecS5[,orders[5]], col = "purple")
# # legend('topright',paste('lam = ', lam_vec), col = c("rosybrown1",'orange','green','blue','purple'), lwd = 2)
# 
# plot(x, -evalS1[orders[6]] * evecS1[,orders[6]] + rep(0.002, n), col = "rosybrown1", type = "l", xlab = "",ylab = "", main = "6th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
# abline(h = 0, lty = 2)
# lines(x, evalS2[orders[6]] * evecS2[,orders[6]] + rep(-0.001, n), col = "orange")
# lines(x,evalS3[orders[6]] * -evecS3[,orders[6]], col = "green")
# lines(x, -evalS4[orders[6]] * evecS4[,orders[6]], col = "blue")
# lines(x, -evalS5[orders[6]] * evecS5[,orders[6]], col = "purple")
# # legend('topright',paste('lam = ', lam_vec), col = c("rosybrown1",'orange','green','blue','purple'), lwd = 2)
# 
# plot(x, evalS1[orders[7]] * evecS1[,orders[7]] + rep(0.002, n), col = "rosybrown1", type = "l", xlab = "",ylab = "", main = "7th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
# abline(h = 0, lty = 2)
# lines(x, evalS2[orders[7]] * evecS2[,orders[7]] + rep(-0.001, n), col = "orange")
# lines(x, -evalS3[orders[7]] * evecS3[,orders[7]], col = "green")
# lines(x, evalS4[orders[7]] * evecS4[,orders[7]], col = "blue")
# lines(x, evalS5[orders[7]] * evecS5[,orders[7]], col = "purple")
# # legend('topright',paste('lam = ', lam_vec), col = c("rosybrown1",'orange','green','blue','purple'), lwd = 2)
# 
# plot(x, -evalS1[orders[8]] * evecS1[,orders[8]] + rep(0.002, n), col = "rosybrown1", type = "l", xlab = "",ylab = "", main = "8th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
# abline(h = 0, lty = 2)
# lines(x, evalS2[orders[8]] * evecS2[,orders[8]] + rep(-0.001, n), col = "orange")
# lines(x, evalS3[orders[8]] * evecS3[,orders[8]], col = "green")
# lines(x, -evalS4[orders[8]] * evecS4[,orders[8]], col = "blue")
# lines(x, -evalS5[orders[8]] * evecS5[,orders[8]], col = "purple")
# # legend('topright',paste('lam = ', lam_vec), col = c("rosybrown1",'orange','green','blue','purple'), lwd = 2)
# 
# plot(x, -evalS1[orders[9]] * evecS1[,orders[9]] + rep(0.002, n), col = "rosybrown1", type = "l", xlab = "",ylab = "", main = "9th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
# abline(h = 0, lty = 2)
# lines(x, evalS2[orders[9]] * evecS2[,orders[9]] + rep(-0.001, n), col = "orange")
# lines(x, evalS3[orders[9]] * evecS3[,orders[9]], col = "green")
# lines(x, -evalS4[orders[9]] * evecS4[,orders[9]], col = "blue")
# lines(x, evalS5[orders[9]] * evecS5[,orders[9]], col = "purple")
# # legend('topright',paste('lam = ', lam_vec), col = c("rosybrown1",'orange','green','blue','purple'), lwd = 2)
# 
# dev.off()








