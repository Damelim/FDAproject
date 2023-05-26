source('G:/My Drive/FunctionalDataAnalysis/Project/kernelridge/kernelridge.R')

load('G:/My Drive/research/Pspline/Psplinedata.rdata')
n = 800  
set.seed(1) ; x = sort(runif(n, 0, 1)) ; x_pred = seq(0.001, 0.999, by = 0.002) ## function evaluation points

train_sigma = train_mediumsigma
test_Y_signal = f8(x_pred)
y = f8(x) + rnorm(n, mean = 0, sd = train_sigma)
mny = mean(y)
ycentered = y - mny


#### gamma = 10 ####

K = generate_kernel_matrix(x, kernel_type = "sigmoid", params = 10)
K_pred = generate_prediction_kernel_matrix(trainx = x, testx = x_pred, kernel_type = "sigmoid", params = 10)

lam_vec = c(0.01, 0.1, 1, 10)

Smoothermat1 = eigenMapMatMult(eigenMapMatMult(K, armaInv(eigenMapMatMult(t(K), K) + lam_vec[1] * diag(1, n))) , t(K) )
Smoothermat2 = eigenMapMatMult(eigenMapMatMult(K, armaInv(eigenMapMatMult(t(K), K) + lam_vec[2] * diag(1, n))) , t(K) )
Smoothermat3 = eigenMapMatMult(eigenMapMatMult(K, armaInv(eigenMapMatMult(t(K), K) + lam_vec[3] * diag(1, n))) , t(K) )
Smoothermat4 = eigenMapMatMult(eigenMapMatMult(K, armaInv(eigenMapMatMult(t(K), K) + lam_vec[4] * diag(1, n))) , t(K) )

eigS1 = eigen(Smoothermat1);eigS2 = eigen(Smoothermat2);eigS3 = eigen(Smoothermat3)  ;eigS4 = eigen(Smoothermat4)

evalS1 = as.numeric(eigS1$values);evalS2 = as.numeric(eigS2$values);evalS3 = as.numeric(eigS3$values)  ;evalS4 = as.numeric(eigS4$values)

evecS1 = eigS1$vectors ; evecS2 = eigS2$vectors ; evecS3 = eigS3$vectors ; evecS4 = eigS4$vectors  
for(i in 1:n){
  for(j in 1:n){
    evecS1[i,j] = as.numeric(evecS1[i,j])
    evecS2[i,j] = as.numeric(evecS2[i,j])
    evecS3[i,j] = as.numeric(evecS3[i,j])
    evecS4[i,j] = as.numeric(evecS4[i,j])
  }
}
ypred_1 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[1] * diag(1, n))) , t(K) )  , ycentered)
ypred_2 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[2] * diag(1, n))) , t(K) )  , ycentered)
ypred_3 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[3] * diag(1, n))) , t(K) )  , ycentered)
ypred_4 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[4] * diag(1, n))) , t(K) )  , ycentered)
# ypred_1 = eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[1] * diag(1, n))) , t(K) )  , y)
# ypred_2 = eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[2] * diag(1, n))) , t(K) )  , y)
# ypred_3 = eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[3] * diag(1, n))) , t(K) )  , y)
# ypred_4 = eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[4] * diag(1, n))) , t(K) )  , y)



par(mfrow = c(1,2))
plot(x, y, xlab = "", ylab ="", cex = 0.4, main = paste("Pspline-out of sample prediction, rank(K) =",Matrix::rankMatrix(K)[1],', n=',n,', sigma=',train_sigma), cex.main = 1)
lines(x_pred, test_Y_signal, col = "black", lwd = 4)
lines(x_pred, ypred_1 , col = "red", lwd = 2)
lines(x_pred, ypred_2, col = "green", lwd = 2)
lines(x_pred, ypred_3, col = "blue", lwd = 2)
lines(x_pred, ypred_4, col = "purple", lwd = 2)
legend('topright',paste('lam = ', lam_vec), col = c('red','green','blue','purple'), lwd = 2)

plot(1:25, seq(-0.2,1.2, length = 25), type = "n", xaxt = "n", xlab = "order", ylab = "eigenvalues", main = paste("Pspline-First 25 eigenvalues of smoother matrix, rank(K) =", Matrix::rankMatrix(K)[1]), cex.main = 1)
abline(h = 1, lty = 2)
abline(h = 0, lty = 2)
points(1:25,evalS1[1:25], type = "o", col = "red", pch = 16)
points(1:25,evalS2[1:25], type = "o", col = "green", pch = 16)
points(1:25,evalS3[1:25], type = "o", col = "blue", pch = 16)
points(1:25,evalS4[1:25], type = "o", col = "purple", pch = 16)
legend('topright',paste('lam = ', lam_vec), col = c('red','green','blue','purple'), lwd = 2)
axis(1, xaxp = c(1,25,24), las = 2)



par(mfrow = c(3,3))

orders = 1:9

plot(x, -evalS1[orders[1]] * evecS1[,orders[1]] + rep(-0.001, n), col = "red", type = "l", xlab = "", ylab = "", main = "1st eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, -evalS2[orders[1]] * evecS2[,orders[1]] + rep(-0.002, n), col = "green")
lines(x, -evalS3[orders[1]] * evecS3[,orders[1]] + rep(0.002, n), col = "blue")
lines(x, -evalS4[orders[1]] * evecS4[,orders[1]], col = "purple")
legend('topright',paste('lam = ', lam_vec), col = c('red','green','blue','purple'), lwd = 2)

plot(x, evalS1[orders[2]] * evecS1[,orders[2]] + rep(0.002, n), col = "red", type = "l", xlab = "", ylab = "", main = "2nd eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, -evalS2[orders[2]] * evecS2[,orders[2]] , col = "green")
lines(x, -evalS3[orders[2]] * evecS3[,orders[2]] , col = "blue")
lines(x, evalS4[orders[2]] * evecS4[,orders[2]], col = "purple")
legend('topright',paste('lam = ', lam_vec), col = c('red','green','blue','purple'), lwd = 2)

plot(x, evalS1[orders[3]] * evecS1[,orders[3]] + rep(0.002, n), col = "red", type = "l", xlab = "", ylab = "", main = "3rd eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, evalS2[orders[3]] * evecS2[,orders[3]] , col = "green")
lines(x, -evalS3[orders[3]] * evecS3[,orders[3]] , col = "blue")
lines(x, -evalS4[orders[3]] * evecS4[,orders[3]], col = "purple")
legend('topright',paste('lam = ', lam_vec), col = c('red','green','blue','purple'), lwd = 2)

plot(x, evalS1[orders[4]] * evecS1[,orders[4]] + rep(0.002, n), col = "red", type = "l", xlab = "",ylab = "", main = "4th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, -evalS2[orders[4]] * evecS2[,orders[4]], col = "green")
lines(x, -evalS3[orders[4]] * evecS3[,orders[4]], col = "blue")
lines(x, evalS4[orders[4]] * evecS4[,orders[4]], col = "purple")
legend('topright',paste('lam = ', lam_vec), col = c('red','green','blue','purple'), lwd = 2)

plot(x, evalS1[orders[5]] * evecS1[,orders[5]] + rep(0.002, n), col = "red", type = "l", xlab = "",ylab = "", main = "5th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, evalS2[orders[5]] * evecS2[,orders[5]], col = "green")
lines(x, -evalS3[orders[5]] * evecS3[,orders[5]], col = "blue")
lines(x, evalS4[orders[5]] * evecS4[,orders[5]], col = "purple")
legend('topright',paste('lam = ', lam_vec), col = c('red','green','blue','purple'), lwd = 2)

plot(x, -evalS1[orders[6]] * evecS1[,orders[6]] + rep(0.002, n), col = "red", type = "l", xlab = "",ylab = "", main = "6th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x,evalS2[orders[6]] * -evecS2[,orders[6]], col = "green")
lines(x, evalS3[orders[6]] * evecS3[,orders[6]], col = "blue")
lines(x, -evalS4[orders[6]] * evecS4[,orders[6]], col = "purple")
legend('topright',paste('lam = ', lam_vec), col = c('red','green','blue','purple'), lwd = 2)

plot(x, evalS1[orders[7]] * evecS1[,orders[7]] + rep(0.002, n), col = "red", type = "l", xlab = "",ylab = "", main = "7th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, evalS2[orders[7]] * evecS2[,orders[7]], col = "green")
lines(x, -evalS3[orders[7]] * evecS3[,orders[7]], col = "blue")
lines(x, evalS4[orders[7]] * evecS4[,orders[7]], col = "purple")
legend('topright',paste('lam = ', lam_vec), col = c('red','green','blue','purple'), lwd = 2)

plot(x, -evalS1[orders[8]] * evecS1[,orders[8]] + rep(0.002, n), col = "red", type = "l", xlab = "",ylab = "", main = "8th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, evalS2[orders[8]] * evecS2[,orders[8]], col = "green")
lines(x, -evalS3[orders[8]] * evecS3[,orders[8]], col = "blue")
lines(x, -evalS4[orders[8]] * evecS4[,orders[8]], col = "purple")
legend('topright',paste('lam = ', lam_vec), col = c('red','green','blue','purple'), lwd = 2)

plot(x, -evalS1[orders[9]] * evecS1[,orders[9]] + rep(0.002, n), col = "red", type = "l", xlab = "",ylab = "", main = "9th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, -evalS2[orders[9]] * evecS2[,orders[9]], col = "green")
lines(x, evalS3[orders[9]] * evecS3[,orders[9]], col = "blue")
lines(x, evalS4[orders[9]] * evecS4[,orders[9]], col = "purple")
legend('topright',paste('lam = ', lam_vec), col = c('red','green','blue','purple'), lwd = 2)













































#### gamma = 100 ####

K = generate_kernel_matrix(x, kernel_type = "sigmoid", params = 100)
K_pred = generate_prediction_kernel_matrix(trainx = x, testx = x_pred, kernel_type = "sigmoid", params = 100)

lam_vec = c(0.01, 0.1, 1, 10)

Smoothermat1 = eigenMapMatMult(eigenMapMatMult(K, armaInv(eigenMapMatMult(t(K), K) + lam_vec[1] * diag(1, n))) , t(K) )
Smoothermat2 = eigenMapMatMult(eigenMapMatMult(K, armaInv(eigenMapMatMult(t(K), K) + lam_vec[2] * diag(1, n))) , t(K) )
Smoothermat3 = eigenMapMatMult(eigenMapMatMult(K, armaInv(eigenMapMatMult(t(K), K) + lam_vec[3] * diag(1, n))) , t(K) )
Smoothermat4 = eigenMapMatMult(eigenMapMatMult(K, armaInv(eigenMapMatMult(t(K), K) + lam_vec[4] * diag(1, n))) , t(K) )

eigS1 = eigen(Smoothermat1);eigS2 = eigen(Smoothermat2);eigS3 = eigen(Smoothermat3)  ;eigS4 = eigen(Smoothermat4)

evalS1 = as.numeric(eigS1$values);evalS2 = as.numeric(eigS2$values);evalS3 = as.numeric(eigS3$values)  ;evalS4 = as.numeric(eigS4$values)

evecS1 = eigS1$vectors ; evecS2 = eigS2$vectors ; evecS3 = eigS3$vectors ; evecS4 = eigS4$vectors  
for(i in 1:n){
  for(j in 1:n){
    evecS1[i,j] = as.numeric(evecS1[i,j])
    evecS2[i,j] = as.numeric(evecS2[i,j])
    evecS3[i,j] = as.numeric(evecS3[i,j])
    evecS4[i,j] = as.numeric(evecS4[i,j])
  }
}
ypred_1 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[1] * diag(1, n))) , t(K) )  , ycentered)
ypred_2 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[2] * diag(1, n))) , t(K) )  , ycentered)
ypred_3 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[3] * diag(1, n))) , t(K) )  , ycentered)
ypred_4 = mny + eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[4] * diag(1, n))) , t(K) )  , ycentered)
# ypred_1 = eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[1] * diag(1, n))) , t(K) )  , y)
# ypred_2 = eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[2] * diag(1, n))) , t(K) )  , y)
# ypred_3 = eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[3] * diag(1, n))) , t(K) )  , y)
# ypred_4 = eigenMapMatMult(  eigenMapMatMult(eigenMapMatMult(K_pred, armaInv(eigenMapMatMult(t(K), K) + lam_vec[4] * diag(1, n))) , t(K) )  , y)



par(mfrow = c(1,2))
plot(x, y, xlab = "", ylab ="", cex = 0.4, main = paste("Sigmoid(2): prediction, rank(K) =",Matrix::rankMatrix(K)[1],', n=',n,', sigma=',train_sigma), cex.main = 1)
lines(x_pred, test_Y_signal, col = "black", lwd = 4)
lines(x_pred, ypred_1 , col = "red", lwd = 2)
lines(x_pred, ypred_2, col = "green", lwd = 2)
lines(x_pred, ypred_3, col = "blue", lwd = 2)
lines(x_pred, ypred_4, col = "purple", lwd = 2)
legend('topright',paste('lam = ', lam_vec), col = c('red','green','blue','purple'), lwd = 2)

plot(1:25, seq(-0.2,1.2, length = 25), type = "n", xaxt = "n", xlab = "order", ylab = "eigenvalues", main = paste("Pspline-First 25 eigenvalues of smoother matrix, rank(K) =", Matrix::rankMatrix(K)[1]), cex.main = 1)
abline(h = 1, lty = 2)
abline(h = 0, lty = 2)
points(1:25,evalS1[1:25], type = "o", col = "red", pch = 16)
points(1:25,evalS2[1:25], type = "o", col = "green", pch = 16)
points(1:25,evalS3[1:25], type = "o", col = "blue", pch = 16)
points(1:25,evalS4[1:25], type = "o", col = "purple", pch = 16)
legend('topright',paste('lam = ', lam_vec), col = c('red','green','blue','purple'), lwd = 2)
axis(1, xaxp = c(1,25,24), las = 2)



par(mfrow = c(3,3))

orders = 1:9

plot(x, -evalS1[orders[1]] * evecS1[,orders[1]] + rep(-0.001, n), col = "red", type = "l", xlab = "", ylab = "", main = "1st eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, -evalS2[orders[1]] * evecS2[,orders[1]] + rep(-0.002, n), col = "green")
lines(x, -evalS3[orders[1]] * evecS3[,orders[1]] + rep(0.002, n), col = "blue")
lines(x, -evalS4[orders[1]] * evecS4[,orders[1]], col = "purple")
legend('topright',paste('lam = ', lam_vec), col = c('red','green','blue','purple'), lwd = 2)

plot(x, evalS1[orders[2]] * evecS1[,orders[2]] + rep(0.002, n), col = "red", type = "l", xlab = "", ylab = "", main = "2nd eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, -evalS2[orders[2]] * evecS2[,orders[2]] , col = "green")
lines(x, -evalS3[orders[2]] * evecS3[,orders[2]] , col = "blue")
lines(x, evalS4[orders[2]] * evecS4[,orders[2]], col = "purple")
legend('topright',paste('lam = ', lam_vec), col = c('red','green','blue','purple'), lwd = 2)

plot(x, evalS1[orders[3]] * evecS1[,orders[3]] + rep(0.002, n), col = "red", type = "l", xlab = "", ylab = "", main = "3rd eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, evalS2[orders[3]] * evecS2[,orders[3]] , col = "green")
lines(x, -evalS3[orders[3]] * evecS3[,orders[3]] , col = "blue")
lines(x, -evalS4[orders[3]] * evecS4[,orders[3]], col = "purple")
legend('topright',paste('lam = ', lam_vec), col = c('red','green','blue','purple'), lwd = 2)

plot(x, evalS1[orders[4]] * evecS1[,orders[4]] + rep(0.002, n), col = "red", type = "l", xlab = "",ylab = "", main = "4th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, -evalS2[orders[4]] * evecS2[,orders[4]], col = "green")
lines(x, -evalS3[orders[4]] * evecS3[,orders[4]], col = "blue")
lines(x, evalS4[orders[4]] * evecS4[,orders[4]], col = "purple")
legend('topright',paste('lam = ', lam_vec), col = c('red','green','blue','purple'), lwd = 2)

plot(x, evalS1[orders[5]] * evecS1[,orders[5]] + rep(0.002, n), col = "red", type = "l", xlab = "",ylab = "", main = "5th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, evalS2[orders[5]] * evecS2[,orders[5]], col = "green")
lines(x, -evalS3[orders[5]] * evecS3[,orders[5]], col = "blue")
lines(x, evalS4[orders[5]] * evecS4[,orders[5]], col = "purple")
legend('topright',paste('lam = ', lam_vec), col = c('red','green','blue','purple'), lwd = 2)

plot(x, -evalS1[orders[6]] * evecS1[,orders[6]] + rep(0.002, n), col = "red", type = "l", xlab = "",ylab = "", main = "6th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x,evalS2[orders[6]] * -evecS2[,orders[6]], col = "green")
lines(x, evalS3[orders[6]] * evecS3[,orders[6]], col = "blue")
lines(x, -evalS4[orders[6]] * evecS4[,orders[6]], col = "purple")
legend('topright',paste('lam = ', lam_vec), col = c('red','green','blue','purple'), lwd = 2)

plot(x, evalS1[orders[7]] * evecS1[,orders[7]] + rep(0.002, n), col = "red", type = "l", xlab = "",ylab = "", main = "7th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, evalS2[orders[7]] * evecS2[,orders[7]], col = "green")
lines(x, -evalS3[orders[7]] * evecS3[,orders[7]], col = "blue")
lines(x, evalS4[orders[7]] * evecS4[,orders[7]], col = "purple")
legend('topright',paste('lam = ', lam_vec), col = c('red','green','blue','purple'), lwd = 2)

plot(x, -evalS1[orders[8]] * evecS1[,orders[8]] + rep(0.002, n), col = "red", type = "l", xlab = "",ylab = "", main = "8th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, evalS2[orders[8]] * evecS2[,orders[8]], col = "green")
lines(x, -evalS3[orders[8]] * evecS3[,orders[8]], col = "blue")
lines(x, -evalS4[orders[8]] * evecS4[,orders[8]], col = "purple")
legend('topright',paste('lam = ', lam_vec), col = c('red','green','blue','purple'), lwd = 2)

plot(x, -evalS1[orders[9]] * evecS1[,orders[9]] + rep(0.002, n), col = "red", type = "l", xlab = "",ylab = "", main = "9th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, -evalS2[orders[9]] * evecS2[,orders[9]], col = "green")
lines(x, evalS3[orders[9]] * evecS3[,orders[9]], col = "blue")
lines(x, evalS4[orders[9]] * evecS4[,orders[9]], col = "purple")
legend('topright',paste('lam = ', lam_vec), col = c('red','green','blue','purple'), lwd = 2)


