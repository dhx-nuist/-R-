# 导入数据
data_x <- read.csv("data_x.csv",header = TRUE,sep = ",")
data_y <- read.csv("data_y.csv",header = TRUE,sep = ",")
data_xy <- cbind(data_x,data_y)
p <- 3


# 计算x、y、xy的协方差矩阵和相关系数矩阵（1）
cov_x <- cov(data_x)
cov_y <- cov(data_y)
cov_xy <- cov(data_xy)
cov_xy <- cov_xy[1:3,4:6]
cor_xy <- cor(data_xy)
cor_xy <- cor_xy[1:3,4:6]

# 计算A矩阵（2）
##计算x矩阵的特征值和特征向量
eigen_x <- eigen(cov_x)
eigen_x_value <- eigen_x$values
eigen_x_vector <- eigen_x$vectors
## 将特征值全部开根号求倒数
eigen_x_value_sqrt <- 1/sqrt(eigen_x_value)
## 计算x矩阵的-1/2次方
nsqrt_x <- matrix(0,ncol = 3,nrow = 3)
for (i in 1:3){
  nsqrt_x <- nsqrt_x+eigen_x_value_sqrt[i]*eigen_x_vector[,i]%*%t(eigen_x_vector[,i])
}
## 计算Y的逆矩阵
n_y <- solve(cov_y)
## 计算A矩阵
A <- nsqrt_x%*%cov_xy%*%n_y%*%t(cov_xy)%*%nsqrt_x
## 计算A矩阵的特征值和特征向量
eigen_A <- eigen(A)
eigen_A_value <- eigen_A$values
eigen_A_vector <- eigen_A$vectors

# 计算B矩阵（3）
## 计算y矩阵的特征值和特征向量
eigen_y <- eigen(cov_y)
eigen_y_value <- eigen_y$values
eigen_y_vector <- eigen_y$vectors
## 将特征值全部开根号求倒数
eigen_y_value_sqrt <- 1/sqrt(eigen_y_value)
## 计算y矩阵的-1/2次方
nsqrt_y <- matrix(0,ncol = 3,nrow = 3)
for (i in 1:3){
  nsqrt_y <- nsqrt_y + eigen_y_value_sqrt[i]*eigen_y_vector[,i]%*%t(eigen_y_vector[,i])
}
## 计算X的逆矩阵
n_x <- solve(cov_x)
## 计算B矩阵
B <- nsqrt_y%*%t(cov_xy)%*%n_x%*%cov_xy%*%nsqrt_y
## 计算B矩阵的特征值和特征向量
eigen_B <- eigen(B)
eigen_B_value <- eigen_B$values
eigen_B_vector <- eigen_B$vectors

# 计算第一组典型变量与原始变量的相关系数（4）
## 计算第一组典型变量
V1 <- t(eigen_A_vector[,1]%*%nsqrt_x%*%t(data_x))
V2 <- t(eigen_A_vector[,2]%*%nsqrt_x%*%t(data_x))
V3 <- t(eigen_A_vector[,3]%*%nsqrt_x%*%t(data_x))
## 计算第一组典型变量与原始变量的相关系数矩阵
V_data <- cbind(V1,V2,V3)
G <- cor(V_data,data_x)
G_1 <- t(G[1,])%*%G[1,]/p
G_2 <- t(G[2,])%*%G[2,]/p
G_3 <- t(G[3,])%*%G[3,]/p
Gs <- G_1+G_2+G_3

# 计算第二组典型变量与原始变量的相关系数（5）
## 计算第二组典型变量
W1 <- t(eigen_B_vector[,1]%*%nsqrt_y%*%t(data_y))
W2 <- t(eigen_B_vector[,2]%*%nsqrt_y%*%t(data_y))
W3 <- t(eigen_B_vector[,3]%*%nsqrt_y%*%t(data_y))
## 计算第二组典型变量与原始变量的相关系数矩阵
W_data <- cbind(W1,W2,W3)
H <- cor(W_data,data_y)
H_1 <- t(H[1,])%*%H[1,]/p
H_2 <- t(H[2,])%*%H[2,]/p
H_3 <- t(H[3,])%*%H[3,]/p
Hs <- H_1+H_2+H_3

# 计算第一组典型变量与第二组变量的相关系数（6）
D <- cor(V_data,data_y)
D_1 <- t(D[1,])%*%D[1,]/p
D_2 <- t(D[2,])%*%D[2,]/p
D_3 <- t(D[3,])%*%D[3,]/p
Ds <- D_1+D_2+D_3

# 计算第二组典型变量与第一组变量的相关系数（7）
P <- cor(W_data,data_x)
P_1 <- t(P[1,])%*%P[1,]/p
P_2 <- t(P[2,])%*%P[2,]/p
P_3 <- t(P[3,])%*%P[3,]/p
Ps <- P_1+P_2+P_3

# 计算冗余度（8）
Rv1 <- G_1*sqrt(eigen_A_value[1])/p
Rv2 <- G_2*sqrt(eigen_A_value[2])/p
Rv3 <- G_3*sqrt(eigen_A_value[3])/p
RW1 <- H_1*sqrt(eigen_B_value[1])/p
RW2 <- H_2*sqrt(eigen_B_value[2])/p
RW3 <- H_3*sqrt(eigen_B_value[3])/p
