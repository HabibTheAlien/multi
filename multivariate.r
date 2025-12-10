data <- read.table("E:/study/7th Semester/4104/Final Codes/skull_data_T6-13.dat")

data
colnames(data)<-c("MaxBreath","BasHeight","BasLength","NasHeight","TunePeriod")

X <- as.matrix(data[,1:4])
X
time <- data[,5]

n1 <- sum(time==1)
n2 <- sum(time==2)
n3 <- sum(time==3)
n <- n1 + n2 + n3
n # Number of observation
ens <- c(n1,n2,n3)
p<-4 # Sample vector এ variable সংখ্যা 
g <- 3 # Number of group


S1<- cov(X[time==1,])
S1
S2<- cov(X[time==2,])
S2
S3<- cov(X[time==3,])
S3

W<- (n1-1)*S1 + (n2-1)*S2 + (n3-1)*S3
S_pooled <- W/(n-g)
S_pooled

# Test equality of the Covariance Matrices or Box's M test

# M formula
M <- sum(ens-1)*log(det(S_pooled))-sum((n1 - 1)*log(det(S1)),
                                       (n2 - 1)*log(det(S2)),
                                       (n3 - 1)*log(det(S3)))
M

# U formula
U <- (sum(1 / (ens - 1)) - 1 / sum(ens - 1)) * ((2 * p^2 + 3 * p - 1) / (6 * (p + 1) * (g - 1)))
U

# C formula
C<-(1-U)*M
C

df<-(p*(p+1)*(g-1))/2
chi_val<-qchisq(0.95,df)  
chi_val

X_bar<-colMeans(X)
X_bar

X1_bar<-colMeans(X[time==1,])
X2_bar<-colMeans(X[time==2,])
X3_bar<-colMeans(X[time==3,])

B<-n1*(X1_bar-X_bar)%*%t(X1_bar-X_bar)+n2*(X2_bar-X_bar)%*%t(X2_bar-X_bar)+n3*(X3_bar-X_bar)%*%t(X3_bar-X_bar)

B

wilks<-det(W)/det(W+B)   
wilks


#test statistic
F_cal<- ((n-p-2)/p)*((1-sqrt(wilks))/sqrt(wilks))
F_cal

# for 1% level of significance 
v1<-2*p
v2<-2*(n-p-2)

F_tab<-qf(p=0.99,v1,v2)
F_tab

#Bartlett's Correction
chi_cal<- (-(n-1-((p+g)/2))*log(wilks))
chi_cal

# for 1% level of significance 
df<-p*(g-1)

chi_tab<-qchisq(0.99,df)
chi_tab







