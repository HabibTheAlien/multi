#-------------------------------------------------------------------------------------
#----------------------------- Hotelling's T Square ---------------------------
#-------------------------------------------------------------------------------------


Y1<-c(15,17,15,13,20,15,15,13,14,17,15,17,15,18,18,15,18,10,18,18,13,16,11,16,16,18)
Y2<-c(24,32,29,10,26,26,26,22,30,30,26,28,29,32,31,26,33,19,30,34,30,16,25,26,23,34)
Y3<-c(14,26,23,16,28,21,22,22,17,27,20,24,24,28,27,21,26,17,29,26,24,16,23,16,21,24)

Y<- data.frame(Y1,Y2,Y3)

mu0<-c(14,25,20)

n<-nrow(Y)
p<-ncol(Y)

Ybar<-colMeans(Y)

S<-var(Y)

T2<-n*t(Ybar- mu0)%*%solve(S)%*%(Ybar-mu0)
T2


Fval<-qf(p=0.95,df1=p,df2=(n-p))
Fval

T2_tab<-((n-1)*p/(n-p))*Fval
T2_tab

#Simultaneous Confidence Interval
wd<-sqrt((n-1)*p/(n-p)*Fval)*sqrt(diag(S)/n)
wd

SCI<-cbind(ybar-wd,ybar+wd)
SCI




#-------------------------------------------------------------------------------------
#----------------------------- Multivariate Regression ---------------------------
#-------------------------------------------------------------------------------------

#=============================================================================================
#======================= Problem:-8 =======================================================
#=============================================================================================
Z0<-c(1,1,1,1,1)
Z1<-c(0,1,2,3,4)
Y1<-c(1,4,3,8,9)
Y2<-c(-1,-1,2,3,2)

Y<-matrix(c(Y1,Y2),ncol=2)
Z<-matrix(c(Z0,Z1),ncol=2)

tZ.Z<-t(Z)%*%Z
inv.tZ.Z<-solve(tZ.Z)

Beta1<-inv.tZ.Z%*%t(Z)%*%Y1
Beta1

Beta2<-inv.tZ.Z%*%t(Z)%*%Y2
Beta2

Beta<-cbind(Beta1,Beta2)
Beta

Yhat<-round(Z%*%Beta)
Yhat

ehat<- Y-Yhat
ehat


# Properties Check

ssr<-t(Yhat)%*%Yhat
ssr

sse<-t(ehat)%*%ehat
sse

sst<-t(Y)%*%Y
sst
ssr+sse







#=============================================================================================
#======================= Example:-4   =======================================================
#=============================================================================================

X0<-rep(1,8)
X1<-c(23,23,30,30,25,25,30,30)
X2<-c(40,42,58,45,55,42,48,51)
Y1<-c(41.5,33.8,27.7,21.7,19.9,15,12.2,19.3)
Y2<-c(45.9,53.3,57.5,58.8,60.6,58,60.6,58)
Y3<-c(11.2,11.2,12.7,16,16.2,22.6,24.5,21.3)

Y<-matrix(c(Y1,Y2,Y3),ncol=3)
X<-matrix(c(X0,X1,X2),ncol=3)

Beta1<-solve((t(X)%*%X))%*%t(X)%*%Y1
Beta1

Beta2<-solve((t(X)%*%X))%*%t(X)%*%Y2
Beta2

Beta3<-solve((t(X)%*%X))%*%t(X)%*%Y3
Beta3

Beta<-cbind(Beta1,Beta2,Beta3)
Beta

# Verify
Beta<-solve((t(X)%*%X))%*%t(X)%*%Y
Beta


Yhat<-round(X%*%Beta)
Yhat

ehat<- Y-Yhat
ehat



# Sum of Squares
ssr<-t(Yhat)%*%Yhat
ssr

sse<-t(ehat)%*%ehat
sse

sst<-t(Y)%*%Y
sst
ssr+sse


# Properties check
t(Yhat)%*%ehat
t(X)%*%ehat

# Solve by Regression model
m<-lm(Y~cbind(X1,X2))
summary(m)
coef(m)
residuals(m)







#-------------------------------------------------------------------------------------
#----------------------------- MANOVA ---------------------------
#-------------------------------------------------------------------------------------







X1bar<-matrix(c(2.066,0.480,0.082, 0.360),ncol=1)
X1bar

X2bar<-matrix(c(2.167,0.596,0.124,0.418),ncol=1)
X2bar

X3bar<-matrix(c(2.273,0.521,0.125,0.383),ncol=1)
X3bar

S1<-matrix(c(
  0.291,-0.001, 0.002, 0.010,
  -0.001, 0.011, 0.000, 0.003,
  0.002, 0.000, 0.001, 0.000,
  0.010, 0.003, 0.000, 0.010
),nrow=4,ncol=4,byrow=TRUE)
S1

S2<-matrix(c(
  0.561, 0.011, 0.001, 0.037,
  0.011, 0.025, 0.004, 0.007,
  0.001, 0.004, 0.005, 0.002,
  0.037, 0.007, 0.002, 0.019
),nrow=4,ncol=4,byrow=TRUE)
S2

S3<-matrix(c(
  0.261, 0.030, 0.003, 0.018,
  0.030, 0.017,-0.000, 0.006,
  0.003,-0.000, 0.004, 0.001,
  0.018, 0.006, 0.001, 0.013
),nrow=4,ncol=4,byrow=TRUE)
S3

n1<-271
n2<-138
n3<-107


n<-n1+n2+n3
n

ens<-c(n1,n2,n3)
g<-3
p<-4

W<-(n1-1)*S1+(n2-1)*S2+(n3-1)*S3
W

Xbar<-(n1*X1bar+n2*X2bar+n3*X3bar)/n
Xbar

B<-n1*(X1bar-Xbar)%*%t(X1bar-Xbar)+n2*(X2bar-Xbar)%*%t(X2bar-Xbar)+n3*(X3bar-Xbar)%*%t(X3bar-Xbar)
B

lambda= det(W)/det(B+W)
lambda

#the test statistics
test<- ((n-p-2)/p)*((1-sqrt(lambda))/sqrt(lambda))
test


qf(0.99,df1=2*p,df2=2*(n-p-2))


tau1<-X1bar-Xbar
tau2<-X2bar-Xbar
tau3<-X3bar-Xbar


tau13_tau33= tau1[3]-tau3[3]
tau13_tau33

alpha<-0.05/(p*g*(g-1))

tval<-qt(alpha,df=(n-g))
res<-(1/n1)+(1/n3)
wd<-tval*sqrt((W[3,3]/(n-g))*res)

SCI<-c(tau13_tau33-wd, tau13_tau33+wd)
SCI





#-------------------------------------------------------------------------------------
#----------------------------- Skull ---------------------------
#-------------------------------------------------------------------------------------



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
W
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



























