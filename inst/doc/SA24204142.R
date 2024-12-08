## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = F, message = F)

## -----------------------------------------------------------------------------
library(microbenchmark)
library(ggplot2)
library(patchwork)
library(MASS)
library(boot)
library(bootstrap)
library(DAAG)
library(lpSolve)
library(Rcpp)
library(SA24204142)

## ----fig.align="center", fig.width = 6, fig.height=6--------------------------
ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
weight <- c(ctl, trt)
group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
lm.D9 <- lm(weight ~ group)
par(mfrow=c(2, 2))
par(mar=c(4, 4, 2.5, 2))
plot(lm.D9)

## -----------------------------------------------------------------------------
set.seed(20240910)
ir = iris[,1:4]
myclust = kmeans(x=ir, centers=3, nstart=100) # 分为3类，迭代100次
labels = myclust$cluster

## -----------------------------------------------------------------------------
targets = t(iris[5])
table(targets, labels)

## -----------------------------------------------------------------------------
set.seed(20240910, kind=NULL)
u <- runif(10000)
x <- log(2*u)*(u<=0.5) - log(2-2*u)*(u>0.5)

## ----fig.align="center", fig.width = 6, fig.height=6--------------------------
hist(x, prob=TRUE, breaks=100)
DensityFunction <- function(x)
  return(0.5*exp(-abs(x)))
curve(DensityFunction, col='red', add=TRUE)

## -----------------------------------------------------------------------------
RayleighGenerator <- function(n, sigma) {
  # n: 生成的随机数数目
  # sigma: 目标分布中的参数sigma
  # return: 一个长度为n的向量
  
  # 1. 生成均匀分布随机数U
  U <- runif(n)
  
  # 2. 输出X=F^{-1}(U)
  X <- sqrt(-2 * sigma^2 * log(1-U))
  
  return(X)
}

## -----------------------------------------------------------------------------
n <- 10000   # 生成随机数数目
sigmas <- c(0.1, 1, 10)

## -----------------------------------------------------------------------------
set.seed(20240920)
X1 <- RayleighGenerator(n, sigmas[1])
X2 <- RayleighGenerator(n, sigmas[2])
X3 <- RayleighGenerator(n, sigmas[3])

# 整成dataframe用于画图
data1 <- data.frame(X = X1)
data2 <- data.frame(X = X2)
data3 <- data.frame(X = X3)

## -----------------------------------------------------------------------------
# 目标分布的密度
TargetDensityGenerator <- function(sigma) {
  target <- function(x) {
    return(
      ( x / (sigma^2) * exp(-x^2/(2*sigma^2)) ) * ( x > 0 )
    )
  }
  return(target)
}

## -----------------------------------------------------------------------------
# 画图函数
myplot <- function(data, sigma) {
  # data: 生成的随机数，必须是dataframe，只有一列，名为"X"
  # sigma: 目标分布的参数
  # return: 一个ggplot图，带有data决定的直方图和sigma决定的目标密度曲线
  
  n <- dim(data)[1]
  Figure <- ggplot(data, aes(x=X)) +
    geom_histogram(aes(y = after_stat(density)),
                   binwidth=3.5*sigma/(n^(1/3)),
                   fill='#74c0fc', 
                   color="white") +   
    # 使用Scott法选取binwidth，经我计算Rayleigh分布的方差就是sigma^2
    stat_function(fun = TargetDensityGenerator(sigma), color = '#fd7e14', linewidth = 1) +
    geom_vline(xintercept=sigma, linetype=2, linewidth=0.5) +
    # Rayleigh分布的众数就是sigma
    labs(x="X", 
         y="Density",
         title=bquote(sigma ==. (sigma))) +   # 动态生成标题
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust=0.5),
          legend.position = "none")
  
  return(Figure)
}

## ----fig.align="center", fig.width = 6, fig.height=2--------------------------
Figure1 <- myplot(data1, sigmas[1])
Figure2 <- myplot(data2, sigmas[2])
Figure3 <- myplot(data3, sigmas[3])
Figure1 + Figure2 + Figure3

## -----------------------------------------------------------------------------
MixNormGenerator <- function(n, p) {
  # n: 生成的随机数数目
  # p: p1的取值
  # return: 一个长度为n的向量
  
  # 1.生成潜变量I
  I <- sample(0:1, size=n, replace=TRUE, prob=c(p, 1-p))
  
  # 2.从两个分量中采样
  X1 <- rnorm(n, 0, 1)
  X2 <- rnorm(n, 3, 1)
  
  # 3.输出X
  X <- X1 * (1-I) + X2 * I
  return(X)
}

## -----------------------------------------------------------------------------
# 这里写一个函数根据p1输出目标密度函数，方便后面画图的时候调用
TargetDensityGenerator <- function(p1) {
  target <- function(x) {
    return(p1*dnorm(x, 0, 1)+(1-p1)*dnorm(x, 3, 1))
  }
}

## -----------------------------------------------------------------------------
set.seed(20240921)
n <- 1000
p1 <- 0.75
X <- MixNormGenerator(n, p1)
data <- data.frame(X=X)

## -----------------------------------------------------------------------------
Figure <- ggplot(data, aes(x=X)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth=3.5*sd(data$X)/(n^(1/3)),
                 fill='#74c0fc', 
                 color="white") + 
  stat_function(fun = TargetDensityGenerator(p1), color = '#fd7e14', linewidth = 1) +
  labs(x="X", 
       y="Density",
       title=bquote(p[1] ==. (p1))) +  
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust=0.5),
        legend.position = "none")

Figure

## -----------------------------------------------------------------------------
n <- 1000
ps <- seq(0.1, 0.9, 0.1)

## -----------------------------------------------------------------------------
set.seed(20240921)

X1 <- MixNormGenerator(n, ps[1])
X2 <- MixNormGenerator(n, ps[2])
X3 <- MixNormGenerator(n, ps[3])
X4 <- MixNormGenerator(n, ps[4])
X5 <- MixNormGenerator(n, ps[5])
X6 <- MixNormGenerator(n, ps[6])
X7 <- MixNormGenerator(n, ps[7])
X8 <- MixNormGenerator(n, ps[8])
X9 <- MixNormGenerator(n, ps[9])

data1 <- data.frame(X=X1)
data2 <- data.frame(X=X2)
data3 <- data.frame(X=X3)
data4 <- data.frame(X=X4)
data5 <- data.frame(X=X5)
data6 <- data.frame(X=X6)
data7 <- data.frame(X=X7)
data8 <- data.frame(X=X8)
data9 <- data.frame(X=X9)

## -----------------------------------------------------------------------------
myplot <- function(data, p) {
  n <- nrow(data)
  Figure <- ggplot(data, aes(x=X)) +
    geom_histogram(aes(y = after_stat(density)),
                   binwidth=3.5*sd(data$X)/(n^(1/3)),
                   fill='#74c0fc',
                   color="white") + 
    labs(x="X", 
         y="Density",
         title=bquote(p[1] ==. (p))) +   # 动态生成标题
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust=0.5),
          legend.position = "none")
  return(Figure)
}

## ----fig.align="center", fig.width = 6, fig.height=6--------------------------
Figure1 <- myplot(data1, ps[1])
Figure2 <- myplot(data2, ps[2])
Figure3 <- myplot(data3, ps[3])
Figure4 <- myplot(data4, ps[4])
Figure5 <- myplot(data5, ps[5])
Figure6 <- myplot(data6, ps[6])
Figure7 <- myplot(data7, ps[7])
Figure8 <- myplot(data8, ps[8])
Figure9 <- myplot(data9, ps[9])

(Figure1+Figure2+Figure3)/(Figure4+Figure5+Figure6)/(Figure7+Figure8+Figure9)

## -----------------------------------------------------------------------------
CPGPLocalGenerator <- function(n, lambda, shape, rate, t){
  # n: 生成的随机数数目
  # lambda: Poisson过程参数
  # shape, rate: Gamma分布参数
  # t: 时间参数
  # return: 一个长度为n的向量
  
  # remark: 这个函数只能生成X(t)的样本，不能生成一条轨道
  
  # 1.生成N
  N <- rpois(n, lambda*t)
  
  # 2.输出X
  X <- rgamma(n, shape=N*shape, rate=rate)
  return(X)
}

## -----------------------------------------------------------------------------
CPGPTraceGenerator <- function(n, step, lambda, shape, rate) {
  # n: 生成的样本点数，算上t=0的点，实际上输出的是n+1个点
  # step: 步长
  # lambda: Poisson过程参数
  # shape, rate: Gamma分布参数
  # return: 一个(n+1)*2的dataframe，列名为"t"和"X"
  
  # 生成t
  t <- seq(from=0, by=step, length.out=n+1)
  
  # 生成X
  DeltaN <- rpois(n, lambda*step)
  DeltaX <- rgamma(n, shape=shape*DeltaN, rate=rate)
  X <- c(0, cumsum(DeltaX))

  # 生成dataframe
  data <- data.frame(t=t, X=X)
  
  return(data)
}

## -----------------------------------------------------------------------------
set.seed(20240921)

lambda <- 1
shape <- 1
rate <- 1
n <- 10000
ts <- c(1, 3, 10)

X1 <- CPGPLocalGenerator(n, lambda, shape, rate, ts[1])
X2 <- CPGPLocalGenerator(n, lambda, shape, rate, ts[2])
X3 <- CPGPLocalGenerator(n, lambda, shape, rate, ts[3])

data1 <- data.frame(X=X1)
data2 <- data.frame(X=X2)
data3 <- data.frame(X=X3)

## -----------------------------------------------------------------------------
TargetDensityGenerator <- function(t) {
  TargetDensity <- function(x) {
    result <- 0
    for (i in 1:20) {
      result <- result + t^i * x^{i-1} * exp(-x-t) / (factorial(i) * factorial(i-1))
    }
    return(result)
  }
  return(TargetDensity)
}

myplot <- function(data, t) {
  n <- nrow(data)
  Figure <- ggplot(data, aes(x=X)) +
    geom_histogram(aes(y = after_stat(density)), bins=50, fill='#74c0fc', color="white") + 
    stat_function(fun = TargetDensityGenerator(t), color = '#fd7e14', linewidth = 0.7) +
    labs(x="X(t)", 
         y="Density",
         title=bquote(t ==. (t))) +  
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust=0.5),
          legend.position = "none")
  return(Figure)
}

## ----fig.align="center", fig.width = 6, fig.height=2--------------------------
Figure1 <- myplot(data1, ts[1])
Figure2 <- myplot(data2, ts[2])
Figure3 <- myplot(data3, ts[3])

Figure1+Figure2+Figure3

## -----------------------------------------------------------------------------
freq <- c(sum(X1==0)/n, sum(X2==0)/n, sum(X3==0)/n)
prob <- exp(-ts)
data.frame(t=ts, freq=freq, prob=prob)

## -----------------------------------------------------------------------------
set.seed(20240921)

lambda <- 1
shape <- 1
rate <- 1
n <- 1000
step <- 0.1

data1 <- CPGPTraceGenerator(n, step, lambda, shape, rate)
data2 <- CPGPTraceGenerator(n, step, lambda, shape, rate)
data3 <- CPGPTraceGenerator(n, step, lambda, shape, rate)
data4 <- CPGPTraceGenerator(n, step, lambda, shape, rate)
data5 <- CPGPTraceGenerator(n, step, lambda, shape, rate)

data <- data.frame(t=data1$t, X1=data1$X, X2=data2$X, X3=data3$X, X4=data4$X, X5=data5$X)

## -----------------------------------------------------------------------------
Figure <- ggplot(data=data, aes(x=t)) +
  geom_line(aes(y=X1), color='#1F77B4FF') +
  geom_line(aes(y=X2), color='#FF7F0EFF') +
  geom_line(aes(y=X3), color='#2CA02CFF') +
  geom_line(aes(y=X4), color='#D62728FF') +
  geom_line(aes(y=X5), color='#9467BDFF') +
  stat_function(fun=function(x) {x}, col='black', linetype=2) +
  labs(x="t", 
       y="X",
       title='Trajectory') +   
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = "none")

Figure

## -----------------------------------------------------------------------------
# 此函数用于生成不同参数下的理论值和估计值
Estimator <- function(n, lambdas, shapes, rates, t) {
  # n: 估计时所用样本数
  # lambdas: Poisson过程参数，一个任意长度的向量
  # shapes, rates: Gamma分布参数，两个都是任意长度的向量
  # t: 时间参数，只能输入一个正数
  # return: 一个dataframe，各列分别为lambda, shape, rate, 
  #         theoretical_mean, estimated_mean, theoretical_variance, estimated_variance
  
  estimated_mean <- c()
  estimated_variance <- c()
  theoretical_mean <- c()
  theoretical_variance <- c()
  
  # 首先将lambdas, shapes, rates的取值的各种组合生成出来
  paras <- expand.grid(lambdas, shapes, rates)
  colnames(paras) <- c('lambda', 'shape', 'rate')
  
  # 生成估计值和理论值
  for (i in 1:nrow(paras)) {
    Tmean <- paras[i, ]$lambda * t * paras[i, ]$shape / paras[i, ]$rate
    Tvariance <- paras[i, ]$lambda * t * (paras[i, ]$shape + paras[i, ]$shape^2) / (paras[i, ]$rate^2)
    theoretical_mean <- c(theoretical_mean, Tmean)
    theoretical_variance <- c(theoretical_variance, Tvariance)
    
    X <- CPGPLocalGenerator(n, paras[i, ]$lambda, paras[i, ]$shape, paras[i, ]$rate, t)
    estimated_mean <- c(estimated_mean, mean(X))
    estimated_variance <- c(estimated_variance, var(X))
  }
  
  # 生成dataframe
  data <- paras
  data$theoretical_mean <- theoretical_mean
  data$estimated_mean <- estimated_mean
  data$theoretical_variance <- theoretical_variance
  data$eatimated_variance <- estimated_variance
  
  return(data)
}

## -----------------------------------------------------------------------------
set.seed(20240922)

n <- 10000
lambdas <- c(0.1, 1, 10)
shapes <- c(1, 2)
rates <- c(1, 2)
t <- 10

data <- Estimator(n, lambdas, shapes, rates, t)

print(data)

## -----------------------------------------------------------------------------
pbetaEstimator <- function(x, n) {
  # x: 目标F(x)中的x
  # n: 估计所使用的随机数数目
  # return: 输出估计
  
  # 0.排除x小于零或者大于一的情况
  if (x < 0) {
    return(0)
  }
  
  if (x > 1) {
    return(1)
  }
  
  # 1.生成U
  u <- runif(n, 0, x)
  
  # 2.输出估计
  return(30*x*mean(u^2*(1-u)^2))
}

## -----------------------------------------------------------------------------
set.seed(20240927)
n <- 1000
estimates <- c()
contrast <- c()
for (i in 1:9) {
  estimates[i] <- pbetaEstimator(0.1*i, n)
  contrast[i] <- pbeta(0.1*i, 3, 3)
}
data.frame(x=seq(0.1, 0.9, 0.1), estimate=estimates, contrast=contrast, error=estimates-contrast)

## -----------------------------------------------------------------------------
inverseF <- function(x, sigma) {
  return(sqrt(-2 * sigma^2 * log(1-x)))
}

## -----------------------------------------------------------------------------
n <- 1000 # 随机数数目的一半
sigma <- 1

## -----------------------------------------------------------------------------
set.seed(20240928)
u1 <- runif(n)
u2 <- runif(n)
v1 <- inverseF(u1, sigma)
v2 <- inverseF(u2, sigma)
variance1 <- var((v1+v2)/2)

## -----------------------------------------------------------------------------
set.seed(20240928)
u1 <- runif(n)
u2 <- 1 - u1
v1 <- inverseF(u1, sigma)
v2 <- inverseF(u2, sigma)
variance2 <- var((v1+v2)/2)

## -----------------------------------------------------------------------------
(variance1-variance2)/variance1

## -----------------------------------------------------------------------------
g <- function(x) {
  return(x^2*exp(-x^2/2)/sqrt(2*pi))
}

curve(g, xlim=c(1,6))

## ----fig.align="center", fig.width = 4, fig.height=4--------------------------
Figure <- ggplot() +
  stat_function(fun = function(x) {3*g(x)}, aes(color="3g"), linewidth = 1.3) +
  stat_function(fun = function(x) {dgamma(x-1, shape=2, rate=2)}, aes(color="f1"), linewidth = 1.3) +
  stat_function(fun = function(x) {dexp(x-1, rate=1)}, aes(color="f2"), linewidth = 1.3) +
  xlim(1, 6) + 
  labs(x="X", 
       y="y",
       title="Importance Functions VS Target",
       color="Functions") + 
  scale_color_manual(values=c("3g"="black", "f1"="#FF410D99", "f2"="#6EE2FF99")) +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5))

Figure

## -----------------------------------------------------------------------------
ISEstimator <- function(n, type) {
  # n: 估计所使用的样本数
  # type: 只能取1和2，表示取f1还是f2
  # return: 估计值

  if (type == 1) {
    x <- rgamma(n, shape=2, rate=2)
    result <- mean(g(x+1)/dgamma(x, shape=2, rate=2))
  }
  if (type == 2) {
    x <- rexp(n, 1)
    result <- mean(g(x+1)/dexp(x, 1))
  }
  
  return(result)
}

## -----------------------------------------------------------------------------
set.seed(20240928)
n <- 1000
target <- exp(-0.5)/sqrt(2*pi) + pnorm(-1) # 根据标准正态分布的对称性，其大于1的概率等于其小于-1的概率
data.frame(f1=ISEstimator(n, 1), f2=ISEstimator(n, 2), target=target)

## -----------------------------------------------------------------------------
m <- 1000
set.seed(20240928)
result1 <- c()
result2 <- c()
for (i in 1:m) {
  result1 <- c(result1, ISEstimator(n, 1))
  result2 <- c(result2, ISEstimator(n, 2))
}
variance1 <- var(result1)
variance2 <- var(result2)
data.frame(variance1=variance1, variance2=variance2)

## -----------------------------------------------------------------------------
QuickSort <- function(v) {
  # 1. 如果长度小于等于一则停止
  if (length(v) <= 1) {
    return(v)
  }
  x <- v[1]
  
  # 2. 分成两个新向量
  v1 <- v[v < x]
  v2 <- v[v > x]
  
  # 3. 对v1和v2重复上述步骤
  result <- c(QuickSort(v1), x, QuickSort(v2))
  
  return(result)
}

## -----------------------------------------------------------------------------
set.seed(20240928)
ns <- c(1e4, 2e4, 4e4, 6e4, 8e4)
an <- c()
for (n in ns) {
  ts <- c()
  for (i in 1:100) {
    v <- sample(1:n)
    
    t1 <- Sys.time()
    v <- QuickSort(v)
    t2 <- Sys.time()
    
    ts <- c(ts, t2-t1)
  }
  an <- c(an, mean(ts))
}

data <- data.frame(tn=ns*log(ns), an=an)
data 
# 我这里设置了seed但在不同的实验中a_n的值还是有波动，我怀疑是计算机自身造成的

## ----fig.align="center", fig.width = 6, fig.height=6--------------------------
Figure <- ggplot(data=data, aes(x=tn, y=an)) +
  geom_point() +
  geom_smooth(method="lm", color="#FF410D99") +
  labs(x=expression(t[n]), 
         y=expression(a[n]),
         title='Scatter Plot and Regression Line') + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust=0.5),
        legend.position = "none")

Figure

## -----------------------------------------------------------------------------
Skewness <- function(x) {
  # x: 样本
  # return: 样本对应的样本偏度系数
  
  a <- mean((x-mean(x))^3)
  b <- var(x)^1.5
  
  return(a/b)
}

## -----------------------------------------------------------------------------
MC <- function(q, n, m) {
  # q: q分位数，这里q可以为向量
  # n: \sqrt{b_1}中的n
  # m: 重复次数
  # return: 一个长度和q相同的向量，每一项为对应的q分位数
  
  bs <- replicate(m, expr = {
    x <- rnorm(n)
    Skewness(x)
  })
  
  return(quantile(bs, q))
}

## -----------------------------------------------------------------------------
set.seed(20241002)
n <- 100
m <- 10000
qs <- c(0.025, 0.05, 0.95, 0.975)

result1 <- MC(qs, n, m)
result1

## -----------------------------------------------------------------------------
SD <- function(q, n, m) {
  # q: q分位数
  # n: \sqrt{b_1}中的n
  # m: 估计分位数时一次MC实验中重复m次
  # return: 使用正态近似\sqrt{b_1}分布后的理论标准差
  
  variance <- 6 * (n-2) / ((n+1) * (n+3))
  f <- dnorm(qnorm(q, 0, sqrt(variance)), 0, sqrt(variance))
  return(q*(1-q)/(m*f^2))
}

## -----------------------------------------------------------------------------
SD(qs, n, m)

## -----------------------------------------------------------------------------
result2 <- qnorm(qs, 0, sqrt(6/n))
data.frame(Estimate=result1, LargeSample=result2)

## -----------------------------------------------------------------------------
result3 <- qnorm(qs, 0, sqrt(6*(n-2)/((n+1)*(n+3))))
data.frame(Estimate=result1, LargeSample=result2, NewLargeSample=result3)

## -----------------------------------------------------------------------------
PowerEstimator <- function(method, Sigma, alpha, n, m) {
  # method: 检验方法，只能输入"Pearson"、“Spearman”和"Kendall"
  # Sigma: 二元正态的协方差矩阵
  # alpha: 置信水平
  # n: 一次检验所使用的样本数
  # m: 一次MC实验中重复的次数
  # return: 对应方法的Power的估计值
  
  tmp <- replicate(m, expr = {
    # 1. 生成n个样本
    samples <- mvrnorm(n, c(0, 0), Sigma)
    x <- samples[, 1]
    y <- samples[, 2]
    
    # 2. 进行检验
    result <- cor.test(x, y, method=method)
    
    # 3. 是否拒绝
    result$p.value <= alpha
  })
  
  return(mean(tmp))
}

## -----------------------------------------------------------------------------
set.seed(20241002)

n <- 100
m <- 1000
alpha <- 0.01
sigmas <- c(0.2, 0.3, 0.4, 0.5)

for (sigma in sigmas) {
  
}

PPowers <- sapply(sigmas, function(x) PowerEstimator("pearson", matrix(c(1,x,x,1),2,2), alpha, n, m))

SPowers <- sapply(sigmas, function(x) PowerEstimator("spearman", matrix(c(1,x,x,1),2,2), alpha, n, m))

KPowers <- sapply(sigmas, function(x) PowerEstimator("kendall", matrix(c(1,x,x,1),2,2), alpha, n, m))

data.frame(sigma=sigmas, Pearson=PPowers, Spearman=SPowers, Kendall=KPowers)

## -----------------------------------------------------------------------------
PowerEstimator <- function(method, alpha, n, m) {
  # method: 检验方法，只能输入"Pearson"、“Spearman”和"Kendall"
  # alpha: 置信水平
  # n: 一次检验所使用的样本数
  # m: 一次MC实验中重复的次数
  # return: 对应方法的Power的估计值
  
  tmp <- replicate(m, expr = {
    # 1. 生成n个样本
    x <- rnorm(n)
    y <- as.numeric(x>0) + rnorm(n)
    
    # 2. 进行检验
    result <- cor.test(x, y, method=method)
    
    # 3. 是否拒绝
    result$p.value <= alpha
  })
  
  return(mean(tmp))
}

## -----------------------------------------------------------------------------
set.seed(20241002)

n <- 100
m <- 1000
alpha <- 0.01

PPower <- PowerEstimator("pearson", alpha, n, m)
SPower <- PowerEstimator("spearman", alpha, n, m)
KPower <- PowerEstimator("kendall", alpha, n, m)

data.frame(Pearson=PPower, Spearman=SPower, Kendall=KPower)

## -----------------------------------------------------------------------------
alpha <- 0.1
m <- 10000
N <- 1000
N0 <- 950

## -----------------------------------------------------------------------------
set.seed(20241016)

result1 <- replicate(m, expr={
  ps <- c(runif(N0), rbeta(N-N0, 0.1, 1))
  ps.Bon <- p.adjust(ps, method="bonferroni")
  ps.BH <- p.adjust(ps, method="fdr")
  
  V.Bon <- sum(ps.Bon[1:N0]<alpha)
  S.Bon <- sum(ps.Bon[(N0+1):N]<alpha)
  R.Bon <- V.Bon + S.Bon
  
  V.BH <- sum(ps.BH[1:N0]<alpha)
  S.BH <- sum(ps.BH[(N0+1):N]<alpha)
  R.BH <- V.BH + S.BH
  
  FWER.Bon <- as.numeric(V.Bon>0)
  FDR.Bon <- V.Bon / R.Bon
  TPR.Bon <- S.Bon / (N - N0)
  
  FWER.BH <- as.numeric(V.BH>0)
  FDR.BH <- V.BH / R.BH
  TPR.BH <- S.BH / (N - N0)
  
  c(FWER.Bon, FDR.Bon, TPR.Bon, FWER.BH, FDR.BH, TPR.BH)
})

result2 <- rowMeans(result1)
result3 <- data.frame(Bonferroni=result2[1:3], BH=result2[4:6])
rownames(result3) <- c("FWER", "FDR", "TPR")
result3

## -----------------------------------------------------------------------------
lambda.hat <- 1/mean(aircondit$hours)
lambda.hat

## -----------------------------------------------------------------------------
set.seed(20241016)
MLE <- function(x, i) 1/mean(x[i])
result <- boot::boot(data=boot::aircondit$hours, statistic=MLE, R=1e4)
round(c(original=result$t0, bias=mean(result$t)-result$t0, se=sd(result$t)), 3)

## -----------------------------------------------------------------------------
set.seed(20241016)
boot.mean <- function(x, i) mean(x[i])
boot.obj <- boot::boot(data=boot::aircondit$hours, statistic=boot.mean, R=1e4)
CI <- boot::boot.ci(boot.obj, conf=0.95, type=c("norm", "basic", "perc", "bca"))
CI

## -----------------------------------------------------------------------------
data <- bootstrap::scor

Sigma.hat <- cov(data)
lambdas.hat <- eigen(Sigma.hat)$values
theta.hat <- lambdas.hat[1]/sum(lambdas.hat)

thetas.jack <- c()
for (i in 1:dim(data)[1]) {
  x <- data[-i, ]
  Sigma.jack <- cov(x)
  lambdas.jack <- eigen(Sigma.jack)$values
  theta.jack <- lambdas.jack[1]/sum(lambdas.jack)
  thetas.jack <- c(thetas.jack, theta.jack)
}

bias.jack <- (dim(data)[1]-1)*(mean(thetas.jack)-theta.hat)
se.jack <- sqrt((dim(data)[1]-1)*mean((thetas.jack-theta.hat)^2))

round(c(original=theta.hat, bias.jack=bias.jack, se.jack=se.jack), 3)

## -----------------------------------------------------------------------------
data <- DAAG::ironslag
L1 <- lm(magnetic ~ chemical, data)
L2 <- lm(magnetic ~ chemical + I(chemical^2), data)
L3 <- lm(log(magnetic) ~ chemical, data)
L4 <- lm(magnetic ~ chemical + I(chemical^2) + I(chemical^3), data)

n <- length(data$magnetic)
e1 <- e2 <- e3 <- e4 <- numeric(n)

for (k in 1:n) {
  y <- data$magnetic[-k]
  x <- data$chemical[-k]
  
  J1 <- lm(y ~ x)
  yhat1 <- J1$coef[1] + J1$coef[2] * data$chemical[k]
  e1[k] <- data$magnetic[k] - yhat1
  
  J2 <- lm(y ~ x + I(x^2))
  yhat2 <- J2$coef[1] + J2$coef[2] * data$chemical[k] + J2$coef[3] * data$chemical[k] ^ 2
  e2[k] <- data$magnetic[k] - yhat2
  
  J3 <- lm(log(y) ~ x)
  logyhat3 <- J3$coef[1] + J3$coef[2] * data$chemical[k]
  yhat3 <- exp(logyhat3)
  e3[k] <- data$magnetic[k] - yhat3
  
  J4 <- lm(y ~ x + I(x^2) + I(x^3))
  yhat4 <- J4$coef[1] + J4$coef[2] * data$chemical[k] + J4$coef[3] * data$chemical[k]^2 + J4$coef[4] * data$chemical[k]^3
  e4[k] <- data$magnetic[k] - yhat4
}

c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))

## -----------------------------------------------------------------------------
summary(L1)

## -----------------------------------------------------------------------------
summary(L2)

## -----------------------------------------------------------------------------
summary(L3)

## -----------------------------------------------------------------------------
summary(L4)

## -----------------------------------------------------------------------------
x <- sort(as.vector(chickwts$weight[chickwts$feed=="soybean"]))
y <- sort(as.vector(chickwts$weight[chickwts$feed=="linseed"]))

## -----------------------------------------------------------------------------
CMtest <- function(x, y) {
  n <- length(x)
  m <- length(y)
  z <- c(x, y)
  
  Fn <- ecdf(x)
  Gm <- ecdf(y)
  
  W <- m*n / (m+n)^2 * sum((Fn(z)-Gm(z))^2)
  
  return(W)
}

## -----------------------------------------------------------------------------
set.seed(20241025)

n <- length(x)
m <- length(y)

R <- 999
z <- c(x, y)

W0 <- CMtest(x, y)
W <- numeric(R)
for (i in 1:R) {
  k <- sample(1:(n+m), size=n, replace=FALSE)
  x1 <- z[k]
  y1 <- z[-k]
  W[i] <- CMtest(x1, y1)
}
p <- mean(c(W0, W) >= W0)
p

## -----------------------------------------------------------------------------
hist(W, main="", freq=FALSE, breaks = "scott")
points(W0, 0, cex=1, pch=16)

## -----------------------------------------------------------------------------
SampleGenerator <- function(n, sigma) {
  # n: 生成的样本数
  # sigma: 正态分布的参数
  # return: x和y
  
  Sigma <- matrix(c(1, sigma, sigma, 1), 2, 2)
  samples <- mvrnorm(n, c(0,0), Sigma)
  return(samples)
}


## -----------------------------------------------------------------------------
SpearmanTest <- function(x, y, R) {
  # x, y: 数据
  # R: permutation重复次数
  # return: 近似ASL和p值
  
  test0 <- cor(x, y, method="spearman")
  tests <- replicate(R, expr = {
    x1 <- sample(x, replace=FALSE)
    cor(x1, y, method="spearman")
  })
  ASL <- mean(c(abs(tests), abs(test0)) >= abs(test0))
  p <- cor.test(x, y, method="spearman")$p.value
  
  return(c(ASL=ASL, p.value=p))
}

## -----------------------------------------------------------------------------
set.seed(20241027)
n <- 25
sigmas <- c(0.1, 0.3, 0.5)
R <- 999

ASLs <- c()
p.values <- c()
for (sigma in sigmas) {
  samples <- SampleGenerator(n, sigma)
  x <- samples[, 1]
  y <- samples[, 2]
  result <- SpearmanTest(x, y, R)
  ASLs <- c(ASLs, result[1])
  p.values <- c(p.values, result[2])
}

data.frame(sigma=sigmas, ASL=ASLs, p.value=p.values)

## -----------------------------------------------------------------------------
CauchyGenerator <- function(N, theta, eta, x0, sigma) {
  # n: 生成的样本数
  # theta, eta: Cauchy分布的参数
  # x0: 初始值
  # sigma: 提议分布的参数
  # return: 长度为n的向量
  
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N) # 用于和接受概率比较
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= dcauchy(y, eta, theta)/dcauchy(x[i-1], eta, theta)) {
      x[i] <- y
    }
    else{
      x[i] <- x[i-1]
    }
  }
  return(x)
}

## -----------------------------------------------------------------------------
set.seed(20241031)

N <- 5000
theta <- 1
eta <- 0
sigma <- 4
x0 <- 30

x <- CauchyGenerator(N, theta, eta, x0, sigma)

## -----------------------------------------------------------------------------
plot(x, type="l")

## -----------------------------------------------------------------------------
hist(x[1001:5000], breaks="scott", proba=TRUE, ylim=c(0, 0.3), main="", xlab="x")
curve(dcauchy, col="red", add=TRUE)

## -----------------------------------------------------------------------------
qs <- seq(0.1, 0.9, 0.1)
result1 <- qcauchy(qs)
result2 <- quantile(x[1001:5000], probs=qs)
data.frame(q=qs, Cauchy=result1, Observation=result2, Differece=result2-result1)

## -----------------------------------------------------------------------------
GibbsSampler <- function(N, z0, a, b, n) {
  # N: 生成样本数
  # z0: 初始值，一个二维向量
  # a, b, n: 目标分布的参数
  # return: Nx2的矩阵
  
  z <- matrix(nrow=N, ncol=2)
  z[1, ] <- z0
  x <- z0[1]
  y <- z0[2]
  for (i in 2:N) {
    x <- rbinom(1, n, y)
    y <- rbeta(1, x+a, n-x+b)
    z[i, ] <- c(x, y)
  }
  return(z)
}  

## -----------------------------------------------------------------------------
set.seed(20241101)
N <- 1000
z0 <- c(0, 0)
a <- 1
b <- 1
n <- 10
z <- GibbsSampler(N, z0, a, b, n)

## -----------------------------------------------------------------------------
plot(z[,1], type="l", ylab="x", ylim=c(-1, 11))
plot(z[,2], type="l", ylab="y", ylim=c(-1, 2))

## -----------------------------------------------------------------------------
GR <- function(psi) {
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi)
  B <- n * var(psi.means)
  psi.w <- apply(psi, 1, "var")
  W <- mean(psi.w)
  v.hat <- W * (n-1)/n + B/n
  r.hat <- v.hat / W
  return(r.hat)
}

## -----------------------------------------------------------------------------
set.seed(20241101)

N <- 8000
burn <- 1000
theta <- 1
eta <- 0
sigma <- 4
x0 <- c(-10, -5, 5, 10)

X <- matrix(nrow=length(x0), ncol=N)
for (i in 1:length(x0)) {
  X[i, ] <- CauchyGenerator(N, theta, eta, x0[i], sigma)
}

psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi)) {
  psi[i, ] <- psi[i, ] / (1:ncol(psi))
}
rhat <- rep(0, N)
for (j in (burn+1):N) {
  rhat[j] <- GR(psi[,1:j])
}

## ----fig.align="center", fig.width = 6, fig.height=6--------------------------
par(mfrow=c(2, 2))
par(mar=c(4, 4, 2, 2))
for (i in 1:length(x0)) {
  plot(psi[i, (burn+1):N], type="l", xlab=i, ylab=bquote(psi))
}

## -----------------------------------------------------------------------------
plot(rhat[(burn+1):N], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
N0 <- which(rhat[(burn+1):N]<1.2)[1]
N0

## ----fig.align="center", fig.width = 6, fig.height=6--------------------------
par(mfrow=c(2, 2))
par(mar=c(4, 4, 2, 2))
for (i in 1:length(x0)) {
  plot(X[i, 1:N0], type="l", xlab=i, ylab="x")
}

## ----fig.align="center", fig.width = 6, fig.height=6--------------------------
set.seed(20241102)

N <- 2000
burn <- 500
x0 <- c(0, 0, 0, 0)
y0 <- c(0.2, 0.4, 0.6, 0.8)
a <- 1
b <- 1
n <- 10

X <- matrix(nrow=length(x0), ncol=N)
Y <- matrix(nrow=length(x0), ncol=N)
for (i in 1:length(x0)) {
  z <- GibbsSampler(N, c(x0[i], y0[i]), a, b, n)
  X[i, ] <- as.numeric(z[, 1])
  Y[i, ] <- as.numeric(z[, 2])
}

psi <- t(apply(Y, 1, cumsum))
for (i in 1:nrow(psi)) {
  psi[i, ] <- psi[i, ] / (1:ncol(psi))
}
rhat <- rep(0, N)
for (j in (burn+1):N) {
  rhat[j] <- GR(psi[,1:j])
}

par(mfrow=c(2,2))
par(mar=c(4, 4, 2, 2))
for (i in 1:length(x0)) {
  plot(psi[i, (burn+1):N], type="l", xlab=i, ylab=bquote(psi))
}

## -----------------------------------------------------------------------------
par(mfrow=c(1,1))
plot(rhat[(burn+1):N], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
N0 <- which(rhat[(burn+1):N]<1.2)[1]
N0

## ----fig.align="center", fig.width = 6, fig.height=6--------------------------
par(mfrow=c(2,2))
par(mar=c(4, 4, 2, 2))
for (i in 1:length(x0)) {
  plot(X[i, 1:N0], type="l", xlab=i, ylab="")
  lines(Y[i, 1:N0], col="red")
}

## -----------------------------------------------------------------------------
kthTerm <- function(d, a, k) {
  return((-1)^k*exp((k+1)*log(sum(a^2))+lgamma((d+1)/2)+lgamma(k+1.5)-lgamma(k+1)-k*log(2)-lgamma(k+d/2+1)))
}

## -----------------------------------------------------------------------------
ks <- c(2, 5, 10, 100, 101, 10000, 10001)
ds <- c(2, 5, 10, 100, 101, 10000, 10001)

set.seed(20241107)
result <- c()
for (k in ks) {
  for (d in ds) {
    a <- rnorm(d) # 随机取一个a
    result <- c(result, kthTerm(d, a, k))
  }
}
result

## -----------------------------------------------------------------------------
SumComputer <- function(d, a, error) {
  result <- kthTerm(d, a, 0)
  k <- 1
  term <- kthTerm(d, a, k)
  while (term > error) {
    result <- result + term
    k <- k + 1
    term <- kthTerm(d, a, k)
  }
  return(result)
}

## -----------------------------------------------------------------------------
set.seed(20241107)
error <- 1e-5
result <- c()
for (d in ds) {
  a <- rnorm(d) # 随机取一个a
  result <- c(result, SumComputer(d, a, error))
}
result

## -----------------------------------------------------------------------------
d <- 2
a <- c(1, 2)
error <- 1e-5
SumComputer(d, a, error)

## -----------------------------------------------------------------------------
# 被积函数
f <- function(x, k) {
  return((1+x^2/k)^(-k/2))
}

# 计算c_k
ck <- function(a, k) {
  return(sqrt(a^2*k/(k+1-a^2)))
}

## -----------------------------------------------------------------------------
EquationSolver <- function(k, interval) {
  target <- function(a) {
    c1 <- ck(a, k-1)
    c2 <- ck(a, k)
    I1 <- integrate(f, lower=0, upper=c1, k=k-1)$value
    I2 <- integrate(f, lower=0, upper=c2, k=k)$value
    return(lgamma(k/2)-0.5*log(k-1)-lgamma((k-1)/2)+log(I1)-lgamma((k+1)/2)+0.5*log(k)+lgamma(k/2)-log(I2))
  }
  solution <- uniroot(target, interval)
  return(solution$root)
}

## -----------------------------------------------------------------------------
ks <- c(4:25, 100, 500, 1000)
result1 <- c()
for (k in ks) {
  result1 <- c(result1, EquationSolver(k, c(0.01, sqrt(k)-0.01)))
}
data.frame(k=ks, solution=result1)

## -----------------------------------------------------------------------------
result2 <- c()
for (k in ks) {
  solution <- uniroot(f=function(a) {
    pt(ck(a, k-1), k-1, lower.tail=FALSE)-pt(ck(a, k), k, lower.tail=FALSE)
    }, 
    c(0.01, sqrt(k)-0.01))
  result2 <- c(result2, solution$root)
}
data.frame(k=ks, sqrt.k=sqrt(ks), Exercise11.4=result2, Exercise11.5=result1)

## -----------------------------------------------------------------------------
Ys <- c(0.54, 0.48, 0.33, 0.43, 1.00, 1.00, 0.91, 1.00, 0.21, 0.85)

## -----------------------------------------------------------------------------
EM <- function(Ys, N, lambda0) {
  lambda <- lambda0
  m <- sum(Ys<1)
  n <- length(Ys)
  for (i in 1:N) {
    lambda <- mean(Ys) + (1 - m/n)*lambda
  }
  return(lambda)
}

## -----------------------------------------------------------------------------
lambda1 <- EM(Ys, 100, 0.1)
lambda1

## -----------------------------------------------------------------------------
lambda2 <- sum(Ys)/sum(Ys<1)
lambda2

## -----------------------------------------------------------------------------
f.obj <- c(4, 2, 9)
f.con <- matrix(c(2, 1, 1, -1, 1, 3), nrow=2, ncol=3)
f.dir <- rep("<=", 2)
f.rhs <- c(2, 3)
res.lp <- lp("min", f.obj, f.con, f.dir, f.rhs)
res.lp$solution

## ----eval=FALSE---------------------------------------------------------------
#  formulas <- list(
#    mpg ~ disp,
#    mpg ~ I(1 / disp),
#    mpg ~ disp + wt,
#    mpg ~ I(1 / disp) + wt
#  )

## -----------------------------------------------------------------------------
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)

## -----------------------------------------------------------------------------
lms21 <- list()
for (i in 1:4) {
  lms21[[i]] <- lm(formulas[[i]], mtcars)
}
lms21

## -----------------------------------------------------------------------------
lms22 <- lapply(formulas, lm, data=mtcars)
lms22

## ----eval=FALSE---------------------------------------------------------------
#  bootstraps <- lapply(1:10, function(i) {
#    rows <- sample(1:nrow(mtcars), rep = TRUE)
#    mtcars[rows, ]
#  })

## -----------------------------------------------------------------------------
set.seed(20241114)
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})

## -----------------------------------------------------------------------------
lms31 <- list()
for (i in 1:10) {
  lms31[[i]] <- lm(mpg ~ disp, bootstraps[[i]])
}
lms31

## -----------------------------------------------------------------------------
lms32 <- lapply(bootstraps, lm, formula=mpg ~ disp)
lms32

## -----------------------------------------------------------------------------
set.seed(20241116)
lms33 <- list()
for (i in 1:10) {
  rows <- sample(1:nrow(mtcars), rep=TRUE)
  lms33[[i]] <- lm(formula=mpg~disp ,data=mtcars[rows,])
}
lms33

## ----eval=FALSE---------------------------------------------------------------
#  rsq <- function(mod) summary(mod)$r.squared

## -----------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared

## -----------------------------------------------------------------------------
rs2 <- lapply(lms22, rsq) # Question2的模型
rs3 <- lapply(lms32, rsq) # Question3的模型
as.numeric(rs2) 
as.numeric(rs3)

## ----eval=FALSE---------------------------------------------------------------
#  trials <- replicate(
#    100,
#    t.test(rpois(10, 10), rpois(7, 10)),
#    simplify = FALSE
#  )

## -----------------------------------------------------------------------------
set.seed(20241116)
trials <- replicate(
  100,
  t.test(rpois(10, 10), rpois(7, 10)), 
  simplify=FALSE
)

sapply(trials, function(trial) trial$p.value)

## -----------------------------------------------------------------------------
sapply(trials, '[[', "p.value")

## -----------------------------------------------------------------------------
Newlappy <- function(fun, ..., FUN.VALUE) {
  tmp <- Map(fun, ...)
  result <- vapply(tmp, FUN=identity, FUN.VALUE=FUN.VALUE)
  return(result)
}

## -----------------------------------------------------------------------------
a <- list(c(1, 2, 3, 4, 5), c(1, 1, 1))
b <- list(c(5, 4, 3, 2, 1), c(2, 2, 2))
Newlappy(function(x,y) mean(x+y), a, b, FUN.VALUE=numeric(1))

## -----------------------------------------------------------------------------
NewChisqTest <- function(x, y) {
  data <- table(x, y)  
  N <- length(x)
  row_nums <- rowSums(data)  
  col_nums <- colSums(data)  
  
  Es <- matrix(0, nrow = nrow(data), ncol = ncol(data))  
  for (i in 1:nrow(data)) {  
    for (j in 1:ncol(data)) {  
      Es[i, j] <- (row_nums[i] * col_nums[j]) / N  
    }  
  }  
  
  data <- as.numeric(data)
  Es <- as.numeric(Es)
  result <- sum((data - Es)^2 / Es)  

  return(result)  
}

## -----------------------------------------------------------------------------
set.seed(20241116)
x <- sample(1:5, 50, replace=TRUE)
y <- sample(1:2, 50, replace=TRUE)

test1 <- as.numeric(chisq.test(table(x, y))$statistic)
test2 <- NewChisqTest(x, y)

time1 <- replicate(
  1000,
  expr= {
    start1 <- Sys.time()
    test1 <- as.numeric(chisq.test(table(x, y))$statistic)
    end1 <- Sys.time()
    end1 - start1
  })

time2 <- replicate(
  1000,
  expr= {
    start2 <- Sys.time()
    test2 <- NewChisqTest(x, y)
    end2 <- Sys.time()
    end2 - start2
  })

data.frame(fun=c("chisq.test", "New"), value=c(test1, test2), time=c(mean(time1), mean(time2)))

## -----------------------------------------------------------------------------
NewTable <- function(x, y) {
  x.new <- unique(x)
  y.new <- unique(y)
  result <- matrix(0, nrow=length(x.new), ncol=length(y.new))
  rownames(result) <- x.new
  colnames(result) <- y.new
  
  for (i in 1:length(x)) {
    row_index <- match(x[i], x.new)
    col_index <- match(y[i], y.new)
    result[row_index, col_index] <- result[row_index, col_index]+1
  }
  
  return(result)
}

## -----------------------------------------------------------------------------
set.seed(20241116)
x <- sample(1:5, 50, replace=TRUE)
y <- sample(1:2, 50, replace=TRUE)

result1 <- table(x, y)
result2 <- NewTable(x, y)

time1 <- replicate(
  1000,
  expr= {
    start1 <- Sys.time()
    test1 <- table(x, y)
    end1 <- Sys.time()
    end1 - start1
  })

time2 <- replicate(
  1000,
  expr= {
    start2 <- Sys.time()
    test2 <- NewTable(x, y)
    end2 <- Sys.time()
    end2 - start2
  })

## -----------------------------------------------------------------------------
result1
mean(time1)

## -----------------------------------------------------------------------------
result2
mean(time2)

## -----------------------------------------------------------------------------
NewNewChisqTest <- function(x, y) {
  data <- NewTable(x, y)  
  N <- length(x)
  row_nums <- rowSums(data)  
  col_nums <- colSums(data)  
  
  Es <- matrix(0, nrow = nrow(data), ncol = ncol(data))  
  for (i in 1:nrow(data)) {  
    for (j in 1:ncol(data)) {  
      Es[i, j] <- (row_nums[i] * col_nums[j]) / N  
    }  
  }  
  
  data <- as.numeric(data)
  Es <- as.numeric(Es)
  result <- sum((data - Es)^2 / Es)  

  return(result)  
}

## -----------------------------------------------------------------------------
set.seed(20241116)
x <- sample(1:5, 50, replace=TRUE)
y <- sample(1:2, 50, replace=TRUE)

test1 <- as.numeric(chisq.test(table(x, y))$statistic)
test2 <- NewChisqTest(x, y)
test3 <- NewNewChisqTest(x, y)


time1 <- replicate(
  1000,
  expr= {
    start1 <- Sys.time()
    test1 <- as.numeric(chisq.test(table(x, y))$statistic)
    end1 <- Sys.time()
    end1 - start1
  })

time2 <- replicate(
  1000,
  expr= {
    start2 <- Sys.time()
    test2 <- NewChisqTest(x, y)
    end2 <- Sys.time()
    end2 - start2
  })

time3 <- replicate(
  1000,
  expr= {
    start3 <- Sys.time()
    test3 <- NewNewChisqTest(x, y)
    end3 <- Sys.time()
    end3 - start3
  })

data.frame(fun=c("chisq.test", "New", "NewNew"), value=c(test1, test2, test3), time=c(mean(time1), mean(time2), mean(time3)))

## -----------------------------------------------------------------------------
GibbsR <- function(N, a, b, n) {
  # N: 生成样本数
  # a, b, n: 目标分布的参数
  # return: Nx2的矩阵
  
  z <- matrix(nrow=N, ncol=2)
  z[1, ] <- c(0, 0)
  x <- 0
  y <- 0.2
  for (i in 2:N) {
    x <- rbinom(1, n, y)
    y <- rbeta(1, x+a, n-x+b)
    z[i, ] <- c(x, y)
  }
  return(z)
}  

## -----------------------------------------------------------------------------
set.seed(20241123)
N <- 1000
burn <- 500
a <- 1
b <- 1
n <- 10
result1 <- GibbsR(N, a, b, n)
result2 <- GibbsC(N, a, b, n)

## -----------------------------------------------------------------------------
qqplot(result1[(burn+1):N,1], result2[(burn+1):N,1], xlab="GibbsR", ylab="GibbsC")
abline(0, 1, col="red")
qqplot(result1[(burn+1):N,2], result2[(burn+1):N,2], xlab="GibbsR", ylab="GibbsC")
abline(0, 1, col="red")

## -----------------------------------------------------------------------------
ts <- microbenchmark(gibbR=GibbsR(N, a, b, n),
                     gibbC=GibbsC(N, a, b, n))
summary(ts)[, c(1, 3, 5, 6)]

