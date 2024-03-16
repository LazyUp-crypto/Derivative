[TOC]



# Longstaff-Schwartz Monte Carlo 

## 前言

近年来，我国衍生品市场，尤其是场外衍生品市场规模迅速扩张，衍生品定价能力愈发成为市场参与者竞争的关键所在，同时，伴随着资管新规下刚性兑付的打破，权益风险居高不下、固收类产品收益率低，投资者对绝对收益产品的偏好越来越强烈，在这样的背景下，衍生品作为风险对冲工具的重要性愈发凸显。因而本文选择关注市场上流动性较强、占比较大的类美式期权具体的定价方法，其中Longstaff&Schwartz[^1]在《 Valuing American Options by Simulation: A Simple Least-Squares Approach》(2001)中提出的基于回归的蒙特拉洛方法是本文主要讨论的对象（后文简称LSMC），使用一方法对美式看跌期权进行定价，并基于美式期权与欧式期权、百慕大期权价格的关系对定价准确性进行验证。同时，结合《American Options and the lsm Algorithm: Quasi-Random Sequences and Brownian Bridges》(2005)[^2]这一文章，分别基于布朗桥和准随机序列对LSMC算法进行优化。在上述基础上，探讨LSMC算法在维度诅咒、基函数选取以及收敛速度上的局限性及对应解决方案。

### 期权类型

市场上的主流期权类型为欧式期权、美式期权以及奇异期权，欧式期权是指只能在到期日选择是否行权，美式期权在到期日之前任意时刻可以行权，奇异期权又包括时间依赖、路径依赖、多资产期权等多种类型，美式期权也具有时间、路径依赖性，下表为比较具有代表性的奇异期权类型。

| 特性     | 期权类型                                               |
| -------- | ------------------------------------------------------ |
| 时间依赖 | 百慕大期权                                             |
| 路径依赖 | 亚式、障碍、回溯、阶梯、呐喊                           |
| 多资产   | 彩虹期权、双币期权Quanto option、篮子期权Basket option |

### 期权定价方法

经典的欧式期权定价一般采用BSM公式，蒙特卡洛模拟，有限差分法求解PDE三种方式，其中BSM公式提供了欧式期权的封闭解，因此在学界和业界都被广泛应用。但随着中国商品期权上市步伐的加快，中国的期权品种中大多数都是美式期权，且多种衍生产品都具有类美式期权的性质，即具有路径依赖性，由于美式期权在到期日之前任何一个交易日都可以行权，其定价复杂度远大于欧式期权，这类期权并不存在数值解且具有路径依赖性，必须求助于数值计算方式。期权定价的数值方法主要包括两大类：一类是确定性方法，比如有限差分和树方法；一类是随机方法，也即Monte Carlo方法。目前业界使用最为广泛的仍然是Monte Carlo方法，因其对于路径依赖的各类奇异期权，多资产期权等具有较强的通用性。理论上来说，随机模拟方法效率和精度很低，但是Monte Carlo算法中模拟路径部分相互独立，因而可以通过并行计算、大幅度提高随即次数，达到所需精度。另外，基于其速度问题，国内外多位学者提出了相应的方差缩减技术，如对偶变量法、控制变量法等。

### LSMC算法介绍

美式期权可以在任何时间行权，因而其定价需要考虑最优停时问题（`Optimal Stopping`)，核心是比较 t时刻行权价值(strike value)和存续价值（continuation value)。传统Monte Carlo方法在欧式期权、路径依赖期权中都得到了广泛应用，但其一般次序是按照时间演化次序，知道了现在某个时间点上的S才可以去求解下一个时间点 $\ t + \Delta_{t}\\$上的价格，无法在每一个节点比较行权价值和存续价值，因而无法对美式期权定价。基于此，`Longstaff&Schwartz`提出了一种最小二乘蒙特卡洛（简称LSMC）算法，也叫做Regression-Based Monte Carlo，其执行策略如下：

- 正向模拟标的路径，假设模拟m条路径，总时间为N
  $$
  (X^{m}_n)^{N}_{n=0}
  $$

- 在到期时刻求得payoff
  $$
  V^{m}_N = g(N,x^{m}_N)
  $$
  
- 向前返回一个时间节点，设置基函数如下形式，常见的基函数形式为多项式函数，指数函数以及对数函数等
  $$
  q_n(x) = \Sigma^{K}_{k=0}\beta_{k}x^{k}
  $$

- 利用未来的折现价值以及最小二乘回归法确定参数$\ \beta_{k}$
  $$
  \Sigma^{M}_{m=1}(e^{-rt_n} V^{m}_{n+1}-q_n(x^{m}_{n}))^2
  $$

- 根据所求基函数估算当前节点的存续价值，比较当前节点的存续价值与立即行权价值，取二者较大的一点作为当前节点价格，也即
  $$
  V^{m}_{n} = 
  \begin{equation}
  \label{eq6}
  [x_{i}]=\left\{
  \begin{aligned}
  g(n,x^{m}_{n}) & , & q_{n}(x^{m}_{n})<g(n,x^{m}_{n})
  \\e^{-rn\frac{T}{N}} V^{m}_{n+1}& , & \mu_{a}(x_{i})< \mu_{b}(x_{i}).
  \end{aligned}
  \right.
  \end{equation}
  $$
  
- 重复上述步骤直至达到初始节点，得到期权价值。

这一方法可以应用于美式期权，百慕大期权等带有最优停时问题的期权定价。

## LSMC算法应用

基于R语言实现美式看跌期权定价以及相应希腊字母计算，并基于欧式期权、百慕大期权与美式期权价格之间的关系检验结果准确性。

### 美式看跌期权定价

```R
# Longstaff-Schwartz算法对美式看跌期权进行定价
LSMprice <- function(S0, K, r, T, sigma, steps, n, degree) {
    # steps 是指时间步数，n是指模拟次数，degree是指多项式的阶数
    dt <- T / steps
    t <- seq(0, T, by = dt)
    t <- matrix(rep(t, each = n), nrow = n)
    W <- matrix(rnorm(n * (steps + 1)), nrow = n, ncol = steps + 1) # 随机过程，生成股票价格路径的随机部分
    S <- matrix(0, nrow = n, ncol = steps + 1)
    S[, 1] <- S0
    for (i in 2:(steps + 1)) {
        S[, i] <- S[, i - 1] * exp((r - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * W[, i - 1])
    }

    payoff <- pmax(K - S, 0)
    V <- payoff[, steps + 1]
    discount_factor <- exp(-r * dt)
    for (i in (steps:1)) {
        V <- discount_factor * V
        itm <- (payoff[, i] > 0)
        X <- S[itm, i]
        Y <- V[itm]
        fit <- lm(Y ~ poly(X, degree, raw = TRUE)) 
        continuation_value <- predict(fit, newdata = data.frame(X = X))
        exercise <- payoff[itm, i]
        V[itm] <- ifelse(exercise > continuation_value, exercise, V[itm])
    }
    return(mean(V) * discount_factor)
}
S0 <- 36
K <- 40
r <- 0.06
T <- 1 
sigma <- 0.2
steps <- 50
n <- 10000
degree <- 5
raw_put_price <- LSMprice(S0, K, r, T, sigma, steps, n, degree)
print(paste("美式看跌期权价格：", raw_put_price))
#输出结果为
[1] "使用Longstaff-Schwartz算法得美式看跌期权的价格为： 4.48577263169445"
```

### 希腊字母计算

在期权中一般使用希腊字母度量期权价值的敏感性，计算期权的风险敞口，实际对冲中，期权希腊字母甚至可能比期权定价更为重要。

常用的希腊字母如下：
$$
\Delta = \frac{\alpha{V}}{\alpha{S}};
\gamma = \frac{\alpha{\Delta}}{\alpha{S}};
\theta =  \frac{\alpha{V}}{\alpha{t}};
\rho = \frac{\alpha{V}}{\alpha{r}};
vega =  \frac{\alpha{V}}{\alpha{\sigma}};
$$
欧式期权希腊字母的计算可以通过对封闭数值解求导得到，其他类型期权的希腊字母的计算则一般基于有限差分法，有限差分法又包括前向差分、中间差分和后向差分法，本文选择采用有限差分法中的中间差分法计算美式看跌期权的希腊字母，具体代码如下：

```R
# 使用有限差分-中间差分计算美式看跌期权Greek值
delta <- function() {
    diff <- S0 * 0.01
    myCall_1 <- LSMprice(S0 + diff, K, r, T, sigma, steps, n, degree)
    myCall_2 <- LSMprice(S0 - diff, K, r, T, sigma, steps, n, degree)
    return((myCall_1 - myCall_2) / (2 * diff))
}

gamma <- function() {
    diff <- S0 * 0.01
    myCall_1 <- LSMprice(S0 + diff, K, r, T, sigma, steps, n, degree)
    myCall_2 <- LSMprice(S0 - diff, K, r, T, sigma, steps, n, degree)
    myCall_0 <- LSMprice(S0, K, r, T, sigma, steps, n, degree)
    return((myCall_1 - 2 * myCall_0 + myCall_2) / (diff^2))
}

theta <- function() {
    diff <- T * 0.01
    myCall_1 <- LSMprice(S0, K, r, T + diff, sigma, steps, n, degree)
    myCall_2 <- LSMprice(S0, K, r, T - diff, sigma, steps, n, degree)
    return((myCall_2 - myCall_1) / (2 * diff))
}

rho <- function() {
    diff <- r * 0.01
    myCall_1 <- LSMprice(S0, K, r + diff, T, sigma, steps, n, degree)
    myCall_2 <- LSMprice(S0, K, r - diff, T, sigma, steps, n, degree)
    return((myCall_1 - myCall_2) / (2 * diff))
}

vega <- function() {
    diff <- sigma * 0.01
    myCall_1 <- LSMprice(S0, K, r, T, sigma + diff, steps, n, degree)
    myCall_2 <- LSMprice(S0, K, r, T, sigma - diff, steps, n, degree)
    return((myCall_1 - myCall_2) / (2 * diff))}

print(paste("delta值为：", delta()))
print(paste("gamma值为：", gamma()))
print(paste("theta值为：", theta()))
print(paste("rho值为：", rho()))
print(paste("vega值为：", vega()))

#输出结果为
[1] "delta值为： -0.700176124258597"
[1] "gamma值为： 0.0603123888370873"
[1] "theta值为： -2.57863625983163"
[1] "rho值为： -16.2026307501795"
[1] "vega值为： 18.2459511275874"
```

### 检验

#### 欧式VS美式

相同参数下，美式期权价值一般高于欧式期权价值，因此利用BS模型、基于相同参数计算欧式期权价值，将二者进行比较，验证定价准确性。

```R
BSMPrice <- function(S0,K,T,r,sigma,Type){
  d1 <- (log(S0/K)+(r+sigma^2/2)*T)/(sigma*sqrt(T))
  d2 <- d1-sigma*sqrt(T)
  if(Type=="Call"){
    return(S0*pnorm(d1)-K*exp(-r*T)*pnorm(d2))
  }else{
    return(K*exp(-r*T)*pnorm(-d2)-S0*pnorm(-d1))
  }
}
BSMPrice(36,40,0.5,0.06,0.2,"Put")
#输出结果为3.809587，小于美式期权价值
```

#### 百慕大VS美式

百慕大期权的特点是在可以在到期日之前特定的日期行权。

将百慕大期权的行权时间设置为任意时间步，则相同参数下百慕大期权价格和美式期权价格接近。

```R
LSMprice_Bermuda <- function(S0, K, r, T, sigma, steps, n, degree, exercise_dates) {
    # steps 是指时间步数，n是指模拟次数，degree是指多项式的阶数
    dt <- T / steps
    t <- seq(0, T, by = dt)
    t <- matrix(rep(t, each = n), nrow = n)
    W <- matrix(rnorm(n * (steps + 1)), nrow = n, ncol = steps + 1) # 随机过程，生成股票价格路径的随机部分
    S <- matrix(0, nrow = n, ncol = steps + 1)
    S[, 1] <- S0
    for (i in 2:(steps + 1)) {
        S[, i] <- S[, i - 1] * exp((r - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * W[, i - 1])
    }

    payoff <- pmax(K - S, 0)
    V <- payoff[, steps + 1]
    discount_factor <- exp(-r * dt)
    for (i in (steps:1)) {
        V <- discount_factor * V
        # 只在特定的日期进行行权
        if (i %in% exercise_dates) {
            itm <- (payoff[, i] > 0)
            X <- S[itm, i]
            Y <- V[itm]
            fit <- lm(Y ~ poly(X, degree, raw = TRUE)) 
            continuation_value <- predict(fit, newdata = data.frame(X = X))
            exercise <- payoff[itm, i]
            V[itm] <- ifelse(exercise > continuation_value, exercise, V[itm])
        }
    }
    return(mean(V) * discount_factor)
}
#exercise参数设置为1：steps
S0 <- 36
K <- 40
r <- 0.06
T <- 1 
sigma <- 0.2
steps <- 50
n <- 10000
degree <- 5
exercise_dates <- 1:steps
# 计算期权价格
option_price <- LSMprice_Bermuda(S0, K, r, T, sigma, steps, n, degree, exercise_dates)
print(paste("百慕大期权价格: ", option_price))
#输出结果为
[1] "百慕大期权价格:  4.47529796207148"
```

可以看到，相同参数、百慕大期权可行权时间设置为任意步的情况下，百慕大期权价格和美式期权价格非常接近。

## LSMC算法优化

接下来本文将以美式看跌期权定价为例，基于准随机序列和布朗桥优化LSMC算法。

### 基于准随机序列优化

准随机序列（Quasi-Random Sequence）是一种特殊的数列，它在整个定义域内的分布更为均匀，相比于纯随机数列，准随机序列能更好地覆盖整个定义域。这种特性使得基于准随机序列的蒙特卡洛模拟（Quasi-Monte Carlo Simulation）在许多情况下比基于纯随机数列的蒙特卡洛模拟（Pure Monte Carlo Simulation）更为精确。在LSMC（Least Squares Monte Carlo）算法中需要模拟大量的股票价格路径，然后基于这些路径来估计期权的价格。如果使用准随机序列来生成股票价格路径，那么这些路径将更好地覆盖所有可能的股票价格，从而使得估计更为精确。此外，由于准随机序列的分布更为均匀，基于准随机序列的蒙特卡洛模拟通常能够更快地收敛，这意味着可能需要更少的模拟次数就能得到一个相当精确的估计，从而提高了算法的效率。

本文选择使用Sobol序列生成股票价格路径的随机部分，以此优化算法。

```R
library(randtoolbox)
# Longstaff-Schwartz算法对美式看跌期权进行定价

LSMprice_sobol <- function(S0, K, r, T, sigma, steps, n, degree) {
    dt <- T / steps
    t <- seq(0, T, by = dt)
    t <- matrix(rep(t, each = n), nrow = n)

    # 生成Sobol序列
    sobol_seq <- sobol(n = n, d = steps + 1)
    # 使用Sobol序列来生成股票价格路径的随机部分
    W <- qnorm(sobol_seq / max(sobol_seq))

    S <- matrix(0, nrow = n, ncol = steps + 1)
    S[, 1] <- S0
    for (i in 2:(steps + 1)) {
        S[, i] <- S[, i - 1] * exp((r - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * W[, i - 1])
    }

    payoff <- pmax(K - S, 0)
    V <- payoff[, steps + 1]
    discount_factor <- exp(-r * dt)
    for (i in (steps:1)) {
        V <- discount_factor * V
        itm <- (payoff[, i] > 0)
        X <- S[itm, i]
        Y <- V[itm]
        fit <- lm(Y ~ poly(X, degree, raw = TRUE)) 
        continuation_value <- predict(fit, newdata = data.frame(X = X))
        exercise <- payoff[itm, i]
        V[itm] <- ifelse(exercise > continuation_value, exercise, V[itm])
    }
    return(mean(V) * discount_factor)
}

S0 <- 36
K <- 40
r <- 0.06
T <- 1 
sigma <- 0.2
steps <- 50
n <- 10000
degree <- 5

# 使用纯随机序列的LSMC算法
random_seq_price <- LSMprice(S0, K, r, T, sigma, steps, n, degree)

# 使用准随机序列的LSMC算法
sobol_seq_price <- LSMprice_sobol(S0, K, r, T, sigma, steps, n, degree)

# 比较两种方法的期权价格估计值
print(paste("使用纯随机序列的期权价格估计值：", random_seq_price))
print(paste("使用准随机序列的期权价格估计值：", sobol_seq_price))

# 比较两种方法的估计值的方差
random_seq_var <- var(replicate(100, LSMprice(S0, K, r, T, sigma, steps, n, degree)))
sobol_seq_var <- var(replicate(100, LSMprice_sobol(S0, K, r, T, sigma, steps, n, degree)))
print(paste("使用纯随机序列的估计值的方差：", random_seq_var))
print(paste("使用准随机序列的估计值的方差：", sobol_seq_var))

#输出结果如下：
[1] "使用纯随机序列的期权价格估计值： 4.48224409141348"
[1] "使用准随机序列的期权价格估计值： 4.52991200998241"
[1] "使用纯随机序列的估计值的方差： 0.000838052806644685"
[1] "使用准随机序列的估计值的方差： 0.000103834288182832"
```

### 基于布朗桥优化

布朗桥（Brownian Bridge）是一种特殊的随机过程，它描述了一个布朗运动在两个时间点的值已知的情况下，中间时间点的值的分布。在蒙特卡洛模拟中，布朗桥可以用来更精确地模拟股票价格的路径。在LSMC（Least Squares Monte Carlo）算法中需要模拟大量的股票价格路径，然后基于这些路径来估计期权的价格。如果直接模拟每个时间点的股票价格，那么这些价格之间的相关性可能会被忽视，从而影响估计的精度。而如果使用布朗桥来模拟股票价格的路径，那么就可以在模拟的过程中考虑到股票价格之间的相关性。具体来说，首先模拟期初和期末的股票价格，然后使用布朗桥来填充中间的股票价格。这样，模拟出的股票价格路径将更接近真实的股票价格路径，从而使我们的估计更为精确。此外，由于布朗桥只需要模拟两个时间点的股票价格，然后使用确定性的方法来填充中间的股票价格，因此使用布朗桥的模拟通常比直接模拟每个时间点的股票价格更为高效。因此，使用布朗桥可以优化LSMC算法，提高其精度和效率。

```r
LSMprice_BrownianBridge <- function(S0, K, r, T, sigma, steps, n, degree) {
    dt <- T / steps
    t <- seq(0, T, by = dt)
    t <- matrix(rep(t, each = n), nrow = n)
    W <- matrix(rnorm(n * (steps + 1)), nrow = n, ncol = steps + 1) # 随机过程，生成股票价格路径的随机部分
    S <- matrix(0, nrow = n, ncol = steps + 1)
    S[, 1] <- S0
    for (i in 2:(steps + 1)) {
        # 使用布朗桥来生成中间的布朗运动值
        W[, i] <- W[, i - 1] + sqrt(dt) * rnorm(n)
        S[, i] <- S[, i - 1] * exp((r - 0.5 * sigma^2) * dt + sigma * (W[, i] - W[, i - 1]))
    }

    payoff <- pmax(K - S, 0)
    V <- payoff[, steps + 1]
    discount_factor <- exp(-r * dt)
    for (i in (steps:1)) {
        V <- discount_factor * V
        itm <- (payoff[, i] > 0)
        X <- S[itm, i]
        Y <- V[itm]
        fit <- lm(Y ~ poly(X, degree, raw = TRUE))
        continuation_value <- predict(fit, newdata = data.frame(X = X))
        exercise <- payoff[itm, i]
        V[itm] <- ifelse(exercise > continuation_value, exercise, V[itm])
    }
    return(mean(V) * discount_factor)
}
S0 <- 36
K <- 40
r <- 0.06
T <- 1 
sigma <- 0.2
steps <- 50
n <- 10000
degree <- 5

# 使用未优化的LSMC算法
unoptimized_price <- LSMprice(S0, K, r, T, sigma, steps, n, degree)

# 使用基于布朗桥优化的LSMC算法
optimized_price <- LSMprice_BrownianBridge(S0, K, r, T, sigma, steps, n, degree)

# 比较两种方法的期权价格估计值
print(paste("使用未优化的LSMC算法的期权价格估计值：", unoptimized_price))
print(paste("使用基于布朗桥优化的LSMC算法的期权价格估计值：", optimized_price))

# 比较两种方法的估计值的方差
unoptimized_var <- var(replicate(100, LSMprice(S0, K, r, T, sigma, steps, n, degree)))
optimized_var <- var(replicate(100, LSMprice_BrownianBridge(S0, K, r, T, sigma, steps, n, degree)))
print(paste("使用未优化的LSMC算法的估计值的方差：", unoptimized_var))
print(paste("使用基于布朗桥优化的LSMC算法的估计值的方差：", optimized_var))

#输出结果
[1] "使用未优化的LSMC算法的期权价格估计值： 4.49541959426898"
[1] "使用基于布朗桥优化的LSMC算法的期权价格估计值： 4.49491435319167"
[1] "使用未优化的LSMC算法的估计值的方差： 0.000918542300193074"
[1] "使用基于布朗桥优化的LSMC算法的估计值的方差： 0.000772058506211493"
```

可以看到两种优化算法下，美式看跌期权价格估计值的方差都缩小，除了上述两种优化方法之外，最小二乘蒙特卡洛（Least Squares Monte Carlo，LSMC）算法还可以使用其他方法进行路径优化：

1. 长步长控制变量（Longstaff-Schwartz）：LSMC算法中的一个关键步骤是回归分析，用于估计衍生品的价格函数。长步长控制变量方法通过在模拟路径中引入控制变量，将衍生品价格分解为基于已知信息的部分和残差部分。这种方法可以减小残差项的方差，提高价格估计的准确性。
2. 提前停止规则（Early Exercise Rules）：某些衍生品具有可提前行权的特性。在LSMC算法中，可以使用提前停止规则来确定在模拟路径上是否应该提前行权。通过合理选择提前停止规则，可以减少模拟路径数量，从而提高计算效率。
3. 变量选择技术（Variable Selection Techniques）：在LSMC算法中，回归分析的准确性和计算效率受到回归变量的选择影响。变量选择技术可以帮助确定对于估计价格函数最相关的变量，从而减少回归分析的复杂度和噪声的影响。
4. 平滑技术（Smoothing Techniques）：在模拟路径中，可能存在噪声或不连续性。平滑技术可以通过对路径进行插值或平滑处理来减少噪声的影响，使路径更加连续和可靠。

## 拓展

在面对类美式期权定价问题上，神经网络和强化学习的引入有利于解决“维度诅咒”问题，且模型泛化能力更强。

### 维度灾难

美式期权定价的核心在于最优停时问题，LSMC算法是动态规划的逆向算法，作者也提到在更高维的美式期权定价问题上这一算法受限，存在“维度灾难”问题。Pricing of High-Dimensional American Options by Neural Networks (2010)[^3]第一次提出将神经网络应用于LSMC算法，对更高维的美式期权进行定价，Deep Optimal Stopping（2019）[^4]沿用LSMC的动态规划思想，但是通过神经网络对行权决策函数进行拟合，从而估计衍生品价格，同LSMC相比，这一方法可以在复杂的收益计算方法和多资产的情况下进行期权定价，适用范围更广泛。

Pricing and hedging American-style options with deep learning（2019）[^5]同样使用神经网络解决美式期权定价和对冲问题，与Deep Optimal Stopping（2019） 使用神经网络训练停时不同，此研究和LSMC类似，将LSMC使用回归方法估计期权存续价值变成了借助神经网络拟合存续价值函数。另有一些学者将强化学习引入类美式期权定价，Regression Methods for Pricing Complex American-Style Options中提到使用fitted Q-Iteration[^6]方法，Learning Exercise Policies for American Options[^7]使用least-squares policy iteration方法进行美式期权定价，Optimal Stopping via Randomized Neural Networks（2021）[^8]提出了基于强化学习与RNN的期权定价算法，并且分别在Markov与非Markov情况下的测试中展现出SOTA的结果。

### 基函数选取

LSMC算法中基函数类型的选取、不同基函数的阶数选取都是需要权衡的问题，阶数太低可能无法充分拟合数据，阶数太高可能过拟合，需要选取合适的阶数来平衡模型的拟合能力和泛化能力。对数函数的增长速度叫较慢，一般需要比较高的阶数，指数函数增长速度较快，一般选取比较低的阶数。本文选择使用5折交叉验证方法验证LSMC算法中使用多项式、指数、对数基函数评估模型性能，将阶数固定为3阶。

```R
# 选取不同的基函数对美式看跌期权进行定价
LSMprice_poly <- function(S0, K, r, T, sigma, steps, n, degree) {
    # steps 是指时间步数，n是指模拟次数
    dt <- T / steps
    t <- seq(0, T, by = dt)
    t <- matrix(rep(t, each = n), nrow = n)
    W <- matrix(rnorm(n * (steps + 1)), nrow = n, ncol = steps + 1) # 随机过程，生成股票价格路径的随机部分
    S <- matrix(0, nrow = n, ncol = steps + 1)
    S[, 1] <- S0
    for (i in 2:(steps + 1)) {
        S[, i] <- S[, i - 1] * exp((r - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * W[, i - 1])
    }
    payoff <- pmax(K - S, 0)
    V <- payoff[, steps + 1]
    discount_factor <- exp(-r * dt)
    for (i in (steps:1)) {
        V <- discount_factor * V
        itm <- (payoff[, i] > 0)
        X <- S[itm, i]
        Y <- V[itm]
        fit <- lm(Y ~ poly(X, degree, raw = TRUE)) # degree是指多项式的阶数
        continuation_value <- predict(fit, newdata = data.frame(X = X))
        exercise <- payoff[itm, i]
        V[itm] <- ifelse(exercise > continuation_value, exercise, V[itm])
    }
    return(mean(V) * discount_factor)
}

LSMprice_exp <- function(S0, K, r, T, sigma, steps, n, degree) {
    # steps 是指时间步数，n是指模拟次数，degree是指多项式的阶数
    dt <- T / steps
    t <- seq(0, T, by = dt)
    t <- matrix(rep(t, each = n), nrow = n)
    W <- matrix(rnorm(n * (steps + 1)), nrow = n, ncol = steps + 1) # 随机过程，生成股票价格路径的随机部分
    S <- matrix(0, nrow = n, ncol = steps + 1)
    S[, 1] <- S0
    for (i in 2:(steps + 1)) {
        S[, i] <- S[, i - 1] * exp((r - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * W[, i - 1])
    }

    payoff <- pmax(K - S, 0)
    V <- payoff[, steps + 1]
    discount_factor <- exp(-r * dt)
    for (i in (steps:1)) {
        V <- discount_factor * V
        itm <- (payoff[, i] > 0)
        X <- S[itm, i]
        Y <- V[itm]
        E <- exp(1:degree * X) # 指数函数作为基函数
        fit <- lm(Y ~ E - 1) # "- 1"表示不包含截距项
        continuation_value <- predict(fit, newdata = data.frame(E = E))
        exercise <- payoff[itm, i]
        V[itm] <- ifelse(exercise > continuation_value, exercise, V[itm])
    }
    return(mean(V) * discount_factor)
}

# Longstaff-Schwartz算法对美式看跌期权进行定价
LSMprice_log <- function(S0, K, r, T, sigma, steps, n, degree) {
    # steps 是指时间步数，n是指模拟次数，degree是指多项式的阶数
    dt <- T / steps
    t <- seq(0, T, by = dt)
    t <- matrix(rep(t, each = n), nrow = n)
    W <- matrix(rnorm(n * (steps + 1)), nrow = n, ncol = steps + 1) # 随机过程，生成股票价格路径的随机部分
    S <- matrix(0, nrow = n, ncol = steps + 1)
    S[, 1] <- S0
    for (i in 2:(steps + 1)) {
        S[, i] <- S[, i - 1] * exp((r - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * W[, i - 1])
    }

    payoff <- pmax(K - S, 0)
    V <- payoff[, steps + 1]
    discount_factor <- exp(-r * dt)
    for (i in (steps:1)) {
        V <- discount_factor * V
        itm <- (payoff[, i] > 0)
        X <- S[itm, i]
        Y <- V[itm]
        L <- log(1:degree * X) # 对数函数作为基函数
        fit <- lm(Y ~ L - 1) # "- 1"表示不包含截距项
        continuation_value <- predict(fit, newdata = data.frame(L = L))
        exercise <- payoff[itm, i]
        V[itm] <- ifelse(exercise > continuation_value, exercise, V[itm])
    }
    return(mean(V) * discount_factor)
}

S0 <- 36
K <- 40
r <- 0.06
T <- 1 # 这里是指1年
sigma <- 0.2
steps <- 50
n <- 10000
degree <- 3
folds <- 5

# 初始化结果向量
results_poly <- numeric(folds)
results_exp <- numeric(folds)
results_log <- numeric(folds)

# 进行交叉验证
for (i in 1:folds) {
    # 生成训练集和验证集
    train_indices <- sample(1:n, n * 0.8)
    test_indices <- setdiff(1:n, train_indices)
    
    # 使用多项式基函数
    put_price_poly <- LSMprice_poly(S0, K, r, T, sigma, steps, train_indices, degree)
    results_poly[i] <- put_price_poly
    
    # 使用指数基函数
    put_price_exp <- LSMprice_exp(S0, K, r, T, sigma, steps, train_indices, degree)
    results_exp[i] <- put_price_exp
    
    # 使用对数基函数
    put_price_log <- LSMprice_log(S0, K, r, T, sigma, steps, train_indices, degree)
    results_log[i] <- put_price_log
}

mean_poly <- mean(results_poly)
sd_poly <- sd(results_poly)
mean_exp <- mean(results_exp)
sd_exp <- sd(results_exp)
mean_log <- mean(results_log)
sd_log <- sd(results_log)

print(paste("多项式基函数的平均价格为：", mean_poly, "，标准差为：", sd_poly))
print(paste("指数基函数的平均价格为：", mean_exp, "，标准差为：", sd_exp))
print(paste("对数基函数的平均价格为：", mean_log, "，标准差为：", sd_log))

#输出结果为
[1] "多项式基函数的平均价格为： 4.51496492445049 ，标准差为： 0.109239440602219"
[1] "指数基函数的平均价格为： 3.99555595125704 ，标准差为： 0.000789493906819222"
[1] "对数基函数的平均价格为： 4.21967976327182 ，标准差为： 0.0317458072974367"
```

根据结果可知不同基函数选取、不同阶数选取对于定价结果的影响较大，因此在实务定价中需要根据市场及认知对参数进行调整。

### 收敛速度

1. 传统Monte Carlo加速

   ​    传统Monte Carlo模拟方法的一大困境是计算速度问题，伴随着GPU的不断发展，存在着一种叫做“Great Monte Carlo”的GPU并行解决方案，即利用GPU完成nested MC的计算，再利用CPU处理最后一层期望，而LSMC算法本身是通过数学算法将高运算量转化为低运算量，有学者测试指出，基于GPU加速的“Great Monte Carlo”效率不高于基于CPU的LSMC算法效率。

2. LSMC加速

      在上述背景下，尝试综合计算性能优化和数学算法上的优化，但是由于回归的存在，LSMC算法下GPU与CPU信息交换非常频繁，而GPU与CPU之间的数据转换速度较慢，无法显著提升收敛速度。即便如此，伴随着GPU和CPU数据转换速度的不断优化，LSMC和GPU并行计算的结合是提升收敛速度的重要方向。

## 参考文献

[^1]: Longstaff F A, Schwartz E S. Valuing American options by simulation: a simple least-squares approach[J]. The review of financial studies, 2001, 14(1): 113-147.
[^2]: Chaudhary S K. American options and the LSM algorithm: quasi-random sequences and Brownian bridges[J]. Journal of Computational Finance, 2005, 8(4).
[^3]: Kohler M, Krzyżak A, Todorovic N. Pricing of High‐Dimensional American Options by Neural Networks[J]. Mathematical Finance: An International Journal of Mathematics, Statistics and Financial Economics, 2010, 20(3): 383-410
[^4]: Becker S, Cheridito P, Jentzen A. Deep optimal stopping[J]. The Journal of Machine Learning Research, 2019, 20(1): 2712-2736.
[^5]: Becker S, Cheridito P, Jentzen A. Pricing and hedging American-style options with deep learning[J]. Journal of Risk and Financial Management, 2020, 13(7): 158.
[^6]: Tsitsiklis J N, Van Roy B. Regression methods for pricing complex American-style options[J]. IEEE Transactions on Neural Networks, 2001, 12(4): 694-703.
[^7]: Li Y, Szepesvari C, Schuurmans D. Learning exercise policies for american options[C]//Artificial intelligence and statistics. PMLR, 2009: 352-359.
[^8]: Herrera C, Krach F, Ruyssen P, et al. Optimal stopping via randomized neural networks[J]. arXiv preprint arXiv:2104.13669, 2021.↩
