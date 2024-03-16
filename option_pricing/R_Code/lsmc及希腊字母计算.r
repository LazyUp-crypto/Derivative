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
        fit <- lm(Y ~ poly(X, degree, raw = TRUE)) # degree是指多项式的阶数
        continuation_value <- predict(fit, newdata = data.frame(X = X))
        exercise <- payoff[itm, i]
        V[itm] <- ifelse(exercise > continuation_value, exercise, V[itm])
    }
    return(mean(V) * discount_factor)
}

# 使用有限差分-中间差分计算期权Greek值

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
    return((myCall_1 - myCall_2) / (2 * diff))
}

S0 <- 36
K <- 40
r <- 0.06
T <- 1 # 这里是指1年
sigma <- 0.2
steps <- 50
n <- 10000
degree <- 5
put_price <- LSMprice(S0, K, r, T, sigma, steps, n, degree)
print(paste("使用Longstaff-Schwartz算法得美式看跌期权的价格为：", put_price))

print(paste("delta值为：", delta()))
print(paste("gamma值为：", gamma()))
print(paste("theta值为：", theta()))
print(paste("rho值为：", rho()))
print(paste("vega值为：", vega()))