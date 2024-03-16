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
        fit <- lm(Y ~ poly(X, degree, raw = TRUE)) # degree是指多项式的阶数
        continuation_value <- predict(fit, newdata = data.frame(X = X))
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