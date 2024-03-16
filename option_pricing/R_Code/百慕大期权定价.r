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
option_price <- LSMprice_Bermuda(S0, K, r, T, sigma, steps, n, degree, exercise_dates)
print(paste("百慕大期权价格: ", option_price))