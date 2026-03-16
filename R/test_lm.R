test_lm <- function(X, Y, network, alpha, cell_num, n = 100, nfold = 10, Mthread = TRUE, Mcore = 24){

    set.seed(1)
    m <- nrow(X)
    index0 <- sample(cut(seq(m), breaks = nfold, labels = F))

    print("|**************************************************|")
    print("Perform cross-validation on X with true label")
    MSE_test_real <- NULL
    pb1 <- progress_bar$new(total = nfold)
    for (j in 1:nfold){
        c_index <- which(index0 == j)
        X_train <- X[-c_index,]
        Y_train <- Y[-c_index]
        fit <- NULL
        while (is.null(fit$fit)){
            set.seed(123)
            fit <- APML1(X_train, Y_train, family = "gaussian", penalty = "Net", alpha = alpha, Omega = network, nlambda = 100)
        }
        index <- which.min(abs(fit$fit$nzero - cell_num))
        Coefs <- as.numeric(fit$Beta[,index])
        Cell1 <- Coefs[which(Coefs > 0)]
        Cell2 <- Coefs[which(Coefs < 0)]

        X_test <- X[c_index,]
        Y_test <- Y[c_index]
        MSE_test_real[j] <- mean((Y_test - X_test%*%Coefs)^2)
        rm(X_train, Y_train, X_test, Y_test, fit, Coefs, Cell1, Cell2)
        gc()

        #pb1$tick()
        Sys.sleep(1 / 100)
        if (j == nfold) cat("Finished!\n")
    }

    print("|**************************************************|")
    print("Perform cross-validation on X with permutated label")
    permutation_payload_bytes <- estimate_parallel_payload_bytes(X, Y, network)
    permutation_results <- parallel_task_apply(seq_len(n), function(i) {
        set.seed(i + 100)
        mse_test_back <- matrix(0, nfold, 1, dimnames = list(paste0("Testing_", 1:nfold), "MSE"))
        Y2 <- Y[sample(m)]
        for (j in seq_len(nfold)) {
            c_index <- which(index0 == j)
            X_train <- X[-c_index,]
            Y_train <- Y2[-c_index]
            fit <- NULL
            while (is.null(fit$fit)) {
                set.seed(123)
                fit <- APML1(X_train, Y_train, family = "gaussian", penalty = "Net", alpha = alpha, Omega = network, nlambda = 100)
            }
            index <- which.min(abs(fit$fit$nzero - cell_num))
            Coefs <- as.numeric(fit$Beta[,index])

            X_test <- X[c_index,]
            Y_test <- Y2[c_index]
            mse_test_back[j] <- mean((Y_test - X_test%*%Coefs)^2)
            rm(X_train, Y_train, X_test, Y_test, fit, Coefs)
            gc()
        }
        mse_test_back
    }, Mthread = Mthread, Mcore = Mcore, payload_bytes = permutation_payload_bytes)
    MSE_test_back <- vector("list", n)
    pb2 <- progress_bar$new(total = n)
    for (i in seq_len(n)) {
        MSE_test_back[[i]] <- permutation_results[[i]]
        if (i == n) cat("Finished!\n")
    }
    statistic  <- mean(MSE_test_real)
    background <- NULL
    for (i in 1:n){
        background[i] <- mean(MSE_test_back[[i]][,1])
    }
    p <- sum(background < statistic)/n

    print(sprintf("Test statistic = %s", formatC(statistic, format = "f", digits = 3)))
    print(sprintf("Reliability significance test p = %s", formatC(p, format = "f", digits = 3)))

    return(list(statistic = statistic,
                p = p,
                MSE_test_real = MSE_test_real,
                MSE_test_back = MSE_test_back))
}
