#=============================================================================#
#=============================================================================#
#' mProbes feature selection algorithm
#' 
#' Implements the mProbes feature selection algorithm for Random Forests
#' 
#' @param x N x D predictors data frame where N - no. of samples, D - no. of features
#' @param y a vector of factors of length N, the target class (e.g as.factor("A", "A", "B", etc.))
#' @param nRepeat no. of times features are permuted (this is the sample size used when comparing importance
#' score for permuted vs real features)
#' @param ... arguments passed to the Random Forest classifier (e.g ntree, sampsize, etc.)
#' @return A list with the following components:
#' \item{impMetric}{2D x nRepeat matrix of variable importance measures for each predictor (permuted and not) 
#' for every repeat (Note: the permuted variables have the suffix "Perm")}
#' \item{FWER}{a numeric vector of length D with the family wise error rate computed for every feature}
#' @references \href{https://academic.oup.com/bioinformatics/article/28/13/1766/234473/Statistical-interpretation-of-machine-learning}{Huynh-Thu VA et al. Bioinformatics 2012}
#' @export
#' @examples 
#' bWant <- iris$Species %in% c("versicolor", "virginica")
#' x <- iris[bWant, 1:4]
#' y <- droplevels(as.factor(iris$Species[bWant]))
#'
#' out <- mProbes(x, y, 100)
#' 
mProbes <- function(x, y, nRepeat=100, ...)
{
    # Initialise some variables
    N <- dim(x)[1] # no. of samples
    D <- dim(x)[2] # no. of features
    impMetric <- matrix(NA, nrow=D*2, ncol=nRepeat)
    rownames(impMetric) <- c(colnames(x), paste0("Perm", colnames(x)))

    # Permute variables and fit forest
    for (iRepeat in 1:nRepeat) 
    {
        # Initialise set of artificial/permuted variables
        xTrainPerm <- matrix(NA, nrow=N, ncol=D) 
        colnames(xTrainPerm) <- paste0("Perm", colnames(x))
        for (iFeat in 1:D) # permute every feature
        {
            xTrainPerm[, iFeat] <- x[sample.int(N), iFeat]    
        }
        xTrainFull <- cbind(x, xTrainPerm) # xTrain + xTrainPerm
        fit <- randomForest::randomForest(y=y, x=xTrainFull, importance=TRUE, ...)
        impMetric[, iRepeat] <- randomForest::importance(fit, type=1, scale=F) 
    }
    
    # Compute FWER
    impMetricRandom <- impMetric[(D+1):(D*2), ]
    FWER <- rep(NA, D)
    for (iFeat in seq(D))
    {
        # column-wise comparison in R
        # http://stackoverflow.com/questions/16469053/r-matrix-vector-comparison
        FWER[iFeat] <- sum(apply(impMetricRandom > impMetric[iFeat, ], 2, any))/nRepeat
    }

    # Return impMetric and FWER
    return(list(impMetric=impMetric, FWER=FWER))
}
#=============================================================================#
#=============================================================================#
#' mProbes feature selection algorithm parallelised 
#' 
#' Same as \code{\link{mProbes}} but parallelised
#' 
#' @param x N x D predictors data frame where N - no. of samples, D - no. of features
#' @param y a vector of factors of length N, the target class (e.g as.factor("A", "A", "B", etc.))
#' @param nRepeat no. of times features are permuted (this is the sample size used when comparing importance
#' score for permuted vs real features)
#' @param ... arguments passed to the Random Forest classifier (e.g ntree, sampsize, etc.)
#' @param nThread no. of threads to use to run it in parallel (default: parallel::detected cores - 1)
#' @return A list with the following components:
#' \item{impMetric}{2D x nRepeat matrix of variable importance measures for each predictor (permuted and not) 
#' for every repeat (Note: the permuted variables have the suffix "Perm")}
#' \item{FWER}{a numeric vector of length D with the family wise error rate computed for every feature}
#' @references \href{https://doi.org/10.1093/bioinformatics/bts238}{Huynh-Thu VA et al. Bioinformatics 2012}
#' @export
#' @examples 
#' bWant <- iris$Species %in% c("versicolor", "virginica")
#' x <- iris[bWant, 1:4]
#' y <- droplevels(as.factor(iris$Species[bWant]))
#'
#' out <- mProbesParallel(x, y, 100, 4)
#' 
mProbesParallel <- function(x, y, nRepeat=100, nThread=parallel::detectCores()-1, ...)
{
    # Initialise some variables
    N <- dim(x)[1] # no. of samples
    D <- dim(x)[2] # no. of features
    
    # Initiate cluster
    cl <- makeCluster(nThread)
    registerDoParallel(cl)
    
    # Permute variables and fit forest
    impMetric <- foreach(i=icount(nRepeat), .packages = "tcltk", .combine=cbind, 
                         .export=c("randomForest", "importance")) %dopar% 
                         {
                             # Create progress bar first time round
                             if (!exists("pb")) 
                             {
                                 pb <- tkProgressBar(sprintf("Iteration %d/%d", i, nRepeat), min=1, max=nRepeat)
                             }
                             # Set progress bar
                             setTkProgressBar(pb, i, sprintf("Iteration %d/%d", i, nRepeat))
                             
                             # Create artificial/permuted variables
                             xTrainPerm <- apply(x, 2, function(x) {x[sample.int(length(x))]})
                             colnames(xTrainPerm) <- paste0("Perm", colnames(x))
                             
                             # Append artificial variables to original data and fit RF
                             xTrainFull <- cbind(x, xTrainPerm) # xTrain + xTrainPerm
                             fit <- randomForest(y=y, x=xTrainFull, importance=TRUE, ...)
                             importance(fit, type=1, scale=F)
                         }
    
    # Shut connection
    stopCluster(cl)
    
    # Name impMetric
    rownames(impMetric) <- c(colnames(x), paste0("Perm", colnames(x)))
    
    # Compute FWER
    impMetricRandom <- impMetric[(D+1):(D*2), ]
    FWER <- rep(NA, D)
    for (iFeat in seq(D))
    {
        # column-wise comparison in R
        # http://stackoverflow.com/questions/16469053/r-matrix-vector-comparison
        FWER[iFeat] <- sum(apply(impMetricRandom > impMetric[iFeat, ], 2, any))/nRepeat
    }
    
    # Return impMetric and FWER
    return(list(impMetric=impMetric, FWER=FWER))
}
#=============================================================================#
#=============================================================================#
#' Feature selection algorithm for Random Forests (RF)
#' 
#' Implements a modified mProbes/xRF feature selection algorithm within a
#' cross-validation loop
#'
#' @param x N x D predictors data frame where N - no. of samples, D - no. of features
#' @param y a vector of factors of length N, the target class (e.g as.factor("A", "A", "B", etc.))
#' @param nRepeat no. of times features are permuted (this is the sample size used when comparing importance
#' score for permuted vs real features)
#' @param kFold no. of cross-validation folds (default: 5)
#' @param rKeep percentage of predictors to ignore after first RF fit (0>rKeep<=1) (default: 0.3)
#' @param bParallel whether to use mProbes() or mProbesParallel() (default: True)
#' @param nThread no. of threads to use when running in parallel (default: parallel::detected cores - 1)
#' @param nSeed seed for cross-validation folds (default: 1983)
#' @param pCutOff Bonferonni corrected adjusted p-value cutoff when comparing importance scores of 
#' permuted vs real predictors (default: 0.05)
#' @param ... arguments passed to the Random Forest classifier (e.g ntree, sampsize, etc.)
#' @return A list with the following components:
#' \item{rKeepPredictors}{rKeep\% predictors kept after first RF fit}
#' \item{topPredictors }{top predictors (Bonferroni corrected p-values<pCutOff) of each fold}
#' \item{pValues}{Bonferroni corrected pValues for rKeepPredictors}
#' \item{ROC}{receiver operating characteristic curve, ROCR object}
#' \item{auc}{area under the ROC curve}
#' \item{confMatrix}{confusion matrix on test data}
#' \item{iiFolds}{indices of the cross-validation folds}
#' @references \href{https://doi.org/10.1093/bioinformatics/bts238}{Huynh-Thu VA et al. Bioinformatics 2012}
#' @references \href{http://dx.doi.org/10.1155/2015/471371}{Nguyen et al. The Scientific World Journal 2015}
#' @export
#' @examples 
#' bWant <- iris$Species %in% c("versicolor", "virginica")
#' x <- iris[bWant, 1:4]
#' y <- droplevels(as.factor(iris$Species[bWant]))
#'
#' out <- featselectRF(x, y, nodesize=3, ntree=1001)
#' 
featselectRF <- function(x, y, nRepeat=100, kFold=5, rKeep=0.3, bParallel=TRUE, 
                         nThread=parallel::detectCores()-1, nSeed=1983, pCutOff=0.05, ...)
{
    # Initialise some variables
    rKeepPredictors <- list() # predictors kept after first RF fit
    topPredictors <- list() # top predictors of each fold
    pValues <- list() # pValues for rKeepPredictors  
    confMatrix <- list() # confusion matrix on test data
    AUC <- list() # area under the curve for test data (Random Forests fit)
    ROC <- list() # ROCR object for receiver operating characteristic curve
    yTestAll <- c() # all testing y data 
    yPredAll <- c() # all predictions of y data
    
    # Cross-validation loop
    set.seed(nSeed) # set seed for reproducibility
    folds <- caret::createFolds(y, k=kFold)  
    for (i in 1:kFold)
    {
        #=====================================#
        # Prepare train and test data
        #=====================================#
        des <- paste0("Fold", i) # description of what fold we're in
        writeLines(paste0(des, "/", kFold))
        # Training data set
        yTrain <- y[-folds[[i]]]
        xTrain <- x[-folds[[i]], ]
        # Testing data set
        yTest <- y[folds[[i]]]
        xTest <- x[folds[[i]], ]
        
        #=====================================#
        # Fit random forest to all data
        #=====================================#
        writeLines("a) Fitting random forest to training data (all predictors)...")
        sampSize <- round(0.85*min(table(yTrain)))
        fit <- randomForest::randomForest(x=xTrain, y=yTrain, importance=TRUE, 
                                          sampsize=rep(sampSize, length(table(yTrain))), ...)
        imp <- randomForest::importance(fit, type=1, scale=F)
        
        # Keep top rKeep predictors
        NKeep <- ceiling(rKeep*dim(xTrain)[2])
        iiSort <- order(imp, decreasing=T) 
        xTrain <- xTrain[, iiSort[1:NKeep]]
        rKeepPredictors[[des]] <- colnames(xTrain)
        
        #=====================================#
        # mProbes on rKeep predictors
        #=====================================#
        writeLines("b) Running mProbes feature selection algorithm on the kept predictors...")
        if (bParallel) # Run on a single or multiple cores
        {
            fit <- mProbesParallel(x=xTrain, y=yTrain, nRepeat=nRepeat, 
                                   sampsize=rep(sampSize, length(table(yTrain))), ...)  
        } else {
            fit <- mProbes(x=xTrain, y=yTrain, nRepeat=nRepeat, 
                           sampsize=rep(sampSize, length(table(yTrain))), ...) 
        }
        impMetric <- fit$impMetric
        
        #=====================================#
        # Compute statistical significance between
        # permuted and real importance scores
        #=====================================#
        writeLines("c) Compute statistical significance between real and random features...")
        D <- dim(xTrain)[2] # no. of features
        impMetricRandom <- impMetric[(D+1):(D*2), ]
        impMetricReal <- impMetric[1:D, ]
        # get vector of maximum importance score of random feature from each repeat
        maxPerm <- apply(impMetricRandom, 2, max)
        # calculate statistical significance scores using Wilcoxon rank-sum test
        pValues[[des]] <- apply(impMetricReal, 1, function(x) wilcox.test(x, maxPerm, paired=T, alternative="greater")$p.value)
        pValues[[des]] <- pValues[[des]]*D # Bonferroni correction                
        topPredictors[[des]] <- colnames(xTrain)[which(pValues[[des]]<pCutOff)] 
        
        #=====================================#
        # Evaluate final model on test data
        #=====================================#
        writeLines("d) Fit final model to top predictors...")
        D <- length(topPredictors[[i]]) # no. of features left
        if (D >= 1) 
        {
            # Set training data set
            xTrain <- xTrain[, topPredictors[[des]]]
            xTest <- xTest[, topPredictors[[des]]]
            
            if (D==1)
            {
                # Fit logistic regression model rather than RF
                df <- data.frame(x=xTrain, y=yTrain)
                fit <- glm(y ~ ., data=df, family=binomial(link='logit'))
                dfTest <- data.frame(x=xTest)
                yPred <- predict(fit, newdata=dfTest, type="response")
                # Test confusion matrix
                predClass <- ifelse(yPred>0.5, levels(yTrain)[2], levels(yTrain)[1])
                confusionMatrix <- table(yTest, predClass)
                classError <- 1 - diag(prop.table(confusionMatrix, 1)) # Compute misclassification rates
                confMatrix[[des]] <- cbind(confusionMatrix, classError)
            } else
            {
                # Fit random forest
                sampSize <- round(0.85*min(table(yTrain)))
                fit <- randomForest::randomForest(x=xTrain, y=yTrain, xtest=xTest, ytest=yTest, 
                                                  sampsize=rep(sampSize, length(table(yTrain))), ...)
                yPred <- fit$test$votes[, 2]
                confMatrix[[des]] <- fit$test$confusion # test confusion matrix
            }
            # Compute area under curve
            predObj <- ROCR::prediction(yPred, yTest)
            AUC[[des]] <- round(unlist(slot(ROCR::performance(predObj, "auc"), "y.values")), digits=3)
            ROC[[des]] <- ROCR::performance(predObj, "tpr", "fpr")
            
            # Save yTest/yPred all
            yTestAll <- c(yTestAll, yTest)
            yPredAll <- c(yPredAll, yPred)
        } else 
        {
            # No predictors chosen
            AUC[[des]] <- NA
            confMatrix[[des]] <- NA
            ROC[[des]] <- NA
        }
    }
    
    #=====================================#
    # Compute average ROC/AUC/confusion matrix
    #=====================================#
    predObj <- ROCR::prediction(yPredAll, yTestAll) 
    AUC[["Average"]] <- round(unlist(slot(ROCR::performance(predObj, "auc"), "y.values")), digits=3)
    ROC[["Average"]] <-  ROCR::performance(predObj, "tpr","fpr") 
    confMatrix[["Average"]] <- Reduce('+', confMatrix)
    confMatrix[["Average"]][, 3] <- 1 - diag(prop.table(confMatrix[["Average"]][, 1:2], 1))
        
    #=====================================#
    # Combine results in a list and return
    #=====================================#
    out <- list(rKeepPredictors=rKeepPredictors, topPredictors=topPredictors, 
                pValues=pValues, auc=AUC, ROC=ROC,
                confMatrix=confMatrix, iiFolds=folds)
    
    return(out)
}