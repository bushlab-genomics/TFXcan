suppressMessages(library(glmnet))
suppressMessages(library(methods))
suppressMessages(library(dplyr))
suppressMessages((library(reshape2)))
suppressMessages((library(parallel)))

argv <- commandArgs(trailingOnly = TRUE)
output_prefix <- argv[1]
n_k_folds <- as.numeric(argv[2])
alpha <- as.numeric(argv[3])
expression_RDS <- argv[4]
geno_file <- argv[5]
snp_annot_rds <- argv[6]
gene_annot_rds <- argv[7]
chrom <- as.numeric(argv[8])

tfscores_file <- paste0(output_prefix, "_results/", output_prefix, "_variant_tfscores.csv")
out_dir <- paste0(output_prefix,"_results/")
covariance_out <- paste0(out_dir,output_prefix,"_",chrom, "_",n_k_folds,"_",alpha, "_tfxcan_covariances.txt")
expression <- readRDS(expression_RDS)
groupings <- data.frame(c(row.names(expression)),c(row.names(expression)),
                        c(sample(1:n_k_folds, length(row.names(expression)),replace = T)),
                        c(sample(1:5, length(row.names(expression)), replace =T)))

colnames(groupings) <- c("V1","V2","grouping","fivefold")
rownames(groupings) <- groupings$V1
generate_fold_ids <- function(n_samples, n_folds=10) {
  n <- ceiling(n_samples / n_folds)
  fold_ids <- rep(1:n_folds, n)
  sample(fold_ids[1:n_samples])
}

calc_R2 <- function(y, y_pred) {
  tss <- sum(y**2)
  rss <- sum((y - y_pred)**2)
  1 - rss/tss
}
weight_mult <- function(mat,vec){
  multi_mat <- matrix(data = NA, ncol = ncol(mat), nrow = nrow(mat))
  for (i in 1:nrow(multi_mat)){
    multi_mat[i,] <- vec*mat[i,]
  }
  return(multi_mat)
}

nested_cv_elastic_net_perf <- function(x, y,n_samples, n_train_test_folds, n_k_folds, alpha) {
  rmse_folds <- rep(0, n_train_test_folds)
  R2_folds <- rep(0, n_train_test_folds)
  corr_folds <- rep(0, n_train_test_folds)
  zscore_folds <- rep(0, n_train_test_folds)
  pval_folds <- rep(0, n_train_test_folds)
  # get the preset fivefold ids:
  train_test_fold_ids <- groupings[rownames(y),'fivefold']
  y=as.vector(y)
  for (test_fold in 1:n_train_test_folds) {
    train_idxs <- which(train_test_fold_ids != test_fold)
    test_idxs <- which(train_test_fold_ids == test_fold)
    x_train <- x[train_idxs, ]
    y_train <- y[train_idxs]
    x_test <- x[test_idxs, ]
    y_test <- y[test_idxs]
    cv_fold_ids <- generate_fold_ids(length(y_train), n_k_folds)
    y_pred <- tryCatch({
      # Fit model with training data.
      fit <- cv.glmnet(x_train, y_train,nfolds = n_k_folds, alpha = alpha, type.measure='mse', foldid = cv_fold_ids)
      # Predict test data using model that had minimal mean-squared error in cross validation.
      predict(fit, x_test, s = 'lambda.min')},
      # if the elastic-net model did not converge, predict the mean of the y_train (same as all non-intercept coef=0)
      error = function(cond) rep(mean(y_train), length(y_test)))
    R2_folds[test_fold] <- calc_R2(y_test, y_pred)
    res <- summary(lm(y_test~y_pred))
    rmse_folds[test_fold] <-sqrt(1-res$r.squared)*(res$sigma)
    
    corr_folds[test_fold] <- ifelse(sd(y_pred) != 0, cor(y_pred, y_test), 0)
    zscore_folds[test_fold] <- atanh(corr_folds[test_fold])*sqrt(length(y_test) - 3) # Fisher transformation
    pval_folds[test_fold] <- ifelse(sd(y_pred) != 0, cor.test(y_pred, y_test)$p.value, runif(1))
  }
  rmse_avg <- mean(rmse_folds)
  R2_avg <- mean(R2_folds)
  R2_sd <- sd(R2_folds)
  rho_avg <- mean(corr_folds)
  rho_se <- sd(corr_folds)
  rho_avg_squared <- rho_avg**2
  zscore_est <- sum(zscore_folds) / sqrt(n_train_test_folds)
  zscore_pval <- 2*pnorm(abs(zscore_est), lower.tail = FALSE)
  pval_est <- pchisq(-2 * sum(log(pval_folds)), 2*n_train_test_folds, lower.tail = F)
  list(rmse_avg=rmse_avg,R2_avg=R2_avg, R2_sd=R2_sd, pval_est=pval_est, rho_avg=rho_avg, rho_se=rho_se, rho_zscore=zscore_est, rho_avg_squared=rho_avg_squared, zscore_pval=zscore_pval)
}
class(expression) <- 'numeric'
genotype <- read.table(geno_file, header = TRUE, row.names = 'SNP_ID', stringsAsFactors = FALSE, check.names = F)
variant_weights <- read.csv(tfscores_file, stringsAsFactors = F)
variant_weights <- variant_weights[variant_weights$cube_root_score_scaled != 0,]
gene_annot <- readRDS(gene_annot_rds)

sample_intersection <- intersect(row.names(expression), colnames(genotype))
snp_annot <- readRDS(snp_annot_rds)


rownames(snp_annot) <- snp_annot$SNP_ID
gene_intersection <- intersect(unique(gene_annot$gene_id), unique(variant_weights$hgnc_symbol))

gene_intersection1 <- intersect(gene_intersection, colnames(expression))
snp_intersection <- intersect(rownames(genotype), unique(variant_weights[,2]))
variant_weights <- variant_weights[variant_weights[,2] %in% snp_intersection,]
genotype <- genotype[snp_intersection,]
expression <- expression[sample_intersection,]
genotype <-  genotype[,sample_intersection]
expression <- expression[, gene_intersection1]
exp_samples <- rownames(expression)
exp_genes <- colnames(expression)
n_samples <- length(exp_samples)
n_genes <- length(exp_genes)
seed <- sample(1:100000,1)
set.seed(seed)
genotype <- t(genotype)
n_train_test_folds <- 5
row.names(gene_annot) <- gene_annot$gene_id
write_covariance <- function(gene, cisgenos, model_SNP_IDs, model_varIDs, covariance_out) {
  model_geno <- cisgenos[,as.character(model_varIDs), drop=FALSE]
  geno_cov <- cov(model_geno)
  cov_df <- data.frame(gene=character(),SNP_ID1=character(),SNP_ID2=character(), covariance=double())
  for (i in 1:length(model_SNP_IDs)) {
    for (j in i:length(model_SNP_IDs)) {
      cov_df <- rbind(cov_df, data.frame(gene=gene,SNP_ID1=model_SNP_IDs[i], SNP_ID2=model_SNP_IDs[j], covariance=geno_cov[i,j]))
    }
  }
  write.table(cov_df, file = covariance_out, append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = " ")
}
resultscol <- c('gene_id','alpha', 'n_snps_in_window', 'n_snps_in_model', 'lambda_min_mse',
                'adjust_R2','rmse', 'test_R2_avg', 'test_R2_sd', 'cv_R2_avg', 'cv_R2_sd', 'in_sample_R2',
                'nested_cv_fisher_pval', 'rho_avg', 'rho_se', 'rho_zscore', 'rho_avg_squared', 'zscore_pval',
                'cv_rho_avg', 'cv_rho_se', 'cv_rho_avg_squared', 'cv_zscore_est', 'cv_zscore_pval', 'cv_pval_est')


resultsarray <- array(0, c(length(exp_genes), 24))
dimnames(resultsarray)[[1]] <- gene_intersection1
dimnames(resultsarray)[[2]] <- resultscol
out_name <- paste0(out_dir,output_prefix,"_",chrom, "_",n_k_folds,"_",alpha,"_tfxcan_results.txt")
write(resultscol, file = out_name, ncolumns = 24, sep = "\t")
weightcol <- c("gene","SNP_ID","ref","alt","beta","alpha")
workingweight <- paste0(out_dir, output_prefix,"_",chrom, "_",n_k_folds,"_",alpha,"_tfxcan_weights.txt")
write(weightcol, file = workingweight, ncol = 6, sep = "\t")
log_df <- data.frame(chrom, length(gene_intersection1), seed, alpha)
colnames(log_df) <- c('chr', 'n_genes', 'seed_for_cv', 'alpha')
write.table(log_df, file = paste0(out_dir,output_prefix,'_chr',chrom,'_tfxcan_model_log.txt'),
            quote = FALSE, row.names = FALSE, sep = "\t")
#The block below runs the tfxcan model, many components of it are derived from the epixcan scripts
for (i in 1:length(gene_intersection1)){
  cat(i, "/", length(gene_intersection1), "\n")  
  gene <- gene_intersection1[i]
  gene_name <- gene
  weights_gene_df <- variant_weights[variant_weights$hgnc_symbol == gene,]
  cis_snps <- weights_gene_df[,2]
  cisgenos <- genotype[,cis_snps]
  cisgenos[is.na(cisgenos)] <- 0
  cis_weights <- as.vector(weights_gene_df[,"cube_root_score_scaled"])
  resultsarray[gene,] <- tryCatch({
    if (ncol(cisgenos) < 2) {
      c(gene, rep(23, NA))
      #bestbetas  <- data.frame()
    } else {
      cisgenos_weights <- weight_mult(cisgenos,cis_weights)
      rownames(cisgenos_weights) <- rownames(cisgenos)
      colnames(cisgenos_weights)  <- colnames(cisgenos)
      exppheno <- expression[,gene]
      exppheno <- scale(exppheno, center = TRUE, scale = TRUE)
      exppheno[is.na(exppheno)] <- 0
      perf_measures <- nested_cv_elastic_net_perf(as.matrix(cisgenos_weights), exppheno,n_samples, n_train_test_folds, n_k_folds, alpha)
      R2_avg <- perf_measures$R2_avg
      rmse_avg <- perf_measures$rmse_avg
      R2_sd <- perf_measures$R2_sd
      pval_est <- perf_measures$pval_est
      rho_avg <- perf_measures$rho_avg
      rho_se <- perf_measures$rho_se
      rho_zscore <- perf_measures$rho_zscore
      rho_avg_squared <- perf_measures$rho_avg_squared
      zscore_pval <- perf_measures$zscore_pval
      # get the preset 10-fold CV grouping ids:
      cv_fold_ids <- groupings[rownames(exppheno),'grouping']
      fit <- tryCatch(cv.glmnet(as.matrix(cisgenos),as.vector(exppheno),
                                nfolds = n_k_folds, alpha = 0.5,
                                type.measure='mse',
                                foldid = cv_fold_ids, keep = TRUE),
                      error = function(cond) {message('Error'); message(geterrmessage()); list()})
      if (length(fit) > 0) {
        fit.df <- data.frame(fit$cvm, fit$lambda, 1:length(fit$cvm))
        cv_R2_folds <- rep(0, n_k_folds)
        cv_corr_folds <- rep(0, n_k_folds)
        cv_zscore_folds <- rep(0, n_k_folds)
        cv_pval_folds <- rep(0, n_k_folds)
        best_lam_ind <- which.min(fit$cvm)
        for (j in 1:n_k_folds) {
          fold_idxs <- which(cv_fold_ids == j)
          adj_expr_fold_pred <- fit$fit.preval[fold_idxs, best_lam_ind]
          cv_R2_folds[j] <- calc_R2(exppheno[fold_idxs], adj_expr_fold_pred)
          cv_corr_folds[j] <- ifelse(sd(adj_expr_fold_pred) != 0, cor(adj_expr_fold_pred, exppheno[fold_idxs]), 0)
          cv_zscore_folds[j] <- atanh(cv_corr_folds[j])*sqrt(length(exppheno[fold_idxs]) - 3) # Fisher transformation
          cv_pval_folds[j] <- ifelse(sd(adj_expr_fold_pred) != 0, cor.test(adj_expr_fold_pred, exppheno[fold_idxs])$p.value, runif(1))
        }
        cv_R2_avg <- mean(cv_R2_folds)
        cv_R2_sd <- sd(cv_R2_folds)
        adj_expr_pred <- predict(fit, as.matrix(cisgenos), s = 'lambda.min')
        training_R2 <- calc_R2(exppheno, adj_expr_pred)
        cv_rho_avg <- mean(cv_corr_folds)
        cv_rho_se <- sd(cv_corr_folds)
        cv_rho_avg_squared <- cv_rho_avg**2
        # Stouffer's method for combining z scores.
        cv_zscore_est <- sum(cv_zscore_folds) / sqrt(n_k_folds)
        cv_zscore_pval <- 2*pnorm(abs(cv_zscore_est), lower.tail = FALSE)
        cv_pval_est <- pchisq(-2 * sum(log(cv_pval_folds)), 2*n_k_folds, lower.tail = F)
        if (fit$nzero[best_lam_ind] > 0) {
          n_fit=fit$nzero[best_lam_ind]
          # Discuss with Yungil: use R2 of correlation between predicted and observed expression values for all samples:
          y_all=predict(fit, as.matrix(cisgenos), s = 'lambda.min')
          weights <- fit$glmnet.fit$beta[which(fit$glmnet.fit$beta[,best_lam_ind] != 0), best_lam_ind]
          weighted_snps <- names(fit$glmnet.fit$beta[,best_lam_ind])[which(fit$glmnet.fit$beta[,best_lam_ind] != 0)]
          # correlation R2
          corr_R=ifelse(sd(y_all) != 0, cor(y_all, exppheno),0)
          corr_R2=corr_R**2
          # adjusted correlation R2
          adj_R2=1-(1-corr_R2)*(n_samples-1)/(n_samples-1-n_fit)
          bestweightlist <- weighted_snps
          bestweightinfo <- snp_annot[bestweightlist,]
          colnames(bestweightinfo)<-c('chr', 'pos', 'refAllele', 'effectAllele','SNP_ID')
          weighttable <- bestweightinfo
          weighttable$weight <- weights
          weightfile <- cbind(gene,as.character(weighttable[,"SNP_ID"]),as.character(weighttable[,"refAllele"]),
                              as.character(weighttable[,"effectAllele"]),weighttable[,"weight"], alpha)
          write(t(weightfile), file = workingweight, ncolumns = 6, append = TRUE, sep = "\t")
          #weightfile$gene=gene 
          write_covariance(gene, cisgenos, weighttable[,"SNP_ID"], weighttable[,"SNP_ID"],covariance_out)
          c(gene,alpha, ncol(cisgenos), fit$nzero[best_lam_ind], fit$lambda[best_lam_ind],adj_R2,rmse_avg, R2_avg, R2_sd, cv_R2_avg, cv_R2_sd, training_R2, pval_est,
            rho_avg, rho_se, rho_zscore, rho_avg_squared, zscore_pval, cv_rho_avg, cv_rho_se, cv_rho_avg_squared, cv_zscore_est, cv_zscore_pval, cv_pval_est)

        }else { # there are no snps selected regarding the gene
          c(gene, alpha, ncol(cisgenos), 0, fit$lambda[best_lam_ind],0,rmse_avg, R2_avg, R2_sd,
            cv_R2_avg, cv_R2_sd, training_R2, pval_est, rho_avg, rho_se, rho_zscore, rho_avg_squared, zscore_pval,
            cv_rho_avg, cv_rho_se, cv_rho_avg_squared, cv_zscore_est, cv_zscore_pval, cv_pval_est)
        }
      }else {# there is no model
        c(gene,rep(NA, 23))
      }
    }
  }, error = function(e) c(gene,rep(NA, 23)) )
  write(resultsarray[gene,], file = out_name, ncolumns = 24, append = TRUE, sep = "\t")
}



