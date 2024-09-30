# For the bootstrap method, p-values are calculated based on the proportion of bootstrap samples that are less than or greater than zero, providing a two-sided p-value.
# For the normal approximation method, p-values are calculated using the standard normal distribution based on the z-scores of the estimates.
#
#

ccmm <- function(y, M, tr, x=NULL, w=NULL, method.est.cov = "bootstrap", n.boot=2000, sig.level=0.05, tol=1e-6, max.iter=5000){
  if(!is.matrix(M)) stop("M must be in matrix format!!!", call.=FALSE);
  if(method.est.cov!="bootstrap" & method.est.cov!="normal")
    stop('Parameter (method.est.cov) must be either "bootstrap" or "normal.approx".', call.=FALSE);
  k <- ncol(M);
  rslt.A <- est.comp.param.boot(M, tr, x, w, method.est.cov, n.boot);
  if(is.null(rslt.A)) stop("Estimation of compositional parameters failed.", call.=FALSE)
  
  rslt.B <- est.debias.B(y, M, tr, x, tol, max.iter);
  
  if(method.est.cov=="bootstrap"){
    boot.debias.B <- mvrnorm(n=n.boot,mu=rslt.B$debias.B, Sigma=rslt.B$cov.debias.B);
    DE <- mean(boot.debias.B[,k+1]);
    DE.CI <- quantile(boot.debias.B[,k+1], probs=c(sig.level/2, 1-sig.level/2), names=TRUE);
    p.DE <- 2 * min(mean(boot.debias.B[,k+1] <= 0), mean(boot.debias.B[,k+1] >= 0))
    
    rslt.A[rslt.A <= 0] <- 1e-8
    log_k_rslt_A <- log(k * rslt.A)
    # Check for NaNs and Infs
    if(any(is.na(log_k_rslt_A)) || any(is.infinite(log_k_rslt_A))){
      stop("Non-positive values encountered in rslt.A after adjustment. Cannot compute logarithm.")
    }
    
    boot.IDEs <- log_k_rslt_A * boot.debias.B[,1:k];
    IDEs <- colMeans(log_k_rslt_A) * colMeans(boot.debias.B[,1:k]);
    names(IDEs) <- paste("C", 1:k, sep="");
    IDE.CIs <- apply(boot.IDEs, 2, function(x) quantile(x, probs=c(sig.level/2, 1-sig.level/2), names=TRUE));
    colnames(IDE.CIs) <- paste("C", 1:k, sep="");
    IDE.pvalues <- sapply(1:k, function(j) 2 * min(mean(boot.IDEs[,j] <= 0), mean(boot.IDEs[,j] >= 0)))
    names(IDE.pvalues) <- paste("C", 1:k, sep = "")
    
    TIDE <- sum(IDEs);
    boot.TIDE <- rowSums(log(rslt.A) * boot.debias.B[,1:k]);
    TIDE.CI <- quantile(boot.TIDE, probs=c(sig.level/2, 1-sig.level/2), names=TRUE);
    p.TIDE <- 2 * min(mean(boot.TIDE <= 0), mean(boot.TIDE >= 0))
    
    return(list(DE=DE, DE.CI=DE.CI, DE.pvalue = p.DE, 
                TIDE=TIDE, TIDE.CI=TIDE.CI, TIDE.pvalue = p.TIDE, 
                IDEs=IDEs, IDE.CIs=IDE.CIs, IDE.pvalues = IDE.pvalues));
  } else{
    DE <- rslt.B$debias.B[k+1];
    Var.DE <- rslt.B$cov.debias.B[k+1,k+1];
    z.DE <- DE / sqrt(Var.DE)
    p.DE <- 2 * (1 - pnorm(abs(z.DE)))
    
    # Ensure rslt.A$E.ln.kA is finite
    if(any(is.na(rslt.A$E.ln.kA)) || any(is.infinite(rslt.A$E.ln.kA))){
      stop("Invalid values encountered in rslt.A$E.ln.kA. Cannot proceed with calculations.")
    }
    
    IDEs <- rslt.A$E.ln.kA * rslt.B$debias.B[1:k];
    names(IDEs) <- paste("C", 1:k, sep="");
    Var.IDEs <- numeric(k);
    for(i in 1:k){
      Var.IDEs[i] <- idvar_j(rslt.A$E.ln.kA[i], rslt.B$debias.B[i], rslt.A$Var.ln.kA[i,i], rslt.B$cov.debias.B[i,i]);
    }
    names(Var.IDEs) <- paste("C", 1:k, sep="");
    IDE.zscores <- IDEs / sqrt(Var.IDEs)
    IDE.pvalues <- 2 * (1 - pnorm(abs(IDE.zscores)))
    names(IDE.pvalues) <- paste("C", 1:k, sep = "")
    
    TIDE <- sum(IDEs);
    Var.TIDE <- tvar(rslt.A$E.ln.kA, rslt.B$debias.B[1:k], rslt.A$Var.ln.kA, rslt.B$cov.debias.B[1:k,1:k]);
    z.TIDE <- TIDE / sqrt(Var.TIDE)
    p.TIDE <- 2 * (1 - pnorm(abs(z.TIDE)))
    
    return(list(DE=DE, Var.DE=Var.DE, DE.pvalue = p.DE, 
                TIDE=TIDE, Var.TIDE=Var.TIDE, TIDE.pvalue = p.TIDE, 
                IDEs=IDEs, Var.IDEs=Var.IDEs, IDE.pvalues = IDE.pvalues));
  }
}

convert2comp <- function(x){
  return(x/sum(x));
}

get_D <- function(k, u){
  D_of_u <- matrix(-u, nrow=(k-1), ncol=(k-1));
  diag(D_of_u) <- (k-1)*u;
  return(D_of_u);
}

build_mat_Ds <- function(n.comp, n.covar){
  n.rc <- (n.comp-1)*(n.covar+2);
  mat_Ds <- matrix(0, nrow=n.rc, ncol=n.rc);
  return(mat_Ds);
}

est.comp.param <- function(M, tr, x, w){
  n <- nrow(M); k <- ncol(M);
  if(is.null(x)){
    mat_Ds <- build_mat_Ds(k, 0);
    sum.wtr.sq <- sum(w*tr^2);
    sum.wtr <- sum(w*tr);
    mat_Ds[1:(k-1), 1:(k-1)] <- get_D(k, sum.wtr.sq);
    mat_Ds[k:(2*(k-1)), 1:(k-1)] <- get_D(k, sum.wtr);
    mat_Ds[k:(2*(k-1)), k:(2*(k-1))] <- get_D(k, sum(w));
    mf <- k*t(log(M)) %*% (w*tr) - sum(log(M)*w*tr);
    mg <- k*t(log(M)) %*% w - sum(log(M)*w);
    vec_ms <- c(mf[-k], mg[-k]);
  } else{
    n.x <- ncol(x);
    mat_Ds <- build_mat_Ds(k, n.x);
    sum.wtr.sq <- sum(w*tr^2);
    sum.wtr <- sum(w*tr);
    mat_Ds[1:(k-1), 1:(k-1)] <- get_D(k, sum.wtr.sq);
    mat_Ds[k:(2*(k-1)), 1:(k-1)] <- get_D(k, sum.wtr);
    mat_Ds[k:(2*(k-1)), k:(2*(k-1))] <- get_D(k, sum(w));
    sum.wxtr <- colSums(w*x*tr);
    sum.wx <- colSums(w*x);
    for(i in 3:(n.x+2)){
      mat_Ds[((i-1)*(k-1)+1):(i*(k-1)), 1:(k-1)] <- get_D(k, sum.wxtr[i-2]);
      mat_Ds[((i-1)*(k-1)+1):(i*(k-1)), k:(2*(k-1))] <- get_D(k, sum.wx[i-2]);
      sum.wxx <- colSums(w*x*x[,(i-2)]);
      for(j in i:(n.x+2)){
        mat_Ds[((j-1)*(k-1)+1):(j*(k-1)), ((i-1)*(k-1)+1):(i*(k-1))] <- get_D(k, sum.wxx[j-2]);
      }
    }
    mf <- k*t(log(M)) %*% (w*tr) - sum(log(M)*w*tr);
    mg <- k*t(log(M)) %*% w - sum(log(M)*w);
    m_psy <- k*t(log(M)) %*% (w*x) - tcrossprod(rep(1,k), colSums(t(log(M)*w) %*% x));
    vec_ms <- c(mf[-k], mg[-k], m_psy[-k,]);
  }
  
  mat_Ds[upper.tri(mat_Ds)] <- t(mat_Ds)[upper.tri(mat_Ds)];
  
  # Check if mat_Ds is singular
  if (det(mat_Ds) == 0 || is.nan(det(mat_Ds))) {
    warning("Matrix mat_Ds is singular or nearly singular. Adding small value to diagonal for regularization.")
    mat_Ds <- mat_Ds + diag(1e-6, nrow(mat_Ds))
  }
  
  # Solve the system
  solution <- try(solve(mat_Ds, vec_ms), silent = TRUE)
  if(inherits(solution, "try-error")){
    warning("Could not solve mat_Ds. Returning NULL.")
    return(NULL)
  }
  
  ests.k_minus_1 <- matrix(exp(solution), nrow=(k-1));
  
  # Handle Inf or NaN values
  if(any(is.infinite(ests.k_minus_1)) || any(is.nan(ests.k_minus_1))){
    warning("Invalid values in ests.k_minus_1. Returning NULL.")
    return(NULL)
  }
  ests_k <- 1/(colSums(ests.k_minus_1)+1);
  ests.k_minus_1 <- ests.k_minus_1 * tcrossprod(rep(1,k-1), ests_k);
  comp.param <- rbind(ests.k_minus_1, ests_k);
  rownames(comp.param) <- 1:k;
  
  if(is.null(x)){
    colnames(comp.param) <- c("a", "m0");
  } else{
    colnames(comp.param) <- c("a", "m0", paste("h", 1:n.x, sep=""));
  }
  return(comp.param);
}

est.comp.param.boot <- function(M, tr, x, w, method.est.cov, n.boot){
  n <- nrow(M); k <- ncol(M);
  if(is.null(w)){
    w <- rep(1, n);
  } else{
    if(!is.vector(w)) w <- as.vector(w);
  }
  if(!is.vector(tr)) tr <- as.vector(tr);
  
  valid_boot <- 0
  max_iter <- n.boot * 2  # Allow up to double the required iterations
  boot.A <- matrix(0, nrow=n.boot, ncol=k);
  iter <- 1
  while(valid_boot < n.boot && iter <= max_iter){
    indx <- sample(1:n, n, replace=TRUE);
    boot.M <- M[indx,];
    boot.tr <- tr[indx];
    if(is.null(x)){
      boot.params <- est.comp.param(boot.M, boot.tr, x, w);
    } else{
      boot.x <- x[indx,];
      if(!is.matrix(boot.x)) boot.x <- as.matrix(boot.x);
      boot.params <- est.comp.param(boot.M, boot.tr, boot.x, w);
    }
    if(!is.null(boot.params)){
      boot.A[valid_boot + 1,] <- boot.params[,"a"];
      valid_boot <- valid_boot + 1
    }
    iter <- iter + 1
  }
  if(valid_boot < n.boot){
    warning("Not enough valid bootstrap samples were obtained.")
    boot.A <- boot.A[1:valid_boot, , drop=FALSE]
  }
  
  if(method.est.cov=="bootstrap"){ ### Return bootstrap samples of parameter A
    return(boot.A);
  } else{ ### Return E(log(kA)) and Var(log(kA))
    boot.ln.kA <- log(k*boot.A);
    E.boot.ln.kA <- colMeans(boot.ln.kA);
    ### Estimate variance-covariance of bootstrap log(ka)
    c.boot.ln.kA <- sweep(boot.ln.kA, 2, E.boot.ln.kA, FUN="-"); #center boot.ln.kA
    n.row.H.mat <- k*(k+1)/2;
    H.mat <- matrix(0, n.row.H.mat, n.row.H.mat);
    diag(H.mat)[1:k] <- 1;
    diag(H.mat)[(k+1):n.row.H.mat] <- -2;
    cb <- combn(1:k,2);
    b <- rep(0,n.row.H.mat);
    for(i in 1:k){
      b[i] <- s.h(c.boot.ln.kA[,i], n.boot);
    }
    for(i in (k+1):n.row.H.mat){
      H.mat[i, cb[,i-k]] = 1;
      b[i] <- s.h(c.boot.ln.kA[,cb[1,i-k]]-c.boot.ln.kA[,cb[2,i-k]], n.boot)
    }
    E.s <- solve(H.mat, b^2);
    E.Var.boot.ln.kA <- matrix(0, k, k);
    diag(E.Var.boot.ln.kA) <- E.s[1:k];
    E.Var.boot.ln.kA[lower.tri(E.Var.boot.ln.kA)] <- E.s[(k+1):n.row.H.mat];
    E.Var.boot.ln.kA <- t(E.Var.boot.ln.kA);
    E.Var.boot.ln.kA[lower.tri(E.Var.boot.ln.kA)] <- t(E.Var.boot.ln.kA)[lower.tri(E.Var.boot.ln.kA)];
    return(list(E.ln.kA=E.boot.ln.kA, Var.ln.kA=E.Var.boot.ln.kA));
  }
}

s.h <- function(x, n.boot){
  s <- floor(0.05*n.boot) + 1;
  l <- n.boot - floor(0.05*n.boot);
  q.boot.h <- quantile(x, probs=seq(s,l)/n.boot);
  s.boot.h <- 1/(0.847*n.boot)*sum(qnorm(seq(s,l)/n.boot)*q.boot.h)+0.103*(q.boot.h[length(q.boot.h)]-q.boot.h[1]);
  return(s.boot.h);
}

### Estimate the initial value of lambda for scaled lasso
get.init.lambda <- function(n.sample, n.feature, tol){
  f <- function(x, p) {(qnorm(1-x/p))^4 + 2*((qnorm(1-x/p))^2) - x;}
  k <- uniroot(f, lower=0, upper=n.feature-1, tol=tol, p=n.feature)$root;
  lambda_0 <- sqrt(2/n.sample)*qnorm(1-k/n.feature);
  return(lambda_0);
}

### Lasso with a linear constraint
lasso_constr <- function(y, x, contr2, lam, tol, max.iter){
  n <- nrow(x); p <- ncol(x); k <- dim(contr2)[1];
  gramC <- crossprod(contr2); gramX <- crossprod(x);
  diagC <- diag(gramC); diagX <- diag(gramX);
  dediagC <- gramC - diag(diagC); dediagX <- gramX - diag(diagX);
  covXY <- crossprod(x, y);
  
  mu <- 1; bet <- rep(1,p)/p; bet0 <- rep(0,p); iter <- 0;
  if (sum(abs(contr2))==0){ #No constraint
    term0 <- (covXY-dediagX%*%bet)/n;
    term2 <- diagX/n;
    while (sum(abs(bet-bet0))>tol & iter<max.iter){
      bet0 <- bet;
      for(j in 1:p){
        term1 <- sign(term0[j])*max(0, abs(term0[j])-lam);
        bet[j] <- term1/term2[j];
        dif <- bet[j] - bet0[j];
        term0 <- term0 - dediagX[,j]*dif/n;
      }
      iter <- iter+1;
    }
  } else{ #With constraint
    ksi <- rep(0,k); ksi0 <- rep(1,k);
    term0 <- (covXY-dediagX%*%bet)/n - mu*(t(contr2)%*%ksi + dediagC%*%bet);
    term2 <- diagX/n + mu*diagC;
    while (sum(abs(ksi-ksi0))>tol && iter<max.iter){
      ksi0 <- ksi; iter2 <- 0; bet0 <- bet0 + 1;
      while (sum(abs(bet-bet0))>tol && iter2<1000){
        bet0 <- bet;
        for(j in 1:p){
          term1 <- sign(term0[j])*max(0, abs(term0[j])-lam);
          bet[j] <- term1/term2[j];
          dif <- bet[j] - bet0[j];
          term0 <- term0 - dediagX[,j]*dif/n - dediagC[,j]*dif*mu;
        }
        iter2 <- iter2 + 1;
      }
      dif2 <- contr2%*%bet;
      ksi <- ksi + dif2; 
      term0 <- term0 - mu*t(contr2)%*%dif2;
      iter <- iter + 1;
    }
  }
  return(bet);
}

### Scaled lasso for compositional data
ccmm_slr <- function(y, Q, contr, tol, max.iter){
  n <- nrow(Q); p <- ncol(Q);
  mn.y <- mean(y); mn.Q <- colMeans(Q);
  cent.Q <- (Q-tcrossprod(rep(1, n), mn.Q));
  contr2 <- t(as.matrix(contr));
  cent.y <- y-mn.y;
  lam0 <- get.init.lambda(n, p, tol);
  sigma <- 1; sigma_2 <- 1; sigma_s <- 0.5; iter <- 1;
  while (abs(sigma-sigma_s)>0.01 & iter<100){
    iter <- iter + 1;
    sigma <- (sigma_s + sigma_2)/2;
    lam <- sigma*lam0;
    bet2 <- lasso_constr(cent.y, cent.Q, contr2, lam, tol, max.iter);
    s <- sum(abs(bet2)>0.001);
    s <- min(s, n-1);
    sigma_s <- base::norm(cent.y-cent.Q%*%bet2, type="2")/sqrt(n-s-1);
    sigma_2 <- sigma;
  }
  if(iter==100) print("Not converge!");
  sigma <- sigma_s;
  bet <- bet2;
  intercp <- mn.y - mn.Q%*%bet;
  return(list(beta=bet, intercept=intercp, sigma=sigma, lambda0=lam0));
}

est.debias.B <- function(y, M, tr, x, tol, max.iter){
  if(!is.vector(y)) y <- as.vector(y);
  if(!is.vector(tr)) tr <- as.vector(tr);
  n <- length(y); k <- ncol(M);
  if(is.null(x)){
    Z <- cbind(log(M), tr);
    n.trx <- 1;
  } else{
    Z <- cbind(log(M), tr, x);
    if(!is.matrix(x)) x <- as.matrix(x);
    n.trx <- 1 + ncol(x);
  }
  n.vrs <- k + n.trx;
  contr <- c(rep(1/sqrt(k), k), rep(0, n.trx));
  Z.til <- Z %*% (diag(n.vrs)-tcrossprod(contr));
  est.param <- ccmm_slr(y, Z.til, contr, tol, max.iter);
  cent.y <- scale(y, scale=FALSE);
  cent.Z.til <- scale(Z.til, scale=FALSE);
  gam0 <- est.param$lambda0/3; Sig <- crossprod(cent.Z.til)/n;
  Sig2 <- Sig - diag(diag(Sig));
  Q <- diag(n.vrs) - tcrossprod(contr);
  M.til <- matrix(0, n.vrs, n.vrs);
  for(i in 1:n.vrs){
    gam <- gam0/2;
    while(gam<0.5){
      gam <- gam*2;
      mi <- rep(1,n.vrs);
      mi0 <- rep(0,n.vrs);
      iter <- 1;
      while(sum(abs(mi-mi0))>tol & iter<max.iter){
        mi0 <- mi;
        for(j in 1:n.vrs){
          v <- -Sig2[j,]%*%mi+Q[j,i];
          mi[j] <- sign(v)*max(0, abs(v)-gam)/Sig[j,j];
        }
        iter <- iter + 1;
      }
      if(iter<max.iter) break;
    }
    M.til[i,] <- mi;
  }
  M.til <- Q%*%M.til;
  debias.B <- est.param$beta + M.til%*%t(cent.Z.til)%*%(cent.y-cent.Z.til%*%est.param$beta-(est.param$intercept*rep(1,n)))/n;
  cov.debias.B <- est.param$sigma^2*M.til%*%Sig%*%t(M.til)/n;
  return(list(debias.B=t(debias.B), cov.debias.B=cov.debias.B));
}

tvar <- function(Ea, Eb, VCa, VCb){
  k <- length(Ea);
  term1 <- sum(Ea^2 * diag(VCb) + Eb^2 * diag(VCa));
  term2 <- VCa * tcrossprod(Eb) + VCb * tcrossprod(Ea);
  term2 <- sum(term2[upper.tri(term2)]);
  return(term1+2*term2);
}

idvar_j <- function(Ea, Eb, Va, Vb){
  return(Ea^2*Vb + Eb^2*Va + Va*Vb);
}

### Add a check whether x is in a matrix form
ccmm.sensitivity <- function(rh, y, M, tr, x=NULL, w=NULL){
  n <- nrow(M); k <- ncol(M)
  if(is.null(w)){
    w <- rep(1, n);
  } else{
    if(!is.vector(w)) w <- as.vector(w);
  }
  if(length(rh)==1){
    v.rho <- rep(rh, k-1)
  } else{
    v.rho <- rh
  }
  if(!is.null(x)) x <- as.matrix(x)
  rvcs <- com.var.cov.sen(y, M, tr, x, w, k)
  rbs <- est.tide.sen(v.rho, rvcs, k)
  return(rvcs$alt_a%*%rbs)
}

est.tide.sen <- function(v.rho, rvcs, k, n.iter=1000, tol.bs=1e-6, init.var.U2=1e-6){
  iter <- 1
  del.bs <- 1
  bs.rho <- rep(0, k-1);
  var.U2 <- var.U2.new <- init.var.U2
  while(del.bs > tol.bs & iter < n.iter){
    bs.rho.new <- est.bs.rho(v.rho, var.U2.new, rvcs$var_alt_U1, rvcs$var_U0, rvcs$cov_U01)
    var.U2.new <- est.var.U2(v.rho, bs.rho.new, rvcs$var_alt_U1, rvcs$var_U0)
    del.bs <- max(c(abs(bs.rho.new-bs.rho), abs(var.U2.new-var.U2)))
    bs.rho <- bs.rho.new
    var.U2 <- var.U2.new
    iter <- iter + 1
  }
  if(iter>n.iter) stop("Not converging!!!", call.=FALSE)
  return(bs.rho.new)
}

com.var.cov.sen <- function(y, M, tr, x, w, k){
  if(is.null(x)){
    U0 <- resid(lm(y~tr))
  } else{
    U0 <- resid(lm(y~tr+x))
  }
  var_U0 <- var(U0)
  comp_params <- est.comp.param(M,tr,x,w)
  alt_M <- sweep(log(M[,-k]), 1, log(M[,k]), "-")
  alt_cparams <- log(sweep(comp_params[-k,],2,comp_params[k,],"/"))
  alt_U1 <- alt_M - (cbind(tr,1,x) %*% t(alt_cparams))
  S_alt_U1 <- cov(alt_U1)
  cov_U01 <- cov(U0, alt_U1)
  return(list(alt_a=alt_cparams[,"a"], var_U0=var_U0, var_alt_U1=S_alt_U1, cov_U01=cov_U01))
}

est.var.U2 <- function(v.rho, v.b, S_alt.U1, var.U0){
  term1 <- t(v.b)%*%S_alt.U1%*%v.b
  term2 <- sqrt(diag(S_alt.U1))%*%(v.b*v.rho)
  if(term1 > (term2^2+var.U0)) stop("Variance must be positive!!! Try with a weaker correlation.", call.=FALSE)
  vr.U2 <- sqrt(term2^2+var.U0-term1)-term2
  if(vr.U2>1e6) stop("Variance approaches to infinity!!! Try with a weaker correlation.", call.=FALSE)
  if(vr.U2<0) stop("Variance must be positive!!! Try with a weaker correlation.", call.=FALSE)
  return(vr.U2)
}

est.bs.rho <- function(v.rho, var.U2, S_alt.U1, var.U0, cov.U01){
  rhs.eq <- cov.U01 - v.rho*var.U2*sqrt(diag(S_alt.U1))
  bs.rh <- solve(S_alt.U1, as.numeric(rhs.eq))
  return(bs.rh)
}

rho.range <- function(y,M,tr,x=NULL,w=NULL){
  rho.L <- -1
  rho.U <- 1
  rho.L <- find.rho.range(rho.L, 0.1, y, M, tr, x=NULL, w=NULL)
  rho.L <- find.rho.range(rho.L, 0.01, y, M, tr, x=NULL, w=NULL)
  rho.U <- find.rho.range(rho.U, -0.1, y, M, tr, x=NULL, w=NULL)
  rho.U <- find.rho.range(rho.U, -0.01, y, M, tr, x=NULL, w=NULL)
  return(c(rho.L+0.02, rho.U-0.02))
}

find.rho.range <- function(rho, stp, y, M, tr, x=NULL, w=NULL){
  err.m <- try(ccmm.sensitivity(rho, y, M, tr, x, w), silent=TRUE)
  while(inherits(err.m, "try-error")){
    rho.prev <- rho
    rho <- rho + stp
    err.m <- try(ccmm.sensitivity(rho, y, M, tr, x, w), silent=TRUE)
  }
  return(rho.prev)
}

ccmm.sa <- function(y, M, tr, x=NULL, w=NULL, stp=0.01){
  range_rho <- rho.range(y, M, tr, x, w)
  x.rho <- seq(range_rho[1], range_rho[2], stp)
  y.tide <- NULL
  for(i in 1:length(x.rho)){
    y.tide[i] <- ccmm.sensitivity(x.rho[i], y, M, tr, x, w)
  }
  return(cbind(rho=x.rho, TIDE=y.tide))
}

tide.ci.zero.rho <- function(y, M, tr, x=NULL, w=NULL, n.boot=2000){
  n <- nrow(M); k <- ncol(M)
  if(is.null(w)){
    w <- rep(1, n);
  } else{
    if(!is.vector(w)) w <- as.vector(w);
  }
  tide_sen <- rep(0, n.boot)
  if(is.null(x)){
    for(i in 1:n.boot){
      indx_sen <- sample(1:n, n, replace=TRUE);
      bts.y <- y[indx_sen]
      bts.M <- M[indx_sen,];
      bts.tr <- tr[indx_sen];
      rvcs <- com.var.cov.sen(bts.y, bts.M, bts.tr, x, w, k)
      rbs <- est.tide.sen(0, rvcs, k)
      tide_sen[i] <- rvcs$alt_a%*%rbs
    }
  } else{
    for(i in 1:n.boot){
      indx_sen <- sample(1:n, n, replace=TRUE);
      bts.y <- y[indx_sen]
      bts.M <- M[indx_sen,];
      bts.tr <- tr[indx_sen];
      bts.x <- as.matrix(x[indx_sen,]);
      rvcs <- com.var.cov.sen(bts.y, bts.M, bts.tr, bts.x, w, k)
      rbs <- est.tide.sen(0, rvcs, k)
      tide_sen[i] <- rvcs$alt_a%*%rbs
    }
  }
  return(tide_sen)
}