.lass<-function(z, lambda){
  test <- z-lambda
  u <- ifelse( as.numeric(test)  <= 0, 0, abs(as.numeric(test)))
  thetahat<-sign(z)*u
  return(thetahat)
}

.group <- function(x){
  group <- sqrt(sum(x^2))
  return(group)
}

.cont_loss <- function(x, y, beta, tau, methods ){

  n <- nrow(x)
  p <- ncol(x)

  beta <- matrix(beta, nrow = p,ncol = 1)
  for( j in 1:p){
    sd <-  sd(x[,j])
    sd <- ifelse(sd == 0, 1, sd)
    x[,j] <- x[,j]/sd
  }




  u <- y- x%*%beta

  if("regression" %in% methods){
    score <- -t(x)%*%u/n
    loglike <- 0.5*t(u)%*%u/n
  }

  if("Huber regression" %in% methods){
    id <- abs(u) > tau
    u[id] <- sign(u[id])*tau
    score <- -t(x)%*%u/n
    loglike <- 0.5*t(u[id == FALSE])%*%(u[id == FALSE])/n + tau*sum(abs(u[id == TRUE]))/n - tau^2/2*sum(id)/n
  }

  if("adaptive Huber" %in% methods){

    tau <- tau*sd(u)*sqrt(n/log(p))/10
    id <- abs(u) > tau
    u[id] <- sign(u[id])*tau
    score <- -t(x)%*%u/n
    loglike <- 0.5*t(u[id == FALSE])%*%(u[id == FALSE])/n + tau*sum(abs(u[id == TRUE]))/n - tau^2/2*sum(id)/n

  }

  return(list(score, loglike, tau, u))
}


.binary_loss <- function(x, y, beta ){

  n <- nrow(x)
  p <- ncol(x)
  resp <- unique(y)
  y <- y==resp[1]

  for( j in 1:p){
    sd <-  sd(x[,j])
    sd <- ifelse(sd == 0, 1, sd)
    x[,j] <- x[,j]/sd
  }

  beta <- matrix(beta, nrow = p,ncol = 1)

    fit <- exp(- x%*%beta)
    mu <- 1/(1+fit)
    den <- c(log(mu[y]), log(1-mu[y]))
    den[which(den == -Inf)] <- -700
    loglike <- sum(-den)/n
    u <- y- mu
    score <-  -t(x)%*%(y- mu)/(n)

  return(list(score, loglike, u))
}




.Main_Algo <- function(x, y, Kn, p, n, lambda, beta = 0.1, tau = 1.45 , gamma = 1,
                      methods = "adaptive Huber" ,  import_w = 1,
                      tol = 0.05, max_iter = 100, Complete = TRUE){


  K1 <- Kn[1]
  K2 <- Kn[2]
  K <- K1+K2


  if( length(import_w ) == 1 ) {
    import_w <-  rep(import_w, K)
  }

  adjustments <-  matrix(0, p, K)
  if (Complete == FALSE) {

    for(k in 1:K){
      for (j in 1:p){
        missing <- is.na(x[[k]][,j])
        if (sum(missing)<=0.1*n[k]){
          x[[k]][missing,j] <- mean(x[[k]][missing==FALSE,j] )
        }
        else{
          x[[k]][,j] <- 0
        }
        adjustments[j, k] = sum(missing)/n[k]
      }
    }
  }
  adj_w <- rowSums(adjustments)
  adj_w <- K/(K-adj_w)





  if (is.null(dim(beta)) ==TRUE & length(beta) == 1) {
    beta <- matrix(beta, nrow = p, ncol = K)
  }


  diff <- distance  <- NULL
  j <- 0
  J <- max_iter
  t_old <- 1
  beta_old <- beta

  while( j < J) {

    score <- likelihood <-  NULL

    if(K1 !=0 ){
      for (k in 1:K1){
        Square <-  .cont_loss(x[[k]],y[[k]],beta_old[,k],tau[k], methods[k])
        score <- cbind(score, Square[[1]])
        likelihood <- c(likelihood, Square[[2]])
      }

    }

    if(K2 != 0 ){
      for (k in (K1+1):K){
        Square <-  .binary_loss(x[[k]],y[[k]],beta_old[,k] )
        score <- cbind(score, Square[[1]])
        likelihood <- c(likelihood, Square[[2]])
      }
    }


    score <- score*import_w
    obj  <-   (likelihood)%*%import_w + lambda* sum(apply(beta_old, 1, .group))
    obj1 <- obj+1

    #eigenvl <- numeric(K)
    #for(k in 1:K){
    #  eigenvl[k] <-  max(eigen(t(x[[k]])%*%x[[k]])$values)
    #}

      #max(abs(unlist(x)))


    i <- 0
    gamma <- gamma
    while (obj <= obj1 ){

      L2 <- sqrt(rowSums((beta_old - score/gamma)^2)/adj_w)
      beta <- adj_w^0.5*c(.lass(L2,as.numeric(lambda/gamma))/L2)*(beta_old - score/gamma)

      if( sum( beta[,1]==0) == p)
        warning( "p == 0, penalty parameter is too large")


      score1  <- likelihood1 <- NULL
      if(K1 !=0 ){
        for (k in 1:K1){
          Square <-  .cont_loss(x[[k]],y[[k]],beta[,k],tau[k], methods[k])
          score1 <- cbind(score1, Square[[1]])
          likelihood1 <- c(likelihood1, Square[[2]])
        }

      }
      if(K2 != 0 ){
        for (k in (K1+1):K){
          Square <-  .binary_loss(x[[k]],y[[k]],beta[,k] )
          score1 <- cbind(score1, Square[[1]])
          likelihood1 <- c(likelihood1, Square[[2]])
        }
      }

      score1 <- score1*import_w
      obj1 <-  likelihood1%*%import_w + sum((score1)*(beta - beta_old))
      + gamma*sum((beta - beta_old)^2)/2 + lambda* sum(apply(beta, 1, .group))

      i <- i + 1
      gamma <- gamma*(i*2)

      if (i > 10)
        break
    }

    distance <- beta-beta_old
    diff <- max(abs(distance))
    j <- j +1
    if (diff <= tol) break

    t_new <- (1 + sqrt(1 + 4 * t_old * t_old)) / 2
    beta_old <- beta + (t_old-1)*distance/t_new
    t_old <- t_new

  }


  Selected_List <- which( rowSums(beta^2)!=0)

  Task_weights <- import_w

  Task_type <- numeric(K)

  name1 <- name2 <- NULL
  if(K1 !=0 ){

  name1 <- c(paste(paste(rep("Reg",K1),1:K1)))
  Task_type[1:K1] <- methods

  }

  if(K2 !=0 ){
    name2  <- c(paste(rep("Class",(K2)),(1+K1):K))
    Task_type[(1+K1):K] <- "logit"
  }

  colnames(beta) <- c(name1, name2)
  rownames(beta) <- paste("Features",1:p)

  beta[which(adjustments==1)] <- NA

  Output <- list("beta" = beta,
                 "Task_type" = Task_type,
                 "Task_weights" = Task_weights,
                 "Selected_List" = Selected_List,
                 "Max_iteration" = j)

  return(Output)

}


.info_fitting <- function(x, y, beta,Kn, n, p, tau, methods, link, import_w = 1, Complete = TRUE,
                         alpha = 1, gamma = 1, tol=0.05){



  K1 <- Kn[1]
  K2 <- Kn[2]
  K <- K1+K2

  com.loglike <- numeric(K)
  u <- matrix(1, ncol=K, nrow = max(n))

  if (Complete == FALSE) {

    for(k in 1:K){
      for (j in 1:p){
        missing <- is.na(x[[k]][,j])
        if (sum(missing)<=0.1*n[k]){
          x[[k]][missing,j] <- mean(x[[k]][missing==FALSE,j] )
        }
        else{
          x[[k]][,j] <- 0
        }
      }
    }
  }


   if( length(import_w ) == 1 ) {
    import_w <-  rep(import_w, K)
  }

  feature <- which(beta[,1] != 0 )
  p1 <- length(feature)

  if (p1 == 0) {
    warning( "p == 0, penalty value is too large")
    info_criterion = NA
    Results = list("Bayesian_Information" = info_criterion)
    return(Results )
  }

  if (p1 >= min(n)) {
    warning("p >= n, models cannot converge")
    info_criterion = NA
    Results = list("Bayesian_Information" = info_criterion)
    return(Results )
  }

  MSE <- numeric(K1)
  MAE <- numeric(K1)


  if(K1 !=0 ){

    for (k in 1:K1){

      method <- methods[k]
      if("regression" %in% method){

        design <- x[[k]][,feature]
        model <-  glm(y[[k]]~design-1)
        com.loglike[k] <- (model$aic - 2*p1)/n[k]
        #linear_pred <- rowSums( t(as.matrix(model$coefficients*t(design) )))
        u[1:n[k],k] <- y[[k]]- model$fitted.values
        MSE[k] <- sum((u[1:n[k],k])^2)/n[k]
        MAE[k] <- sum(abs(u[1:n[k],k]))/n[k]

      }

      else {
        design <- x[[k]][,feature]
        sel_beta <- beta[feature,k]

        J <- 100
        j <- 0
        diff <- 1
        while( j < J) {

          sel_beta_old <- sel_beta
          Square <-  .cont_loss(design,y[[k]], sel_beta,tau[k], method)
          score <-  Square[[1]]
          obj <-  Square[[2]]
          obj1 <- obj+1


          i <- 0

          while (obj <= obj1 ){

            sel_beta <-  sel_beta_old  - gamma^{-1}*score
            Square <-  .cont_loss(design,y[[k]], sel_beta,tau[k], method)
            obj1 <-  Square[[2]]
            i <- i + 1
            gamma <- gamma*(i*2)

            if (i > 500) break

          }

          distance <- sel_beta- sel_beta_old
          diff <- max(abs(distance))
          j <- j +1
          if (diff <= tol) break


        }
        u[1:n[k],k] <-  Square[[4]]
        com.loglike[k] <- obj1


        MSE[k] <- sum((u[1:n[k],k])^2)/n[k]
        MAE[k] <- sum(abs(u[1:n[k],k]))/n[k]




      }

    }

  }



  AUC <- numeric(K2)
  DEV <- numeric(K2)

  if(K2 !=0 ){



    for (k in (K1+1):K){

      design <- x[[k]][,feature]


        model <-  glm(y[[k]]~design-1, family = binomial(link = "logit"))
        linear_pred <- rowSums( t(as.matrix(model$coefficients*t(design) )))

        fit <-  exp(linear_pred)
        phi<-1
        DEV[k-K1] <- model$deviance
        AUC[k-K1] <-  auc(y[[k]], model$fitted)[1]

      com.loglike[k] <- (model$aic - 2*p1)/n[k]

      u[1:n[k],k] <- (y[[k]]- model$fitted.values)*phi
    }


  }

  ## Penalty
  h.matrix <- j.matrix <- 0
  ind.score <-  matrix(0, nrow = max(n),ncol = 1)

  for( k in 1:K){

    design <- x[[k]][,feature]
    h.matrix <- bdiag(h.matrix, design )
    #ind.score <- bdiag(ind.score, design*u[1:n[k],k])

    if( length(unique(n))==1){
    ind.score <- cbind(ind.score,design*u[,k])
    }else{
    j.matrix <- bdiag(j.matrix, var(design*u[1:n[k],k]))
    }
  }
  h.matrix <- as.matrix(h.matrix[-1,-1])
  ind.score <- ind.score[ ,-1]
  if( length(unique(n))==1){
  j.matrix <- cov(as.matrix(ind.score))
  }else{
    j.matrix <- as.matrix(j.matrix[-1,-1])
  }
  h.matrix <- t(h.matrix)%*%h.matrix
  df <- alpha*sum(diag(j.matrix%*%solve(h.matrix)))*log(p1)



  info_criterion <- 2*com.loglike%*%import_w + df
  degree_of_freedom <- df/alpha
  compo_likelihood <- com.loglike%*%import_w

  Bayesian_Information <- c(compo_likelihood, degree_of_freedom,info_criterion)
  names(Bayesian_Information ) <- c("Compo_likelihood", "Degree of freedom", "Info criterion")

  Reg_Error <- Class_Perform <-NULL
  if(K1 !=0 ){
  Reg_Error <- rbind(MSE, MAE)
  rownames(Reg_Error) <- c("MSE", "MAE")
  colnames(Reg_Error) <- c(paste(rep("Reg",K1),1:K1))
  }

  if(K2 !=0 ){
  Class_Perform <- rbind(AUC, DEV)
  rownames(Class_Perform) <- c("AUC", "DEV")
  colnames(Class_Perform) <-c(paste(rep("Class",(K2)),(1+K1):(K1+K2)))
  }

    Residuals <- u


    Results <- list("Bayesian_Information" = Bayesian_Information,
                  "Reg_Error" = Reg_Error,
                  "Class_Perform" = Class_Perform,
                  "Residuals" = Residuals)

  return(Results)


}

