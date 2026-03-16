

#################################
#####  Logistic Regression  #####
#################################

LogL0=function(x, y, Omega=NULL, alpha=1.0, lambda=NULL, nlambda=100, rlambda=NULL, nfolds=1, foldid=NULL, inzero=TRUE, wbeta=rep(1,ncol(x)), sgn=rep(1,ncol(x)), isd=FALSE, keep.beta=FALSE, thresh=1e-6, maxit=1e+5, threshP=1e-5, Mthread=FALSE, Mcore=1L) {
  # alpha=alphai; lambda=NULL; nlambda=nlambdai; rlambda=NULL; inzero=T; isd=FALSE; keep.beta=TRUE; thresh=1e-8; maxit=1e+5; threshP=1e-5; wbeta=rep(1,ncol(x)); Omega=NULL

  N0=nrow(x); p=ncol(x)

  ### Adaptive
  aPen=ifelse(all(wbeta>0), TRUE, FALSE)

  ### Lambda path
  if (is.null(lambda)) {
    ilambda=1
    if (is.null(rlambda)) {
      rlambda=ifelse(N0>p, 0.0001, 0.01)
    }
    lambda=(rlambda)^(c(0:(nlambda-1))/(nlambda-1))
  } else {
    ilambda=0
    nlambda=length(lambda)
  }


  if (is.null(Omega)) {
    penalty=ifelse(alpha==1, "Lasso", "Enet")
    adaptive=ifelse(any(wbeta!=1), TRUE, FALSE)
  } else {
    penalty=ifelse(alpha==1, "Lasso", "Net")
    adaptive=c(ifelse(any(wbeta!=1), TRUE, FALSE),ifelse(any(sgn!=1), TRUE, FALSE))
    W=prepare_net_graph(Omega, sgn, add_intercept=TRUE)
  }


  #####  Run  #####
  x1=cbind(1.0,x); storage.mode(y)="double"; p1=p+1;
  wbeta1=c(0.0, wbeta); wbetai=rep(1.0,p1); wbetai[1]=0.0

  ##
  x1i=x1; p1i=p1
  wbeta1i=wbeta1; wbetaii=wbetai

  indexi=which(apply(x1i,2,sd)==0)[-1]
  if (length(indexi)>0) {
    x1i=x1i[,-indexi]; p1i=p1i-length(indexi)
    wbeta1i=wbeta1i[-indexi]; wbetaii=wbetaii[-indexi]
  }

  out=switch(penalty,
             "Net"=NetLogC(x1i, y, alpha, lambda, nlambda, ilambda, wbeta1i, wbetaii, W$Omega, W$loc, W$nadj, p1i, N0, thresh, maxit, threshP),
             EnetLogC(x1i, y, alpha, lambda, nlambda, ilambda, wbeta1i, wbetaii, p1i, N0, thresh, maxit, threshP)
  )
  nlambdai=out$nlambda ## number of lambdas
  if (nlambdai==0)
    return(NULL)
  lambdai=out$lambda[1:nlambdai]

  if (length(indexi)>0) {
    betai=matrix(0,nrow=p1,ncol=nlambdai)
    betai[-indexi,]=out$Beta
    out$Beta=betai

    betai=matrix(0,nrow=p1,ncol=nlambdai)
    betai[-indexi,]=out$BetaSTD
    out$BetaSTD=betai
  }

  out$Beta=Matrix(out$Beta[, 1:nlambdai,drop=F], sparse=TRUE)
  out$BetaSTD=Matrix(out$BetaSTD[, 1:nlambdai], sparse=TRUE)
  out$nzero=apply(out$Beta!=0, 2, sum)
  out$flag=out$flag[1:nlambdai]
  out$LL=out$LL[1:nlambdai]
  out$nzero=apply(out$Beta!=0,2,sum)

  if (nfolds==1 & is.null(foldid)) {
    fit=data.frame(lambda=lambdai, pDev=(out$ll0-out$LL)/out$ll0, nzero=pmax(0,out$nzero-1))
    if (!isd) {
      return(list(Beta=out$Beta, fit=fit, penalty=penalty, adaptive=adaptive, flag=out$flag))
    } else {
      return(list(Beta=out$BetaSTD, fit=fit, penalty=penalty, adaptive=adaptive, flag=out$flag))
    }
  } else {

    ###  Split data for cross-validation
    if (is.null(foldid)) {
      foldid=sample(rep(seq(nfolds), length=N0))
    } else {
      nfolds=max(foldid)
    }
    tb=table(foldid)
    N0i=numeric(nfolds); Nf=numeric(nfolds)
    for (i in 1:nfolds) {
      N0i[i]=sum(tb[-i]); Nf[i]=tb[i]
    }
    weighti=as.vector(tapply(rep(1,N0), foldid, sum))
    fold_index=seq_len(nfolds)

    compute_trim_log_rows <- function(Betai, BetaSTDi, numi, numi2) {
      fold_rows <- parallel_task_apply(fold_index, function(i) {
        temid=foldid==i; numj=min(Betao[i], numi)
        Betaj=Betai[, i]; BetaSTDj=BetaSTDi[, i]

        BetaSTDjj=BetaSTDj
        BetaSTDjj[wbeta1==0]=max(abs(BetaSTDj))+1
        temo=rank(-abs(BetaSTDjj), ties.method="min")

        temo=data.frame(temo[which(temo<=numj)], which(temo<=numj))
        temo=temo[order(temo[, 1]), ]

        x1j=x1[temid, ,drop=F]
        cvTrimLogC(Betaj[temo[, 2]], numj, numi2, temo[, 2]-1, x1j, y[temid], Nf[i], threshP)
      }, Mthread = Mthread, Mcore = Mcore)

      do.call(rbind, fold_rows)
    }


    #####  Cross-validation estimates  #####
    # ypred=matrix(0,nrow=N0,ncol=nlambdai)

    fold_results <- parallel_task_apply(fold_index, function(i) {
      temid=(foldid==i)
      cv_row=rep(NA_real_, nlambdai)

      if (any(y[!temid]==0) & any(y[!temid]==1)) {
        x1i=x1[!temid, ,drop=F]; x1j=x1[temid, ,drop=F]; p1i=p1
        wbeta1i=wbeta1; wbetaii=wbetai

        indexi=which(apply(x1i,2,sd)==0)[-1]
        if (length(indexi)>0) {
          x1i=x1i[,-indexi,drop=F]; x1j=x1j[,-indexi,drop=F]
          p1i=p1i-length(indexi)
          wbeta1i=wbeta1i[-indexi]; wbetaii=wbetaii[-indexi]
        }

        outi_i=switch(penalty,
                      "Net"=cvNetLogC(x1i, y[!temid], alpha, lambdai, nlambdai, wbeta1i, wbetaii, W$Omega, W$loc, W$nadj, p1i, N0i[i],thresh, maxit, x1j, y[temid], Nf[i], threshP),
                      cvEnetLogC(x1i, y[!temid], alpha, lambdai, nlambdai, wbeta1i, wbetaii, p1i, N0i[i],thresh, maxit, x1j, y[temid], Nf[i], threshP)
        )

        if (length(indexi)>0) {
          betai=matrix(0,nrow=p1,ncol=outi_i$nlambda)
          betai[-indexi,]=outi_i$Beta
          outi_i$Beta=betai

          betai=matrix(0,nrow=p1,ncol=outi_i$nlambda)
          betai[-indexi,]=outi_i$BetaSTD
          outi_i$BetaSTD=betai
        }

        if (outi_i$nlambda > 0) {
          cv_row[1:outi_i$nlambda]=2*(0-outi_i$LLF)*Nf[i]
        }
      } else {
        outi_i=list()
        outi_i$Beta=matrix(0,nrow=p1,ncol=nlambdai)
        outi_i$BetaSTD=matrix(0,nrow=p1,ncol=nlambdai)
      }

      list(outi=outi_i, cv=cv_row)
    }, Mthread = Mthread, Mcore = Mcore)
    outi=lapply(fold_results, function(res) res$outi)
    cvRSS=do.call(rbind, lapply(fold_results, function(res) res$cv))

    cvRSS=cvRSS[, 1:nlambdai,drop=F]
    cvraw=cvRSS/weighti; nfoldi=apply(!is.na(cvraw), 2, sum); #rm(cvRSS) #
    cvm=apply(cvraw, 2, weighted.mean, w=weighti, na.rm=TRUE)
    cvse=sqrt(apply(sweep(cvraw, 2, cvm, "-")^2, 2, weighted.mean, w=weighti, na.rm=TRUE)/(nfoldi-1))

    indexi=which.min(cvm)
    indexij=which(cvm<=(cvm[indexi]+cvse[indexi]))[1]
    temi=rep("", nlambdai)
    temi[indexi]="*";#temi[indexij]=ifelse(temi[indexij]=="", "*", "***")
    #temCV=data.frame(lambda=lambdai, cvm=cvm, cvse=cvse, nzero=out$nzero, index=temi,stringsAsFactors=FALSE)
    temCV=data.frame(lambda=lambdai, pDev=(out$ll0-out$LL)/out$ll0, cvm=cvm, cvse=cvse, nzero=pmax(out$nzero-1,0), index=temi, stringsAsFactors=FALSE)

    if (!inzero) {
      rm(outi)
      if (!keep.beta) {
        # lambda.1se=lambdai[indexij]
        return(list(Beta=out$Beta[, indexi], fit=temCV, lambda.min=lambdai[indexi], penalty=penalty, adaptive=adaptive, flag=out$flag))
      } else {
        return(list(Beta=out$Beta, fit=temCV, lambda.min=lambdai[indexi], penalty=penalty, adaptive=adaptive, flag=out$flag))
      }
    }

    #####  Cross-validation hard threshold  #####
    il0=indexi; cvm=list(); cv.min=rep(NA, nlambdai)
    repeat {
      numi=out$nzero[il0]
      Betai=matrix(sapply(outi, function(x){x$Beta[, il0,drop=F]}), nrow=p1)
      BetaSTDi=matrix(sapply(outi, function(x){x$BetaSTD[, il0,drop=F]}), nrow=p1)

      Betao=apply(Betai!=0, 2, sum)
      numi2=pmax(min(max(Betao), numi),1)

      cvRSS=compute_trim_log_rows(Betai, BetaSTDi, numi, numi2)

      cvraw=cvRSS/weighti; nfoldi=apply(!is.na(cvraw), 2, sum); #rm(cvRSS) #
      cvm[[il0]]=apply(cvraw, 2, weighted.mean, w=weighti, na.rm=TRUE)
      temi=cvm[[il0]]

      cv.min[il0]=ifelse(length(temi)>sum(wbeta1==0),min(temi[-c(1:sum(wbeta1==0))]),temi[sum(wbeta1==0)])

      il1=c(il0-1, il0+1)
      for (j in 1:2) {
        if (il1[j]>=1 & il1[j]<=nlambdai) {
          if (is.na(cv.min[il1[j]])) {
            numi=out$nzero[il1[j]]
            Betai=matrix(sapply(outi, function(x){x$Beta[, il1[j],drop=F]}), nrow=p1)
            BetaSTDi=matrix(sapply(outi, function(x){x$BetaSTD[, il1[j],drop=F]}), nrow=p1)

            Betao=apply(Betai!=0, 2, sum)
            numi2=pmax(min(max(Betao), numi),1)

            cvRSS=compute_trim_log_rows(Betai, BetaSTDi, numi, numi2)

            cvraw=cvRSS/weighti;nfoldi=apply(!is.na(cvraw), 2, sum); #rm(cvRSS)
            cvm[[il1[j]]]=apply(cvraw, 2, weighted.mean, w=weighti, na.rm=TRUE)
            temi=cvm[[il1[j]]]

            cv.min[il1[j]]=ifelse(length(temi)>sum(wbeta1==0),min(temi[-c(1:sum(wbeta1==0))]),temi[sum(wbeta1==0)])
          }
        } else {
          break
        }
      }
      # if (il1[j]==1 | il1[j]==nlambdai)
      #   break
      if (il0==which.min(cv.min)) {
        break
      } else {
        il0=which.min(cv.min)
      }
    }
    index0=which.min(cv.min)

    Beta0=out$Beta[,index0]
    BetaSTD0=out$BetaSTD[,index0]

    temi=cvm[[index0]]
    cuti=ifelse(length(temi)>sum(wbeta1==0),which.min(temi[-c(1:sum(wbeta1==0))])+sum(wbeta1==0),sum(wbeta1==0))
    # cuti=which.min(cvm[[index0]])

    Beta0j=out$BetaSTD[,index0]
    Beta0j[which(wbeta1==0)]=max(abs(Beta0j))+1

    Beta0[abs(Beta0j)<=sort(abs(Beta0j),TRUE)[cuti+1]]=0
    BetaSTD0[abs(Beta0j)<=sort(abs(Beta0j),TRUE)[cuti+1]]=0

    temCV0=data.frame(lambda=lambdai[index0],cvm=cv.min[index0],nzero=cuti-1)

    if (!keep.beta) {
      if (!isd) {
        return(list(Beta=out$Beta[, indexi], Beta0=Beta0, fit=temCV, fit0=temCV0, lambda.min=lambdai[indexi], lambda.opt=lambdai[index0], penalty=penalty, adaptive=adaptive, flag=out$flag))
      } else {
        return(list(Beta=out$BetaSTD[, indexi], Beta0=BetaSTD0, fit=temCV, fit0=temCV0, lambda.min=lambdai[indexi], lambda.opt=lambdai[index0], penalty=penalty, adaptive=adaptive, flag=out$flag))
      }
    } else {
      if (!isd) {
        return(list(Beta=out$Beta, Beta0=Beta0, fit=temCV, fit0=temCV0, lambda.min=lambdai[indexi], lambda.opt=lambdai[index0], penalty=penalty, adaptive=adaptive, flag=out$flag))
      } else {
        return(list(Beta=out$BetaSTD, Beta0=BetaSTD0, fit=temCV, fit0=temCV0, lambda.min=lambdai[indexi], lambda.opt=lambdai[index0], penalty=penalty, adaptive=adaptive, flag=out$flag))
      }
    }

  } # folder
}


