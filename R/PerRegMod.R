

####""""""""""the following function gives A Kronecker product B
A_x_B=function(A,B){
  l_A=dim(A)[1]
  C_A=dim(A)[2]
  l_B=dim(B)[1]
  C_B=dim(B)[2]
  C=matrix(0,l_A*l_B,C_A*C_B)
  for (i in 1:l_A) {
    for (j in 1:C_A) {
      for (k in 1:l_B) {
        for (s in 1:C_B) {
          C[(i-1)*l_B+k ,(j-1)*C_B+s]=A[i,j]*B[k,s]

        }

      }

    }

  }
  return(C)

}
#######""""""""""""""""""LSE

####""""""" x a list of all independent variables
####""""""" y a dependent variable
###"""""""" s the perioid of the regression
###"""""""" p the number of independent variables

#########"""""""""""""""estimation of mu_s and beta_s
LSE_Reg_per=function(x,y,s){
  p=length(x)
  n=length(y)
  if(n%%s==0){
    m=n/s
  }else{message(paste(n," does not divide ", s) ) }
  Y=matrix(0,1,n)
  i=1
  while (i<=s) {
    for (r in 0:m-1) {
      a=(i-1)*m+1
      b=i*m
      j=a:b
      Y[1,j[r+1]]=y[r*s+i]
    }
    i=i+1
  }
  Y=t(Y)
  II=diag(rep(1,s) )
  l=matrix(1,m,1)
  D=A_x_B(II,l)

  ###"""""""""x
  t=list()
  length(t)=p
  for (i in 1:p) {
    t[[i]]=matrix(0,n,s+1)
  }
  for (k in 1:p) {
    i=2
    while (i<=s+1) {
      for (r in 0:m-1) {
        a=(i-2)*m+1
        b=(i-1)*m
        j=a:b
        t[[k]][j[r+1],i]=x[[k]][r*s+i-1]
      }
      i=i+1
    }
  }
  for (k in 1:p) {
    t[[k]]=t[[k]][,-1]
  }

  ###################""""""""""
  if(p==1){
  X=cbind(D,t[[1]])
  }else{
    X=cbind(D,t[[1]])
    for (k in 2:p) {
      X=cbind(X,t[[k]])
    }
  }
  beta=solve(t(X)%*%X)%*%t(X)%*%Y
 return(list(beta=beta,X=X,Y=Y))
}
#####""""""""""""""""""estimation of sigma_s
sd_estimation_for_each_s=function(x,y,s,beta_hat){
  p=length(x)
  n=length(y)
  if(n%%s==0){
    m=n/s
  }else{message(paste(n," does not divide ", s) ) }

  ##########"""""""sum of x
  xx=rep(0,n)
  for (k in 1:p) {
    for (i in 1:s) {
      q=seq(i,n,s)
      xx[q]=xx[q]+beta_hat[i+s*k,1]*x[[k]][q]
    }
  }

  ###################""""
  z=rep(0,n)
  for (k in 1:p) {
    for (i in 1:s) {
      q=seq(i,n,s)
      z[q]=y[q]-xx[q]
    }
  }
  for (i in 1:s) {
    q=seq(i,n,s)
    z[q]=z[q]-beta_hat[i,1]
  }
  sigma=matrix(0,s,1)
  for (i in 1:s) {
    q=seq(i,n,s)
    sigma[i,1]=sum(z[q]^2 )/(m-p-1)
  }
  return(sd=sqrt(sigma) )
}

###########""""""""""""""""""""""""
lm_per=function(x,y,s){
  p=length(x)
  n=length(y)
  m=n/s
  beta_hat=LSE_Reg_per(x,y,s=s)$beta
  sd=sd_estimation_for_each_s(x,y,s,beta_hat)

  xx=rep(0,n)
  for (k in 1:p) {
    for (i in 1:s) {
      q=seq(i,n,s)
      xx[q]=xx[q]+beta_hat[i+s*k,1]*x[[k]][q]
    }
  }
  ######""
  residuals=rep(0,n)
  for (k in 1:p) {
    for (i in 1:s) {
      q=seq(i,n,s)
      residuals[q]=y[q]-xx[q]
    }
  }
  for (i in 1:s) {
    q=seq(i,n,s)
    residuals[q]=residuals[q]-beta_hat[i,1]
  }
  if(p==1){
    ceof=cbind("intercept","x[[ 1 ]]")
  }else{
    ceof=cbind("intercept","x[[ 1 ]]")
    for (k in 2:p) {
      ceof=cbind(ceof,paste("[[",k,"]]") )
    }
  }
  beta_hat_mat=matrix(0,s,p+1)
  for (i in 1:s) {
    for (k in 1:p) {
      beta_hat_mat[i,k+1]=beta_hat[s+(k-1)*s+i ,1]
    }

  }
  for (i in 1:s) {
    beta_hat_mat[i,1]=beta_hat[i ,1]
  }
  newdata=as.data.frame(beta_hat_mat)
  names(newdata)=ceof
  if(s==2){
    s_name=rbind("s= 1","s= 2")
  }else{
    s_name=rbind("s= 1","s= 2")
    for (i in 3:s) {
      s_name=rbind(s_name,paste("s=",i) )
    }
  }
  row.names(newdata)=s_name
  newdata=data.frame(newdata,sd )
  RMSE=sqrt(mean(residuals^2) )
  #####"""""""""""R^2
  y_hat=rep(0,n)
  for (i in 1:s) {
    q=seq(i,n,s)
    y_hat[q]=beta_hat[i]
  }
  for (j in 1:p) {
    for (i in 1:s) {
      q=seq(i,n,s)
      y_hat[q]=y_hat[q]+beta_hat[i+j*s]*x[[j]][q]
    }
  }
  sst=rep(0,s)
  for (i in 1:s) {
    q=seq(i,n,s)
    sst[i]=sum((y[q]-mean(y[q]))^2)
  }
  sst_f=sum(sst)
  R_squared_for_per=1-sum( (y-y_hat)^2)/sst_f
  R_squared_adjusted_for_per=1-(sum( (y-y_hat)^2)/(m-p-1) )/(sst_f/(m-1) )
  verbose=TRUE
  if(verbose)cat("Residuals:\n")
  if(verbose)print(quantile(residuals, probs = c(0, 0.25, 0.5, 0.75, 1)),digits = 5)
  if(verbose)cat("\n")
  if(verbose)cat("Coefficients:\n")
  if(verbose)print(newdata, digits = 5)
  if(verbose)cat("\n")
  if(verbose)cat(paste("R-squared: ", round(R_squared_for_per,6),"Adjusted R-squared: ",round(R_squared_adjusted_for_per,6) ))
  if(verbose)cat("\n")
  if(verbose)cat("Root mean square error: ", RMSE)

}

########""""""""""
#n=20
#s=2
#x1=rnorm(n,0,1)
#x2=rnorm(n,0,2)
#x3=rnorm(n,0,3)
#x4=rnorm(n,0,2.7)
#y=rnorm(n,0,2.5)
#x=list(x1,x2,x3,x4)
#beta_hat=LSE_R_p_p_multiple_p_3(x,y,s=s)$beta
#sd_estimation_for_each_s(x,y,s,beta_hat)
#lm_per_AE(x,y,s)

########""""""""""""""""""""""Adaptive method
########""""""""" define central sequence, Delta, for periodic regression model
#x=list(x1,x2)
DELTA=function(x,phi,s,e,sigma){
  p=length(x)
  n=length(phi)
  m=n/s
  z=rep(0,n)
  for (i in 1:s) {
    for (r in 0:m-1) {
      z[s*r+i]=e[s*r+i]/sigma[i]
    }

  }
  delta=rep(0,p*s+s+s)
  a=rep(0,s)
  for (i in 1:s) {
    a[i]=0
    for (r in 1:m) {
      a[i]=a[i]+(1/sigma[i])*phi[s*(r-1)+i]
    }

  }
  for (i in 1:s) {
    delta[i]=a[i]
  }

  for (i in 1:s) {
    delta[i+s]=0
    for (r in 1:m) {
      delta[i+s]=delta[i+s]+1/(2*sigma[i]^2)*(phi[s*(r-1)+i]*z[s*(r-1)+i]-1)
    }
  }

  ######
  M=list()
  length(M)=p
  for (i in 1:p) {
    M[[i]]=rep(0,s)
  }
  k=list()
  length(k)=p
  for (i in 1:p) {
    k[[i]]=rep(0,s)
  }

  for (j in 1:p) {
    i=1
    while(i<=s) {
      p_1=(i-1)*m+1
      p_2=i*m
      M[[j]][i]=(1/m)*sum(x[[j]][p_1:p_2]^2)
      k[[j]][i]=(M[[j]][i])^(-1/2)
      i=i+1
    }
  }

  for (j in 1:p) {
    for (i in 1:s) {
      for (r in 1:m) {
        delta[2*s+(j-1)*s+i]=delta[2*s+(j-1)*s+i]+(1/sigma[i])*phi[s*(r-1)+i]*x[[j]][s*(r-1)+i]*k[[j]][i]
      }

    }
  }

  delta=(n^(-1/2))*delta
  return(delta)
}
#DELTA(x,phi_n(z) ,s=s,e=z,sigma = sigmma)
#####################################"""""""""""""gamma
GAMMA=function(x,phi,s,z,sigma){
  p=length(x)
  n=length(x)
  m=n/s
  ########## information of Fisher I
  I_s=rep(0,s)
  for (i in 1:s) {
    I_s[i]=0
    for (r in 1:m) {
      I_s[i]=I_s[i]+(1/m)*phi[s*(r-1)+i]^2

    }

  }
  I=sum(I_s/s)

  ################## estimate J
  J_s=rep(0,s)
  for (i in 1:s) {
    J_s[i]=0
    for (r in 1:m) {
      J_s[i]=J_s[i]+(1/m)*phi[s*(r-1)+i]^2*z[s*(r-1)+i]^2

    }

  }
  J=sum(J_s/s)

  ########################### estimate K
  kk_s=rep(0,s)
  for (i in 1:s) {
    kk_s[i]=0
    for (r in 1:m) {
      kk_s[i]=kk_s[i]+(1/m)*phi[s*(r-1)+i]^2*z[s*(r-1)+i]

    }

  }
  kk=sum(kk_s/s)
  #################################""" Gamma
  gamma=matrix(0,p*s+2*s,p*s+2*s)
  for (i in 1:s) {
    gamma[i,i]=(I/s)*(1/sigma[i]^2)
  }
  i=1
  while(i<=s) {
    gamma[i+s,i+s]=(J-1)/(4*s*(sigma[i])^4)
    i=i+1
  }
  ###########"
  for (j in 1:p) {
    i=1
    while(i<=s) {
      gamma[2*s+(j-1)*s+i,2*s+(j-1)*s+i]=(I/s)*(1/sigma[i]^2)
      i=i+1
    }
  }
  #####################"""""""
  for (i in 1:s) {
    gamma[i,i+s]=kk/(2*s*sigma[i]^3)
    gamma[i+s,i]=kk/(2*s*sigma[i]^3)

  }
  ########"""""""""""""""""""""""""
  M=list()
  length(M)=p
  for (j in 1:p) {
    M[[j]]=rep(0,s)
  }
  k=list()
  length(k)=p
  for (j in 1:p) {
    k[[j]]=rep(0,s)
  }
  ########"""""""""""""""
  for (j in 1:p) {
    i=1
    while(i<=s) {
      p_1=(i-1)*m+1
      p_2=i*m
      M[[j]][i]=(1/m)*sum(x[[j]][p_1:p_2]^2)
      k[[j]][i]=(M[[j]][i])^(-1/2)
      i=i+1
    }
  }
  ######"""""""""
  for (j in 1:p) {
    for (l in 1:p) {
      if(j!=l){
        for (i in 1:s) {
          p_1=(i-1)*m+1
          p_2=i*m
          gamma[2*s+(j-1)*s+i,2*s+(l-1)*s+i]=k[[j]][i]*k[[l]][i]*sum(x[[j]][p_1:p_2]*x[[l]][p_1:p_2]) /(m*sigma[i]^2)

        }

      }

    }

  }




  return(gamma)
}
#GAMMA_multi(x,phi_n(z),s=s,z,sigma = sigmma)
######""""""""""""""""""""""""phi_n for kernel gaussain

#phi_n<- function(z) { return(z)}
phi_n<- function(x) {
  b_n=0.002
  l=length(x)
  uu=0
  for (i in 1:l) {
    uu=uu+(x-x[i])*exp(-(x-x[i])^2/(2*b_n^2) )
  }
  vv=0
  for (i in 1:l) {
    vv=vv+exp(-(x-x[i])^2/(2*b_n^2) )
  }

  u=1/b_n^2*uu/vv
  return(u)}


estimate_para_adaptive_method=function(n,s,y,x){
  p=length(x)
  m=n/s
  ########""" verify assumption B
  x_bar=list()
  length(x_bar)=p
  for (k in 1:p) {
    x_bar[[k]]=matrix(0,1,s)

  }
  #########"""""""""""
  for (k in 1:p) {
    i=1
    while(i<=s){
      a=(m-1)*s+i
      j=seq(i,a,s)
      for (r in 1:m) {
        x_bar[[k]][1,i]=x_bar[[k]][1,i]+x[[k]][j[r]]
      }
      #############""""""""""""""""""""""""""""""""
      i=i+1
    }
  }
  for (k in 1:p) {
    for (i in 1:s) {
      x_bar[[k]][1,i]=(1/m)*x_bar[[k]][1,i]
    }

  }

  ###########estimating sigma
  beta_l=LSE_Reg_per(x,y,s=s)$beta
  ##########"""""""sum of beat*x
  xx=rep(0,n)
  for (k in 1:p) {
    for (i in 1:s) {
      q=seq(i,n,s)
      xx[q]=xx[q]+beta_l[i+s*k,1]*x[[k]][q]
    }
  }


  eee=rep(0,n)
  sigmma=rep(0,s)
  for (i in 1:s) {
    for (r in 1:m) {
      eee[i+(r-1)*s] =(y[i+(r-1)*s]-beta_l[i,1]-xx[i+(r-1)*s])^2
    }
  }
  for (i in 1:s) {
    q=seq(i,n,s)
    sigmma[i]=1/(m-p-1)*sum(eee[q] )

  }

  #######"""""""""""
  ee2=rep(0,n)
  for (i in 1:s) {
    q=seq(i,n,s)
    ee2[q] =y[q]-beta_l[i,1]-xx[q]
  }
  z=rep(0,n)
  for (i in 1:s) {
    for (r in 0:(m-1) ) {
      z[s*r+i]=ee2[s*r+i]/sigmma[i]
    }

  }

  DALTA1=DELTA(x,phi_n(z) ,s=s,e=ee2,sigma = sigmma)
  GAMMA1=GAMMA(x,phi_n(z),s=s,z,sigma = sigmma)
  ##########""""

  BTEA2=matrix(0,p*s+2*s,1)
  for (i in 1:s) {
    BTEA2[i,1]=beta_l[i,1]
    BTEA2[i+s,1]=sigmma[i]
  }

  for (k in 1:p) {
    for (i in 1:s) {
      BTEA2[i+(k-1)*s +2*s,1]=beta_l[i+k*s,1]
    }
  }

  #""""""""""
  result=BTEA2+1/sqrt(n)*solve(GAMMA1)%*%DALTA1
  #return(list(beta_ad=result,beta_lse=BTEA2))
  return(list(beta_ad=result))

}
#para_ad=estimate_para_adaptive_method(n=length(y),s,y,x)$beta_ad
#########""""""""""""""""""""""""

lm_per_AE=function(x,y,s){
  p=length(x)
  n=length(y)
  m=n/s
  beta_hat=estimate_para_adaptive_method(n,s,y,x)$beta_ad
  sd_p=matrix(0,s,1)
  for (i in 1:s) {
    sd_p[i,1]=beta_hat[s+i,1]
  }

  xx=rep(0,n)
  for (k in 1:p) {
    for (i in 1:s) {
      q=seq(i,n,s)
      xx[q]=xx[q]+beta_hat[i+s*(k+1),1]*x[[k]][q]
    }
  }
  ######""
  residuals=rep(0,n)
  for (k in 1:p) {
    for (i in 1:s) {
      q=seq(i,n,s)
      residuals[q]=y[q]-xx[q]
    }
  }
  for (i in 1:s) {
    q=seq(i,n,s)
    residuals[q]=residuals[q]-beta_hat[i,1]
  }
  if(p==1){
    ceof=cbind("intercept","x[[ 1 ]]")
  }else{
    ceof=cbind("intercept","x[[ 1 ]]")
    for (k in 2:p) {
      ceof=cbind(ceof,paste("[[",k,"]]") )
    }
  }
  beta_hat_mat=matrix(0,s,p+1)
  for (i in 1:s) {
    for (k in 1:p) {
      beta_hat_mat[i,k+1]=beta_hat[s+k*s+i ,1]
    }

  }
  for (i in 1:s) {
    beta_hat_mat[i,1]=beta_hat[i ,1]
  }
  newdata=as.data.frame(beta_hat_mat)
  names(newdata)=ceof
  if(s==2){
    s_name=rbind("s= 1","s= 2")
  }else{
    s_name=rbind("s= 1","s= 2")
    for (i in 3:s) {
      s_name=rbind(s_name,paste("s=",i) )
    }
  }
  row.names(newdata)=s_name
  newdata=data.frame(newdata,SD=sd_p )
  RMSE=sqrt(mean(residuals^2) )
  #####"""""""""""
  #y_hat=rep(0,n)
  #for (i in 1:s) {
  #  q=seq(i,n,s)
  # y_hat[q]=beta_hat[i]
  #}
  #for (j in 1:p) {
  # for (i in 1:s) {
  #  q=seq(i,n,s)
  # y_hat[q]=y_hat[q]+beta_hat[i+(j+1)*s]*x[[j]][q]
  #}
  #}
  #sst=rep(0,s)
  #for (i in 1:s) {
  # q=seq(i,n,s)
  #sst[i]=sum((y[q]-mean(y[q]))^2)
  #}
  #sst_f=sum(sst)
  #R_squared_for_per=1-sum( (y-y_hat)^2)/sst_f
  #R_squared_adjusted_for_per=1-(sum( (y-y_hat)^2)/(m-p-1) )/(sst_f/(m-1) )

  verbose=TRUE

  #return(list(cat("Residuals:\n"),print(quantile(residuals, probs = c(0, 0.25, 0.5, 0.75, 1)), digits = 5),cat("\n"),cat("Coefficients:\n"),print(newdata, digits = 5),cat("\n"),cat("Root mean square error: ", RMSE)))
  if(verbose)cat("Residuals:\n")
  if(verbose)print(quantile(residuals, probs = c(0, 0.25, 0.5, 0.75, 1)),digits = 5)
  if(verbose)cat("\n")
  if(verbose)cat("Coefficients:\n")
  if(verbose)print(newdata, digits = 5)
  if(verbose)cat("\n")
  #if(verbose)cat(paste("R-squared: ", round(R_squared_for_per,6),"Adjusted R-squared: ",round(R_squared_adjusted_for_per,6) ))
  #if(verbose)cat("\n")
  if(verbose)cat("Root mean square error: ", RMSE)

}

#######""""""""""""""""""#########"""""""""""""""""""""""" statistic test

###x is a list
pseudo_gaussian_test=function(x,z,s){
  phi_gaussian=function(x){x}
  psi_gaussian=function(x){x^2-1}
  p=length(x)
  n=length(x[[1]])
  m=n/s
  ###""""""""""""""""""""""""""""""""""""
  X_bar=list()
  length(X_bar)=p
  for (i in 1:p) {
    X_bar[[i]]=matrix(0,1,s)

  }
  for (i in 1:p) {
    for (j in 1:s) {
      q=seq(j,n,s)
      X_bar[[i]][j] =mean(x[[i]][q] )

    }

  }
  for (i in 1:p) {
    for (j in 1:s) {
      q=seq(j,n,s)
      x[[i]][q]=x[[i]][q]-X_bar[[i]][j]
    }
  }

  ###""""""""""""""""""""""""""
  delta__2=matrix(0,s-1,1)
  for ( i in 1:(s-1) ) {
    q=seq(i,n,s)
    d=seq(s,n,s)
    delta__2[i,1]=n^(-0.5)*sum(phi_gaussian(z[q])-phi_gaussian(z[d]) )

  }

  delta__4=matrix(0,s-1,1)
  for ( i in 1:(s-1) ) {
    q=seq(i,n,s)
    d=seq(s,n,s)
    delta__4[i,1]=n^(-0.5)/2*sum(psi_gaussian(z[q])-psi_gaussian(z[d]) )
  }
  ##########""""""""""""""""""Delta 4 sans sigma
  kk=MM=list()
  for (i in 1:s ) {
    kk[[i]]=MM[[i]]=matrix(0,p,p)

  }
  for (i in 1:s) {
    q=seq(i,n,s)
    for (j in 1:p) {
      for (l in 1:p) {
        MM[[i]][j,l]=mean(x[[j]][q]*x[[l]][q])
      }
    }

  }
  for (i in 1:s) {
    kk[[i]]=solve(sqrtm(MM[[i]]))##### M1^-0.5
  }



  delta__6=matrix(0,p*(s-1) ,1)

  for (j in 1:(s-1) ) {
    q=seq(j,n,s)
    d=seq(s,n,s)
    for (i in 1:p) {
      delta__6[(j-1)*p+i,1]=n^(-0.5)*sum(x[[i]][q]*phi_gaussian(z[q] )-x[[i]][d]*phi_gaussian(z[d] )  )
    }
  }

  ##########"""""""""
  delta__5=matrix(0,p,1)
  for ( i in 1:p) {
    delta__5[i,1]=n^(-0.5)*sum(x[[i]]*phi_gaussian(z) )

  }


  ##############""""""""""""""""""""""""""""""
  I=1
  A1=matrix(0,s-1,s-1)
  for (i in 1:(s-1) ) {
    for (j in 1:(s-1) ) {
      if(i==j){
        A1[i,j]=2
      }else{A1[i,j]=1  }

    }

  }
  Gamma___22=I/s*A1

  ##################"""""""""
  m_4=mean(z^4)
  m_3=mean(z^3)



  Gamma___44=(m_4-1)/(4*s)*A1
  Gamma___24=m_3/(2*s)*A1

  ######""""""""""""""""""""""""""""""##############""""""""""""""""""""""""""""""Gamma 44 sans sigma
  Gamma___66=matrix(0,p*(s-1),p*(s-1) )
  for (j in 1:(s-1) ) {
    q=seq(j,n,s)
    d=seq(s,n,s)
    for (i in 1:p) {
      for (l in 1:p) {
        if(i==l){
          Gamma___66[(j-1)*p+i,(j-1)*p+l]=sum(x[[i]][q]^2+x[[i]][d]^2 )

        }else{
          Gamma___66[(j-1)*p+i,(j-1)*p+l]=sum(x[[i]][q]*x[[l]][q] +x[[i]][d]*x[[l]][d] )

        }
      }

    }

  }

  for (j in 1:(s-1)) {
    for (k in 1:(s-1) ) {
      d=seq(s,n,s)
      for (i in 1:p) {
        for (l in 1:p) {
          if(j!=k){
            Gamma___66[(j-1)*p+i,(k-1)*p+l]=sum(x[[i]][d]*x[[l]][d])
          }

        }

      }

    }

  }


  Gamma___66=I/n*Gamma___66
  #####################""""""""""""
  Gamma___55=matrix(0,p,p)

  for (i in 1:p) {
    for (k in 1:p) {
      if(i==k){
        Gamma___55[i,k]=sum(x[[i]]^2)
      }else{
        Gamma___55[i,k]=sum(x[[i]]*x[[k]])

      }

    }

  }
  Gamma___55=I/n*Gamma___55
  #####################""""""""""""""""""""""
  Gamma___56=matrix(0,p,p*(s-1))
  for (j in 1:(s-1) ) {
    q=seq(j,n,s)
    d=seq(s,n,s)
    for (i in 1:p) {
      for (l in 1:p) {
        Gamma___56[i,(j-1)*p+l]=sum(x[[i]][q]*x[[l]][q]-x[[i]][d]*x[[l]][d])

      }
    }

  }

  Gamma___56=I/n*Gamma___56


  #####################""""""""""""""""""""""Gamma global
  Gamma_auxi=matrix(0,s-1,p*(s-1))
  l1=cbind(Gamma___22,Gamma___24)
  l2=cbind(Gamma___24,Gamma___44)
  #l3=cbind(t(Gamma_auxi),t(Gamma_auxi),Gamma___66)

  Gamma_final=rbind(l1,l2)

  #######""""""" delta global

  delta=cbind(t(delta__2),t(delta__4))

  #""""""""""
  delta2=delta__6-t(Gamma___56)%*%solve(Gamma___55)%*%delta__5

  ########"""" Gamma**


  gamma2=Gamma___66-t(Gamma___56)%*%solve(Gamma___55)%*%Gamma___56



  T=delta%*%solve(Gamma_final)%*%t(delta)+t(delta2)%*%solve(gamma2)%*%delta2

  return(T)
}

check_periodicity=function(x,y,s){
  p=length(x)
  n=length(x[[1]] )
  Y=matrix(0,n,1)
  for (i in 1:n) {
    Y[i,1]=y[i]
  }


  X=matrix(0,n,p+1)
  for (i in 1:n) {
    X[i,1]=1
  }
  for (k in 1:p) {
    for (i in 1:n) {
      X[i,k+1]=x[[k]][i]
    }
  }
  esti_beta=solve(t(X)%*%X)%*%t(X)%*%Y
  #####"""""" sum of beta*x
  xx=rep(0,n)
  for (i in 1:n) {
    for (k in 1:p) {
      xx[i]=xx[i]+esti_beta[k+1,1]*x[[k]][i]
    }
  }
  xx=xx+esti_beta[1,1]
  ###########""""
  z=y-xx
  test=pseudo_gaussian_test(x,z,s)
  df=(s-1)*(p+2)
  p_value=1-pchisq(as.numeric(test) , df=df )
  verbose=TRUE
  if(verbose)cat("Chi-squared test for detecting periodicity of coefficients")
  if(verbose)cat("\n")
  if(verbose)cat(paste("Chi-squared statistic:",round(test,5)," on df=",df,", p-value: ",round(p_value,7) ))
}

