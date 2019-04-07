test = function(){  # use by copy & paste

    source("hg.R")

    source("pfaffian-FB.r")
    source("Gradient_mle.r")
    alpha0 = c(0.1, 0.2, 0.3, 0.4, 0, 0, 0, 0)
    alpha1 = c(1, 2, 3, 4, sqrt(2), sqrt(3), sqrt(4), sqrt(5))
    G0 = G.FB(alpha0, method="lpinv", withvol=TRUE)
    G1 = G.FB(alpha1, method="MC", withvol=TRUE)
    G1a = hg(alpha0, G0, alpha1, dG.fun.FB)
    alpha2 = c(0.1, 0.2, 0.3, 0.4, 0.1, 0.2, 0.3, 0.4)
    H0 = G.FB(alpha2, method="MC", ns=c(2,3,1,4), withvol=TRUE)
    H1 = G.FB(alpha1, method="MC", ns=c(2,3,1,4), withvol=TRUE)
    H1a = hg(alpha2, H0, alpha1, dG.fun.FB, fn.params=list(ns=c(2,3,1,4)))

    #source("/Users/alfredkume/Dropbox/Tomonari/R/150429ver/pfaffian-FB.R")
    #source("pfaffian-FB.R")
    source("hg.R")
    alpha0 = c(0.1, 0.3, 0.1, 0.1)
    alpha1 = c(10, 30, 10, 10)
    ns = c(1,2)  # multiplicities
    H0 = G.FB(alpha0, ns=ns, withvol=TRUE)   # power series
    H1 = G.FB(alpha1, ns=ns, method="MC", withvol=TRUE)   # MC
    H1a = hg(alpha0, H0, alpha1, dG.fun.FB, fn.params=list(ns=ns))$G   # HG
    n=length(alpha0[1:4])+1
    v=(2*pi^(n/2)/gamma(n/2))

    saddleaprox.FB.revised(sort(c(0,alpha0[1:4]))+1,dub=1,order=2)*exp(1)/v
    hg.Bingham(-alpha0[1:4])$G

    n=2
    l=sort(abs(rnorm(n,sd=2)),decreasing = FALSE)*3
    alpha0=c(l,rep(0,length(l)))

    G0=G.FB(c(alpha0[1:n],rep(0,n)), method="lpinv", withvol=TRUE);G0[1]
    #G0=G.FB(c(0,alpha0[1:4],rep(0,5)), method="lpinv", withvol=TRUE);G0
    SPA=saddleaprox.FB.revised(sort(alpha0[1:n])+1,M=rep(1,n)*0,dub=1,order=2)*exp(1)
    G0[1]/SPA
    p=3
    alpha1=c(sort(abs(rnorm(p,sd=10))),abs(rnorm(p,sd=1)))
    ####
    # alpha1=c(0.01 * (1:p), abs(rnorm(p)))
    # alpha1[1:5]=alpha1[1:p]-alpha1[1]
    # alpha1[-(1:5)]=sort(alpha1[-(1:p)])
    ####

    #alpha1[3]=alpha1[2]+.1
    #alpha1[2]=alpha1[1]+.002
    alpha1[1:p]=alpha1[1:p]-alpha1[1]

    alpha1[-(1:p)]=sort(alpha1[-(1:p)])

    alpha0=c(alpha1[1:p],rep(0,p))
    G0=G.FB(alpha0, method="lpinv", withvol=TRUE);G0[1]
    SPA0=saddleaprox.FB.revised(sort(alpha1[1:p])+1,M=alpha1[(p+1):(2*p)]*0,dub=1,order=3)*exp(1)
    SPA0/G0[1]

    G1a = hg(alpha0, G0, alpha1, dG.fun.FB);G1a$G[1]# note that the mean parameter here is twice that in saddlepoint approx
    SPA1=saddleaprox.FB.revised(sort(alpha1[1:p])+1,M=alpha1[(p+1):(2*p)]/2,dub=1,order=3)*exp(1)
    SPA1/G1a$G[1]
}

#Test:
# source("hg.R")
# source("Laplace2.r")
# L = 1:4
# gamma = sqrt(2:5)
# hg.FB(L, gamma)
# G.FB(c(L,gamma), method="MC", withvol=TRUE)

#may break down for high dimension:
# for(i in 2:30) show(G.FB(c(0.01*1:i,rep(0,i)),method="lpinv",withvol=TRUE)[1])

# for(i in 2:20) show(integlpinvR(c(0.01*1:i))-G.FB(c(0.01*1:i,rep(0,i)),method="lpinv",withvol=TRUE)[1])

# for(i in 2:20) show(G.FB(c(0.01*1:i,rep(0,i)),method="MC",withvol=TRUE)[1])

integlpinvR=function(L)
    #for even p
{
    Y=0
    p1=floor(length(L)/2)
    for (i in 1:p1)
    {
        f=function(x) {return(Partlast(L,i,x))}
        f=Vectorize(f)
        #print(integrateR(f,0,1/sqrt(2)))
        #v=seq(0,1/sqrt(2),length=100);plot(v,f(v));readline()
        Y=(-1)^(i+1)*integrateR(f,0,1/sqrt(2),ord=20)$v+Y
    }
    Y
    #f=function(x) {return(Part(L,1,x))}
    if((length(L)-2*p1)>0)
    {
        f1=function(x) {return(Partinf(L,x))}
        f1=Vectorize(f1)
        Y=Y+(-1)^p1*integrate(f1,0,Inf)$v
    }
    Y=Y*2*pi^(length(L)/2-1)
    Y
}
Loglikelihood=function(theta,gamma, O=diag(rep(1,length(gamma))),A,B,n=1,method="hg")
    #Calculates the loglikelohood where A and B are observed second and first moments A=sum X_i #X_i^t and B=sum X_i.
    #O is the orthogonal component of the covariance matrix
{
    p=length(theta)
    ths=sort(theta,index=T)
    #ths$ix=1:p
    alpha1=c(theta[ths$ix],gamma[ths$ix])
    alpha0=c(theta[ths$ix],gamma*0)
    #

    #G0[2]=0

    #if(method=="MC")
    #G0 = G.FB(alpha0, method="MC", withvol=TRUE)



    if(method=="MC")
        hgout=G.FB(alpha1, method="MC", withvol=TRUE)

    if(method=="hgold")
    {
        G0 = G.FB(alpha0, method="lpinv", withvol=TRUE)
        #hgout=hg.log(alpha0, G0, alpha1, dG.fun.FB)$G
        #hgout[2]=0#to make sure that the smallest theta is 0
        #hgout[2:(p+1)]=hgout[2:(p+1)][sort(ths$ix,index=T)$ix]
        #hgout[p+(2:(p+1))]=hgout[p+(2:(p+1))][sort(ths$ix,index=T)$ix]

        hgout=hgm.FB.2(alpha1)

    }

    if(method=="hg")
    {
        hgout=hgm.FB.2(alpha1)
    }


    if(method=="SPA")
        hgout=c(saddleaprox.FB.revised(sort(alpha1[1:p])+1,M=alpha1[(p+1):(2*p)]/2,dub=1,order=3)*exp(1))

    if(method=="lpinv")
        hgout=c(lpinv(alpha1))

    nc=hgout[1]
    l=-n*log(nc)-sum(diag(A%*%t(O)%*%diag(theta)%*%O+gamma%*%t(B)%*%t(O)))

    grad=-n*hgout[-1]/nc+c(-diag(O%*%A%*%t(O)),O%*%B)
    # grad[1]=0#seems to be very important
    #grad[1:p]=-grad[1:p]

    y=list(log=l,grad=grad)
    y
}

# gg=1:4;th=1:4;A=gg%*%t(th);B=th;
# X=rnorm(40);
# X=array(X,c(4,10))
# for (i in 1:10)
# X[i,]=X[i,]/norm(X[i,])

# A=X%*%t(X); B=apply(X,1,mean)*10


# Loglikelihood(th,gg,A=A,B=B,n=60)

#dth=0.001;Loglikelihood(th,gg,A=A,B=B,n=60)$log-Loglikelihood(th+c(dth,0,0),gg,A=A,B=B,n=60)$log/dth

Grad_mle_update=function(th,gg,A=A,B=B,O=diag(rep(1,length(gg))),n=1)
{
    p=length(th)
    current=Loglikelihood(th,gg,A=A,B=B,O=O,n=n)
    delta=1
    current$grad[1]=0
    newth=th+delta*current$grad[1:(p)]
    #delta=ifelse(min(newth)<0,-min(th/current$grad[1:(p)]),delta)
    newth=th+delta*current$grad[1:(p)]
    newth[1]=0# might need that for guaranteeing convergence
    newgg=gg+delta*current$grad[(p+1):(2*p)]
    new=Loglikelihood(newth,newgg,A=A,B=B,n=n,O=O)
    while(new$log<current$log&delta>1e-3)
    {
        delta=delta/2
        newth[2:p]=th[2:p]+delta*current$grad[2:(p)]
        newgg=gg+delta*current$grad[(p+1):(2*p)]
        new=Loglikelihood(newth,newgg,A=A,B=B,n=n,O=O)
        print(delta,new$log<current$log)
        #readline()
    }
    y=list(th=newth,gg=newgg,new,delt=delta)
    y
}

Grad_mle_update_optim=function(th,gg,A=A,B=B,O=diag(rep(1,length(gg))),n=1)
{
    p=length(th)
    current=Loglikelihood(th,gg,A=A,B=B,O=O,n=n)
    delta=1
    newth=th+delta*current$grad[1:(p)]
    newth[1]=0# might need that for guaranteeing convergence
    newgg=gg+delta*current$grad[(p+1):(2*p)]
    f=function(dl)
    {
        newth[2:p]=th[2:p]+dl*current$grad[2:(p)]
        newgg=gg+dl*current$grad[(p+1):(2*p)]
        x=-Loglikelihood(newth,newgg,A=A,B=B,n=n,O=O)$log
        x
    }
    dmax=min(c(2,abs(1/current$grad[2:p]*th[-1])))

    delt_opt=optimize(f,c(0,dmax))
    delta=delt_opt$min
    newth[2:p]=th[2:p]+delta*current$grad[2:(p)]
    newgg=gg+delta*current$grad[(p+1):(2*p)]
    new=Loglikelihood(newth,newgg,A=A,B=B,n=n,O=O)
    y=list(th=newth,gg=newgg,new,delt=delta)
    y
}

Grad_mle_update_optim_sqrt=function(th,gg,A=A,B=B,O=diag(rep(1,length(gg))),n=1,dmax=1)
    #since th>0 and gg>0 we can write th=x^2 and gg=y^2 and optimize wrt x and y
{
    p=length(th)
    current=Loglikelihood(th,gg,A=A,B=B,O=O,n=n)
    x=sqrt(th);y=sqrt(gg)
    delta=10
    newx=x*(1+delta*current$grad[1:(p)])
    newy=y*(1+delta*current$grad[(p+1):(2*p)])
    # newth=th+delta*current$grad[1:(p)]
    newx[1]=0;current$grad[1]=0# might need that for guaranteeing convergence
    # newgg=gg+delta*current$grad[(p+1):(2*p)]
    f=function(dl)
    {
        x=sqrt(th);y=sqrt(gg)
        newx=x*(1+dl*current$grad[1:(p)])
        newy=y*(1+dl*current$grad[(p+1):(2*p)])

        z=-Loglikelihood(newx^2,newy^2,A=A,B=B,n=n,O=O)$log
        z
    }
    # dmax=min(c(2,abs(1/current$grad[2:p]*th[-1])))
    # dmax=10
    delt_opt=optimize(f,c(0,dmax))
    delta=delt_opt$min
    newx=x*(1+delta*current$grad[1:(p)])
    newy=y*(1+delta*current$grad[(p+1):(2*p)])

    newth=newx^2
    newgg=newy^2
    new=Loglikelihood(newth,newgg,A=A,B=B,n=n,O=O)
    y=list(th=newth,gg=newgg,new,delt=delta)
    y
}

library(Matrix)
Grad_mle_orth_update=function(th,gg,A=A,B=B,O=diag(rep(1,length(gg))),n=1,tol=1e-3)
{
    AA=diag(th)%*%O%*%A%*%t(O)-O%*%A%*%t(O)%*%diag(th)+gg%*%t(B)%*%t(O)
    vhat=AA-t(AA)
    cc=-sum(diag(A%*%t(O)%*%diag(th)%*%O+gg%*%t(B)%*%t(O)))
    current=Loglikelihood(th,gg,A=A,B=B,n=n,O=O)
    delta=1*sign(sum(diag(AA%*%vhat)))
    newO=as.matrix(expm(vhat*delta))%*%O
    newcc=-sum(diag(A%*%t(newO)%*%diag(th)%*%newO+gg%*%t(B)%*%t(newO)))
    new=Loglikelihood(th,gg,A=A,B=B,n=n,O=newO)
    # print(c(delta,current$log-new$log,cc-newcc))


    f0=function(t0)
    {
        Y=-Loglikelihood(th,gg,A=A,B=B,n=n,O=as.matrix(expm(vhat*t0))%*%O)$log
        Y
    }
    out0=optimise(f0,100*c(-1,1))
    newO1=as.matrix(expm(vhat*out0$min))%*%O
    print(c(out0$min,sum(vhat^2)))
    new1=Loglikelihood(th,gg,A=A,B=B,n=n,O=newO1)
    if(new$log>new1$log)
    {
        newO=newO1
    }
    else{
        print("too small vhat?")
    }
    # while(new$log-current$log<tol)
    # {
    # delta=delta/2
    # newO=as.matrix(expm(vhat*delta))%*%O
    # new=Loglikelihood(th,gg,A=A,B=B,n=n,O=newO)
    # print(c(delta,current$log,new$log))
    # }
    # print(c(delta,current$log,new$log))
    y=list(O=newO,new1,AA=AA)
    y
}
# newO=diag(rep(1,length(gg)))
# O=Grad_mle_orth_update(th,gg,A=A,B=B,O=newO,n=60)$O

# O=Grad_mle_orth_update(th,gg,A=A,B=B,O=O,n=60)$O

optimal_orth=function(th,gg,A=A,B=B,O=diag(rep(1,length(gg))),n=1,tol=1e-2)
    #provides the  optimal orthogonal matrix using the gradient method as in the paper.
{
    GG=Grad_mle_orth_update(th,gg,A=A,B=B,O=O,n=n)
    newO=GG$O
    a=norm(O-newO)
    k=1
    while(a>tol)
    {
        Oold=newO
        GG=Grad_mle_orth_update(th,gg,A=A,B=B,O=Oold,n=n)
        newO=GG$O
        a=norm(Oold-newO)
        a=norm(GG$AA-t(GG$AA))
        print(a)
        #readline()
        k=k+1
    }
    list(O=GG$O,AA=GG$AA)
}
# #OO=optimal_orth(th,gg,A=A,B=B,O=diag(rep(1,length(th))),n=60,tol=1e-2)
# gg%*%t(B)%*%t(OO$O)
# OO$O%*%A%*%t(OO$O)
# th
optimisation=function(th,gg,A=A,B=B,n=1,O=diag(rep(1,length(gg))),orth="no",tol=1e-3,iters=200)
    #This optimises the likelihood for a fixed orthogonal matrix O
{
    ll=Grad_mle_update(th,gg,A=A,B=B,O=O,n=1)
    llnew=Grad_mle_update(ll$th,ll$gg,A=A,B=B,n=n)
    a=norm(ll$th-llnew$th)
    b=norm(ll$gg-llnew$gg)
    print(c(a,b,llnew[[3]]$log))
    newO=O
    iter=1
    c=norm(ll$""$grad-llnew$""$grad)
    while(max(c,a+b)>tol&&iter<iters)
    {
        a=norm(ll$th-llnew$th)
        b=norm(ll$gg-llnew$gg)
        c=norm(llnew$""$grad)
        ll=llnew
        if(orth=="yes")
        {
            newO=Grad_mle_orth_update(ll$th,ll$gg,A=A,B=B,O=newO,n=n)$O
            llnew=Grad_mle_update_optim_sqrt(ll$th,ll$gg,A=A,B=B,O=newO,n=n)
        }
        else
        {
            llnew=Grad_mle_update_optim_sqrt(ll$th,ll$gg,A=A,B=B,n=n)
            # print(c(a,b,llnew[[3]]$log,llnew[[3]]$log-ll[[3]]$log))
        }
        print(c("iteration"))
        print(c(iter,a,b,c,llnew[[3]]$log,llnew[[3]]$log-ll[[3]]$log))
        iter=iter+1
    }
    llnew$O=newO
    llnew
}
# th[1]=0
# Grad_mle_update=Grad_mle_update_optim
# outo=optimisation(th,gg,A=A,B=B,tol=1e-3)
# outo=optimisation(outo$th,outo$gg,A=A,B=B,tol=1e-3)
# out1o=optimisation(th,gg,A=A,B=B,n=1,tol=1e-2,orth="yes")
# out2o=optimisation(out1o$th,out1o$gg,A=A,B=B,n=60,O=out1o$O,tol=1e-2,orth="yes")
# out1o=optimisation(th,gg+3,A=A,B=B,n=60,tol=1e-3)
###################
saddleaprox.FB.revised<-function(L,M=L*0,dub=3, order=3)
    #calculates the normalising constant of Fisher-Bingham distribution in pre-shape space there are three methods as described
    #L is the vector of positive values
    #M is the vecor of mu<-i's
    #dub is the number at which each entry of L is doubled
{
    #L<-sort(rep(t(L),dub))
    L<-rep(t(L),dub)
    a<-prod(L/pi)^(-1/2)
    Y<-0

    KM<-function(t)
    {
        #Y<-sum(-1/2*log(1-t/L)+L*M^2/(1-t/L)-1*L*M^2)
        Y<-sum(-1/2*log(1-t/L)+M^2/(1-t/L)/L)
        Y
    }
    KM1<-function(t)
    {
        Y<-sum(1/2*1/(L-t)+M^2/(L-t)^2)
        Y
    }
    KM2<-function(t)
    {
        Y<-sum(1/2*1/(L-t)^2+2*M^2/(L-t)^3)
        Y
    }
    KM3<-function(t)
    {
        Y<-sum(1/(L-t)^3+6*M^2/(L-t)^4)
        Y
    }
    KM4<-function(t)
    {
        Y<-sum(3/(L-t)^4+24*M^2/(L-t)^5)
        Y
    }
    sol<-function(KM1,y)
    {
        #Y<-optimize(loc<-function(t) {abs(KM1(t)-1)},c(min(L)-length(L)/(2*y)-.Machine$double.eps^1,min(L)), tol = .Machine$double.eps^2)$min
        Y<-optimize(loc<-function(t) {abs(KM1(t)-1)},c(min(L)-length(L)/(4*y)-sqrt(length(L)^2/4+length(L)*max(L)^2*max(M)^2),min(L)), tol = .Machine$double.eps^2)$min
        Y
    }
    ##
    that<-sol(KM1,1)
    #cat("that=",that,exp(KM(that))*a,KM2(that));
    Y<-2*a/sqrt(2*pi*KM2(that))*exp(KM(that)-that)
    #Y<-2*a/sqrt(2*pi*KM2(that))*exp(KM(that)-that)*exp(sum(L*M))
    #Y<-1/sqrt(2*pi*KM2(that))*exp(KM(that)-that)
    if (order==3)
    {
        rho3sq<-KM3(that)^2/KM2(that)^3
        rho4<-KM4(that)/KM2(that)^(4/2)
        Rhat<-3/24 *rho4-5/24*rho3sq
        Y<-Y*exp(Rhat)
    }
    if (order==2)
    {
        rho3sq<-KM3(that)^2/KM2(that)^3
        rho4<-KM4(that)/KM2(that)^(4/2)
        Rhat<-3/24 *rho4-5/24*rho3sq
        Y<-Y*(1+Rhat)
    }
    Y
}

########
#some support functions here
supfun=1
if(supfun==1)
{
    convert=function(p)
        # vector of length 6 which contains the entries of a 3x3 symmetric matrix(output)
        # with offdiagonal entries values halved
        # this function is used to use the data from Seietal paper as in lower triagonal entries. A as in our paper can be adopted as 2convert(S1..S6)-its(diagonal)
    {
        A=array(0,c(3,3))
        A[1,1]=p[1]
        A[1,2]=p[2]/2; A[2,1]=A[1,2]
        A[1,3]=p[3]/2;A[3,1]=A[1,3]
        A[2,2]=p[4]
        A[2,3]=p[5]/2; A[3,2]=A[2,3]
        A[3,3]=p[6]
        A
    }


    my.ode.hg.log=function(tau, GL, params)
        # this is applying the log transform and working with the  log(G) vector whose derivative is dGL below
    {
        th = params$th
        dG.fun = params$dG.fun
        v = params$v
        fn.params = params$fn.params
        th = th + tau * v
        G=GL*exp(GL[1])
        G[1]=exp(GL[1])
        if(is.null(fn.params))  dG = dG.fun(th, G)
        else dG = dG.fun(th, G, fn.params)
        dG=v %*% dG
        dGL=dG/G[1]-G*dG[1]/G[1]^2
        dGL[1]=dGL[1]+G[1]*dG[1]/G[1]^2
        G.rhs=dGL
        list(G.rhs)
    }
    hg.log=function(th0, G0, th1, dG.fun, times=seq(0,1,length=101), fn.params=NULL, show.trace=FALSE)
    {
        params = list(th = th0, dG.fun = dG.fun.FB, v = th1-th0, fn.params = fn.params)
        G0L=G0/G0[1]; G0L[1]=log(G0[1]);

        rk.res = rk(G0L, times, my.ode.hg.log, params,rtol=.Machine$double.eps^0.25, atol =.Machine$double.eps^0.25)
        if(show.trace) trace = rk.res
        else trace=NULL
        list(GL = rk.res[nrow(rk.res), 1+(1:length(G0))], trace = trace)

        GL = rk.res[nrow(rk.res), 1+(1:length(G0))]
        G=GL
        G[1]=exp(GL[1])
        G[-1]=G[-1]*G[1]
        list(G=G,trace=trace)
    }

    hg=function(th0, G0, th1, dG.fun, times=seq(0,1,length=101), fn.params=NULL, show.trace=FALSE)
    {
        params = list(th = th0, dG.fun = dG.fun, v = th1-th0, fn.params = fn.params)
        rk.res = rk(G0, times, my.ode.hg, params,rtol=.Machine$double.eps^0.25, atol =.Machine$double.eps^0.25)
        if(show.trace) trace = rk.res
        else trace=NULL
        list(G = rk.res[nrow(rk.res), 1+(1:length(G0))], trace = trace)
    }


    # Loglikelihood(th,gg,A=A,B=B,n=1)
    Numeical_optim=function(th,gg,A=A,B=B,n=1)
    {
        p=length(gg)
        vec=c(th,gg)
        fgrad=function(x)
        {
            p=(length(x)+1)/2
            Y=Loglikelihood(c(0,x[1:(p-1)]),x[p-1+(1:p)],A=A/n,B=B/n,n=1)$grad[-1]
            as.vector(Y)
        }
        flik=function(x)
        {
            p=(length(x)+1)/2
            Y=-Loglikelihood(c(0,x[1:(p-1)]),x[p-1+(1:p)],A=A/n,B=B/n,n=1)$log
            Y
        }
        f=function(x)
        {
            fgrad=function(y)
            {
                Loglikelihood(c(0,y[1:(p-1)]),y[p-1+(1:p)],A=A/n,B=B/n,n=1)$grad[-1]
            }
            Y=-Loglikelihood(c(x[1:(p)]),x[p+(1:p)],A=A/n,B=B/n,O=O,n=1)$log
            # attr(Y,"gradient")=Vectorize(fgrad(x))
            Y
        }
        # optim(vec,flik,fgrad, method = "L-BFGS-B",control = list(maxit = 20))
        out=nlm(f,vec,iterlim = 20)
        thhat=c(0,out$est[1:(p-1)])
        gghat=out$est[-(1:(p-1))]
        z=list(th=thhat, gg=gghat,Loglikelihood(thhat,gghat,A=A/n,B=B/n,n=1))
        z
    }

    numer.opt.F.B=function(inth,ingg,A=A,B=B,O=diag(rep(1,dim(A)[1])),Bing=FALSE,iter=199,print=0,Kent="",meth="hg")
    {
        # This function optimises numerically the likelihood fnction when p=3. when the data are given as in the paper of Seietal
        #Bingh=TRUE reduces to the Bingham distribution
        # If Kent is numeric and Bing=FALSE(default) we can run the Kent distribution with the Kent=numeric eigenvector present in the model.
        likelihood.total.opt=function(par,p=3)
        {
            n=1
            #th=abs(c(0,par[1:2]))
            th=c(0,par[1:2])
            gg=par[3:5]
            if (is.numeric(Kent)==TRUE)
            {
                gg[-Kent]=gg[-Kent]*0
                th[2]=-th[3]
                th=th+abs(min(th))
            }

            if(Bing==TRUE)
                gg=abs(par[3:5])*0


            vhat=array(0,c(p,p))
            vhat[2,1]=par[6]
            vhat[3,1]=par[7]
            vhat[3,2]=par[8]
            VV=vhat-t(vhat)
            O=as.matrix(expm(VV))
            # cc=sign(par[3:5]);Y=-Loglikelihood(th,gg*cc,A=A1/n,B=B/n,O=diag(cc)%*%O,n=1)$log
            Y=-Loglikelihood(th,gg,A=A/n,B=B/n,O=O,n=1,method=meth)$log
            #Y=-Loglikelihood(th,gg,A=A/n,B=B/n,O=O,n=1)$log
            Y
        }
        orth=function(x)
            # generates the orth matrix given its variable in exp proj mapping
        {
            vhat=array(0,c(3,3))
            vhat[2,1]=x[1]
            vhat[3,1]=x[2]
            vhat[3,2]=x[3]
            VV=vhat-t(vhat)
            O=as.matrix(expm(VV))
            O
        }
        ort=logm(O)
        par=c(inth[2:3]-inth[1],ingg,ort[2,1],ort[3,1],ort[3,2])
        opt=nlm(likelihood.total.opt,par,iterlim=iter,print.level=print)
        if (is.numeric(Kent)==TRUE)
        { opt$est[1]=-opt$est[2]
        opt$est[4:5]=opt$est[4:5]*0
        }

        opt
        z=list(th=c(0,opt$est[1:2]),gg=opt$est[3:5],O=orth(opt$est[6:8]),out=opt)
        z
    }

}
#########

