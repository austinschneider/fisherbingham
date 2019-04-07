##########################################
# hgm-FB.R
# ver. 2017-04-27
##########################################
#
# Description
#
#  Computing the normalizing constant for Fisher-Bingham distribution on S_{p-1} via
#  the holonomic gradient method (HGM). Degenerate cases are available.
#  The initial value is obtained by power series or Monte Carlo.
#
#
# Usage
#
#  hgm.FB(alpha, ns=rep(1,length(alpha)/2), alpha0=NULL, G0=NULL, withvol=TRUE)
#  hgm.FB.2(alpha, ns=rep(1,length(alpha)/2), withvol=TRUE)
#
#  # supplementary functions:
#
#  dG.fun.FB(alpha, G, fn.params=NULL)  # the Pfaffian
#  G.FB(alpha, method="power", withvol=TRUE)  # initial value by power series
#  G.FB(alpha, method="MC", withvol=TRUE)  # initial value by Monte Carlo
#  SPA(alpha, ns=rep(1,length(alpha)/2), withvol=TRUE)  # saddle point approximation
#
#
# Arguments
#
#  alpha: (2*p)-dim vector consisting of theta and xi (= gamma in the paper)
#  ns: p-dim vector denoting the multiplicities (default is the all-one vector)
#
#
# Details
#
#  The function hgm.FB computes a (2*p+1)-dim vector G, that consists of
#  the normalising constant of the Fisher-Bingham distribution and its derivatives.
#  The path connecting the starting point alpha0 and the target alpha is the straight line,
#  where alpha0 is proportional to alpha with small norm if it is not specified.
#  If the argument withvol is TRUE, the total base measure is the volume of the sphere.
#  (1 if FALSE)
#
#  The function hgm.FB.2 computes the G vector by HGM along a suitable path,
#  which is proposed by Koyama et al. (2013) and Koyama and Takemura (2015).
#
#  The function dG.fun.FB computes the right hand side of the Pfaffian eq.
#
#  The function G.FB computes the initial value of G by power series or Monte Carlo.
#
#
# Examples (see more in test())
#
#  source("hgm-FB.R") # this file
#
#  alpha1 = c(1, 2, 3, 4, 7, 6, 8, 5)  # parameter
#  ns1 = c(2,2,1,3)  # multiplicities
#  hgm.FB(alpha1, ns=ns1)  # HGM
#  hgm.FB.2(alpha1, ns=ns1)  # HGM via a suitable path (the same result)
#
#
#  For Saddle point approximation the follwing lines need to run
#  SPA(alpha1, ns=ns1)
#
#
#  G.FB(alpha1, ns=ns1, method="MC")  # direct computation by Monte Carlo
#  # G.FB(alpha1, ns=ns1, method="power")  # power series (do not work)
#
#  alpha2 = c(0.067, 0.294, 0.311, 0.405, 0.663, 0.664, 0.667, 0.712, 0.784, 0.933, 0.070, 0.321, 0.288, 0.367, 0.051, 0.338, 0.798, 0.968, 0.506, 0.590)  # near singular
#  ns2 = rep(1, 10)
#  hgm.FB(alpha2, ns=ns2)  # HGM (does not work)
#  hgm.FB.2(alpha2, ns=ns2)  # HGM via a suitable path
#  G.FB(alpha2, ns=ns2, method="MC", withvol=TRUE)  # direct computation by Monte Carlo
#  SPA(alpha2, ns2)  # saddle point approximation
#
##########################################

library(deSolve)

### Pfaffian for Fisher-Bingham distribution
dG.fun.FB = function(alpha, G, fn.params=list(ns=rep(1,length(alpha)/2), s=1)){
	# s: scale
	d = length(alpha)
	p = d/2
	r = length(G) # r should be equal to d+1
	th = alpha[1:p]
	xi = alpha[p+(1:p)]
	gam = xi^2/4
	Gth = G[1+(1:p)]
	Gxi = G[1+p+(1:p)]
	dG = array(0,c(d, r))
	ns = fn.params$ns
	s = fn.params$s

	dG[,1] = G[1+(1:d)]
	# derivative of dC/dth[i] with respect to th[j]
	for(i in 1:p){
		dG[i,1+i] = -s*Gth[i]
		for(j in 1:p){
			if(j != i){
				a1 = - ( ns[j]/2/(th[j]-th[i]) + xi[j]^2/4/(th[j]-th[i])^2 )
				a2 = - ( ns[i]/2/(th[i]-th[j]) + xi[i]^2/4/(th[i]-th[j])^2 )
				a3 = - ( ns[j]*xi[i]/4/(th[j]-th[i])^2 + xi[i]*xi[j]^2/4/(th[j]-th[i])^3 )
				a4 = - ( ns[i]*xi[j]/4/(th[i]-th[j])^2 + xi[i]^2*xi[j]/4/(th[i]-th[j])^3 )
				dG[j,1+i] = a1*Gth[i] + a2*Gth[j] + a3*Gxi[i] + a4*Gxi[j]
				dG[i,1+i] = dG[i,1+i] - dG[j,1+i]
			}
		}
	}
	# derivative of dC/dxi[i] with respect to th[j]
	# derivative of dC/dth[j] with respect to xi[i]
	for(i in 1:p){
		dG[i,1+p+i] = -s*Gxi[i]
		for(j in 1:p){
			if(j != i){
				b2 = xi[i]/2/(th[i]-th[j])
				b3 = - ( ns[j]/2/(th[j]-th[i]) + xi[j]^2/4/(th[j]-th[i])^2 )
				b4 = xi[i]*xi[j]/4/(th[i]-th[j])^2
				dG[j,1+p+i] = b2*Gth[j] + b3*Gxi[i] + b4*Gxi[j]
				dG[p+i,1+j] = dG[j,1+p+i]
				dG[i,1+p+i] = dG[i,1+p+i] - dG[j,1+p+i]
			}
		}
		dG[p+i,1+i] = dG[i,1+p+i]
	}
	# derivative of dC/dxi[i] with respect to xi[j]
	for(i in 1:p){
		for(j in 1:p){
			if(j != i){
				c3 = xi[j]/2/(th[j]-th[i])
				c4 = -xi[i]/2/(th[j]-th[i])
				dG[p+j,1+p+i] = c3*Gxi[i] + c4*Gxi[j]
			}
		}
		if(ns[i] == 1 || abs(xi[i]) < 1e-10)  dG[p+i,1+p+i] = -Gth[i]  # cheat
		#  	if(ns[i] == 1)  dG[p+i,1+p+i] = -Gth[i]
		else  dG[p+i,1+p+i] = -Gth[i] - (ns[i]-1)/xi[i]*Gxi[i]  # singular if xi[i] = 0
	}
	dG
}

### Initial value for Pfaffian system (by power series expansion)
# parametrization: exp(th*x^2 + xi*x)
C.FB.power = function(alpha, v=rep(0,length(alpha)), d=rep(1,length(alpha)/2), Ctol=1e-6, alphatol=1e-10, Nmax=ceiling(10^(10/length(alpha)))){
	p = length(alpha)/2  # not (p-1)
	alpha.sum = sum(abs(alpha[1:p])) + sum(abs(alpha[(p+1):(2*p)]))
	N0 = max(ceiling(alpha.sum), 1)
	if(N0 > Nmax) stop("too large value: alpha.")
	for(N in N0:Nmax){ # ToDo: modify
		logep = N*log(alpha.sum) - lfactorial(N) + log((N+1)/(N+1-alpha.sum))
		if(logep < log(Ctol)) break
	}
	if(N > Nmax) stop("too large value: alpha.")
	f = function(k){
		kn = length(k)
		if(kn == 2*p){
			k1 = k[1:p]
			k2 = k[(p+1):(2*p)]
			v1 = v[1:p]
			v2 = v[(p+1):(2*p)]
			if(any((k2+v2) %% 2 == 1)) return(0)
			w = which(k > 0)
			a = prod( alpha[w]^(k[w]) )
			b1 = sum( - lfactorial(k1) - lfactorial(k2) + lgamma(k1 + v1 + (k2 + v2)/2 + (d/2)) + lgamma((k2+v2)/2 + (1/2)) - lgamma((k2+v2)/2 + (d/2)))
			b2 = - lgamma(sum(k1 + v1 + (k2 + v2)/2 + (d/2)))
			b3 = lgamma(sum(d)/2) - p*lgamma(1/2)
			return( a * exp(b1+b2+b3) )
		}else{ # recursive part
			knxt = kn + 1
			if(abs(alpha[knxt]) < alphatol) return( f(c(k,0)) ) # speed-up
			a = 0
			imax = N - sum(k)
			for(i in 0:imax) a = a + f(c(k,i))
			return(a)
		}
	}
	ret = f(c())
    ret
}

G.FB.power = function(alpha, ns=rep(1,length(alpha)/2), Ctol=1e-6, alphatol=1e-10){
	# not restricted to Bingham
	p = length(alpha)/2
	th = alpha[1:p]
	xi = alpha[(p+1):(2*p)]
	C = C.FB.power(c(-th, xi), d=ns, Ctol=Ctol, alphatol=alphatol)  # note: parameterization
	dC = numeric(2*p)
	for(i in 1:p){
		e = rep(0,2*p); e[i] = 1
		dC[i] = -C.FB.power(c(-th, xi), v=e, d=ns, Ctol=Ctol, alphatol=alphatol)  # note: parametrerization
		e = rep(0,2*p); e[p+i] = 1
		dC[p+i] = C.FB.power(c(-th, xi), v=e, d=ns, Ctol=Ctol, alphatol=alphatol)  # note: parametrerization
	}
	c(C,dC)
}


### Initical value for Pfaffian system (by Monte Carlo)
rsphere = function(N, p){
	x = matrix(rnorm(N*p), N, p)
	r = sqrt(apply(x^2, 1, sum))
	x / r
}

G.FB.MC = function(alpha, ns=rep(1,length(alpha)/2), N=1e6, t=NULL){
	p = length(alpha)/2
	th = alpha[1:p]
	xi = alpha[p+(1:p)]
	G = numeric(2*p+1)
	if(is.null(t)) t = rsphere(N, sum(ns))
	idx = c(0, cumsum(ns))
	t2 = t1 = matrix(0, N, p)
	for(i in 1:p){
		t2[,i] = rowSums(t[,(idx[i]+1):idx[i+1],drop=FALSE]^2)
		t1[,i] = t[,idx[i]+1]
	}
	itrd = exp(t2 %*% (-th) + t1 %*% xi)  # integrand
	G[1] = mean(itrd)
	for(i in 1:p){
		G[1+i] = -mean(itrd * t2[,i])
		G[1+p+i] = mean(itrd * t1[,i])
	}
	return( G )
}


### Initial value (wrapper)
G.FB = function(alpha, ns=rep(1,length(alpha)/2), method="power", withvol=TRUE){
	dsum = sum(ns)
	v0 = 2 * pi^(dsum/2)/gamma(dsum/2)
	v = ifelse(withvol, v0, 1)
	if(method == "power"){
		return( v * G.FB.power(alpha, ns=ns) )
	}
	if(method == "MC"){
		return( v * G.FB.MC(alpha, ns=ns) )
	}
	stop("method not found")
}

### ODE
my.ode.hg = function(tau, G, params){
	th = params$th
	dG.fun = params$dG.fun
	v = params$v
	fn.params = params$fn.params
	th = th + tau * v
	if(is.null(fn.params))  dG = dG.fun(th, G)
	else dG = dG.fun(th, G, fn.params)
	G.rhs = v %*% dG
	list(G.rhs)
}

### hg main
hg = function(th0, G0, th1, dG.fun, times=seq(0,1,length=101), fn.params=NULL, show.trace=FALSE){
	params = list(th = th0, dG.fun = dG.fun, v = th1-th0, fn.params = fn.params)
	rk.res = rk(G0, times, my.ode.hg, params)
	if(show.trace) trace = rk.res
	else trace=NULL
	list(G = rk.res[nrow(rk.res), 1+(1:length(G0))], trace = trace)
}

### ODE (for square-root transformation)
my.ode.hg.mod = function(tau, G, params){
	th0 = params$th0
	th1 = params$th1
	dG.fun = params$dG.fun

	p = length(th0) / 2
	th = v = numeric(2*p)
	w = 1:p
	th[w] = th0[w] + tau * (th1[w] - th0[w])
	v[w] = th1[w] - th0[w]
	w = (p+1):(2*p)
	th[w] = sign(th1[w]) * sqrt(th0[w]^2 + tau * (th1[w]^2 - th0[w]^2))
	wz = (th[w] == 0)
	w = w[!wz]
	v[w] = (th1[w]^2 - th0[w]^2) / 2 / th[w]

	fn.params = params$fn.params
	if(is.null(fn.params))  dG = dG.fun(th, G)
	else dG = dG.fun(th, G, fn.params)
	G.rhs = v %*% dG
	list(G.rhs)
}

### hg main (for square-root transformation)
hg.mod = function(th0, G0, th1, dG.fun, times=seq(0,1,length=101), fn.params=NULL, show.trace=FALSE){
	params = list(th0 = th0, th1 = th1, dG.fun = dG.fun, fn.params = fn.params)
	rk.res = rk(G0, times, my.ode.hg.mod, params)
	if(show.trace) trace = rk.res
	else trace=NULL
	list(G = rk.res[nrow(rk.res), 1+(1:length(G0))], trace = trace)
}

### evaluating FB normalising constant by HGM
hgm.FB = function(alpha, ns=rep(1,length(alpha)/2), alpha0=NULL, G0=NULL, withvol=TRUE){
	p = length(alpha) / 2
	r = sum(abs(alpha))
	N = max(r, 1)^2 * 10
	if(is.null(alpha0))  alpha0 = alpha[1:(2*p)] / N
	if(is.null(G0)) G0 = G.FB(alpha0, ns=ns, method="power", withvol=withvol)
	fn.params = list(ns=ns, s=1)
    #print("G0")
    #print(G0)
    #print("alpha0")
    #print(alpha0)
    #print("ns")
    #print(ns)
	as.vector(hg(alpha0, G0, alpha, dG.fun.FB, fn.params=fn.params)$G)
}

### evaluating FB normalising constant by HGM (via square-root transformation)
hgm.FB.2 = function(alpha, ns=rep(1,length(alpha)/2), withvol=TRUE){
	p = length(alpha) / 2
	r = sum(abs(alpha))
	N = max(r, 1)^2 * 10
	alpha0 = c(alpha[1:p] / N, alpha[(p+1):(2*p)] / sqrt(N))
	G0 = G.FB(alpha0, ns=ns, method="power", withvol=withvol)
	fn.params = list(ns=ns, s=1)
	res = hg.mod(alpha0, G0, alpha, dG.fun.FB, fn.params=fn.params)$G
	as.vector(res)
}

####SPA calculation
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

    #cat("L=", L, "\n")
    #cat("M=", M, "\n")
    #cat("a=", a, "\n")
    #cat("Y=", Y, "\n")

	KM<-function(t)
	{
		#Y<-sum(-1/2*log(1-t/L)+L*M^2/(1-t/L)-1*L*M^2)
		Y<-sum(-1/2*log(1-t/L)+M^2/(1-t/L)/L)
        #cat("t=", t, ", KM=", Y, "\n")
		Y
	}
	KM1<-function(t)
	{
		Y<-sum(1/2*1/(L-t)+M^2/(L-t)^2)
        #cat("t=", t, ", KM1=", Y, "\n")
		Y
	}
	KM2<-function(t)
	{
		Y<-sum(1/2*1/(L-t)^2+2*M^2/(L-t)^3)
        #cat("t=", t, ", KM2=", Y, "\n")
		Y
	}
	KM3<-function(t)
	{
		Y<-sum(1/(L-t)^3+6*M^2/(L-t)^4)
        #cat("t=", t, ", KM3=", Y, "\n")
		Y
	}
	KM4<-function(t)
	{
		Y<-sum(3/(L-t)^4+24*M^2/(L-t)^5)
        #cat("t=", t, ", KM4=", Y, "\n")
		Y
	}
	sol<-function(KM1,y)
	{
		#Y<-optimize(loc<-function(t) {abs(KM1(t)-1)},c(min(L)-length(L)/(2*y)-.Machine$double.eps^1,min(L)), tol = .Machine$double.eps^2)$min
        #cat("min(L)=", min(L), "\n")
        #cat("-length(L)/(4*y)=", -length(L)/(4*y),"\n")
        #cat("length(L)^2/4=", length(L)^2/4, "\n")
        #cat("length(L)*max(L)^2*max(M)^2=", length(L)*max(L)^2*max(M)^2, "\n")
        #fmin<-min(L)-length(L)/(4*y)-sqrt(length(L)^2/4+length(L)*max(L)^2*max(M)^2)
        #fmax<-min(L)
        #cat("fmin=", fmin, "\n")
        #cat("fmax=", fmax, "\n")
        #loc<-function(t) {abs(KM1(t)-1)}
        #cat("f(fmin)=", loc(fmin), "\n")
        #cat("f(fmax)=", loc(fmax), "\n")
        #cat("####")
		Y<-optimize(loc<-function(t) {abs(KM1(t)-1)},c(min(L)-length(L)/(4*y)-sqrt(length(L)^2/4+length(L)*max(L)^2*max(M)^2),min(L)), tol = .Machine$double.eps^2)$min
		Y
	}
	##
	that<-sol(KM1,1)
    #cat("that=", that, "\n")
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

SPA=function(alpha,ns=rep(1,length(alpha)/2), withvol=TRUE)
{
	p=length(alpha)/2
	theta=rep(alpha[1:(p)],ns)
    #cat("aaahhh=", alpha[-(1:(p))], "\n")
    #cat("aaahhh=", rep(alpha[-(1:(p))], ns), "\n")
	mu=rep(alpha[-(1:(p))]/sqrt(ns),ns)/2

	dsum = sum(ns)
	v0 = 2 * pi^(dsum/2)/gamma(dsum/2)
	coef = ifelse(withvol, 1, 1/v0)

	#cat("theta=", theta, "\n")
	#cat("mu=", mu, "\n")
	#cat("coef=", coef, "\n")

	return(saddleaprox.FB.revised(theta+1,mu,dub=1,order=3)*exp(1) / coef)
}

Loglikelihood=function(theta,gamma, O=diag(rep(1,length(gamma))),A,B,n=1,method="hg")
    #Calculates the loglikelohood where A and B are observed second and first moments A=sum X_i #X_i^t and B=sum X_i.
    #O is the orthogonal component of the covariance matrix
{
    print(theta)
    print(gamma)
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

    print("O")
    print(O)
    print("A")
    print(A)
    print("B")
    print(B)

    print("A%*%t(O)")
    print(A%*%t(O))

    l=-n*log(nc)-sum(diag(A%*%t(O)%*%diag(theta)%*%O+gamma%*%t(B)%*%t(O)))

    grad=-n*hgout[-1]/nc+c(-diag(O%*%A%*%t(O)),O%*%B)
    # grad[1]=0#seems to be very important
    #grad[1:p]=-grad[1:p]

    y=list(log=l,grad=grad)
    y
}

test = function(){
	alpha1 = c(1, 2, 3, 4, 7, 6, 8, 5)  # parameter
	ns1 = c(2,2,1,3)  # multiplicities
	cat("alpha=", alpha1, "\n")
	cat("ns=", ns1, "\n")
	G1 = hgm.FB(alpha1, ns=ns1) # HGM
	cat("HGM: ", G1, "\n")
	G2 = hgm.FB.2(alpha1, ns=ns1) # HGM via a suitable path (the same result)
	cat("HGM2:", G2, "\n")
	G3 = G.FB(alpha1, ns=ns1, method="MC") # direct computation by Monte Carlo
	cat("MC:  ", G3, "\n")
	# G.FB(alpha1, ns=ns1, method="power") # power series (do not work)
	G4 = SPA(alpha1, ns1)
	cat("SPA: ", G4, "\n")

    x=c(0.3119,0.0292,0.0707,0.3605,0.0462,0.3276, -0.0063,-0.0054,-0.0762)
    B=c(-0.0063,-0.0054,-0.0762)
    A=array(0,c(3,3))
    A[1,1]=0.3119
    A[1,2]=0.0292/2; A[2,1]=A[1,2]
    A[1,3]=0.0707/2;A[3,1]=A[1,3]
    A[2,2]=0.3605
    A[2,3]=0.0462/2; A[3,2]=A[2,3]
    A[3,3]=0.3276
    Ast=A
    Bst=-B
    A1st=2*Ast-diag(diag(Ast))
    alpha = c(0,1,2,1,2,3)
    p = length(alpha) / 2
    print(alpha)
    ll = Loglikelihood(alpha[1:p],alpha[p+(1:p)],A=A1st,B=Bst)
    print(ll$grad)

	cat("\nNear singular case:\n")
	alpha2 = c(0.067, 0.294, 0.311, 0.405, 0.663, 0.664, 0.667, 0.712, 0.784, 0.933, 0.070, 0.321, 0.288, 0.367, 0.051, 0.338, 0.798, 0.968, 0.506, 0.590)  # near singular
	ns2 = rep(1, 10)
	cat("alpha=", alpha2, "\n")
	cat("ns=", ns2, "\n")
	G1 = hgm.FB(alpha2, ns=ns2)  # HGM (does not work)
	cat("HGM: ", G1, "\n")
	G2 = hgm.FB.2(alpha2, ns=ns2)  # HGM via a suitable path
	cat("HGM2:", G2, "\n")
	G3 = G.FB(alpha2, ns=ns2, method="MC")  # direct computation by Monte Carlo
	cat("MC:  ", G3, "\n")
	G4 = SPA(alpha2, ns2)
	cat("SPA: ", G4, "\n")
}

test()
