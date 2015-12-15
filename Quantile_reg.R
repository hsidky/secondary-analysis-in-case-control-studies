library(quantreg)
# ---------------------- WEE secondary quantile regression -------------------------- # 
#                                                                                     #
#  This code serves as the function to conduct secondary analysis for the continuous  #
#  secondary traits using quantile regression in genetic case-control studies. 	      #
#  Note: the quantile regression package "quantreg" is required. Users can use the    #
#  following codes to obtain "quantreg" package.        			      #
#  > install.packages("quantreg")                                                     #
#  > library(quantreg)                  		                              #
#                                                                                     #
#  *function* 				  					      #			  
#	-secondary_rq(formula, D, tau, data, pd_pop, iter=10, boot=0): obtain the     #
#	point estimate, bootstrap SE, Wald test statistics and p-values of the        #
#	associations between the continous secondary trait $Y$ and the variants $X$   #
#	adjusting for the covariates under quantile regression.                       #
#                                                                                     #
#  *input*									      #
#	-data: dataset with real observation.					      #
#	-formula: the secondary trait given SNPs and covariates.  e.g. y~x+z	      #
# 	-D: primary disease (case-control status). 				      #
# 	-tau: the quantile level to be estimated. Multiple taus can be chosen         # 	
#	-pd_pop: the population disease prevelance of primary disease.		      #
#	-iter: number of generating pseudo observations. (iter=10 by default)	      #
#	-boot: number of bootstrape samples. (boot=0 by default)		      #
#                                                                                     #
#  *output*									      #
#  -if boot=0, output the point estimates					      #
#  -if boot!=0, output the point estimates, bootstrap SE, Wald test and p-values      #   
#                                                                                     #
#  *Example:*                                                   		      #
#  - create the sample dataset                                                        #                                                                       
#   x=cbind(rbinom(3000, 2, 0.3),rbinom(3000, 2, 0.3))                                #
#   y1=rnorm(3000)                                                                    #
#   D1=c(rep(0, 2000), rep(1, 1000))                                                  #
#   dat=data.frame(cbind(x, y1, D1))                                                  #
#   pd=0.1                                                                            #
#   colnames(dat)[1:2]=c("x1", "x2")                                                  #
#  - run the proposed approach based on quantile regression                           #
#   secondary_rq(y1~x1, D1, tau=tau=0.5, dat, pd)                                     #
#   secondary_rq(y1~x1+x2, D1, tau=c(0.25, 0.5), dat, pd)                             #
#   secondary_rq(y1~x1+x2, D1, tau=c(0.25, 0.5), dat, pd, boot=100)                   #
# ----------------------------------------------------------------------------------- # 


secondary_rq <-function(formula, D, data, pd_pop, tau, iter=10, boot=0) {
	mf<-model.frame(formula, data=data)
	y<-model.response(mf, "numeric")
	namesx=all.vars(formula)[-1]
	xx<-model.matrix(attr(mf, "terms"), data=mf)
	x=matrix(xx[,-1], nrow=nrow(data), byrow=F); colnames(x)=namesx; p=dim(x)[2]
	temp.data=data.frame(cbind(D, y, x))
	f=update(formula, y  ~ . )
    n1=sum(D==1);n0=sum(D==0)      
	ltau=length(tau);flag=seq(1, ltau)    
	
	# compute the weight p(D|X)
	gamma=coef(glm(D~x, family="binomial"))
	PO = function(gamma0){
		gamma[1]=gamma0
  		(mean( exp(xx%*%gamma)/(1+exp(xx%*%gamma)))  -pd_pop)^2
	}; gamma[1]= suppressWarnings(optim(gamma[1],PO)$par)
	temp.data$estpx = exp(xx%*%gamma)/(1+exp(xx%*%gamma))
    
	# estimate py in cases and controls separately
  
	pyD1=rq(f,  tau=-1, data=temp.data[which(temp.data$D==1),])
	pyD0=rq(f,  tau=-1, data=temp.data[which(temp.data$D==0),])
	xcase=xx[which(D==1),];xcontrol=xx[which(D==0),]
	pred1=	as.matrix(xcase) %*% pyD0$sol[4:length(pyD0$sol[,1]), ] # pseudo control
	pred0=  as.matrix(xcontrol) %*% pyD1$sol[4:length(pyD0$sol[,1]), ]

	step1=apply(cbind(pred1[,1], pred1), 1, function(x) rearrange(stepfun(pyD0$sol[1,], x)))
	step0=apply(cbind(pred0[,1], pred0), 1, function(x) rearrange(stepfun(pyD1$sol[1,], x)))

	
	# generate pseudo observations for $iter$ times and get the averaged $iter$ estimates as coefficient
	pseudo=NULL
	for (iiter in 1:iter) {
		pseudo1=unlist(lapply(step1, function(x) x(runif(1))))
		pseudo0=unlist(lapply(step0, function(x) x(runif(1))))
		data1=cbind( rep(0, n1), pseudo1, temp.data[which(D==1), c(namesx, "estpx")] ) 
		data0=cbind( rep(1, n0), pseudo0, temp.data[which(D==0), c(namesx, "estpx")] )
		colnames(data1)[1:2]=c("D",  "y")
		colnames(data0)[1:2]=c("D", "y")
		alldat=rbind(temp.data, data1, data0)
		alldat$estpx[which(alldat$D==0)] = 1- alldat$estpx[which(alldat$D==0)]
		pseudo =rbind(pseudo, rq(f , tau=tau, data= alldat, weight=estpx)$coef)		
	}
	# the point estmate
	if (ltau ==1) {
		coef=list(apply(pseudo, 2, mean)[-1])
	} else {  
		coef_est=do.call("rbind", lapply(  (seq(0:p)%%(p+1))[-1]  , function(x) apply(pseudo[which(seq(nrow(pseudo))%%(p+1)== x),], 2, mean)))
	coef=split(coef_est, seq(ncol(coef_est)))   } 
	



	# bootstrap SE	
	if(boot==0){ 	names(coef)=paste("tau=", tau)
		for (iflag in flag) {names(coef[[iflag]])=namesx}; coef
		} else {
	sample_cases = temp.data[which(temp.data$D==1),]
	sample_controls= temp.data[which(temp.data$D==0),]
	
	bootcoef=vector("list", ltau)
	for (iboot in 1: boot ) {

		boot_cases_sample = sample_cases[sample(n1, n1, replace=T),]
		boot_controls_sample = sample_controls[sample(n0, n0, replace=T),]

		bootsample= rbind(boot_cases_sample, boot_controls_sample)
		bootmf<-model.frame(f, data=bootsample)
		bootxx<-model.matrix(attr(bootmf, "terms"), data=bootmf)
		xcase= bootxx[which(bootsample$D==1), ]
		xcontrol= bootxx[which(bootsample$D==0), ]

		# compute p(D|X)
		gamma=coef(glm(update(formula, D  ~ . ), family="binomial", data=bootsample))
		PO = function(gamma0){
			gamma[1]=gamma0
  			(mean( exp(bootxx%*%gamma)/(1+exp(bootxx%*%gamma)))  -pd_pop)^2
		}; gamma[1]= suppressWarnings(optim(gamma[1],PO)$par)
		bootsample$estpx = exp(bootxx%*%gamma)/(1+exp(bootxx%*%gamma))


    		pyD1=rq(f, tau=-1,  data= boot_cases_sample)
		pyD0=rq(f, tau=-1,  data= boot_controls_sample)
		pred1 = as.matrix(xcase) %*% pyD0$sol[4:length(pyD0$sol[,1]), ]
		pred0 = as.matrix(xcontrol) %*% pyD1$sol[4:length(pyD0$sol[,1]), ]
		
		step1=apply(cbind(pred1[,1], pred1), 1, function(x) rearrange(stepfun(pyD0$sol[1,], x)))
		step0=apply(cbind(pred0[,1], pred0), 1, function(x) rearrange(stepfun(pyD1$sol[1,], x)))
		
	
		# generate pseudo observations for T(=10) times and get the averaged T estimates as coefficient
		pseudo=NULL
		for (iiter in 1:iter) {

		pseudo1=unlist(lapply(step1, function(x) x(runif(1))))
		pseudo0=unlist(lapply(step0, function(x) x(runif(1))))

		data1=cbind( rep(0, n1), pseudo1, boot_cases_sample[, c(namesx,  "estpx")] ) 
		data0=cbind( rep(1, n0), pseudo0, boot_controls_sample[, c(namesx, "estpx")] )
		colnames(data1)[1:2]=c("D",  "y")
		colnames(data0)[1:2]=c("D", "y")

		alldat=rbind(bootsample, data1, data0)
		alldat$estpx[which(alldat$D==0)] = 1- alldat$estpx[which(alldat$D==0)]

		pseudo=rbind(pseudo, rq(f, tau=tau, data= alldat, weight=estpx)$coef)
		}

	
		if (ltau==1) {
		bootcoef=rbind(unlist(bootcoef), apply(pseudo, 2, mean)[-1])
		} else {  
		temp=t(do.call("rbind", lapply(    (seq(0:p)%%(p+1))[-1]  , function(x) apply(pseudo[which(seq(nrow(pseudo))%%(p+1)== x),], 2, mean)))  ); for (iflag in flag) {
		bootcoef[[iflag]]=rbind(bootcoef[[iflag]], temp[iflag,])}
	}
	

	}
	if (ltau==1) {bootcoef=list(bootcoef)}
	var=lapply(seq(1:ltau), function(x) apply(bootcoef[[x]], 2, var))
	wald= lapply(seq(1:ltau), function(x) coef[[x]]^2/var[[x]])
	pvalue=lapply(seq(1:ltau), function(x) pchisq(wald[[x]], df=1, lower.tail=F))
	
		
	TAB=vector("list", ltau)
	for ( iflag in flag) {
		TAB[[iflag]]<-cbind(Estimate= coef[[iflag]], StdErr=sqrt(var[[iflag]]), Wald=wald[[iflag]], p.value=pvalue[[iflag]])					
		rownames(TAB[[iflag]])= namesx
	}
	names(TAB)=paste('tau=', tau)
 	TAB    
	}	

}




