rm(list=ls())
#Chisq
# ---------------------------- secondary linear regression -------------------------- # 
#                                                                                     #
#  This code serves as the WEE function to conduct secondary analysis for continuous  #
#  secondary traits using linear regression in genetic case-control studies           #
#                                                                                     #
#  *function* 																		  #			  
#	-WEE_cts(formula, D, data, pd_pop,  boot=0): obtain the point estimate,	   	      #
#		bootstrap SE, Wald test statistics and p-values of the associations between   #
#		the	continuous secondary trait $Y$ and the genetic variants $X$	adjusting     #
#		for the covariates under least square regression.                             #
#                                                                                     #
#  *input*																			  #
#	-data: dataset with real observation.											  #
#	-formula: the secondary trait given SNPs and covariates.  e.g. y~x+z				  #
# 	-D: primary disease (case-control status). 										  #
#	-pd_pop: the population disease prevelance of primary disease.					  #
#	-boot: number of bootstrape samples. (boot=0 by default)				  		  #
#                                                                                     #
#  *output*																			  #
#  -if boot=0, output the point estimates											  #
#  -if boot!=0, output the point estimates, bootstrap SE, Chisq test and p-values 	  #
#                                                                                     #
#  *Example*                                                                          #
# 	x=cbind(rbinom(3000, 2, 0.3),rbinom(3000, 2, 0.3))                                #
# 	y1=rnorm(3000)                                                                    #
# 	D1=c(rep(0, 1000), rep(1, 2000))                                                  #
# 	dat=data.frame(cbind(x, y1, D1))                                                  #
# 	pd=0.1                                                                            #
# 	colnames(dat)[1:2]=c("x1", "x2")                                                  #
# 	secondary_cts(y1~x1+x2, D1, dat, pd)                                              #
# 	secondary_cts(y1~x1+x2, D1, dat, pd, boot=100)                                    #                                                                            
# ----------------------------------------------------------------------------------- # 



secondary_cts<-function(formula, D, data, pd_pop,  boot=0) {
	#formula=y1~x1+x2
	mf<-model.frame(formula, data=data)
	y<-model.response(mf, "numeric")
	namesx=all.vars(formula)[-1]
	xx<-model.matrix(attr(mf, "terms"), data=mf)
	x=matrix(xx[,-1], nrow=nrow(data), byrow=F); colnames(x)=namesx
	#D=as.numeric(D1)
	temp.data=data.frame(cbind(D, y, x))
	f=update(formula, y  ~ . )
    n1=sum(D==1);n0=sum(D==0)      
    	
	# compute the weight p(D|X)
	gamma=coef(glm(D~x, family="binomial"))
	PO = function(gamma0){
		gamma[1]=gamma0
  		(mean( exp(xx%*%gamma)/(1+exp(xx%*%gamma)))  -pd_pop)^2
	}
	gamma[1]= suppressWarnings(optim(gamma[1],PO)$par)
	temp.data$estpx = exp(xx%*%gamma)/(1+exp(xx%*%gamma))
    
    # estimate py in cases and controls separately  
    pyD1=lm(f, data=temp.data[which(temp.data$D==1),])
	pyD0=lm(f, data=temp.data[which(temp.data$D==0),])
	py1 = predict(pyD0, newdata=data.frame(temp.data[which(temp.data$D==1),])) # pseudo control
	py0 = predict(pyD1, newdata=data.frame(temp.data[which(temp.data$D==0),]))
	data1=cbind( rep(0, n1), py1, temp.data[which(D==1), c(namesx, "estpx")] ) 
	data0=cbind( rep(1, n0), py0, temp.data[which(D==0), c(namesx, "estpx")] )
	colnames(data1)[1:2]=c("D",  "y")
	colnames(data0)[1:2]=c("D", "y")
	alldat=rbind(temp.data, data1, data0)
	alldat$estpx[which(alldat$D==0)] = 1- alldat$estpx[which(alldat$D==0)]
	# the point estmate
	coef=lm(f ,data= alldat, weight=estpx)$coef[-1]



	# bootstrap SE	
	if(boot==0){list(coefficients =coef)} else 
	{ 
	sample_cases = temp.data[which(temp.data$D==1),]
	sample_controls= temp.data[which(temp.data$D==0),]
	
	bootcoef=NULL	
	for (iboot in 1: boot ) {
	
	boot_cases_sample = sample_cases[sample(n1, n1, replace=T),]
	boot_controls_sample = sample_controls[sample(n0, n0, replace=T),]
	
	bootsample= rbind(boot_cases_sample,boot_controls_sample)
	
	bootmf<-model.frame(f, data=bootsample); 
	bootxx<-model.matrix(attr(bootmf, "terms"), data=bootmf)
	# compute p(D|X)
	gamma=coef(glm(update(formula, D  ~ . ), family="binomial", data=bootsample))
	PO = function(gamma0){
		gamma[1]=gamma0
  		(mean( exp(bootxx%*%gamma)/(1+exp(bootxx%*%gamma)))  -pd_pop)^2
	}
	gamma[1]= suppressWarnings(optim(gamma[1],PO)$par)
	bootsample$estpx = exp(bootxx%*%gamma)/(1+exp(bootxx%*%gamma))
	
		
    	pyD1=lm(f, data=boot_cases_sample)
		pyD0=lm(f, data= boot_controls_sample)
		py1 = predict(pyD0, newdata=data.frame(boot_cases_sample)) # pseudo control
		py0 = predict(pyD1, newdata=data.frame(boot_controls_sample))
		data1=cbind( rep(0, n1), py1, boot_cases_sample[, c(namesx,  "estpx")] ) 
		data0=cbind( rep(1, n0), py0, boot_controls_sample[, c(namesx, "estpx")] )
		colnames(data1)[1:2]=c("D",  "y")
		colnames(data0)[1:2]=c("D", "y")
	
		alldat=rbind(bootsample, data1, data0)
		alldat$estpx[which(alldat$D==0)] = 1- alldat$estpx[which(alldat$D==0)]
		bootcoef=rbind(bootcoef, lm(f ,data= alldat, weight=estpx)$coef[-1])

	}
	

	var=apply(bootcoef, 2, var)
	chisq= coef^2/var
	pvalue=pchisq(chisq, df=1, lower.tail=F)
	
	TAB<-cbind(Estimate=coef, StdErr=sqrt(var), Chisq = chisq, p.value=pvalue)

 	TAB    
	}

}


	x=cbind(rbinom(3000, 2, 0.3),rbinom(3000, 2, 0.3))                                #
	y1=rnorm(3000)                                                                    #
	D1=c(rep(0, 1000), rep(1, 2000))                                                  #
	dat=data.frame(cbind(x, y1, D1))                                                  #
	pd=0.1                                                                            #
	colnames(dat)[1:2]=c("x1", "x2")                                                  #
	secondary_cts(y1~x1+x2, D1, dat, pd)                                              #
	secondary_cts(y1~x1+x2, D1, dat, pd, boot=100) 




