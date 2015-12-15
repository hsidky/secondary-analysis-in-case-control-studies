# --------------------- WEE secondary logistic regression --------------------------- # 
#  Contact: Xiaoyu Song (xs2148@cumc.columbia.edu)                                    #
#                                                                                     #
#  This code serves as the function to conduct the secondary analysis for the binary  #
#  secondary traits using logistic regression in genetic case-control studies         #
#                                                                                     #
#  *function* 									      #									  #
#                                                                                     #			  
#	-WEE_binary(formula, D, data, pd_pop, iter=10, boot=0): obtain the point      #
#	estimate, bootstrap SE, Wald test statistics and p-values of the associations #
#	between the binary secondary trait $Y$ and the genetic variants $X$ adjusting #
#	for the covariates under logistic regression.                                 #
#                                                                                     #
#  *input*									      #
#	-data: dataset with real observation.					      #
#	-formula: the secondary trait given SNPs and covariates.  e.g. y~x+z	      #
# 	-D: primary disease (case-control status). 				      #
#	-pd_pop: the population disease prevelance of primary disease.		      #
#	-iter: number of generating pseudo observations. (iter=10 by default)	      #
#	-boot: number of bootstrape samples. (boot=0 by default)	      	      #
#                                                                                     #
#  *output*								  	      #
#  -if boot=0, output the point estimates					      #
#  -if boot~=0, output the point estimates, bootstrap SE, Wald test and p-values      #                                                                            
#                                                                                     #
#  *example*								  	      #
#  x=cbind(rbinom(3000, 2, 0.3),rbinom(3000, 2, 0.2))				      #
#  y1=rbinom(3000, 1, 0.3)							      #
#  D1=c(rep(0, 1000), rep(1, 2000))						      #
#  dat=data.frame(cbind(x, y1, D1))						      #
#  pd=0.1					  				      #
#  colnames(dat)[1:2]=c("x1", "x2")		  				      #
#  secondary_binary(y1~x1+x2, D1, dat, pd)	  				      #
#  secondary_binary(y1~x1+x2, D1, dat, pd, boot=100)	  			      #
# ----------------------------------------------------------------------------------- # 


secondary_binary <-function(formula, D, data, pd_pop, iter=10, boot=0) {
	mf<-model.frame(formula, data=data)
	y<-model.response(mf, "numeric")
	namesx=all.vars(formula)[-1]
	xx<-model.matrix(attr(mf, "terms"), data=mf)
	x=matrix(xx[,-1], nrow=nrow(data), byrow=F); colnames(x)=namesx
	temp.data=data.frame(cbind(D, y, x))
	f=update(formula, y  ~ . )
    n1=sum(D==1);n0=sum(D==0)      
    	
	# compute the weight p(D|X)
	gamma=coef(glm(D~x, family="binomial"))
	PO = function(gamma0){
		gamma[1]=gamma0
  		(mean( exp(xx%*%gamma)/(1+exp(xx%*%gamma)))  -pd_pop)^2
	}; gamma[1]= suppressWarnings(optim(gamma[1],PO)$par)
	temp.data$estpx = exp(xx%*%gamma)/(1+exp(xx%*%gamma))
    
    	# estimate P(Y=1) in cases and controls separately
    	pyD1=glm(f, family="binomial", data=temp.data[which(temp.data$D==1),])
	pyD0=glm(f, family="binomial", , data=temp.data[which(temp.data$D==0),])
	pred1= predict(pyD0, newdata=data.frame(temp.data[which(temp.data$D==1),])) # generate pseudo control
	pred0= predict(pyD1, newdata=data.frame(temp.data[which(temp.data$D==0),])) # generate pseudo case
	py1=exp(pred1)/(1+exp(pred1)) 
	py0=exp(pred0)/(1+exp(pred0)) 
	
	# generate pseudo observations for $iter$ times and get the averaged $iter$ estimates as coefficient
	pseudo=NULL
	for (iiter in 1:iter) {
		pseudo1=rbinom(n1, 1, py1)
		pseudo0=rbinom(n0, 1, py0)
		data1=cbind( rep(0, n1), pseudo1, temp.data[which(D==1), c(namesx, "estpx")] ) 
		data0=cbind( rep(1, n0), pseudo0, temp.data[which(D==0), c(namesx, "estpx")] )
		colnames(data1)[1:2]=c("D",  "y")
		colnames(data0)[1:2]=c("D", "y")
		alldat=rbind(temp.data, data1, data0)
		alldat$estpx[which(alldat$D==0)] = 1- alldat$estpx[which(alldat$D==0)]
		pseudo=rbind(pseudo, suppressWarnings(glm(f , family="binomial", data= alldat, weight=estpx)$coef[-1]))
	}
	# the point estmate
	coef=apply(pseudo, 2, mean)

	# bootstrap SE	
	if(boot==0){list(coefficients =coef)} else {
	sample_cases = temp.data[which(temp.data$D==1),]
	sample_controls= temp.data[which(temp.data$D==0),]

	bootcoef=NULL	
	for (iboot in 1: boot ) {
		boot_cases_sample = sample_cases[sample(n1, n1, replace=T),]
		boot_controls_sample = sample_controls[sample(n0, n0, replace=T),]
		bootsample= rbind(boot_cases_sample, boot_controls_sample)
		bootmf<-model.frame(f, data=bootsample)
		bootxx<-model.matrix(attr(bootmf, "terms"), data=bootmf)
		# compute the weight p(D|X)
		gamma=coef(glm(update(formula, D  ~ . ), family="binomial", data=bootsample))
		PO = function(gamma0){
			gamma[1]=gamma0
	  		(mean( exp(bootxx%*%gamma)/(1+exp(bootxx%*%gamma)))  -pd_pop)^2
		}; gamma[1]= suppressWarnings(optim(gamma[1],PO)$par)
		bootsample$estpx = exp(bootxx%*%gamma)/(1+exp(bootxx%*%gamma))

	    	pyD1=glm(f, family="binomial", data= boot_cases_sample)
		pyD0=glm(f, family="binomial", data= boot_controls_sample)
		pred1= predict(pyD0, newdata=data.frame(boot_cases_sample)) # pseudo control
		pred0= predict(pyD1, newdata=data.frame(boot_controls_sample))
		py1=exp(pred1)/(1+exp(pred1)) 
		py0=exp(pred0)/(1+exp(pred0)) 
	
		# generate pseudo observations for T(=10) times and get the averaged T estimates as coefficient
		pseudo=NULL
		for (iiter in 1:iter) {
			pseudo1=rbinom(n1, 1, py1)  
			pseudo0=rbinom(n0, 1, py0)
			data1=cbind( rep(0, n1), pseudo1, boot_cases_sample[, c(namesx,  "estpx")] ) 
			data0=cbind( rep(1, n0), pseudo0, boot_controls_sample[, c(namesx, "estpx")] )
			colnames(data1)[1:2]=c("D",  "y")
			colnames(data0)[1:2]=c("D", "y")

			alldat=rbind(bootsample, data1, data0)
			alldat$estpx[which(alldat$D==0)] = 1- alldat$estpx[which(alldat$D==0)]

			pseudo=rbind(pseudo, suppressWarnings(glm(f, family="binomial", data= alldat, weight=estpx)$coef[-1]))
		}
	bootcoef=rbind(bootcoef, apply(pseudo, 2, mean))
	}
	var=apply(bootcoef, 2, var)
	wald= coef^2/var
	pvalue=pchisq(wald, df=1, lower.tail=F)
	
	TAB<-cbind(Estimate=coef, StdErr=sqrt(var), Wald=wald, p.value=pvalue)
 	TAB    
	}
}
