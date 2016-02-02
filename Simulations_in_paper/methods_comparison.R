# This code serves as comparing WEE with exisitng methods. The exisiting methods compared here include
# 	1. classical methods: direct regression in combined case-control sample, case sample, control sample and combined case-control sample adjusting for primary disease
# 	2. IPW approach
# 	3. SPML approach
# Note: SPML approach has its indepedent software SPREG. To compare with it, run SPREG software separately using shell scrpit, and import its results. For details about SPREG, please visit http://dlin.web.unc.edu/software/spreg-2/


result=NULL
for (i in 1:500) {
	print(i)
	
	# get the data
	string=paste("Bbasecase/input", i, ".txt", sep="")
	sample=as.data.frame(t(read.table(string)))
	
	# classical method: direct regression in combined case-control sample
	combined=glm(y~x+z, data=sample, family="binomial")
	
	# classical method: direct regression in case sample
	case=glm(y~x+z, data=sample, family="binomial", subset=D==1)

	# classical method: direct regression in control sample
	control=glm(y~x+z, data=sample, family="binomial", subset=D==0)

	# classical method: direct regression in combined case-control sample adjusting for primary disease status
	stratified=glm(y~x+z+D, data=sample, family="binomial")
	
	# IPW approach under random sample 
	sample$weight[which(sample$D==0)]=(450000/2000)/500000
	sample$weight[which(sample$D==1)]=(50000/2000)/500000
	ft=glm(y~x+z, data=sample, family="binomial", weights=weight)

	# WEE
	wee= WEE_binary(y~x+z, data=sample, D=D, pd=0.1)
	
	# combine the results of 6 approaches
	result=rbind(result, cbind(coef(combined)[2], coef(case)[2],coef(control)[2],coef(stratified)[2], coef(ft)[2], coef(wee)[1]))	

}

(apply(result, 2, mean)-0.2)/0.2
apply(result, 2, sd)

sum((result[,1]-0.2)^2)/500
sum((result[,2]-0.2)^2)/500
sum((result[,3]-0.2)^2)/500
sum((result[,4]-0.2)^2)/500
sum((result[,5]-0.2)^2)/500
sum((result[,6]-0.2)^2)/500

# For SPML approach, we run shell codes using existing SPREG software. Here is a sample shell code

#!/bin/bash
for i in {1..500}
do
   echo "Processing $i"
   ./spreg input${i}.txt output${i}.txt 0.1 1
done

# Codes to combine all the output files in shell
# awk ' NR %2 == 0’ output* > Bbasecase_lin.txt


# We import the result for comparison
dat=read.table("/Lin/Bbasecase_lin.txt")
colnames(dat)<-c("Gene_number", "Estimate", "Std_Error", "Z_stat", "p_value")

(mean(dat[,2])-0.2)/0.2
sd(dat[,2])
sum((dat[,2]-0.2)^2)/500
