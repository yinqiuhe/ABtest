# Running time: 17.61s
#cat("The range of test statistics, from 0 to (default=10): ")
#int<-as.numeric(readLines(con=stdin(),1))
#if (is.na(int)) 
int<-10

B<-10000; pdf<-NULL
for (i in 1:B){
	pdfi<-besselK(x=int*i/B, nu=0)
	pdf<-c(pdf, pdfi)
	#print(i); flush.console()
}
print(int)

myp<-function(cut){
	select<-(int*1:B/B)>cut
	pdf.sub<-pdf[select]
	pval<-sum(pdf.sub)/sum(pdf)
	return(pval)
}

MT_Comp<-function(a, b){
	ab<-a*b
	pp0<-sapply(abs(ab)/sqrt(1), myp)
	pp1<-sapply(abs(ab)/sd(a), myp)
	pp2<-sapply(abs(ab)/sd(b), myp)
	pp.comp<-pp1+pp2-pp0
	return(pp.comp)
}

My_MT_Comp<-function(mydata, n, num_alpha, num_beta){
  AB = computeTest_C(mydata, n, num_alpha, num_beta)
  pp.comp<-myp( abs(AB[1]*AB[2]))
  return(pp.comp)
}



test_mediate<-function(n, mydata_original_tmp, other_covariates, sim_num = 10^5){
  #process data for using mediation package
  my_data1 <- cbind( mydata_original_tmp[,"exposure"], mydata_original_tmp[,"mediator"], other_covariates )
  colnames(my_data1)[c(1,2)] <- c("exposure", "mediator")
  my_data2 <- cbind( mydata_original_tmp[,"mediator"], mydata_original_tmp[,"outcome"],  mydata_original_tmp[,"exposure"], other_covariates)
  colnames(my_data2)[c(1,2,3)] <- c( "mediator", "outcome", "exposure")
  
  #fit linear regression model
  my_data1 = as.data.frame(my_data1)
  colnames(my_data1)[1:5] = c("E", "M", "SEX", "AGE", "X0") #E: exposure, M: mediator, X0: intercept
  my_data2 = as.data.frame(my_data2)
  colnames(my_data2)[1:6] = c("M", "O",  "E", "SEX", "AGE", "X0")  #E: exposure, M: mediator, O: outcome, X0: intercept
  
  fit.em <- lm( M ~ 0 + ., data = my_data1) #+0 because intercept is already included
  fit.mo <- lm( O ~  0 + ., data = my_data2)
  
  med_res <- summary(mediate(fit.em, fit.mo, treat ="E", mediator ="M", boot=TRUE, sims=sim_num))
  #return p-value
  return(med_res$d1.p)
}

