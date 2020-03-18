#library(foreign)
#library(sn)

# Setting the number of students:
studentN <- 25
sampleN <- 1000

# Generating the results dataframe:
results <- data.frame(	id=seq(1:studentN),
						lm1_bmi_coeff=0,
						lm1_bmi_CIl=0,
						lm1_bmi_CIu=0,
						lm1_bmi_p=0,
						lm1_r2=0,
						lm1_adjr2=0,
						lm2_bmi_coeff=0,
						lm2_bmi_CIl=0,
						lm2_bmi_CIu=0,
						lm2_bmi_p=0,
						lm2_r2=0,
						lm2_adjr2=0,
						lg1_bmi_or=0,
						lg1_bmi_CIl=0,
						lg1_bmi_CIu=0,
						lg1_bmi_p=0,
						#lg1_pr2=0,
						lg2_bmi_or=0,
						lg2_bmi_CIl=0,
						lg2_bmi_CIu=0,
						lg2_bmi_p=0,
						lg2_bd_or=0,
						lg2_bd_CIl=0,
						lg2_bd_CIu=0,
						lg2_bd_p=0,
						#lg2_pr2=0,
						lg3_bmi_or=0,
						lg3_bmi_CIl=0,
						lg3_bmi_CIu=0,
						lg3_bmi_p=0,
						lg3_bd_or=0,
						lg3_bd_CIl=0,
						lg3_bd_CIu=0,
						lg3_bd_p=0#,
						#lg3_pr2=0
						)
							
# Looping over the number of students:
for (i in 1:studentN) {

	# Covariance matrix request correlation of 0.15 between BMI and Bone Density:
	corrMat <- matrix(c(1, 0.15, 0.15, 1), nrow=2, ncol=2)

	# Creating crude data:
	smoke <- as.numeric(runif(sampleN) > 0.8)
	bmi_z <- rnorm(sampleN)
	BD_z <- rnorm(sampleN)

	# Introducing correlations
	r <- matrix(c(bmi_z,BD_z), nrow=sampleN, ncol=2)
	corrL <- t(chol(corrMat))
	unc_data <- t(r)
	c_data <- t(corrL) %*% unc_data
	data1 <- as.data.frame(t(c_data))
	data1$smoke <- smoke
	names(data1) <- c('BMI','BD_zscore','Smoker')

	# Introducing smoking as a confounder
	data1$BD_zscore <- round(data1$BD_zscore - (data1$Smoker*rnorm(sampleN,1,0.25)), digits=3)
	data1$BMI <- round((data1$BMI + (data1$Smoker*rnorm(sampleN,1,0.25)))*2 + 22.5, digits=2)

	# Generating a logistic model for broken bone:
	b_0 <- -5
	b_bmi <- 0.1
	b_bd <- -1
	prob <- exp(b_0 + b_bmi*data1$BMI + b_bd*data1$BD_zscore)/ (1 + exp(b_0 + b_bmi*data1$BMI + b_bd*data1$BD_zscore))
	data1$Break <- rbinom(n=sampleN, size=1, prob=prob)
	
	# Running the analyses
	lm_model1 <- lm(BD_zscore ~ BMI, data=data1)
	lm_model2 <- lm(BD_zscore ~ BMI + smoke, data=data1)
	lg_model1 <- glm(Break ~ BMI,family=binomial(link="logit"), data=data1)
	lg_model2 <- glm(Break ~ BMI + BD_zscore,family=binomial(link="logit"), data=data1)
	lg_model3 <- glm(Break ~ BMI + BD_zscore + Smoker,family=binomial(link="logit"), data=data1)

	# Storing the results:
	results[i,]$lm1_bmi_coeff <- lm_model1$coefficients[[2]]
	results[i,]$lm1_bmi_CIl <- lm_model1$coefficients[[2]] - (summary(lm_model1)$coefficients[4]*1.96)
	results[i,]$lm1_bmi_CIu <- lm_model1$coefficients[[2]] + (summary(lm_model1)$coefficients[4]*1.96)
	results[i,]$lm1_bmi_p <- summary(lm_model1)$coefficients[8]
	results[i,]$lm1_r2 <- summary(lm_model1)$r.squared
	results[i,]$lm1_adjr2 <- summary(lm_model1)$adj.r.squared

	results[i,]$lm2_bmi_coeff <- lm_model2$coefficients[[2]]
	results[i,]$lm2_bmi_CIl <- lm_model2$coefficients[[2]] - (summary(lm_model2)$coefficients[5]*1.96)
	results[i,]$lm2_bmi_CIu <- lm_model2$coefficients[[2]] + (summary(lm_model2)$coefficients[5]*1.96)
	results[i,]$lm2_bmi_p <- summary(lm_model2)$coefficients[11]
	results[i,]$lm2_r2 <- summary(lm_model2)$r.squared
	results[i,]$lm2_adjr2 <- summary(lm_model2)$adj.r.squared

	results[i,]$lg1_bmi_or <- exp(summary(lg_model1)$coefficients[2])
	results[i,]$lg1_bmi_CIl <- exp(summary(lg_model1)$coefficients[2] - summary(lg_model1)$coefficients[4])
	results[i,]$lg1_bmi_CIu <- exp(summary(lg_model1)$coefficients[2] + summary(lg_model1)$coefficients[4])
	results[i,]$lg1_bmi_p <- summary(lg_model1)$coefficients[8]

	results[i,]$lg2_bmi_or <- exp(summary(lg_model2)$coefficients[2])
	results[i,]$lg2_bmi_CIl <- exp(summary(lg_model2)$coefficients[2] - summary(lg_model2)$coefficients[5])
	results[i,]$lg2_bmi_CIu <- exp(summary(lg_model2)$coefficients[2] + summary(lg_model2)$coefficients[5])
	results[i,]$lg2_bmi_p <- summary(lg_model2)$coefficients[11]
	results[i,]$lg2_bd_or <- exp(summary(lg_model2)$coefficients[3])
	results[i,]$lg2_bd_CIl <- exp(summary(lg_model2)$coefficients[3] - summary(lg_model2)$coefficients[6])
	results[i,]$lg2_bd_CIu <- exp(summary(lg_model2)$coefficients[3] + summary(lg_model2)$coefficients[6])
	results[i,]$lg2_bd_p <- summary(lg_model2)$coefficients[12]

	results[i,]$lg3_bmi_or <- exp(summary(lg_model3)$coefficients[2])
	results[i,]$lg3_bmi_CIl <- exp(summary(lg_model3)$coefficients[2] - summary(lg_model3)$coefficients[6])
	results[i,]$lg3_bmi_CIu <- exp(summary(lg_model3)$coefficients[2] + summary(lg_model3)$coefficients[6])
	results[i,]$lg3_bmi_p <- summary(lg_model3)$coefficients[14]
	results[i,]$lg3_bd_or <- exp(summary(lg_model3)$coefficients[3])
	results[i,]$lg3_bd_CIl <- exp(summary(lg_model3)$coefficients[3] - summary(lg_model3)$coefficients[7])
	results[i,]$lg3_bd_CIu <- exp(summary(lg_model3)$coefficients[3] + summary(lg_model3)$coefficients[7])
	results[i,]$lg3_bd_p <- summary(lg_model3)$coefficients[15]

	# Writing files:
	datafilename <- paste("PHPS40460_Data_group",i,".csv")
	write.csv(data1, file=datafilename, row.names=F)

}

resultsfilename <- paste("results.csv")
write.csv(results, file=resultsfilename, row.names=F)

