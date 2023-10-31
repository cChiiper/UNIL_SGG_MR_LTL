## MR functions, all functions adapted from TwoSampleMR package

mr_ivw <- function(b_exp, b_out, se_exp, se_out)
{
	if(sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 2)
	return(list(b=NA, se=NA, pval=NA, nsnp=NA))

	ivw.res <- summary(lm(b_out ~ -1 + b_exp, weights = 1/se_out^2))
	b <- ivw.res$coef["b_exp","Estimate"]
	se <- ivw.res$coef["b_exp","Std. Error"]/min(1,ivw.res$sigma) #sigma is the residual standard error
	pval <- 2 * pnorm(abs(b/se), lower.tail=FALSE)
	Q_df <- length(b_exp) - 1
	Q <- ivw.res$sigma^2 * Q_df
	Q_pval <- pchisq(Q, Q_df, lower.tail=FALSE)
	# from formula phi =  Q/DF rearranged to to Q = phi*DF, where phi is sigma^2
	# Q.ivw<-sum((1/(se_out/b_exp)^2)*(b_out/b_exp-ivw.reg.beta)^2)
	return(list(b = b, se = se, pval = pval, nsnp=length(b_exp), Q = Q, Q_df = Q_df, Q_pval = Q_pval))
}

weighted_median <- function(b_iv, weights)
{
	betaIV.order <- b_iv[order(b_iv)]
	weights.order <- weights[order(b_iv)]
	weights.sum <- cumsum(weights.order)-0.5*weights.order
	weights.sum <- weights.sum/sum(weights.order)
	below <- max(which(weights.sum<0.5))
	b = betaIV.order[below] + (betaIV.order[below+1]-betaIV.order[below])*
	(0.5-weights.sum[below])/(weights.sum[below+1]-weights.sum[below])
	return(b)
}

weighted_median_bootstrap <- function(b_exp, b_out, se_exp, se_out, weights, nboot)
{
    b_exp.boot = rnorm(length(b_exp), mean=b_exp, sd=se_exp)
    b_out.boot = rnorm(length(b_out), mean=b_out, sd=se_out)
    betaIV.boot = b_out.boot/b_exp.boot
    med = weighted_median(betaIV.boot, weights)
	return (med)
}

mr_simple_median <- function(b_exp, b_out, se_exp, se_out, nboot)
{
	if(sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 3)
	return(list(b=NA, se=NA, pval=NA, nsnp=NA))

	b_iv <- b_out / b_exp
    weights = rep(1/length(b_exp), length(b_exp))
	b <- weighted_median(b_iv, weights)

    b_boot <- replicate(weighted_median_bootstrap(b_exp, b_out, se_exp, se_out, weights), n = nboot)
    se = sd(b_boot)

	pval <- 2 * pnorm(abs(b/se), lower.tail=FALSE)
	return(list(b=b, se=se, pval=pval, nsnp=length(b_exp)))
}

mr_mode <- function(dat, phi, alpha, penk, nboot, mode_method="all")
{
	b_exp <- dat$beta.exposure
	b_out <- dat$beta.outcome
	se_exp <- dat$se.exposure
	se_out <- dat$se.outcome

	#--------------------------------------#
	#Function to compute the point estimate#
	#--------------------------------------#
	#BetaIV.in: ratio estimates
	#seBetaIV.in: standard errors of ratio estimates

	beta <- function(BetaIV.in, seBetaIV.in, phi)
	{
		#Bandwidth rule - modified Silverman's rule proposed by Bickel (2002)
		s <- 0.9*(min(sd(BetaIV.in), mad(BetaIV.in)))/length(BetaIV.in)^(1/5)

		#Standardised weights
		weights <- seBetaIV.in^-2/sum(seBetaIV.in^-2)

		beta <- NULL

		for(cur_phi in phi)
		{
			#Define the actual bandwidth
			h <- max(0.00000001, s*cur_phi)
			#Compute the smoothed empirical density function
			densityIV <- density(BetaIV.in, weights=weights, bw=h)
			#Extract the point with the highest density as the point estimate 
			beta[length(beta)+1] <- densityIV$x[densityIV$y==max(densityIV$y)]
		}
		return(beta)
	}

	#------------------------------------------#
	#Function to estimate SEs through bootstrap#
	#------------------------------------------#
	#BetaIV.in: ratio estimates
	#seBetaIV.in: standard errors of ratio estimates
	#beta_Mode.in: point causal effect estimates
	boot <- function(BetaIV.in, seBetaIV.in, beta_Mode.in, nboot)
	{
		#Set up a matrix to store the results from each bootstrap iteration
		beta.boot <- matrix(nrow=nboot, ncol=length(beta_Mode.in))

		for(i in 1:nboot) 
		{
			#Re-sample each ratio estimate using SEs derived not assuming NOME
			BetaIV.boot      <- rnorm(length(BetaIV.in), mean=BetaIV.in, sd=seBetaIV.in[,1])
			#Re-sample each ratio estimate using SEs derived under NOME
			BetaIV.boot_NOME <- rnorm(length(BetaIV.in), mean=BetaIV.in, sd=seBetaIV.in[,2])

			#Simple mode, not assuming NOME
			beta.boot[i,1:length(phi)] <- beta(BetaIV.in=BetaIV.boot, seBetaIV.in=rep(1, length(BetaIV)), phi=phi)
			#Weighted mode, not assuming NOME
			beta.boot[i,(length(phi)+1):(2*length(phi))] <- beta(BetaIV.in=BetaIV.boot, seBetaIV.in=seBetaIV.in[,1], phi=phi)
			#Penalised mode, not assuming NOME
			weights <- 1/seBetaIV.in[,1]^2
			penalty <- pchisq(weights * (BetaIV.boot-beta.boot[i,(length(phi)+1):(2*length(phi))])^2, df=1, lower.tail=FALSE)
			pen.weights <- weights*pmin(1, penalty*penk)

			beta.boot[i,(2*length(phi)+1):(3*length(phi))] <- beta(BetaIV.in=BetaIV.boot, seBetaIV.in=sqrt(1/pen.weights), phi=phi)

			#Simple mode, assuming NOME
			beta.boot[i,(3*length(phi)+1):(4*length(phi))] <- beta(BetaIV.in=BetaIV.boot_NOME, seBetaIV.in=rep(1, length(BetaIV)), phi=phi)
			#Weighted mode, assuming NOME
			beta.boot[i,(4*length(phi)+1):(5*length(phi))] <- beta(BetaIV.in=BetaIV.boot_NOME, seBetaIV.in=seBetaIV.in[,2], phi=phi)
		}
		return(beta.boot)
	}

	#Ratio estimates
	BetaIV   <- b_out/b_exp

	#SEs of ratio estimates
	seBetaIV <- cbind(sqrt((se_out^2)/(b_exp^2) + ((b_out^2)*(se_exp^2))/(b_exp^4)), #SEs NOT assuming NOME
	se_out/abs(b_exp)) #SEs ASSUMING NOME

	#Point causal effect estimate using the simple mode
	beta_SimpleMode <- beta(BetaIV.in=BetaIV, seBetaIV.in=rep(1, length(BetaIV)), phi=phi)

	#Point causal effect estimate using the weighted mode (not asusming NOME)
	beta_WeightedMode <- beta(BetaIV.in=BetaIV, seBetaIV.in=seBetaIV[,1], phi=phi)
	weights <- 1/seBetaIV[,1]^2
	penalty <- pchisq(weights * (BetaIV-beta_WeightedMode)^2, df=1, lower.tail=FALSE)
	pen.weights <- weights*pmin(1, penalty*penk) # penalized 

	beta_PenalisedMode <- beta(BetaIV.in=BetaIV, seBetaIV.in=sqrt(1/pen.weights), phi=phi)

	#Point causal effect estimate using the weighted mode (asusming NOME)
	beta_WeightedMode_NOME <- beta(BetaIV.in=BetaIV, seBetaIV.in=seBetaIV[,2], phi=phi)

	#Combine all point effect estimates in a single vector
	beta_Mode <- rep(c(beta_SimpleMode, beta_WeightedMode, beta_PenalisedMode, beta_SimpleMode, beta_WeightedMode_NOME))

	#Compute SEs, confidence intervals and P-value
	beta_Mode.boot <- boot(BetaIV.in=BetaIV, seBetaIV.in=seBetaIV, beta_Mode.in=beta_Mode, nboot=nboot)
	se_Mode <- apply(beta_Mode.boot, 2, mad)

	CIlow_Mode <- beta_Mode-qnorm(1-alpha/2)*se_Mode
	CIupp_Mode <- beta_Mode+qnorm(1-alpha/2)*se_Mode

	P_Mode <- pt(abs(beta_Mode/se_Mode), df=length(b_exp)-1, lower.tail=F)*2

	#Vector to indicate the method referring to each row
	Method <- rep(c('Simple mode', 'Weighted mode', 'Penalised mode', 'Simple mode (NOME)', 'Weighted mode (NOME)'), each=length(phi))

	#Return a data frame containing the results
	id.exposure <- ifelse("id.exposure" %in% names(dat), dat$id.exposure[1], "")
	id.outcome <- ifelse("id.outcome" %in% names(dat), dat$id.outcome[1], "")
	Results <- data.frame(
		id.exposure = id.exposure, 
		id.outcome = id.outcome, 
		method = Method, 
		nsnp = length(b_exp), 
		b = beta_Mode, 
		se = se_Mode, 
		ci_low = CIlow_Mode, 
		ci_upp = CIupp_Mode, 
		pval = P_Mode, 
		stringsAsFactors=FALSE
	)

	if(mode_method == "all")
	{
		return(Results)
	} else {
		stopifnot(all(mode_method %in% Results$method))
		i <- which(Results$method == mode_method)
		return(list(b = Results$b[i], se = Results$se[i], pval=Results$pval[i], nsnp=length(b_exp)))
	}

	return(Results)
}

mr_simple_mode <- function(b_exp, b_out, se_exp, se_out, phi = 1, alpha = 0.05, penk = 20, nboot = 1000) 
{
	index <- !is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)
	if(sum(index) < 3)
	return(list(b=NA, se=NA, pval=NA, nsnp=NA))

	b_exp <- b_exp[index]
	b_out <- b_out[index]
	se_exp <- se_exp[index]
	se_out <- se_out[index]

	return(mr_mode(data.frame(beta.exposure=b_exp, beta.outcome=b_out, se.exposure=se_exp, se.outcome=se_out), phi, alpha, penk, nboot, mode_method="Simple mode"))
}