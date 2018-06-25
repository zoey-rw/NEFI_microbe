#This function takes a vector of fixed effects
#writes a linear beta regression model for JAGS, and priors.
#saves as a file.

write_jags_beta <- function(fixed.effects,file.path){
  #text conenction for file
  fileConn <- file(file.path)
  
  #setup priors block
  prior.lines <- list()
  for(i in 1:length(fixed.effects)){
    parm <- paste0('m',i)
    prior.line <- paste0(parm,'~ dnorm(0, .001)')
    prior.lines[[i]] <- prior.line
  }
  prior.lines <- do.call('rbind',prior.lines)
  #add intercept and ucnertainty prior
  prior.lines <- rbind('m0~ dnorm(0, .001)',prior.lines)
  prior.lines <- rbind(prior.lines,'tau ~ dgamma(.1,.1)')
  
  #setup model block
  #write model lines
  pre.mod <- 
    rbind('for(i in 1:N){',
          'y[i] ~ dbeta(p[i], q[i])',
          'p[i] <- mu[i] * tau',
          'q[i] <- (1 - mu[i]) * tau'
    )
  
  #linear statement
  linear.out <- list()
  for(i in 1:length(fixed.effects)){
    linear.out[[i]] <- paste0('m',i,' * x',i,'[i]')
  }
  linear <- do.call('rbind',linear.out)
  linear <- paste0(paste0(linear, sep = ' + '), collapse = '')
  linear <- paste(linear, 'm0')
  linear <- paste('logit(mu[i]) <-',linear)
  
  #wrap up model lines
  model.lines <- rbind(pre.mod,
                       linear,
                       '}')
  
  #put it all together
  all.lines <- rbind('model{',
                     prior.lines,
                     model.lines,
                     '}')
  
  #write to file connection.
  writeLines(all.lines,fileConn)
}