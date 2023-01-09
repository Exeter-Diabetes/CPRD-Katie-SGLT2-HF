#Run HTE models and extract coefs

#Define subsets of interest to calc HTE 
#overall
overall <- train
#sglt2.best
sglt2.best <- train %>% dplyr::filter(bestdrug=="SGLT2")
#dpp4.best
dpp4.best <- train %>% dplyr::filter(bestdrug=="GLP1")
#Subgroups by predicted treatment difference
#df1 <- train %>% dplyr::filter(hba1c_diff<= -10) 
df2 <- train %>% dplyr::filter(hba1c_diff<= -5) 
df3 <- train %>% dplyr::filter(hba1c_diff> -5 & hba1c_diff<= -3) 
df4 <- train %>% dplyr::filter(hba1c_diff> -3 & hba1c_diff <= 0) 
df5 <- train %>% dplyr::filter(hba1c_diff> 0 & hba1c_diff < 3) 
df6 <- train %>% dplyr::filter(hba1c_diff>= 3 & hba1c_diff < 5) 
df7 <- train %>% dplyr::filter(hba1c_diff>= 5) 
#df8 <- train %>% dplyr::filter(hba1c_diff>= 10)


#Run HTE models and extract coefs for each subset

hte.model.coefs <- function(x,nmodels) {
  mnumber = c(1:nmodels)
  models <- as.list(1:nmodels)
  nobs <- vector()
  coef <- vector()
  lower <- vector()
  upper <- vector()
  pvalue <- vector()
  data <- x
  
  for(i in mnumber) {
    models[[i]] <- lm(as.formula(f[[i]]),data=data)
    nobs <- append(nobs,nobs(models[[i]]))
    coef <- append(coef,models[[i]]$coefficients[2])
    confint_all <- confint(models[[i]], levels=0.95)
    lower <- append(lower,confint_all[2,1])
    upper <- append(upper,confint_all[2,2])
    pvalue <- append(pvalue,summary(models[[i]])$coefficients[2,4])
  }
  
  datasetname = c(deparse(substitute(x)),deparse(substitute(x)),deparse(substitute(x)))
  x <- data.frame(datasetname,modelname,cbind(nobs,coef,lower,upper,pvalue))
  rownames(x) <- c()
  return(x)
}  
#define model formula
formula1 <- "posthba1cfinal~drugclass"
formula2 <- "posthba1cfinal~drugclass+prehba1c+ncurrtx+drugline+rcs(hba1cmonth,3)"
formula2 <- "posthba1cfinal~drugclass+prehba1c"
formula3 <- "posthba1cfinal~drugclass+rcs(prehba1c,3)+ncurrtx+drugline+rcs(hba1cmonth,3)+
rcs(preegfr,3)+
rcs(prealtlog,3)+
rcs(agetx,3)+
rcs(prebmi,3)+
sex+
rcs(t2dmduration,3)+
rcs(pretotalcholesterol,3)"

# rcs(prehdl,3)+
#   rcs(prebilirubinlog,3)+
#   rcs(prealbumin_blood,3)+

# formula3 <- "posthba1cfinal~drugclass+rcs(prehba1c,3)+ncurrtx+drugline+rcs(hba1cmonth,3)+
# rcs(preegfr,3)+
# rcs(prealtlog,3)+
# rcs(agetx,3)+
# rcs(prebmi,3)+
# sex+
# rcs(t2dmduration,3)+
# rcs(prehdl,3)+
# rcs(prebilirubinlog,3)+
# rcs(prealbumin_blood,3)+
# rcs(pretotalcholesterol,3)"

modelname <- c("unadj","simple","full")
f <- as.list(c(formula1,formula2,formula3))

dflist <- list(overall,sglt2.best,dpp4.best,df2,df3,df4,df5,df6,df7)
res.list <- lapply(dflist, function(df) {
  hte.model.coefs(df,3)
})

res.list <- bind_rows(res.list, .id="column_label")
res.lab <- c(rep("overall",3),rep("sglt2.best",3),rep("glp1.best",3),rep("sglt.best5",3),rep("sglt.best3",3),rep("sglt.best0-3",3),
             rep("glp1.best0-3",3),rep("glp1.best3",3),rep("glp1.best5",3))
res.list <- data.frame(res.lab,res.list) %>% dplyr::select(-column_label,-datasetname)
res.list %>% dplyr::filter(modelname=="full") #train adjusted list

#Calibration plot HTE

#Unadjusted obs vs pred (may be different treatment line etc.)

#Define tenths
train <- train %>% mutate(hba1c_diff.q = ntile(hba1c_diff, 10))    

#define dataset with predicted values
t1 <- ddply(train, "hba1c_diff.q", dplyr::summarise,
            N    = length(hba1c_diff),
            hba1c_diff.pred = mean(hba1c_diff))

#check some patients actually prescribed both drugs in each tenth
# ddply(train, c("hba1c_diff.q","drugclass"), dplyr::summarise,
#       N    = length(posthba1c_train),
#       posthba1c_train.m = mean(posthba1c_train),
#       se = sd(posthba1c_train)/sqrt((length(posthba1c_train))))

#obs vs pred, by decile of predicted treatment difference
#For Formula 1-3
mnumber = c(1:10)
models  <- as.list(1:10)

hba1c_diff.obs.unadj <- vector()
lower.unadj <- vector()
upper.unadj <- vector()
hba1c_diff.obs.sim <- vector()
lower.sim <- vector()
upper.sim <- vector() 
hba1c_diff.obs.adj <- vector()
lower.adj <- vector()
upper.adj <- vector() 

#Unadj
for(i in mnumber) {
  models[[i]] <- lm(as.formula(formula1),data=train,subset=hba1c_diff.q==i)
  hba1c_diff.obs.unadj <- append(hba1c_diff.obs.unadj,models[[i]]$coefficients[2])
  confint_all <- confint(models[[i]], levels=0.95)
  lower.unadj <- append(lower.unadj,confint_all[2,1])
  upper.unadj <- append(upper.unadj,confint_all[2,2])
}
#Simple 
for(i in mnumber) {
  models[[i]] <- lm(as.formula(formula2),data=train,subset=hba1c_diff.q==i)
  hba1c_diff.obs.sim <- append(hba1c_diff.obs.sim,models[[i]]$coefficients[2])
  confint_all <- confint(models[[i]], levels=0.95)
  lower.sim <- append(lower.sim,confint_all[2,1])
  upper.sim <- append(upper.sim,confint_all[2,2])
}
#Full
for(i in mnumber) {
  models[[i]] <- lm(as.formula(formula3),data=train,subset=hba1c_diff.q==i)
  hba1c_diff.obs.adj <- append(hba1c_diff.obs.adj,models[[i]]$coefficients[2])
  confint_all <- confint(models[[i]], levels=0.95)
  lower.adj <- append(lower.adj,confint_all[2,1])
  upper.adj <- append(upper.adj,confint_all[2,2])
}

#train data.frame  
t1 <- t1[1:10,]
t1 <- data.frame(t1,cbind(hba1c_diff.obs.unadj,lower.unadj,upper.unadj,
                          hba1c_diff.obs.sim,lower.sim,upper.sim,
                          hba1c_diff.obs.adj,lower.adj,upper.adj))

#Function to output HTE by subgroup
hte_plot <- function(data,pred,obs,obslowerci,obsupperci) {
  
  #ymin <- min(data$lci); ymax <- max(data$uci);yminr  <- 2*round(ymin/2);  ymaxr <- 2*round(ymax/2)
  ymin  <- -10;  ymax <- 10
  
  ggplot(data=data,aes_string(x=pred,y=obs)) +
    geom_point(alpha=1) + theme_bw() +
    geom_errorbar(aes_string(ymin=obslowerci, ymax=obsupperci), colour="black", width=.1) +
    ylab("Observed HbA1c difference (mmol/mol)") + xlab("Predicted HbA1c difference (mmol/mol)") +
    scale_x_continuous(limits=c(ymin,ymax),breaks=c(seq(ymin,ymax,by=2))) +
    scale_y_continuous(limits=c(ymin,ymax),breaks=c(seq(ymin,ymax,by=2))) +
    # scale_x_continuous(limits=c(ymin,ymax),breaks=c(seq(yminr,ymaxr,by=2))) +
    # scale_y_continuous(limits=c(ymin,ymax),breaks=c(seq(yminr,ymaxr,by=2))) +
    theme_base() + geom_abline(intercept=0,slope=1, color="red", lwd=0.75) + ggtitle("") +
    geom_vline(xintercept=0, linetype="dashed", color = "grey60") + geom_hline(yintercept=0, linetype="dashed", color = "grey60") 
}

#simple adj
plotdata <- t1 %>% dplyr::mutate(obs=hba1c_diff.obs.unadj,lci=lower.unadj,uci=upper.unadj)
hte_plot(plotdata,"hba1c_diff.pred","obs","lci","uci") 
#simple adj
plotdata <- t1 %>% dplyr::mutate(obs=hba1c_diff.obs.sim,lci=lower.sim,uci=upper.sim)
hte_plot(plotdata,"hba1c_diff.pred","obs","lci","uci") 
#splie adj
plotdata <- t1 %>% dplyr::mutate(obs=hba1c_diff.obs.adj,lci=lower.adj,uci=upper.adj)
hte_plot(plotdata,"hba1c_diff.pred","obs","lci","uci")  