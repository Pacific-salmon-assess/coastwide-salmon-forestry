d_chm_cpd=read.csv(here('stan models','outs','fits','posterior','bh_chm_cpd.csv'),check.names=F)

d_alpha_chm_eca=data.frame(d_chm_eca[,grepl('alpha_t',colnames(d_chm_eca))])
dalpha0_chm_eca=apply(d_alpha_chm_eca,1,median)

x_n_chm_eca=seq(min(ch20r$sqrt.ECA.std),max(ch20r$sqrt.ECA.std),length.out=100)
chm_eca_n=(x_n_chm_eca*sd(ch20r$sqrt.ECA)+mean(ch20r$sqrt.ECA))^2


exp(dalpha0_chm_eca[1]+d_chm_cpd$b_for[1]*x_n_chm_eca)/exp(dalpha0_chm_eca[1]+d_chm_cpd$b_for[1]*x_n_chm_eca[1])
exp(dalpha0_chm_eca[1]+d_chm_cpd$b_for[1]*x_n_chm_eca-log(1+(exp(dalpha0_chm_eca[1])/d_chm_cpd$`Rk[100]`[1])*100))/exp(dalpha0_chm_eca[1]+d_chm_cpd$b_for[1]*x_n_chm_eca[1]-log(1+(exp(dalpha0_chm_eca[1])/d_chm_cpd$`Rk[100]`[1])*100))
exp(dalpha0_chm_eca[1]+d_chm_cpd$b_for[1]*x_n_chm_eca-log(1+(exp(dalpha0_chm_eca[1])/d_chm_cpd$`Rk[100]`[1])*7000))/exp(dalpha0_chm_eca[1]+d_chm_cpd$b_for[1]*x_n_chm_eca[1]-log(1+(exp(dalpha0_chm_eca[1])/d_chm_cpd$`Rk[100]`[1])*7000))

alpha_t=d_chm_cpd[,grepl('alpha_t',colnames(d_chm_cpd))==T]

ch20r$tt=as.numeric(factor(ch20r$BroodYear))
for(j in 1:length(unique(ch20r$River))){
  x=subset(ch20r,River==unique(ch20r$River)[j])
  alpha_x=alpha_t[,x$tt]
  #logRS -> RS -> RS*S -> pred R
  x$pred_Reca=exp(eqn)*S
  x$pred_R0=exp(eqn)*S
  
  exp(dalpha0_chm_eca[1]+d_chm_cpd$b_for[1]*x_n_chm_eca-log(1+(exp(dalpha0_chm_eca[1])/d_chm_cpd$`Rk`[1])*7000))/exp(dalpha0_chm_eca[1]+d_chm_cpd$b_for[1]*x_n_chm_eca[1]-log(1+(exp(dalpha0_chm_eca[1])/d_chm_cpd$`Rk[100]`[1])*7000))
  
}


#list - length J for stocks... each matrix has N years columns, n iterations rows - name columns by brood year
#to_vector(matrix)-> predicted recruits; rep(x$broodyear,4000)
#long format... 4000 x N year rows... 2 columns...: predicted recrutis & broodyear

#do.call(rbind, list)
#order(broodyear)

#rlist 

#diff. in predicted recruits by brood year... matrix(4000 rows, n columns for stocks represented in that year)...apply(,1,sum) = total loss of recruits from forestry (4000 estimates - take median and confidence intervals)

#from that gather all columns together that are in the same year across stocks


#predicted Recruits - at each brood year, across all stocks in that year
