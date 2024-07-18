d_chm_cpd=read.csv(here('stan models','outs','fits','posterior','bh_chm_cpd.csv'),check.names=F)

d_alpha_chm_eca=data.frame(d_chm_eca[,grepl('alpha_t',colnames(d_chm_eca))])
dalpha0_chm_eca=apply(d_alpha_chm_eca,1,median)

x_n_chm_eca=seq(min(ch20r$sqrt.ECA.std),max(ch20r$sqrt.ECA.std),length.out=100)

exp(dalpha0_chm_eca[1]+d_chm_cpd$b_for[1]*x_n_chm_eca)/exp(dalpha0_chm_eca[1]+d_chm_cpd$b_for[1]*x_n_chm_eca[1])
exp(dalpha0_chm_eca[1]+d_chm_cpd$b_for[1]*x_n_chm_eca-log(1+(exp(dalpha0_chm_eca[1])/d_chm_cpd$`Rk[100]`[1])*100))/exp(dalpha0_chm_eca[1]+d_chm_cpd$b_for[1]*x_n_chm_eca[1]-log(1+(exp(dalpha0_chm_eca[1])/d_chm_cpd$`Rk[100]`[1])*100))
exp(dalpha0_chm_eca[1]+d_chm_cpd$b_for[1]*x_n_chm_eca-log(1+(exp(dalpha0_chm_eca[1])/d_chm_cpd$`Rk[100]`[1])*7000))/exp(dalpha0_chm_eca[1]+d_chm_cpd$b_for[1]*x_n_chm_eca[1]-log(1+(exp(dalpha0_chm_eca[1])/d_chm_cpd$`Rk[100]`[1])*7000))


