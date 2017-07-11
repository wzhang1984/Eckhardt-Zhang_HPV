

d=read.delim('../data/glmBinomial_anovaChisq/df4glm_HNSC_CESC.txt',header=T,row.names=1)
og_tsg_df=read.delim('/cellar/data/users/wzhang1984/forNBS/oncogene_tsg.txt',header=F,row.names=1)

out=as.data.frame(matrix(ncol=8))
i=1
for (gene in colnames(d)[1:(length(colnames(d))-3)]) {
    HPV_p_vector=d[which(d$HPV==1),gene]
    HPV_n_vector=d[which(d$HPV==0),gene]
    rate_p=sum(HPV_p_vector,na.rm=T)/sum(1-is.na(HPV_p_vector))
    rate_n=sum(HPV_n_vector,na.rm=T)/sum(1-is.na(HPV_n_vector))
    if (rate_p==0 || rate_n==0) {next}
    fit=glm(d[,gene]~d[,'disease']+d[,'Chicago']+d[,'HPV'],family='binomial')
    fit0=glm(d[,gene]~d[,'disease']+d[,'Chicago'],family='binomial')
    anovaChisq=anova(fit0,fit,test='Chisq')
    Deviance=anovaChisq$Deviance[length(anovaChisq$Deviance)]
    pval=anovaChisq[,'Pr(>Chi)'][length(anovaChisq[,'Pr(>Chi)'])]
    score=0
    if (fit$coefficients[4]<0) {score=Deviance} else {score=0}
    mut_rate=sum(d[,gene],na.rm=T)/sum(1-is.na(d[,gene]))
    mut_type='loss'
    if (og_tsg_df[gene,'V2'] %in% c('Oncogene','Amplification_Oncogene')) {mut_type='gain'}
    out[i,1]=gene
    out[i,2]=mut_type
    out[i,3]=pval
    out[i,4]=rate_p
    out[i,5]=rate_n
    out[i,6]=mut_rate
    out[i,7]=-fit$coefficients[4]
    out[i,8]=score
    i=i+1
    print(c(gene,score))
}

write.table(out,'../data/glmBinomial_anovaChisq/glmBinomialAnovaChisq_HNSC_CESC_2cov_2.txt',sep="\t",row.name=F, quote =F,col.name=F)
