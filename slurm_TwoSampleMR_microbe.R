library(TwoSampleMR)
library(data.table)
library(R.utils)
##Utilized the selected IVs of gut microbiome(stroke for reverse MR).
##The outcome data should be fixed as the following column.
##SNP     A1      A2      freq    b       se      p
#rs4951859       C       G       0.1968  -0.0282 0.0132  0.03232
#rs58276399      T       C       0.8364  0.033   0.0127  0.009287
#rs58276399      T       C       0.8364  0.033   0.0127  0.009287
#rs141242758     T       C       0.8434  0.0289  0.0113  0.01099
#rs79010578      A       T       0.1607  -0.0244 0.0117  0.03738
#rs28527770      T       C       0.8299  0.0161  0.0113  0.1536
expose_path<-'/share/home/microbiome_GCTA/'
outcome_path<-'/share/home/stroke_EAS/'
results_path<-'/share/home/Gut_microbiome_stroke_EAS/'
IV_path<-"/share/home/microbiome_EAS/IV_selected_1e_5_001/"
IVfile<-list.files(IV_path,pattern = "iv$")
for(expos in IVfile){
  name=sub(".iv","",expose)
  cc <- paste0(IV_path,expos)   
  IV_exposure<-read.table(cc,header=T,sep="\t",stringsAsFactors = F)#读取IV文件
  IV_exposure<-subset(IV_exposure,Chr!="Chr")
  gwasresults_exposure<-fread(paste(expos_path,name,
                                     '_GCTA.assoc.raw',sep = ''),header=T,sep="\t",stringsAsFactors = F)
  gwasresults_exposure<-subset(gwasresults_exposure,SNP %in% IV_exposure[,2])
  gwasresults_exposure$p<-as.numeric(gwasresults_exposure$p)
  gwasresults_exposure$phenotype<-name
  exp_dat <- format_data(
    gwasresults_exposure,
    phenotype_col = "phenotype",
    type='exposure',
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col ="A1",
    other_allele_col = "A2",
    eaf_col = "freq",
    pval_col = "p")
  outcomefile<-list.files(outcome_path,pattern = ".summary.gwas.format.gz")
  for(outfile in outcomefile){
    outname=sub(".summary.gwas.format.gz","",outfile)
    #outname=sub("summary_file_","",outname)
    if(name==outname){next}else{
      outcc <- paste(outcome_path,outname,'.summary.gwas.format.gz',sep="")
      gwasresults_outcome<-fread(outcc, header=TRUE, stringsAsFactors=FALSE, sep="\t")
      gwasresults_outcome$phenotype<-outname
      gwasresults_outcome<-subset(gwasresults_outcome,p_value!='NA')
      snps<-subset(gwasresults_outcome,gwasresults_outcome$SNP %in% exp_dat$SNP)
      if(length(snps$SNP)==0){next}else{
        out_dat <- format_data(
          dat=gwasresults_outcome,
          type = "outcome",
          snps = exp_dat$SNP,
          header = TRUE,
          phenotype_col = "phenotype",
          snp_col = "SNP",
          beta_col = "b",
          se_col = "se",
          effect_allele_col = "A1",
          other_allele_col = "A2",
          eaf_col = "freq",
          pval_col = "p")
        
        #harmonise_data
        mydata <- harmonise_data(
          exposure_dat=exp_dat,
          outcome_dat=out_dat,
          action= 1)
        write.table(mydata,file = paste(results_path,'Gut_microbiome_stroke_EAS_1e_5_001_mydata.txt',sep=''),append=T,row.names=F,quote = F,sep = ',')
        steiger<-mr_steiger(p_exp=mydata$pval.exposure,p_out=mydata$pval.outcome,gwasresults_exposure$N[1],gwasresults_outcome$N[1],r_exp = mydata$beta.exposure,r_out = mydata$beta.outcome)
        s<-data.frame(steiger[1:12])
        s$exposure<-name
        s$outcome<-outname
        write.table(s,file=paste(results_path,'Gut_microbiome_stroke_EAS_1e_5_001_steiger.csv',sep = ''),quote = F,append = T,sep=',',row.names = F)
        #run mr
        res <- mr(mydata)
        if(length(res)==0){next}else{
          write.table(res,file=paste(results_path,'Gut_microbiome_stroke_EAS_1e_5_001_MR.csv',sep = ''),quote = F,append = T,sep=',',row.names = F)
          if(res$nsnp[1]>1){
            het <- mr_heterogeneity(mydata)
            colnames(het)[5]<-'method_het'
            if(length(het$Q_pval)==2&het$Q_pval[1]>het$Q_pval[2]){het<-het[2,]}else{het<-het[1,]}
            res1<-merge(res,het,all.x = T,by=c('id.exposure','id.outcome','outcome','exposure'))
            pleio <- mr_pleiotropy_test(mydata)
            colnames(pleio)[6:7]<-c('pleio_se','pleio_p')
            res1<-merge(res1,pleio,all.x = T,by=c('id.exposure','id.outcome','outcome','exposure'))
            write.table(res1,file=paste(results_path,'Gut_microbiome_stroke_EAS_1e_5_001_snesitivity_analysis.csv',sep = ''),quote = F,append = T,sep = ',',row.names = F)
            pdf(file = paste(results_path,'/',name,'_',outname,'.pdf',sep=''))
            single <- mr_leaveoneout(mydata)
            print(mr_leaveoneout_plot(single))
            print(mr_scatter_plot(res,mydata))
            res_single <- mr_singlesnp(mydata)
            print(mr_forest_plot(res_single))
            print(mr_funnel_plot(res_single))
            dev.off()
          }else{write.table(res,file=paste(results_path,'Gut_microbiome_stroke_EAS_1e_6_001_Wald_Ratio.csv',sep = ''),quote = F,append = T,sep = ',',row.names = F)}}}}}}
