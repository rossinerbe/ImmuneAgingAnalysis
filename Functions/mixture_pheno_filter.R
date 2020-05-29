#finds pheno samples which match desired phenotype inputs
#uses TCGA_ID to take only these samples from the mixture_data file
#args
#mixture_data - a result from a CIBERSORT run on RNAseq data
#pheno - phenotype data downloaded from cBioPortal
#age_range - range of ages to filter data by (e.g. c(70,80))
#cancer_types - vector of cancer_types to filter by using TCGA acronyms (e.g. c(BRCA, GBM))
#returns filtered mixture_data with attached phenotype data
mixture_pheno_filter = function(mixture_data, pheno, age_range = NULL,
                                cancer_types = NULL) {
  
  #if an age range is given, select phenotype data indices within the input range at diagnosis
  if(!is.null(age_range)) {
    age_matches = which(pheno$`Diagnosis Age`> age_range[1] & pheno$`Diagnosis Age` < age_range[2])
  }
  
  #if a list of cancer types is input
  if(!is.null(cancer_types)) {
    cancer_matches = c()
    #iterate through list of cancer types finding matches for each within phenotype data
    for(i in 1:length(cancer_types)) {
      tmpinds = which(pheno$`TCGA PanCanAtlas Cancer Type Acronym` == cancer_types[i])
      cancer_matches = c(cancer_matches, tmpinds)
      
      #find the caseID for phenotype matches
      if(!is.null(age_range)) {
        matches_either = c(age_matches, cancer_matches)
        matches_both = matches_either[duplicated(matches_either)]
        match_IDs = pheno$`Patient ID`[matches_both]
      }
      else{match_IDs = pheno$`Patient ID`[cancer_matches]}
    }
  }
  else {match_IDs = pheno$`Patient ID`[age_matches]}
  
  #find intersect of caseIDs between phenotype data and CIBERSORT data
  mixture_IDs = rownames(mixture_data)
  mixture_sub_IDs = intersect(match_IDs, mixture_IDs)
  #subset CIBERSORT data by CaseID
  mixture_sub = mixture_data[which(mixture_IDs %in% mixture_sub_IDs),]
  
  #add phenotype data to cibersort data frame
  mixture_pheno = data.frame()
  for(i in 1:nrow(mixture_sub)){
    matched_pheno_ind = which(pheno$`Patient ID` == rownames(mixture_sub)[i])
    if(length(matched_pheno_ind) == 0) {next}
    if(length(matched_pheno_ind) > 1) {next}
    cp_row = as.data.frame(c(mixture_sub[i,], pheno[matched_pheno_ind,]))
    mixture_pheno = rbind(mixture_pheno, cp_row)
  }
  
  
  #if NA for any survival data, remove rows
  if(sum(is.na(mixture_pheno$Disease.specific.Survival.status)) != 0) {
    mixture_pheno = mixture_pheno[-which(is.na(mixture_pheno$Disease.specific.Survival.status)),]
  }
  if(sum(is.na(mixture_pheno$Months.of.disease.specific.survival)) != 0) {
    mixture_pheno = mixture_pheno[-which(is.na(mixture_pheno$Months.of.disease.specific.survival)),]
  }
  
  
  #return merged subset of CIBERSORT and phenotype data
  return(mixture_pheno)
}
