coxph_model_mixture = function(mixture_data, time, event, covars = NULL, checkInteraction = FALSE, order = TRUE, startcolumn = 1, endcolumn = 22, strat = FALSE) {
  
  #necesssary libraries
  require(survival)
  require(broom)
  
  #create data frame for summary data from coxph model
  totalcolumns = endcolumn - startcolumn + 1
  immune_cells = colnames(mixture_data[,startcolumn:endcolumn])
  model_pvalue <- data.frame(matrix(ncol = 6, nrow = totalcolumns))
  colnames(model_pvalue) <- c("Cell_Type", "estimate", "std.error", "statistic", "p.value", "exp")
  model_pvalue[,1] <- colnames(mixture_data[,startcolumn:endcolumn])
  #if no covars
  if(is.null(covars)) {
    #iteratively run coxph model over all immune cell types
    for(i in 1:totalcolumns) {
      temp_form <- as.formula(paste("Surv(", time, ", event = ", event, ") ~", immune_cells[i]))
      temp_model <- coxph(temp_form, data = mixture_data)
      temp_model <- tidy(temp_model)
      #collect summary data
      model_pvalue[i, 2:5] <- temp_model[1, 2:5]
      model_pvalue[i, 6] = exp(model_pvalue[i,2])
    }
  }
  #if covars are given
  else {
    if(checkInteraction == TRUE) {
    #iteratively run coxph model over all immune cell types
    for(i in 1:totalcolumns) {
      temp_form <- as.formula(paste("Surv(", time, ", event = ", event, ") ~", immune_cells[i], "*", covars))
      temp_model <- coxph(temp_form, data = mixture_data)
      temp_model <- tidy(temp_model)
      if(temp_model[3,5] < 0.05) {
        print(paste("Interaction term significant for:", immune_cells[i]), quote = F)
      }
      else {
        temp_form <- as.formula(paste("Surv(", time, ", event = ", event, ") ~", immune_cells[i], "+", covars))
        temp_model <- coxph(temp_form, data = mixture_data)
        temp_model <- tidy(temp_model)
      }
      #collect summary data
      model_pvalue[i, 2:5] <- temp_model[1, 2:5]
      model_pvalue[i, 6] = exp(model_pvalue[i,2])
    }
    }
    else {
      for(i in 1:totalcolumns) {
        if(strat == TRUE) {
          temp_form <- as.formula(paste("Surv(", time, ", event = ", event, ") ~", immune_cells[i], "+", "strata(", covars, ")"))
        }
        else{
          temp_form <- as.formula(paste("Surv(", time, ", event = ", event, ") ~", immune_cells[i], "+", covars))
        }
        temp_model <- coxph(temp_form, data = mixture_data)
        temp_model <- tidy(temp_model)
        #collect summary data
        model_pvalue[i, 2:5] <- temp_model[1, 2:5]
        model_pvalue[i, 6] = exp(model_pvalue[i,2])
      }
    }
  }
  
  #adjust p-values by bonferroni and benjamini-hochberg
  bonf_adj_significant<- model_pvalue$Cell_Type[model_pvalue$p.value < 0.05/totalcolumns]
  bh_adjusted_sig <- model_pvalue$Cell_Type[p.adjust(model_pvalue$p.value, method = "BH") < 0.05]
  #ordered list of immune cells by p-value
  ordered_model <- model_pvalue[order(model_pvalue$p.value),]
  
  if(order == TRUE) {
    return(list("Bonferroni Significant" = bonf_adj_significant, "BH FDR Significant" = bh_adjusted_sig, "Ordered Cell Types" = ordered_model))
  }
  else {
    return(list("Bonferroni Significant" = bonf_adj_significant, "BH FDR Significant" = bh_adjusted_sig, "Cell Types" = model_pvalue))
  }
}
