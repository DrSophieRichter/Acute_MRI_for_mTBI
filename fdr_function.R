fdr_maxp <- function(my_p_values, threshold) 
{
  #Collect all p-values and count how many there are
  p <- my_p_values
  n <- length(p)
  
  #sort and rank p-values and collect this info in a dataframe fd (false discovery)
  p <- sort(p)
  fd <- as.data.frame(p)
  fd$i <-order(fd$p)
  
  #Set maximum false discovery rate 
  d <- threshold
  
  #Calculate d*i/n
  fd$din <- (d*fd$i/n) %>% round(5)
  
  #Decide whether p <= din
  fd$significant <- ifelse(fd$p <= fd$din, "yes", "no")
  
  #Extract highest significant p-value
  yes <- fd %>% filter(significant == "yes") 
  maxp <- ifelse(nrow(yes) == 0, 0, max(yes$p))
  
  return(maxp)
}