

qreg_cdf <- function(model_list, x, m)
{
  #Evaluates F_m(x) where F_m is the CDF of data length m implied by quantile regression
  #model_list needs to be sorted in ascending order of quantiles
  
  # For values of x matching an existing quantile, we return that quantile
  # For values of x within two quantiles we return a linear interpolation
  # For values of x more extreme than all quantiles we assign prob 0 or 1 
  
  pred_input = data.frame(n = m)
  
  if(length(model_list) == 1)
  {
    stop("We require more quantile models to form a sensible CDF")
  } 
  
  if(length(model_list) == 2)
  {
    model_l = model_list[[1]]
    model_r = model_list[[2]]
    
    pred_l = predict(model_l, newdata = pred_input)[[1]]
    pred_r = predict(model_r, newdata = pred_input)[[1]]
  
    if(pred_l == x){return(model_l$tau)}
    else if(pred_r == x){return(model_r$tau)} 
    else if(x < pred_l){return(0)}
    else if(x > pred_r){return(1)}
    else 
    {
      weight = (x - pred_l) / (pred_r - pred_l)
      prob = (1 - weight) * model_l$tau + weight * model_r$tau
      return(prob)
    }
  }
  
  else
  {
    index = floor(mean(c(1, length(model_list)))) 
    model = model_list[[index]]
    pred = predict(model, newdata = pred_input)[[1]]

    if(pred == x){return(model$tau)}
    else if(x < pred){return(qreg_cdf(model_list[1 : index], x, m))}
    else if(x > pred){return(qreg_cdf(model_list[index : length(model_list)], x, m))}
  }
  
  return(p_val)
}

