library(quantreg)
library(splines2)

mass_train <- function(train_data, population, quantile)
{
  ##fits a quantile reg model to each set of knots in the new population
  models = list()
  
  for(i in 1 : length(population))
  {
    models[[i]] = rq(response ~ bSpline(n, knots = population[[i]]), tau = quantile, data = train_data)
  }
  
  return(models)
}

model_error <- function(model, test_data, quantile)
{
  #mean quantile absolute error
  y_pred = predict(model, newdata = test_data) 
  y_true = test_data$response
  y_diff = y_true - y_pred
  
  qae = quantile * sum(y_diff[y_diff>0])
  qae = qae - (1 - quantile) * sum(y_diff[y_diff<0])
  
  return(qae)
}

mass_fitness <- function(mass_models, train_data, test_data, quantile)
{
  #Computes the Pseudo-Rsq statistic from Koenker and Machado (1999) for all of the models trained
  #this is 1 - the ratio of the quantile error for each model divided by the quantile error of the constant model
  
   fitness = c()
  
   constant_model = rq(response ~ 1, tau = quantile, data = train_data)
   constant_error = model_error(constant_model, test_data, quantile)  
    
   for(i in 1 : length(mass_models))
   {
     temp_fit = 1 - (model_error(mass_models[[i]], test_data, quantile) / constant_error)
     fitness = c(fitness, temp_fit)
   }
   return(fitness)
}

