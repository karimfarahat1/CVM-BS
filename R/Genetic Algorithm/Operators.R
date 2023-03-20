library(quantreg)
library(splines2)

breed <- function(current_population, fitness_levels, cross_over_rate, mutation_rate)
{
  n = length(current_population)
  num_cross_overs = rbinom(1, n, cross_over_rate)
  num_elite = n - num_cross_overs
  
  #Elite selection
  sorted_fitness = sort(fitness_levels, decreasing = T)
  elite_fitness = sorted_fitness[1 : num_elite]
  elite_knots = current_population(which(fitness_levels) == elite_fitness)
  new_population = elite_knots
  
  #Cross-overs
  prob_vec = fitness_levels / sum(fitness_levels)
  
  for(i in 1 : num_cross_overs)
  {
    parent_indices = sample(x = 1 : n, size = 2, prob = prob_vec)
    
    p1 = current_population[[parent_indices[1]]]
    l1 = length(p1)
    p2 = current_population[[parent_indices[2]]]
    l2 = length(p2)
    
    if(l1>=l2)
    {
      lp = p1
      sp = p2
      lf = fitness_levels[parent_indices[1]]
      sf = fitness_levels[parent_indices[2]]  
    }
    
    else
    {
      lp = p2
      sp = p1
      lf = fitness_levels[parent_indices[2]]
      sf = fitness_levels[parent_indices[1]]
    }
    
    child = cross_over(long_parent = lp, short_parent = sp, long_fitness = lf, short_fitness = sf)
    
    #adding new element to list
    new_population[[length(new_population) + 1]] = child
    
  }
  
  #Mutations
  num_mutations = rbinom(1, n, mutation_rate)
  indices = sample(x = 1 : n, size = num_mutations)
  
  for(i in indices)
  {
    new_population[[i]] = mutate(new_population[[i]])
  }
   
}

cross_over <- function(long_parent, short_parent, long_fitness, short_fitness)
{
  n = length(long_parent)
  m = length(short_parent)
  
  indicators = rbinom(n, 1, long_fitness / (short_fitness + long_fitness))
  
  child = c()
  
  for(i in 1:n)
  {
    indic = indicators[i]
    
    if(i<=m)
    {
      if(indic == 1)
      {
        child = c(child, long_parent[i])  
      }
      
      else
      {
        child = c(child, short_parent[i])
      }
    }
    
    else
    {
      if(indic == 1)
      {
        child = c(child, long_parent[i])
      }
    }
  }
}

mutate <- function(knots, max_n, mut_prob)
{
  indicator = rbinom(1, 1, mut_prob)
  
  if(indicator == 1)
  {
    new_knot  = sample(1 : max_n)
    return(c(knots, new_knot))  
  }
  
  else
  {
    delete_index = sample(1 : length(knots))
    knots = c(knots[1 : delete_index - 1], knots[1 : delete_index + 1])
   
    return(knots)
  }
}

