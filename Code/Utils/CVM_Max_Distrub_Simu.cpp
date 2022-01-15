#include <Rcpp.h>
using namespace Rcpp;
#include <math.h>
#include <set>
#include <iostream>
#include<map>
#include <vector>
#include <limits>
#include<algorithm>
#include<ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>

using namespace std;
using namespace __gnu_pbds;

typedef tree<int,
             null_type,
             less<int>,
             rb_tree_tag,
             tree_order_statistics_node_update>
  order_stat_tree;

int getSum(vector<int> BITree, int index)
{
  int sum = 0; // Initialize result
  
  // index in BITree[] is 1 more than the index in arr[]
  index = index + 1;
  
  // Traverse ancestors of BITree[index]
  while (index>0)
  {
    // Add current element of BITree to sum
    sum += BITree[index];
    
    // Move index to parent node in getSum View
    index -= index & (-index);
  }
  return sum;
}

void updateBIT(vector<int> &BITree, int n, int index, int val)
{
  // index in BITree[] is 1 more than the index in arr[]
  index = index + 1;
  
  // Traverse all ancestors and add 'val'
  while (index <= n)
  {
    // Add 'val' to current node of BI Tree
    BITree[index] += val;
    
    // Update index to that of parent in update View
    index += index & (-index);
  }
}

vector<int> constructBITree(int n)
{
  // Create and initialize BITree[] as 0
  vector<int> BITree(n+1);
  
  // Store the actual values in BITree[] using update()
  for (int i=0; i < n; i++)
  {
    updateBIT(BITree, n, i, i + 1);
  }
  
  return BITree;
}

void rank_transform_int(vector<int> &input)
{
  // create an empty ordered map
  map<int, int> map;
  
  // store (element, index) pair in a map
  for (unsigned int i = 0; i < input.size(); i++)
  {
    map[input[i]] = i;
  }
  
  // keys are stored in sorted order in an ordered map
  
  // rank starts from 1
  int rank = 1;
  
  // iterate through the map and replace each element with its rank
  std::map<int, int>::iterator i;
  
  for(i = map.begin(); i != map.end(); i++)
  {
    input[(*i).second] = rank++;
  }
}

vector<int> rank_transform_double(vector<double> input)
{
  // create an empty ordered map
  map<double, unsigned int> temp_map;
  
  // store (element, index) pair in a map
  for (unsigned int i = 0; i < input.size(); i++) {
    temp_map[input[i]] = i;
  }
  
  // keys are stored in sorted order in an ordered map
  
  // rank starts from 1
  int rank = 1;
  vector<int> output(input.size());
  
  // iterate through the map and replace each element with its rank
  
  std::map<double, unsigned int>::iterator i;

  for(i = temp_map.begin(); i != temp_map.end(); i++)
  {
    output[(*i).second] = rank++;
  }
  
  return output;
}

NumericVector changepoint_scan(const vector<int> time_series)
{
  double cp_index = 2; 
  double cp_cost = 0;
  double data_length = time_series.size(); 
  double right_length = data_length - 1; 
  double cp_mean = (data_length + 1) / (6 * data_length);
  
  NumericVector statistics(data_length - 3);
  
  vector<int> left_fenwick_tree(data_length + 1);
  updateBIT(left_fenwick_tree, data_length, time_series[0], time_series[0]);
  
  vector<int> right_fenwick_tree = constructBITree(data_length + 1);
  updateBIT(right_fenwick_tree, data_length, time_series[0], -1 * time_series[0]);
  
  int left_sum_cost = pow(time_series[0] - 1, 2);
  int right_sum_cost = (data_length - time_series[0]); 
  
  order_stat_tree right_order_tree;
  
  for(int i = 1; i < data_length; i++)
  {
    right_order_tree.insert(time_series[i]);
  }
  
  double cp_variance = 0;
  int local_rank = 0;
  int global_rank = 0;
  
  for(; cp_index != (data_length - 1); cp_index += 1)
  {
    right_length -= 1;
    global_rank = time_series[cp_index - 1];
    
    local_rank = right_order_tree.order_of_key(global_rank);
    right_order_tree.erase(global_rank);
    
    //Computing sum recurrence
    left_sum_cost += pow(local_rank, 2) + cp_index - (global_rank - local_rank) + (cp_index - 1) * (cp_index) - (global_rank - (local_rank + 1)) * (global_rank - local_rank) - 2 *(getSum(left_fenwick_tree, data_length - 1) - getSum(left_fenwick_tree, global_rank - 1));
    right_sum_cost +=  -1 * pow(global_rank - (local_rank + 1), 2) + right_length - (local_rank) + -1 * ((right_length + 1 )* (right_length + 2) - (local_rank + 1) * (local_rank + 2)) + 2 * (getSum(right_fenwick_tree, data_length - 1) - getSum(right_fenwick_tree, global_rank - 1));
    
    //Two-sample CVM
    cp_cost = ((cp_index * left_sum_cost + right_length * right_sum_cost) / (data_length * cp_index * right_length)) - (4 * cp_index * right_length - 1)/(6 * data_length);
    
    //Standardising CVM
    cp_variance = ((1 - 3 / (4 * cp_index)) * pow(data_length, 2) + (1 - cp_index) * data_length - cp_index) * (data_length + 1) / (45 * pow(data_length, 2) * (data_length - cp_index));
    cp_variance = sqrt(cp_variance);

    cp_cost = (cp_cost - cp_mean) / cp_variance;
    
    updateBIT(left_fenwick_tree, data_length, global_rank - 1, global_rank);
    updateBIT(right_fenwick_tree, data_length, global_rank - 1, -1 * (global_rank));
    
    // output
    statistics[cp_index - 2] = cp_cost;
    
  }
  return statistics;
}

// [[Rcpp::export]]

NumericVector CVM_stats(vector<double> input_time_series)
{
  order_stat_tree global_order_tree;
  vector<int> time_series = rank_transform_double(input_time_series);
  
  return changepoint_scan(time_series);
}


/*** R

data_lengths = c(25, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 25000)
num_reps = 50
num_obs = 10000

norm_cvm = array(0, c(length(data_lengths), num_reps, num_obs))

setwd("C:\\Users\\karim\\OneDrive\\Desktop\\Mathematics\\Changepoint Research\\CVM-BS\\Code\\Simulations\\p-values sims")

for(i in 1 : length(data_lengths))
{
  t1 = Sys.time()
  
  for(j in 1 : num_reps)
  {
   for(k in 1 : num_obs)
   {
     
     #distribution params
     mean = runif(1, -5, 5)
     var = runif(1, 0, 10)
     
     #data instance
     data = rnorm(data_lengths[i], mean, var)
    
     #calling algo for the stats 
     output_stats = CVM_stats(data) 
     cvm_max = max(output_stats)
     
     #storing the stat to our estimated distribution
     norm_cvm[i,j,k] = cvm_max
     
   } 
  }
  
  t2 = Sys.time()
  
  print('Data length:')
  print(i)
  print('Iter Time:')
  print(t2-t1)
  
}

norm_file = paste('CVM MAX Normal Distrub.RData')
save(norm_cvm, file = norm_file)

*/
