#include <iostream>
#include <fstream>
#define ARMA_DONT_USE_STD_MUTEX
#include <armadillo>
#include <math.h>
#include <set>
#include <map>
#include <vector>
#include <limits>
#include <algorithm>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>

using namespace __gnu_pbds;

typedef tree<int,
             null_type,
             std::less<int>,
             rb_tree_tag,
             tree_order_statistics_node_update>
  order_stat_tree;


int getSum(std::vector<int> BITree, int index)
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

void updateBIT(std::vector<int> &BITree, int n, int index, int val)
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

std::vector<int> constructBITree(int n)
{
  // Create and initialize BITree[] as 0
  std::vector<int> BITree(n+1);

  // Store the actual values in BITree[] using update()
  for (int i=0; i < n; i++)
  {
    updateBIT(BITree, n, i, i + 1);
  }

  return BITree;
}

void rank_transform_int(std::vector<int> &input)
{
  // create an empty ordered map
  std::map<int, int> map;

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

std::vector<int> rank_transform(std::vector<double> input)
{
  // create an empty ordered map
  std::map<double, unsigned int> temp_map;

  // store (element, index) pair in a map
  for (unsigned int i = 0; i < input.size(); i++) {
    temp_map[input[i]] = i;
  }

  // keys are stored in sorted order in an ordered map

  // rank starts from 1
  int rank = 1;
  std::vector<int> output(input.size());

  // iterate through the map and replace each element with its rank

  std::map<double, unsigned int>::iterator i;

  for(i = temp_map.begin(); i != temp_map.end(); i++)
  {
    output[(*i).second] = rank++;
  }

  return output;
}


double qreg_pred(arma::mat knots, arma::mat beta, double intercept, double m)
{
  //Implements De-Boor's algorithm to efficiently compute predictions using B-Splines

  arma::vec d = arma::zeros(beta.n_rows);
  int deg = 3; // assuming degree 3 splines are being used

  knots = knots.elem(find(knots > 0)); // removing the 0's added to equalise df of all splines
  beta = beta.elem(find(beta != 0)); // removing the 0's added to equalise the df of all splines

  arma::vec pad_knot = arma::zeros(knots.n_rows + 2 * (deg));

    //boundaries
  pad_knot[0] = 20;
  pad_knot[1] = 20;
  pad_knot[2] = 20;
  pad_knot[knots.n_rows + 3] = 3000;
  pad_knot[knots.n_rows + 4] = 3000;
  pad_knot[knots.n_rows + 5] = 3000;

  for(unsigned int i = 0; i < knots.n_rows; i++)
  {
    pad_knot[i + 3] = knots(i,0);
  }

  unsigned int k = std::upper_bound(pad_knot.begin(), pad_knot.end(), m) - pad_knot.begin() - 1; // index of the pair [knot_k, knot_k+1] containing x

  if(k < 2){k = 2;}
  else if(k > pad_knot.size() - 4){k = pad_knot.size() - 4;}

  if(k > 2){d[0] = beta(k - deg, 0);}

  for(int i = 1; i <= deg; i++)
  {
    d[i] = beta(i + k - deg, 0);
  }

  double alpha;

  for(int i = 1; i <= deg; i++)
  {
    for(int j = deg; j >= i; --j)
    {
      alpha = (m - pad_knot(j + k - deg, 0)) / (pad_knot(j + 1 + k - i, 0) - pad_knot(j + k - deg, 0));
      d[j] = (1 - alpha) * d[j - 1] + alpha * d[j];
    }
  }

  return(intercept + d[deg]);
}

double qreg_cdf(arma::mat knots, arma::mat beta, std::vector<double> intercepts, std::vector<double> q, double x, double m)
{

  if(knots.n_rows == 2)
  {
    arma::mat beta_l = beta.row(0);
    arma::mat beta_r = beta.row(1);
    arma::mat knots_l = knots.row(0);
    arma::mat knots_r = knots.row(1);
    double q_l = q[0];
    double q_r = q[1];
    double pred_l = qreg_pred(knots_l, beta_l, intercepts[0], m);
    double pred_r = qreg_pred(knots_r, beta_r, intercepts[1], m);

    if(pred_l == x){return q_l;}
    else if(pred_r == x){return q_r;}
    else if(x < pred_l){return 0;}
    else if(x > pred_r){return 1;}
    else
    {
      double weight = (x - pred_l) / (pred_r - pred_l);
      double prob = (1 - weight) * q_l + weight * q_r;
      return prob;
    }
  }

  else
  {

    int index = (knots.n_rows - 1) / 2;
    double pred = qreg_pred(knots.row(index), beta.row(index), intercepts[index], m);

    if(pred == x){return q[index];}

    else if(x < pred){return qreg_cdf(knots.rows(0, index),
            beta.rows(0, index),
            std::vector<double>(intercepts.begin(), intercepts.begin() + index + 1),
            std::vector<double>(q.begin(), q.begin() + index + 1),
            x,
            m);}

    else if(x > pred){return qreg_cdf(knots.rows(index, knots.n_rows - 1),
            beta.rows(index, beta.n_rows - 1),
            std::vector<double>(intercepts.begin() + index, intercepts.end()),
            std::vector<double>(q.begin() + index, q.end()),
            x,
            m);}
  }
}

std::pair<double,int> max_cost_cp(const std::vector<int> time_series)
{
  double cp_index = 2;
  double cp_cost = 0;
  double data_length = time_series.size();
  double right_length = data_length - 1;
  double cvm_mean = (data_length + 1) / (6 * data_length); // per Anderson's paper

  std::pair<double, int> output;

  //Initialising trees & changepoint cost
  std::vector<int> left_fenwick_tree(data_length + 1);
  updateBIT(left_fenwick_tree, data_length, time_series[0], time_series[0]);

  std::vector<int> right_fenwick_tree = constructBITree(data_length + 1);
  updateBIT(right_fenwick_tree, data_length, time_series[0], -1 * time_series[0]);

  int left_sum_cost = pow(time_series[0] - 1, 2);
  int right_sum_cost = (data_length - time_series[0]);

  order_stat_tree order_tree;

  for(int i = 1; i < data_length; i++)
  {
    order_tree.insert(time_series[i]);
  }

  double cvm_std = 0;
  int local_rank = 0;
  int global_rank = 0;

  for(; cp_index != (data_length - 1); cp_index += 1)
  {
    right_length -= 1;
    global_rank = time_series[cp_index - 1];

    local_rank = order_tree.order_of_key(global_rank);
    order_tree.erase(global_rank);

    //Computing sum recurrence
    left_sum_cost += pow(local_rank, 2) + cp_index - (global_rank - local_rank) + (cp_index - 1) * (cp_index) - (global_rank - (local_rank + 1)) * (global_rank - local_rank) - 2 *(getSum(left_fenwick_tree, data_length - 1) - getSum(left_fenwick_tree, global_rank - 1));
    right_sum_cost +=  -1 * pow(global_rank - (local_rank + 1), 2) + right_length - (local_rank) + -1 * ((right_length + 1 )* (right_length + 2) - (local_rank + 1) * (local_rank + 2)) + 2 * (getSum(right_fenwick_tree, data_length - 1) - getSum(right_fenwick_tree, global_rank - 1));

    //Two-sample CVM
    cp_cost = ((cp_index * left_sum_cost + right_length * right_sum_cost) / (data_length * cp_index * right_length)) - (4 * cp_index * right_length - 1)/(6 * data_length);

    //Standardising the CVM: variance expression per Anderson's paper
    cvm_std = ((1 - 3 / (4 * cp_index)) * pow(data_length, 2) + (1 - cp_index) * data_length - cp_index) * (data_length + 1) / (45 * pow(data_length, 2) * (data_length - cp_index));
    cvm_std = sqrt(cvm_std);

    //studentising the CVM
    cp_cost = (cp_cost - cvm_mean) / cvm_std;

    //storing output
    if(cp_cost > output.first)
    {
      output.first = cp_cost;
      output.second = cp_index;
    }

    //updating trees for next iter
    updateBIT(left_fenwick_tree, data_length, global_rank - 1, global_rank);
    updateBIT(right_fenwick_tree, data_length, global_rank - 1, -1 * (global_rank));

  }

  return output;
}

bool fdr_control(std::map<double, int> changepoints, double fdr, double new_signif)
{
  double num_cps = changepoints.size() + 1;

  double pval;
  int cp_count = 0;

  double denom = num_cps * (log(num_cps) + 0.57721 + 1/(2 * num_cps)); //Benjamini-Yekutieli: harmonic series approx. by log



  for(std::map<double, int>::iterator it = changepoints.begin(); it != changepoints.end(); it++)
  {
    pval = (*it).first;
    cp_count++;
    if(pval > (cp_count/ denom) * fdr){return false;}
  }

  if(new_signif <= (num_cps / denom) *fdr){return true;}
  else{return false;}
}


void CVM_BinSeg(const std::vector<int> time_series,
                std::map<double, int> &changepoints,
                const int changepoint_offset,
                const double fdr,
                const arma::mat knots,
                const arma::mat beta,
                const std::vector<double> intercepts,
                const std::vector<double> q)
{
  std::pair<double, int> cp_vec = max_cost_cp(time_series);

  cp_vec.first = 1 - qreg_cdf(knots, beta, intercepts, q, cp_vec.first, time_series.size());

  if(fdr_control(changepoints, fdr, cp_vec.first))
  {

    cp_vec.second = cp_vec.second + changepoint_offset;
    changepoints.insert(cp_vec);

    if(cp_vec.second - changepoint_offset > 10)
    {
    std::vector<int> left_series = std::vector<int>(time_series.begin(), time_series.begin() + cp_vec.second - changepoint_offset);
    rank_transform_int(left_series);

    CVM_BinSeg(left_series,
               changepoints,
               changepoint_offset,
               fdr,
               knots,
               beta,
               intercepts,
               q);
    }

    if(time_series.size() - (cp_vec.second - changepoint_offset) > 10)
    {
      std::vector<int> right_series = std::vector<int>(time_series.begin() + cp_vec.second - changepoint_offset, time_series.end());
      rank_transform_int(right_series);

      CVM_BinSeg(right_series,
                 changepoints,
                 cp_vec.second,
                 fdr,
                 knots,
                 beta,
                 intercepts,
                 q);
    }
  }
}


int main()
{
    std::vector<double> x;
    double temp;
    std::ifstream fin;
    fin.open("data.txt");

    for(int i = 0; i < 400; i++)
    {
        fin>>temp;
        x.push_back(temp);
    }

    fin.close();

    std::vector<int> time_series = rank_transform(x);

    std::map<double, int> changepoints;
    double fdr = 0.01;

    arma::mat knots;
    knots.load("knots.txt", arma::arma_ascii);
    arma::mat beta;
    beta.load("beta.txt", arma::arma_ascii);

    fin.open("qs.txt");
    std::vector<double> qs;

    for(int i = 0; i < 99; i++)
    {
        fin>>temp;
        qs.push_back(temp);
    }

    fin.close();

    fin.open("intercepts.txt");
    std::vector<double> intercepts;

    for(int i = 0; i < 99; i++)
    {
        fin>>temp;
        intercepts.push_back(temp);
    }

    fin.close();

    CVM_BinSeg(time_series, changepoints, 0, fdr, knots, beta, intercepts, qs);

    for(std::map<double, int>::iterator it = changepoints.begin(); it != changepoints.end(); it++)
    {
        std::cout<<"\n Changepoint Significance: " << (*it).first << " Changepoint Location " << (*it).second;
    }

    return 0;
}
