#include<iostream>
#include<random>
#include<chrono>
#include <ext/pb_ds/assoc_container.hpp> // Common file
#include <ext/pb_ds/tree_policy.hpp> // Including tree_order_statistics_node_update
using namespace std;
using namespace __gnu_pbds;

mt19937_64 rng;
// // initialize the random number generator with time-dependent seed
// uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
seed_seq ss{uint32_t(0xffffffff), uint32_t(1)};

// mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

typedef tree<double, null_type, greater<double>, rb_tree_tag, 
tree_order_statistics_node_update> ordered_set;

vector<double> prob;

double getprob(int target){
	
	int sz = prob.size();
	double dp[sz+1][target+1];
	for(int i = 0; i < sz+1; i++){
		for(int j = 0; j < target+1; j++){
			dp[i][j] = 0;
		}
	}
	dp[0][0] = 1;
	for(int i = 1; i <= sz; i++){
		for(int j = 0; j < target+1; j++){
			if(j != 0) dp[i][j] += prob[i-1]*dp[i-1][j-1];
			dp[i][j] += (1-prob[i-1])*dp[i-1][j];
		}
	}

	double total = 0;
	for(int i = 0; i <= target; i++) total += dp[sz][i];

	return total;
}

double perform_experiment(int n, int k, int sample){

	// cout<<n<<" "<<k<<" "<<sample<<endl;

	double total = 0;	
	prob.clear();

	for(int i = sample+1; i <= n; i++){
		double p = min(1.0,(double)k/i);
		total += 1.0/n*getprob(min((int)prob.size(),k-1));
		prob.push_back(p);
	}

	return total;
}


int main()
{

	int n = 150; 
	
	cout<<"n k comp ratio samplefactor"<<endl;

	for(int k = 1; k <= n; k++){
		double mx = -1, optsamp = -1;
		for(int sample = 0; sample <= n; sample++){
			double ans = perform_experiment(n,k,sample);
			if(ans > mx){
				mx = ans;
				optsamp = (double)sample/n;
			}
		}
		cout<<n<<" "<<k<<" "<<mx<<" "<<optsamp<<endl;
	}

	return 0;
}