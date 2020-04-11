#include<iostream>
#include<random>
#include<chrono>
#include <ext/pb_ds/assoc_container.hpp> // Common file
#include <ext/pb_ds/tree_policy.hpp> // Including tree_order_statistics_node_update
using namespace std;
using namespace __gnu_pbds;

// mt19937_64 rng;
// // initialize the random number generator with time-dependent seed
// uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
// seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};

mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

typedef tree<double, null_type, greater<double>, rb_tree_tag, 
tree_order_statistics_node_update> ordered_set;

double perform_experiment(int n, int k, double factor){

	int samplesize = factor*n;
	double generatednumbers[n];

	// uniform_real_distribution<long double> unif(0, 1); // ready to generate random numbers
	normal_distribution<double> gaussian(0,1.0);

	vector<double> gennumbers;

	for(int i = 0; i < n; i++){
		generatednumbers[i] = abs(gaussian(rng));
		gennumbers.push_back(generatednumbers[i]);
	}

	sort(gennumbers.begin(),gennumbers.end(),greater<double>());

	double opt = 0;
	for(int i = 0; i < k; i++) opt += gennumbers[i];

	int taken = 0;
	double curtotal = 0;

	ordered_set curnumbers;

	for(int i = 0; i < samplesize; i++){
		curnumbers.insert(generatednumbers[i]);
	}

	for(int i = samplesize; i < n; i++){
		if(taken<k){
			if(curnumbers.size()<k || (generatednumbers[i]>*curnumbers.find_by_order(k-1))){
				taken++;
				curtotal += generatednumbers[i];
			}
		}
		curnumbers.insert(generatednumbers[i]);
	}

	return (curtotal/opt);
}


int main()
{

	int iter = 100;
	int n = 1000;
	double factor = 0; 
	double lower = 0.001;
	double upper = 1;
	double stepsize = 0.001;
	
	for(int k = lower*n; k <= upper*n; k += stepsize*n){
		double countsuccess = 0;
		for(int j = 0; j < iter; j++){
			countsuccess += perform_experiment(n,k,factor);
		}
		cout<<n<<", "<<k<<", "<<(double)countsuccess/iter<<endl;
	}



	return 0;
}