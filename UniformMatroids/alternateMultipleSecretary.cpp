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

	double curtotal = 0;

	ordered_set curnumbers;

	for(int i = 0; i < samplesize; i++){
		curnumbers.insert(generatednumbers[i]);
	}

	bool taken[k] = {0};

	for(int i = samplesize; i < n; i++){

		curnumbers.insert(generatednumbers[i]);
		int pos = curnumbers.order_of_key(generatednumbers[i]);

		if(pos < k && taken[pos] == 0){
			curtotal += generatednumbers[i];
			taken[pos] = 1;
		}
		
	}

	return (curtotal/opt);
}


int main()
{

	int iter = 100;
	int n = 100; 
	
	cout<<"n k comp ratio samplefactor"<<endl;

	for(int k = 1; k <= n; k += 1){
		double mx = 0, samplesize = 0;
		for(int g = 1; g <= 100; g++){
			double countsuccess = 0;
			for(int j = 0; j < iter; j++){
				countsuccess += perform_experiment(n,k,0.01*g);
			}
			countsuccess /= iter;
			if(mx<countsuccess){
				mx = countsuccess;
				samplesize = 0.01*g;
			}
		}
		cout<<n<<" "<<k<<" "<<mx<<" "<<samplesize<<endl;
	}

	return 0;
}