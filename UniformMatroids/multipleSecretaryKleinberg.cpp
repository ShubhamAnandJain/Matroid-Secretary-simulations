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

void multiple_secretary(int a, int b, int k, double* generatednums, double &curtotal, int &acttake){
	
	if(b<a) return;
	int len = b-a+1;

	if(k == 0) return;
	if(k == 1){
		int samplesize = exp(-1)*(double)len;
		double takemx = 0;
		for(int i = a; i < a+samplesize; i++){
			takemx = max(takemx,generatednums[i]);
		}
		for(int i = a+samplesize; i < b; i++){
			if(generatednums[i] > takemx){
				curtotal += generatednums[i];
				acttake = 1;
				return;
			}
		}
		return;
	}
	
	binomial_distribution<> takek(len, 0.5);
	int numfirst = takek(rng);
	int done = 0;
	int take = min(numfirst,(int)floor(k/2));
	int didtk = 0;
	multiple_secretary(a,a+numfirst-1,take,generatednums,curtotal,didtk);
	take = didtk;
	acttake += didtk;
	vector<int> chkthres;
	for(int i = a; i < a+numfirst; i++){
		chkthres.push_back(generatednums[i]);
	}
	sort(chkthres.begin(),chkthres.end(),greater<int>());
	for(int i = a+numfirst; i <= b; i++){
		if(take == 0){
			if(acttake < k){
				acttake++;
				curtotal += generatednums[i];
			}
		}
		else if(generatednums[i] > chkthres[take-1] && acttake < k){
			acttake++;
			curtotal += generatednums[i];
		}
	}

}

double perform_experiment(int n, int k){

	double generatednumbers[n];

	// uniform_real_distribution<long double> unif(0, 1); // ready to generate random numbers
	normal_distribution<double> gaussian(0,1.0);

	vector<double> gennumbers;

	for(int i = 0; i < n; i++){
		generatednumbers[i] = abs(gaussian(rng));
		gennumbers.push_back(generatednumbers[i]);
	}

	double opt = 0, curtotal = 0;
	int acttake = 0;

	sort(gennumbers.begin(),gennumbers.end(),greater<double>());
	for(int i = 0; i < k; i++) opt += gennumbers[i];
	
	multiple_secretary(0,n-1,k,generatednumbers,curtotal,acttake);

	assert(opt >= 1e-9);

	return (curtotal/opt);
}


int main()
{

	int iter = 5000;
	int n = 1000; 
	
	cout<<"n k comp ratio"<<endl;

	for(int k = 1; k <= n; k += 1){
		double countsuccess = 0;
		for(int j = 0; j < iter; j++){
			countsuccess += perform_experiment(n,k);
		}
		cout<<n<<" "<<k<<" "<<(double)countsuccess/iter<<endl;
	}

	return 0;
}