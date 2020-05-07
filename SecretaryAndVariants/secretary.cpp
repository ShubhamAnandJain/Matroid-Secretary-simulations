#include<iostream>
#include<random>
#include<chrono>
using namespace std;

mt19937_64 rng;
// initialize the random number generator with time-dependent seed
uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};

bool perform_experiment(int n){

	// Runs one instance of the secretary problem
	// Returns 1 if best is hired, 0 if not

	uniform_real_distribution<long double> unif(0, 1); // ready to generate random numbers
	// normal_distribution<long double> gaussian(0,100);

	long double e = exp(1);

	int threshold = n/e;
	int maxindex = -1;
	long double maxcur = 0;
	long double maxsample = 0;
	bool choose = 0;

	long double generatednumbers[n];

	for(int i = 0; i < n; i++){
		generatednumbers[i] = unif(rng);
		if(generatednumbers[i]>maxcur){
			maxcur = generatednumbers[i];
			maxindex = i;
		}
	}

	for(int i = 0; i < threshold; i++){
		maxsample = max(maxsample, generatednumbers[i]);
	}

	for(int i = threshold; i < n; i++){
		if(generatednumbers[i]>maxsample){
			if(i == maxindex) return 1;
			return 0;
		}
	}

	return 0;
}


int main()
{

	int iter = 1000; 
	int lower = 10000;
	int upper = 100000;
	int stepsize = 10000;

	for(int n = lower; n <= upper; n += stepsize){
		int countsuccess = 0;
		for(int j = 0; j < iter; j++){
			countsuccess += (int)perform_experiment(n);
		}
		cout<<n<<", "<<(double)countsuccess/iter<<endl;
	}



	return 0;
}