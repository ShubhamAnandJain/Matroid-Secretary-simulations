#include<iostream>
#include<random>
#include<chrono>
#include<queue>
using namespace std;

// mt19937_64 rng;
// initialize the random number generator with time-dependent seed
// uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
// seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

bool perform_experiment(int n, double sample, double myu, double lambda){

	// Runs one instance of the stochastic secretary problem
	// Returns 1 if best is hired, 0 if not
	// Takes sample*n as the sample

	// Arrival distribution, taken as poisson with myu as mean

	poisson_distribution<int> distribution1 (myu);

	// Departure distribution, taken as exponential with 1/lambda as mean

	exponential_distribution<double> distribution2 (lambda);

	// Distribution of skills, taken to be uniform

	uniform_real_distribution<double> unif(0, 1); // ready to generate random numbers

	// Actual code follows

	int threshold = n*sample;
	int maxindex = -1;
	int holduntil = n;
	double maxcur = 0;
	double maxsample = 0;
	double maxalgo = 0;
	int indcurr = -1;

	double generatednumbers[n];

	for(int i = 0; i < n; i++){
		generatednumbers[i] = unif(rng);
		if(generatednumbers[i]>maxcur){
			maxcur = generatednumbers[i];
			maxindex = i;
		}
	}

	priority_queue<pair<double,pair<int,int>>> events;

	for(int i = 0; i < n; i++){
		events.push(make_pair(-distribution1(rng),make_pair(1,i)));
	}

	int left = 0;

	while(!events.empty()){
		auto curr = events.top();
		int curindex = curr.second.second;
		events.pop();

		if(curr.second.first == 1){
			if(generatednumbers[curindex]>maxcur){
				maxindex = curindex;
				maxcur = generatednumbers[maxindex];
			}
			double endtime = curr.first - distribution2(rng);
			events.push(make_pair(endtime,make_pair(0,curindex)));
			if(generatednumbers[curindex] > maxalgo){
				maxalgo = generatednumbers[curindex];
				indcurr = curindex;
			}

			continue;
		}

		else if(curr.second.first == 0){
			left++;
			if(left < threshold && curindex == maxindex){
				return 0;
			}
			else if(left < threshold){
				maxsample = max(maxsample, generatednumbers[curindex]);
				continue;
			}
			if(curindex == maxindex){
				return 1;
			}
			if(generatednumbers[curindex] >= maxsample && indcurr == curindex){
				return 0;
			}

		}

		else{
			cout<<"Something went wrong\n";
			return 0;
		}
	}

	return 0;
}


int main()
{

	int iter = 1000; 
	int lower = 1000;
	int upper = 10000;
	int stepsize = 1000;
	double sample = 0.44;

	for(int n = lower; n <= upper; n += stepsize){
		int countsuccess = 0;
		for(int j = 0; j < iter; j++){
			countsuccess += (int)perform_experiment(n,sample,10,1);
		}
		cout<<n<<", "<<(double)countsuccess/iter<<endl;
	}



	return 0;
}