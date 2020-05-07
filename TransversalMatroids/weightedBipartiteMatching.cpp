#include<iostream>
#include<vector>
#include<random>
#include<chrono>
#include<queue>
#include <ext/pb_ds/assoc_container.hpp> // Common file
#include <ext/pb_ds/tree_policy.hpp> // Including tree_order_statistics_node_update
using namespace std;
using namespace __gnu_pbds;

// mt19937_64 rng;
// // initialize the random number generator with time-dependent seed
// uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
// seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};

mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

struct Hungarian{
	vector<vector<double>> costmatrix;
	vector<double> u, v;
	vector<int> p, way;
	int l,r;

	Hungarian(int l, int r): l(l), r(r), u(l+1,0), v(r+1,0), p(r+1,0), way(r+1,0), costmatrix(l+1,vector<double>(r+1,0)){}
	void addEdge(int l, int r, double val){
		costmatrix[l][r] = val;
	}
	double solveAssignmentProblem(){

		for(int i = 0; i <= l; i++) u[i] = 0;
		for(int i = 0; i <= r; i++){
			v[i] = 0;
			p[i] = 0;
			way[i] = 0;
		}

		for(int i = 1; i <= l; i++){
			p[0] = i;
			int j0 = 0;
			vector<double> minv(r+1,2e17+10);
			vector<bool> used(r+1,0);
			do{
				used[j0] = true;
				int i0 = p[j0];
				double delta = 2e17+10;
				int j1 = 0;
				for(int j = 1; j <= r; j++){
					if(!used[j]){
						double curr = costmatrix[i0][j]-u[i0]-v[j];
						if(curr < minv[j]){
							minv[j] = curr;
							way[j] = j0;
						}
						if(minv[j] < delta){
							delta = minv[j];
							j1 = j;
						}
					}
				}
				for(int j = 0; j <= r; j++){
					if(used[j]){
						u[p[j]] += delta;
						v[j] -= delta;
					}
					else minv[j] -= delta;
				}
				j0 = j1;
			}
			while(p[j0] != 0);
			do{
				int j1 = way[j0];
				p[j0] = p[j1];
				j0=j1;
			}
			while(j0!=0);
		}
		return -v[0];
	}
};

double perform_experiment(int n, double p, double factor){

	// Runs Erdos-Renyi random graph with parameters (2n,p)
	// Only considers bipartite graph, left side n nodes and right side n nodes

	int samplesize = factor*n;

	uniform_real_distribution<double> unif(0, 1); // ready to generate random numbers
	normal_distribution<double> gaussian(0,1.0); // ready to generate edge weights

	Hungarian weightedmatch(n,n);

	double edgeweight[n][n];
	vector<int> adj[n]; // adjacency matrix 

	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			edgeweight[i][j] = 0; // initially edge does not exist
			if(unif(rng) <= p){
				edgeweight[i][j] = abs(gaussian(rng));
				adj[i].push_back(j);
			}
		}
	}

	for(int i = 0; i < samplesize; i++){
		for(auto &v: adj[i]){
			weightedmatch.addEdge(v+1,i+1,-edgeweight[i][v]);
		}
	}

	int taken = 0;

	double curtotal = 0, opt = 0;
	bool assignedr[n] = {0};

	for(int i = samplesize; i < n; i++){
		for(auto &v: adj[i]){
			weightedmatch.addEdge(v+1,i+1,-edgeweight[i][v]);
		}
		weightedmatch.solveAssignmentProblem();
		
		int rightmatch = weightedmatch.p[i+1];
		if(rightmatch != 0 && assignedr[rightmatch-1] == 0 && edgeweight[i][rightmatch-1]>1e-9){
			assignedr[rightmatch-1] = 1;
			taken++;
			curtotal += edgeweight[i][rightmatch-1];
		}	
	}


	opt = weightedmatch.solveAssignmentProblem();
	opt = -opt;

	if(abs(opt)<1e-9) return 1;
	return (curtotal/opt);
}


int main()
{

	int iter = 10;
	int n = 50;
	double factor = exp(-1); 
	double lower = 0.1;
	double upper = 0.1;
	double stepsize = 0.1;
	
	for(double p = lower; p <= upper; p += stepsize){
		double countsuccess = 0;
		for(int j = 0; j < iter; j++){
			countsuccess += perform_experiment(n,p,factor);
		}
		cout<<n<<", "<<p<<", "<<(double)countsuccess/iter<<endl;
	}



	return 0;
}