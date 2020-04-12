#include<iostream>
#include<random>
#include<chrono>
#include<set>
#include<vector>
#include<algorithm>
#include<utility>
using namespace std;

// mt19937_64 rng;
// initialize the random number generator with time-dependent seed
// uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
// seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};

mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

struct DisjointSetUnion{

	int n;
	int par[10005], size[10005];
	
	DisjointSetUnion(int n): n(n){
		for(int i = 0; i < n; i++){
			par[i] = i;
			size[i] = 1;
		}
	}	

	int find(int node){
		if(par[node] == node) return node;
		return par[node] = find(par[node]);
	}
	
	bool unite(int node1, int node2){
		node1 = find(node1);
		node2 = find(node2);
		if(node1 == node2) return 0;
		if(size[node1]<size[node2]) swap(node1,node2);
		size[node1] += size[node2];
		par[node2] = node1;
		return 1;
	}
};


bool kruskal(int n, set<pair<double,pair<int,int>>> edges, pair<int,int> target){

	// Input edges should be sorted in decreasing order
	// For this, we put the edge weights as negative of the true weight

	DisjointSetUnion dsu(n);

	for(auto &v : edges){
		if(v.second == target){
			return dsu.unite(v.second.first, v.second.second);
		}
		dsu.unite(v.second.first, v.second.second);
	}
}

double max_weight_kruskal(int n, set<pair<double,pair<int,int>>> edges){

	DisjointSetUnion dsu(n);
	double currtot = 0;

	for(auto &v : edges){
		if(dsu.unite(v.second.first,v.second.second)) currtot += v.first;
	}
	currtot = -currtot;

	return currtot;
}

double perform_experiment(int n, double factor){

	// Runs one instance of the secretary problem
	// Returns 1 if best is hired, 0 if not

	// uniform_real_distribution<double> unif(0, 1); // ready to generate random numbers
	normal_distribution<double> gaussian(0,1.0);

	int m = n*(n-1)/2; // Number of edges

	int samplesize = m*factor;
	set<pair<double,pair<int,int>>> edges; //Maintains a running total of edge weights

	double edgeweights[n][n];
	vector<pair<int,int>> edgeorder;

	for(int i = 0; i < n; i++){
		for(int j = 0; j < i; j++){
			edgeweights[i][j] = abs(gaussian(rng));
			edgeorder.push_back(make_pair(i,j));
			// if(i == 1 && j == 0) edgeweights[i][j] = 1e9; 
		}
	}

	shuffle(edgeorder.begin(),edgeorder.end(),rng);

	for(int i = 0; i < samplesize; i++){
		int node1 = edgeorder[i].first;
		int node2 = edgeorder[i].second;
		edges.insert(make_pair(-edgeweights[node1][node2],edgeorder[i]));
	}

	DisjointSetUnion indset(n);
	double curtotal = 0;

	for(int i = samplesize; i < m; i++){
		// if(i%100 == 0) cout<<i<<endl;
		int node1 = edgeorder[i].first;
		int node2 = edgeorder[i].second;
		edges.insert(make_pair(-edgeweights[node1][node2],edgeorder[i]));
		if(indset.find(node1) == indset.find(node2)) continue;
		bool toadd = kruskal(n, edges, edgeorder[i]);
		if(toadd){
			indset.unite(node1,node2);
			curtotal += edgeweights[node1][node2];
		} 
	}

	double opt = max_weight_kruskal(n,edges);

	if(abs(opt)<1e-9) return 1;
	return (curtotal/opt);
}


int main()
{

	int iter = 100; 
	int lower = 100;
	int upper = 100;
	int stepsize = 10;
	double factor = exp(-1);

	for(int n = lower; n <= upper; n += stepsize){
		double countsuccess = 0;
		for(int j = 0; j < iter; j++){
			countsuccess += perform_experiment(n, factor);
			cout<<"iteration number "<<j<<" current success "<<(double)countsuccess/(j+1)<<endl;
		}
		cout<<n<<", "<<(double)countsuccess/iter<<endl;
	}



	return 0;
}