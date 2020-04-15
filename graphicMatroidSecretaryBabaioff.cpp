#include<iostream>
#include<random>
#include<chrono>
#include<set>
#include<vector>
#include<algorithm>
#include<utility>
using namespace std;

mt19937_64 rng;
// initialize the random number generator with time-dependent seed
// uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
seed_seq ss{uint32_t(0xffffffff), uint32_t(1)}; //For fixed seed
// mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

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


pair<double,pair<int,int>> kruskal(int n, set<pair<double,pair<int,int>>> indsetsample, 
	set<pair<double,pair<int,int>>> indsetchose, pair<double,pair<int,int>> current){

	// Input edges should be sorted in decreasing order
	// For this, we put the edge weights as negative of the true weight

	DisjointSetUnion dsu(n);

	for(auto &v : indsetchose){
		dsu.unite(v.second.first, v.second.second);
	}
	if(dsu.unite(current.second.first,current.second.second) == 0){
		return current;
	}
	for(auto &v : indsetsample){
		if(dsu.unite(v.second.first, v.second.second) == 0){
			if(-v.first < -current.first) return v;
			return current;
		}
	}

	return make_pair(0,make_pair(n,n));
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

set<pair<double,pair<int,int>>> kruskal_get_edges(int n, set<pair<double,pair<int,int>>> edges){

	DisjointSetUnion dsu(n);
	set<pair<double,pair<int,int>>> taken;

	for(auto &v : edges){
		if(dsu.unite(v.second.first,v.second.second)) taken.insert(v);
	}
	return taken;
}

double perform_experiment(int n, double factor){

	uniform_real_distribution<double> unif(0, 1); // ready to generate random numbers
	// normal_distribution<double> gaussian(0,1.0);

	uniform_real_distribution<double> unif1(0.001, 0.002);
	uniform_real_distribution<double> unif2(0.002, 0.003);

	int m = n*(n-1)/2; // Number of edges
	m = 2*(n-2)+1;

	int samplesize = m*factor;
	set<pair<double,pair<int,int>>> edges; //Maintains a running total of edge weights

	double edgeweights[n][n];
	bool insample[n+1][n+1];
	vector<pair<int,int>> edgeorder;
	insample[n][n] = 1;

	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			insample[i][j] = 0;
			edgeweights[i][j] = 0;
		}
	}

	// for(int i = 0; i < n; i++){
	// 	for(int j = 0; j < i; j++){
	// 		insample[i][j] = 0;
	// 		edgeweights[i][j] = abs(unif(rng));
	// 		edgeweights[i][j] = 1e5;
	// 		edgeorder.push_back(make_pair(i,j));
	// 	}
	// }

	for(int i = 0; i < 2; i++){
		for(int j = i+1; j < n; j++){
			if(i == 0 && j == 1) edgeweights[i][j] = n-1;
			else if(i == 0) edgeweights[i][j] = unif1(rng);
			else if(i == 1) edgeweights[i][j] = unif2(rng);
			edgeorder.push_back(make_pair(i,j));
		}
	}

	shuffle(edgeorder.begin(),edgeorder.end(),rng);

	for(int i = 0; i < samplesize; i++){
		int node1 = edgeorder[i].first;
		int node2 = edgeorder[i].second;
		insample[node1][node2] = 1;
		edges.insert(make_pair(-edgeweights[node1][node2],edgeorder[i]));
	}

	set<pair<double,pair<int,int>>> indsetsample = kruskal_get_edges(n,edges);
	set<pair<double,pair<int,int>>> indsetchose;
	double curtotal = 0;
	for(int i = samplesize; i < m; i++){
		int node1 = edgeorder[i].first;
		int node2 = edgeorder[i].second;
		pair<double,pair<int,int>> curr = make_pair(-edgeweights[node1][node2],make_pair(node1,node2));

		edges.insert(curr);
		pair<double,pair<int,int>> checkedge = kruskal(n,indsetsample, indsetchose, curr);
		
		if(insample[checkedge.second.first][checkedge.second.second]){
			curtotal += edgeweights[node1][node2];
			indsetchose.insert(curr);
			if(checkedge.second.first == n) continue;
			indsetsample.erase(checkedge);
		}
	}

	double opt = max_weight_kruskal(n,edges);

	// cout<<curtotal<<" "<<opt<<endl;

	if(abs(opt)<1e-9) return 1;
	return (curtotal/opt);
}


int main()
{

	int iter = 1000; 
	int lower = 50;
	int upper = 50;
	int stepsize = 10;
	double factor = exp(-1);

	for(int n = lower; n <= upper; n += stepsize){
		double countsuccess = 0;
		for(int j = 0; j < iter; j++){
			countsuccess += perform_experiment(n, factor);
			// cout<<"iteration number "<<j<<" current success "<<(double)countsuccess/(j+1)<<endl;
		}
		cout<<n<<", "<<(double)countsuccess/iter<<endl;
	}



	return 0;
}