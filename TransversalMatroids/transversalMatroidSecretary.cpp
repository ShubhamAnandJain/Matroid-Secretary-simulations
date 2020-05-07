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

// Max matching 
//1 indexed Hopcroft-Karp Matching in O(E sqrtV)
struct Hopcroft_Karp{
	static const int inf = 1e9;
	int n;
	vector<int> matchL, matchR, dist;
	vector<vector<int> > g;
	Hopcroft_Karp(int n):n(n),matchL(n+1),matchR(n+1),dist(n+1),g(n+1){}
	void addEdge(int u, int v){
		g[u].push_back(v);
	}
	void remove(int u){
		g[u].clear();
	}
	bool bfs(){
		queue<int> q;
		for(int u=1;u<=n;u++){
			if(!matchL[u]){
				dist[u]=0;
				q.push(u);
			}else dist[u]=inf;
		}
		dist[0]=inf;
		while(!q.empty()){
			int u=q.front();
			q.pop();
			for(auto v:g[u]){
				if(dist[matchR[v]] == inf){
					dist[matchR[v]] = dist[u] + 1;
					q.push(matchR[v]);
				}
			}
		}
		return (dist[0]!=inf);
	}
	bool dfs(int u){
		if(!u) return true;
		for(auto v:g[u]){
			if(dist[matchR[v]] == dist[u]+1 &&dfs(matchR[v])){
				matchL[u]=v;
				matchR[v]=u;
				return true;
			}
		}
		dist[u]=inf;
		return false;
	}
	int max_matching(){

		for(int i = 0; i <= n; i++){
			matchL[i] = 0;
			matchR[i] = 0;
			dist[i] = 0;
		}

		int matching=0;
		while(bfs()){
			for(int u=1;u<=n;u++){
				if(!matchL[u])
					if(dfs(u)) matching++;
			}
		}
		return matching;
	}
};

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
	Hopcroft_Karp chkindset(n+n); // n left nodes + n right nodes

	bool edge[n][n];
	double vertexwt[n] = {0};
	vector<int> adj[n]; // adjacency matrix 

	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			edge[i][j] = 0; // initially edge does not exist
			if(unif(rng) <= p){
				edge[i][j] = 1;
				edge[j][i] = 1;
				adj[i].push_back(j);
			}
		}
		vertexwt[i] = abs(gaussian(rng));
	}

	for(int i = 0; i < samplesize; i++){
		for(auto &v: adj[i]){
			weightedmatch.addEdge(v+1,i+1,-vertexwt[i]);
		}
	}

	int taken = 0;

	double curtotal = 0, opt = 0;

	for(int i = samplesize; i < n; i++){

		for(auto &v: adj[i]){
			chkindset.addEdge(i+1,n+v+1);
			weightedmatch.addEdge(v+1,i+1,-vertexwt[i]);
		}

		int chkind = chkindset.max_matching();
		if(chkind != taken+1){
			chkindset.remove(i+1);
			continue;
		}

		weightedmatch.solveAssignmentProblem();
		
		int rightmatch = weightedmatch.p[i+1];
		if(rightmatch != 0 && edge[i][rightmatch-1] == 1){
			taken++;
			curtotal += vertexwt[i];
		}	
		else{
			chkindset.remove(i+1);
		}
	}


	opt = weightedmatch.solveAssignmentProblem();
	opt = -opt;

	if(abs(opt)<1e-9){
		return 1;	
	} 
	assert(curtotal <= opt+1e-9);
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