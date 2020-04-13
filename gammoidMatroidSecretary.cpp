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

struct Edge{
    int from, to, capacity;
    double cost;
};


struct MinCostFlow{

	//reference: https://cp-algorithms.com/graph/min_cost_flow.html
	// This is a slightly edited version
	// O()

	vector<vector<int>> adj, capacity;
	vector<vector<double>> cost;

	const double INF = 1e9;

	void shortest_paths(int n, int v0, vector<double>& d, vector<int>& p) {
	    d.assign(n, INF);
	    d[v0] = 0;
	    vector<bool> inq(n, false);
	    queue<int> q;
	    q.push(v0);
	    p.assign(n, -1);

	    while (!q.empty()) {
	        int u = q.front();
	        q.pop();
	        inq[u] = false;
	        for (int v : adj[u]) {
	            if (capacity[u][v] > 0 && d[v] > d[u] + cost[u][v]) {
	                d[v] = d[u] + cost[u][v];
	                p[v] = u;
	                if (!inq[v]) {
	                    inq[v] = true;
	                    q.push(v);
	                }
	            }
	        }
	    }
	}

	pair<double,bool> min_cost_flow(int N, vector<Edge> edges, int K, int s, int t, Edge target) {
	    adj.assign(N, vector<int>());
	    cost.assign(N, vector<double>(N, 0));
	    capacity.assign(N, vector<int>(N, 0));
	    for (Edge e : edges) {
	        adj[e.from].push_back(e.to);
	        adj[e.to].push_back(e.from);
	        cost[e.from][e.to] = e.cost;
	        cost[e.to][e.from] = -e.cost;
	        capacity[e.from][e.to] = e.capacity;
	    }

	    int flow = 0;
	    double cost = 0;
	    vector<double> d;
	    vector<int> p;
	    while (flow < K) {
	        shortest_paths(N, s, d, p);
	        if (d[t] == INF)
	            break;

	        // find max flow on that path
	        int f = K - flow;
	        int cur = t;
	        while (cur != s) {
	            f = min(f, capacity[p[cur]][cur]);
	            cur = p[cur];
	        }

	        // apply flow
	        flow += f;
	        cost += f * d[t];
	        cur = t;
	        while (cur != s) {
	            capacity[p[cur]][cur] -= f;
	            capacity[cur][p[cur]] += f;
	            cur = p[cur];
	        }
	    }

	    if (flow < K)
	        return make_pair(1,0);
	    else
	    	if(capacity[target.from][target.to] < 1e-9){
	    		return make_pair(cost,1);
	    	}
	    	else{
	    		return make_pair(cost,0);
	    	}
	}
};	

struct Max_Flow{

	// reference: https://cp-algorithms.com/graph/push-relabel-faster.html
	// O(VE + V^2sqrt(E))

	const int inf = 1000000000;

	int n;
	vector<vector<int>> capacity, flow;
	vector<int> height, excess;

	Max_Flow(int n){
		this->n = n;
		capacity.assign(n,vector<int>(n,0));
	}

	void push(int u, int v)
	{
	    int d = min(excess[u], capacity[u][v] - flow[u][v]);
	    flow[u][v] += d;
	    flow[v][u] -= d;
	    excess[u] -= d;
	    excess[v] += d;
	}

	void relabel(int u)
	{
	    int d = inf;
	    for (int i = 0; i < n; i++) {
	        if (capacity[u][i] - flow[u][i] > 0)
	            d = min(d, height[i]);
	    }
	    if (d < inf)
	        height[u] = d + 1;
	}

	vector<int> find_max_height_vertices(int s, int t) {
	    vector<int> max_height;
	    for (int i = 0; i < n; i++) {
	        if (i != s && i != t && excess[i] > 0) {
	            if (!max_height.empty() && height[i] > height[max_height[0]])
	                max_height.clear();
	            if (max_height.empty() || height[i] == height[max_height[0]])
	                max_height.push_back(i);
	        }
	    }
	    return max_height;
	}

	int max_flow(int s, int t)
	{
	    height.assign(n, 0);
	    height[s] = n;
	    flow.assign(n, vector<int>(n, 0));
	    excess.assign(n, 0);
	    excess[s] = inf;
	    for (int i = 0; i < n; i++) {
	        if (i != s)
	            push(s, i);
	    }

	    vector<int> current;
	    while (!(current = find_max_height_vertices(s, t)).empty()) {
	        for (int i : current) {
	            bool pushed = false;
	            for (int j = 0; j < n && excess[i]; j++) {
	                if (capacity[i][j] - flow[i][j] > 0 && height[i] == height[j] + 1) {
	                    push(i, j);
	                    pushed = true;
	                }
	            }
	            if (!pushed) {
	                relabel(i);
	                break;
	            }
	        }
	    }

	    int max_flow = 0;
	    for (int i = 0; i < n; i++)
	        max_flow += flow[0][i];

	    return max_flow;
	}

};

double perform_experiment(int n, int s, int t, double p, double factor){

	// s is the number of source vertices
	// t is the number of destination vertices
	// p is the edge existence probability
	// n is the number of the vertices in the graph excluding s and t
	// factor is the sampling fraction

	int samplesize = factor*s;

	uniform_real_distribution<double> unif(0, 1); // ready to generate random numbers
	normal_distribution<double> gaussian(0,1.0); // ready to generate edge weights

	double curtotal = 0, opt = 0;

	double vertexwt[s] = {0};

	for(int i = 0; i < s; i++){
		vertexwt[i] = abs(gaussian(rng));
	}

	int edgeweight[s+n+t][s+n+t];

	for(int i = 0; i < s+n+t; i++){
		for(int j = 0; j < s+n+t; j++){
			edgeweight[i][j] = 0;
		}
	}

	int numedges = 0;

	for(int i = 0; i < s; i++){
		for(int j = s; j < s+n+t; j++){
			if(unif(rng) <= p){
				edgeweight[i][j] = 1;
				numedges++;
			}
		}
	}

	// for(int i = 0; i < s+n+t; i++){
	// 	for(int j = 0; j < s+n+t; j++){
	// 		if(edgeweight[i][j]) cout<<"edge "<<i<<" "<<j<<endl;
	// 	}
	// }

	vector<Edge> edges;

	for(int i = s; i < s+n; i++){
		for(int j = i+1; j < s+n+t; j++){
			if(unif(rng) <= p){
				edgeweight[i][j] = 1;
				numedges++;
			}
		}
	}

	Max_Flow chkindset(1+s+n+t+1+numedges);
	MinCostFlow maxweightindset;

	int cnt = 0;

	for(int i = s; i < s+n; i++){
		for(int j = i+1; j < s+n+t; j++){
			if(edgeweight[i][j] == 1){
				cnt++;
				Edge x;
				x.from = i+1;
				x.to = 1+n+s+t+cnt;
				x.capacity = 1;
				x.cost = 0;
				edges.push_back(x);
				x.from = 1+n+s+t+cnt;
				x.to = j+1;
				edges.push_back(x);
				x.from = j+1;
				x.to = i+1;
				edges.push_back(x);
				chkindset.capacity[i+1][1+n+s+t+cnt] = 1;
				chkindset.capacity[1+n+s+t+cnt][j+1] = 1;
				chkindset.capacity[j+1][i+1] = 1;
			}
		}
	}

	for(int i = n+s; i < n+s+t; i++){
		Edge x;
		x.from = i+1;
		x.to = s+n+t+1;
		x.capacity = 1;
		x.cost = 0;
		edges.push_back(x);
		chkindset.capacity[i+1][s+n+t+1] = 1;
	}

	int taken = 0;

	for(int i = 0; i < samplesize; i++){
		Edge x;
		x.from = 0;
		x.to = i+1;
		x.capacity = 1;
		x.cost = -vertexwt[i];
		edges.push_back(x);

		for(int j = s; j < n+s+t; j++){
			if(edgeweight[i][j] == 1){
				cnt++;
				x.from = i+1;
				x.to = 1+n+s+t+cnt;
				x.capacity = 1;
				x.cost = 0;
				edges.push_back(x);
				x.from = 1+n+s+t+cnt;
				x.to = j+1;
				edges.push_back(x);
				x.from = j+1;
				x.to = i+1;
				edges.push_back(x);
			}
		}
	}

	for(int i = samplesize; i < s; i++){
		Edge x, target;
		x.from = 0;
		x.to = i+1;
		target.from = 0;
		target.to = i+1;
		target.capacity = 1;
		target.cost = 0;
		x.capacity = 1;
		x.cost = -vertexwt[i];
		edges.push_back(x);
		chkindset.capacity[0][i+1] = 1;
		int last = cnt;

		for(int j = s; j < n+s+t; j++){
			if(edgeweight[i][j] == 1){
				cnt++;
				x.from = i+1;
				x.to = 1+n+s+t+cnt;
				x.capacity = 1;
				x.cost = 0;
				edges.push_back(x);
				x.from = 1+n+s+t+cnt;
				x.to = j+1;
				edges.push_back(x);
				x.from = j+1;
				x.to = i+1;
				edges.push_back(x);
				chkindset.capacity[i+1][n+s+t+1+cnt] = 1;
				chkindset.capacity[n+s+t+1+cnt][j+1] = 1;
				chkindset.capacity[j+1][i+1] = 1;
			}
		}

		if(chkindset.max_flow(0,s+n+t+1) != taken+1){
			// cout<<"wut lol "<<i<<endl;
			chkindset.capacity[0][i+1] = 0;
			for(int j = s; j < n+s+t; j++){
				if(edgeweight[i][j] == 1){
					last++;
					chkindset.capacity[i+1][n+s+t+1+last] = 0;
					chkindset.capacity[1+n+s+t+last][j+1] = 0;
					chkindset.capacity[j+1][i+1] = 0;
				}
			}
			continue;
		}

		double mx = -1e9;
		bool isinc = 0;

		for(int j = 1; j <= s; j++){
			pair<double,bool> costchk = maxweightindset.min_cost_flow(1+n+s+t+1+numedges,
				edges,j,0,n+s+t+1,target);
			if((mx < -costchk.first) || (mx == -costchk.first && costchk.second == 1)){
				mx = -costchk.first;
				isinc = costchk.second;
			}
		}
		if(isinc){
			taken++;
			curtotal += vertexwt[i];
		}
		else{
			chkindset.capacity[0][i+1] = 0;
			for(int j = s; j < n+s+t; j++){
				if(edgeweight[i][j] == 1){
					last++;
					chkindset.capacity[i+1][1+n+s+t+last] = 0;
					chkindset.capacity[1+n+s+t+last][j+1] = 0;
					chkindset.capacity[j+1][i+1] = 0;
				}
			}
			continue;
		}
	}

	Edge target;
	target.from = 0;
	target.to = 0;

	for(int j = 1; j <= s; j++){
		opt = max(opt,-maxweightindset.min_cost_flow(1+n+s+t+1+numedges,edges,
		j,0,n+s+t+1,target).first);
	}

	if(opt<1e-9){
		return 1;
	}

	// cout<<samplesize<<endl;
	// for(int i = 0; i < s; i++) cout<<i<<" "<<vertexwt[i]<<endl;

	// cout<<"try "<<curtotal<<" "<<opt<<endl;

	return (curtotal/opt);
}


int main()
{

	int iter = 100;
	int n = 10;
	int s = 50;
	int t = 5;
	double factor = exp(-1); 
	double lower = 0.3;
	double upper = 0.3;
	double stepsize = 0.3;
	
	for(double p = lower; p <= upper; p += stepsize){
		double countsuccess = 0;
		for(int j = 0; j < iter; j++){
			countsuccess += perform_experiment(n,s,t,p,factor);
			cout<<"running success iteration "<<j+1<<" "<<countsuccess/(j+1)<<endl;
		}
		cout<<n<<", "<<p<<", "<<(double)countsuccess/iter<<endl;
	}



	return 0;
}