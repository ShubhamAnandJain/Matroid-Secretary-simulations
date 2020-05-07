#include<iostream>
#include<random>
#include<algorithm>
#include<chrono>
#include <ctime>    // For time()
#include <cstdlib>  // For srand() and rand()


using namespace std;

void random_perm(int ranks[], int n){
	// random permutation
	//int ranks[n];
	for (int i=0; i<n; i++) ranks[i]=i+1;

	//for (int i=0; i<n; i++) cout << ranks[i] << ' ';

	for (int j=n; j>1; j--){
		int r = rand()%(j);
		swap(ranks[n-j],ranks[n-j+r]);
	}
}

int comp_ratio(int ranks[], int n, int k, int sampleSize){
	int topk[k]; // stores top k largest elements till now.
	//int isTopk[n]; //denote if an element is in top k till now.
	
	for (int i=0; i<k; i++) topk[i]=0;
	//for (int i=0; i<n; i++) isTopk[i]=0;
	sort(ranks,ranks+sampleSize);	
	int minimum = min(k, sampleSize);
	for (int i=0; i<minimum; i++) topk[i]=ranks[sampleSize-i-1]; 
	
	int comp_ratio=0;
	int curr_size=0;


	for (int i =sampleSize; (i<n && curr_size<k); i++){
		int curr=ranks[i];
		int pos;
		for (pos=0; pos<k; pos++){
			if (topk[pos]<curr) break; 
		}
		for (int j= k-1; j>pos; j--){
			topk[j]=topk[j-1];
		}
		if(pos<k){
			topk[pos]=curr;
			curr_size++;
			if(ranks[i]>n-k) comp_ratio++;
		}
		
	}



	return comp_ratio;
	
}

int main()
{

	int n = 100;
	int iter = n*n;
	int n_by_e = n*exp(-1);
	int n_by_k = 3;
	//int k = 1;
	//int sampleSize=7;

	srand(time(0));  // Initialize random number generator.
	
	//random_perm(ranks, n);
	
	//for (int i=0; i<n; i++) cout << ranks[i] << ' ';
	//cout << '\n';
	
	cout << "n k best_comp_ratio best_sample_factor\n";
	
	for(int k=1; k<n; k=k+1){
		iter = n*n;
		int ranks[n];

		double best_comp_ratio=0;
		double bestSampleSize=0;

		for(int sampleSize=1; sampleSize<=n-1; sampleSize++){
			double avg_comp_ratio=0;
			for(int i=0; i<iter; i++){
				random_perm(ranks, n);
				avg_comp_ratio+=comp_ratio(ranks,n,k,sampleSize);
			}
			avg_comp_ratio = avg_comp_ratio/iter;
			if (avg_comp_ratio >= best_comp_ratio){
				best_comp_ratio = avg_comp_ratio;
				bestSampleSize=sampleSize;
			}
		}
		cout << n << " " << k << " " << best_comp_ratio/k << " " <<  bestSampleSize/n << '\n';
	}
	//double avg_comp_ratio=0;
	//
	//int comp= comp_ratio(ranks,n,k,sampleSize);
	//cout << comp << '\n';

	return 0;
}