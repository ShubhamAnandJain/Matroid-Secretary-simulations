#include<iostream>
#include<random>
#include<algorithm>
#include<chrono>
#include <ctime>    // For time()
#include <cstdlib>  // For srand() and rand()


using namespace std;

int bernoulli(double p){
	// random 0-1 variable with pr p.
	double r = ((double) rand() / (RAND_MAX));
	cout << r << " " << RAND_MAX << "\n";
	if(r<p)
		return 1; 
	else 
		return 0;
}

void sum_k_atATime(double a[], double sum[], int n, int k){
	// sum should be of size  k
	
	if (n==1){
		sum[0]=a[0];	
	} 
	else{
		sum_k_atATime(a, sum, n-1, k);
		int m = min(n,k);
		for (int j=m-1; j>0; j--){
			// sum j+1 at a time.
			sum[j]=sum[j]+a[n-1]*sum[j-1];
		}
		sum[0]=sum[0]+a[n-1];
	}

				
}

double probability(int n, int k, int t){
		
	double p[n-t]; // prob of being in top k
	for (int i=0; i<n-t; i++) p[i] =  (double)k/(i+t+1);
	
	double prod = 1;
	for (int i=0; i<n-t; i++) prod = prod*(1-p[i]); 
	
	double pbyoneminusp[n-t];
	for (int i=0; i<n-t; i++) pbyoneminusp[i] = (double)k/(i+t+1-k);
	
	
	double sum[k-1];
	for (int i=0;i<k-1;i++) sum[i]=0.0;

	sum_k_atATime(pbyoneminusp, sum, n-t, k-1);
		
	double prob= 1;
	for (int i=0;i<k-1;i++) prob = prob + sum[i];
	prob = prob*prod;
	return prob;
}

int main()
{

	int n = 100;
	//int k = 2;
	int t = n*exp(-1); //sample size

	// "probability that the there are less than k elements till (n-1)th element that are in top k till now" 
	srand(time(0));  // Initialize random number generator.
	
	cout << "n \t k \t t \t prob \n";
	for (int k=1;k<=t;k++)
		cout << n << "\t " << k << "\t " << t << "\t " <<  probability(n,k,t) << "\n";
	return 0;
}