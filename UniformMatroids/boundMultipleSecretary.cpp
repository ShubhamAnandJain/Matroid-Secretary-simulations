#include<iostream>
#include<cmath>
using namespace std;

long double fastpow(long double a, int p){
	if(p == 0) return 1;
	long double z = fastpow(a,p/2);
	z = z*z;
	if(p%2) z = a*z;
	return z;
}

int main(){

	cout<<"n k comp ratio sample_size"<<endl;

	for(int n = 500; n <= 2000; n++){
		for(int k = 500; k <= 500; k++){
		long double best = 0;
		long double mul = fastpow(exp(1),k);
		long double xval = -1;

		for(long double x = 0.0002; x < (double)k/n; x += 0.0002){
			long double curr = 1.0/n*(min(n,(int)(k*fastpow(exp(1),(int)((n*x)/k))))-(double)n*x);
			long double mul2 = fastpow(exp(1),n*x);
			for(int l = k; l < min(n,(int)(k*fastpow(exp(1),(n*x)/k))); l++){
				curr -= 1.0/n*mul2*fastpow((double)k/l,k)*fastpow(k/(n*x)*log((l-1)/(k-1)),n*x);
			}
			if(curr>best){
				best = curr;
				xval = x*n;
			}
		}

			for(long double x = (double)k/n; x <= 1; x += 0.002){
				long double curr = 0;
				for(int l = n*x; l < min(n,(int)((n*x-1)*exp(1)+1)); l++){
					curr += 1.0/n*max((long double)0, 1-fastpow(n*x/l,k)*mul*fastpow(log((l-1)/(n*x-1)),k));
				}
				if(curr>best){
					best = curr;
					xval = x*n;
				}
			}
			cout<<n<<" "<<k<<" "<<best<<" "<<xval<<endl;
		}
	}

}