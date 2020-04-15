#include<iostream>
#include<cmath>
using namespace std;

long double fastpow(long double a, int p){
	cout<<a<<" "<<p<<endl;
	if(p == 0) return 1;
	long double z = fastpow(a,p/2);
	z = z*z;
	if(p%2) z = a*z;
	return z;
}

int main(){

	int n = 100;

	for(int k = 1; k <= 100; k++){
		long double best = 0;
		long double mul = fastpow(exp(1),k);
		long double xval = -1;

		for(long double x = 0.02; x <= 1; x += 0.01){
			long double curr = 0;
			for(int l = n*x; l <= n; l++){
				curr += 1.0/n*max((long double)0, 1-mul*fastpow(n*x/l,k)*fastpow(log((l-1)/(n*x-1)),k));
			}
			if(curr>best){
				best = curr;
				xval = x*n;
			}
		}
		cout<<n<<" "<<k<<" "<<best<<" "<<xval<<endl;
	}
	
}