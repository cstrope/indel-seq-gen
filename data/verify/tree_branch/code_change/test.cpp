#include<cmath>
#include<stdio.h>

using namespace std;

int main() {
	for (int i = 0; i < 1000; i++) {
		printf("%d: %0.50f\n", i, exp(-i));
	}	
	return 0;
}
