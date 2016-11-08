#include <stdio.h>
#include <stdlib.h>

double gr(unsigned long long n)
{
	double cf = 0.0;

	while (n--) {
		cf = 1.0 / (1.0 + cf);
	}

	return 1.0 + cf;
}


double pi(unsigned long long n)
{
	double cf = 0.0, a, b;

	for (; n; n--) {
		a = (2 * n - 1) * (2 * n - 1);
		b = 2;
		cf = a / (b + cf);
	}

	return 4 / (1.0);
}


int main(int argc, char *argv[])
{
	if (argc <= 1) {
		return 1;
	}

	unsigned long long n = strtoull(argv[1], NULL, 10);

	printf("R = %.12g\n", gr(n));

	return 0;
}
