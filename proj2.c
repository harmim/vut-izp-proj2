/**
 * @name Projekt 2 - Iterační výpočty
 * @author Dominik Harmim <xharmi00@stud.fit.vutbr.cz>
 */

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define PRINT_ERR(s, ...) fprintf(stderr, s "\n", __VA_ARGS__)
#define RESULT_PRECISION_FOR_OUTPUT 12
#define EPS 1e-8


/**
 * constant string - help
 */
const char usage_string[] =
	"	proj2 --log X N (calculating the natural logarithm of the number X in\n"
	"		N iterations)\n"
	"	proj2 --pow X Y N (calculate exponential function of the number Y with\n"
	"		the basis X in N iterations)";


/**
 * function prototypes
 */
double taylor_log(double x, unsigned int n);
double cfrac_log(double x, unsigned int n);
double taylor_pow(double x, double y, unsigned int n);
double taylorcf_pow(double x, double y, unsigned int n);
double mylog(double x);
double mypow(double x, double y);
float check_spc_argv_of_logarithm(double x);


/**
 * check special value of logarithm argument
 *
 * @param x number from which we calculate the logarithm
 * @return if argument of logaritm is special value, return special result, otherwise -1
 */
float check_spc_argv_of_logarithm(double x)
{
	if (fabs(x) == 0.0) {
		return -INFINITY;
	} else if (x < 0) {
		return NAN;
	} else if (x == 1.0) {
		return 0.0;
	} else if (isinf(x)) {
		return INFINITY;
	} else if (isnan(x)) {
		return NAN;
	}

	return -1;
}


/**
 * calculating the natural logarithm of the number `x` in `n` iterations using Taylor polynomial
 *
 * @param x number from which we calculate the logarithm
 * @param n number of members in polynomial > 0
 * @return logarithm of the number x
 */
double taylor_log(double x, unsigned int n)
{
	double spc_result;
	if ((spc_result = check_spc_argv_of_logarithm(x)) != -1) {
		return spc_result;
	}

	double sum = 0.0, numerator = 1.0, frac;
	unsigned int i;

	if (x < 1.0) {
		// algorithm for (0,1)
		x = 1.0 - x;
		for (i = 1; i <= n; i++) {
			numerator *= x;
			frac = numerator / i;
			if (isinf(frac)) {
				break;
			}
			sum -= frac;
		}
	} else {
		// algorithm for <1,INF)
		for (i = 1; i <= n; i++) {
			numerator *= (x - 1.0) / x;
			frac = numerator / i;
			if (isinf(frac)) {
				break;
			}
			sum += frac;
		}
	}

	return sum;
}


/**
 * calculating the natural logarithm of the number `x` in `n` iterations using continued fractions
 *
 * @param x number from which we calculate the logarithm
 * @param n number of steps in continued fraction > 0
 * @return logarithm of the number x
 */
double cfrac_log(double x, unsigned int n)
{
	double spc_result;
	if ((spc_result = check_spc_argv_of_logarithm(x)) != -1) {
		return spc_result;
	}

	x = (x - 1.0) / (x + 1.0);
	unsigned int coef = 2 * n - 1, i = n;
	double sum = coef, pow2x = x * x, frac;
	while (i > 1) {
		i--;
		coef -= 2;
		frac = i * i * pow2x / sum;
		if (isinf(frac)) {
			break;
		}
		sum = coef - frac;
	}

	return 2.0 * x / sum;
}


/**
 * calculate exponential function of the number `y` with the basis `x` in `n` iterations
 * for the calculation of the natural logarithm uses function taylor_log
 *
 * @see taylor_log
 * @param x number from which we calculate the exponencial function
 * @param y exponent
 * @param n number of members in polynomial > 0
 * @return the base to the exponent power
 */
double taylor_pow(double x, double y, unsigned int n)
{
	if (x <= 0.0) {
		return NAN;
	}

	double log_x = taylor_log(x, n);
	if (isinf(log_x)) {
		return log_x;
	}

	double sum = 1.0, pow_y = 1.0, pow_ln_x = 1.0, denominator = 1.0, frac;
	for (unsigned int i = 1; i <= n; i++) {
		pow_y *= y;
		pow_ln_x *= log_x;
		denominator *= i;
		frac = pow_y * pow_ln_x / denominator;
		if (isinf(frac)) {
			break;
		}
		sum += frac;
	}

	return sum;
}


/**
 * calculate exponential function of the number `y` with the basis `x` in `n` iterations
 * for the calculation of the natural logarithm usees function cfrac_log
 *
 * @see cfrac_log
 * @param x number from which we calculate the exponencial function
 * @param y exponent
 * @param n number of members in polynomial > 0
 * @return the base to the exponent power
 */
double taylorcf_pow(double x, double y, unsigned int n)
{
	if (x <= 0.0) {
		return NAN;
	}

	double log_x = cfrac_log(x, n);
	if (isinf(log_x)) {
		return log_x;
	}

	double sum = 1.0, pow_y = 1.0, pow_ln_x = 1.0, denominator = 1.0, frac;
	for (unsigned int i = 1; i <= n; i++) {
		pow_y *= y;
		pow_ln_x *= log_x;
		denominator *= i;
		frac = pow_y * pow_ln_x / denominator;
		if (isinf(frac)) {
			break;
		}
		sum += frac;
	}

	return sum;
}


/**
 * calculating the natural logarithm of the number `x`
 * elects to accurate calculation type and minimum number of iterations for required accuracy EPS
 *
 * @param x number from which we calculate the logarithm
 * @return logarithm of the number x
 */
double mylog(double x)
{
	double spc_result;
	if ((spc_result = check_spc_argv_of_logarithm(x)) != -1) {
		return spc_result;
	}

	double taylor_result = 0.0, taylor_prev_result, taylor_numerator = 1.0, taylor_frac;
	double cfrac_result = 0.0, cfrac_prev_result;
	unsigned int n = 1;

	do {
		// Taylor
		taylor_prev_result = taylor_result;

		if (x < 1.0) {
			// algorithm for (0,1)
			taylor_numerator *= 1.0 - x;
			taylor_frac = taylor_numerator / n;
			if (isinf(taylor_frac)) {
				break;
			}
			taylor_result -= taylor_frac;
		} else {
			// algorithm for <1,INF)
			taylor_numerator *= (x - 1.0) / x;
			taylor_frac = taylor_numerator / n;
			if (isinf(taylor_frac)) {
				break;
			}
			taylor_result += taylor_frac;
		}

		// Cfrac
		cfrac_prev_result = cfrac_result;
		cfrac_result = cfrac_log(x, n);

		n++;
	} while (fabs(taylor_result - taylor_prev_result) > EPS && fabs(cfrac_result - cfrac_prev_result) > EPS);

	return fabs(taylor_result - taylor_prev_result) <= EPS ? taylor_result : cfrac_result;
}


/**
 * calculate exponential function of the number `y` with the basis `x`
 * elects to accurate calculation type and minimum number of iterations for required accuracy EPS
 *
 * @param x number from which we calculate the exponencial function
 * @param y exponent
 * @return the base to the exponent power
 */
double mypow(double x, double y)
{
	if (x <= 0.0) {
		return NAN;
	}

	double log_x = mylog(x);
	if (isinf(log_x)) {
		return log_x;
	}

	double sum = 1.0, pow_y = 1.0, pow_ln_x = 1.0, denominator = 1.0, frac = 0.0, prev_frac;
	unsigned int i = 1;
	do {
		prev_frac = frac;
		pow_y *= y;
		pow_ln_x *= log_x;
		denominator *= i;
		frac = pow_y * pow_ln_x / denominator;
		if (isinf(frac)) {
			break;
		}
		sum += frac;
		i++;
	} while (fabs(frac - prev_frac) > EPS);

	return sum;
}


/**
 * converts value of argument in string format to double
 * if in value there is non-convertible part, print it to stderr
 *
 * @param value value of argument in string format (from argv)
 * @param arg name of argument that value belongs to
 * @param error set to true if in value there is non-convertible part
 * @return converted double value
 */
double value_of_arg_to_double(const char *value, const char *arg, bool *error)
{
	char *endptr = NULL;
	double number = strtod(value, &endptr);
	if (*endptr) {
		PRINT_ERR("Value of argument %s must be real number but there is non-convertible part: %s.", arg, endptr);
		*error = true;
	}

	return number;
}


/**
 * converts value of argument in string format to unsigned int
 * if in value there is non-convertible part, or minus sign, print it to stderr
 *
 * @param value value of argument in string format (from argv)
 * @param arg name of argument that value belongs to
 * @param error set to true if in value there is non-convertible part or minus sign
 * @return converted unsigned int value
 */
unsigned int value_of_arg_to_unsigned_int(const char *value, const char *arg, bool *error)
{
	if (strchr(value, '-') != NULL) {
		PRINT_ERR("Value of argument %s must be a positive number, given %s.", arg, value);
		*error = true;
		return 0;
	}

	char *endptr = NULL;
	unsigned int number = strtoul(value, &endptr, 10);
	if (*endptr) {
		PRINT_ERR("Value of argument %s must be integer but there is non-convertible part: %s.", arg, endptr);
		*error = true;
	}

	return number;
}


/**
 * print results of logarithm functions to stdout
 *
 * @param x number from which we calculate the logarithm
 * @param n number of iterations for calculation
 */
void print_log_results(double x, unsigned int n)
{
	printf("       log(%g) = %.*g\n", x, RESULT_PRECISION_FOR_OUTPUT, log(x));
	printf(" cfrac_log(%g) = %.*g\n", x, RESULT_PRECISION_FOR_OUTPUT, cfrac_log(x, n));
	printf("taylor_log(%g) = %.*g\n", x, RESULT_PRECISION_FOR_OUTPUT, taylor_log(x, n));
#ifdef DEBUG
	printf("       log(%g) = %.7e\n", x, log(x));
	printf("     mylog(%g) = %.7e\n", x, mylog(x));
#endif
}


/**
 * print results of exponencial functions to stdout
 *
 * @param x number from which we calculate expinencial function
 * @param y basis of exponencial function
 * @param n number of iterations for calculation
 */
void print_pow_results(double x, double y, unsigned int n)
{
	printf("         pow(%g,%g) = %.*g\n", x, y, RESULT_PRECISION_FOR_OUTPUT, pow(x, y));
	printf("  taylor_pow(%g,%g) = %.*g\n", x, y, RESULT_PRECISION_FOR_OUTPUT, taylor_pow(x, y, n));
	printf("taylorcf_pow(%g,%g) = %.*g\n", x, y, RESULT_PRECISION_FOR_OUTPUT, taylorcf_pow(x, y, n));
#ifdef DEBUG
	printf("         pow(%g,%g) = %.7e\n", x, y, pow(x, y));
	printf("       mypow(%g,%g) = %.7e\n", x, y, mypow(x, y));
#endif
}


/**
 * process input arguments
 *
 * @param argc count of program arguments
 * @param argv program arguments
 * @param help true if usage string should be listed, false otherwise
 * @return true if arguments was processed successfully, false otherwise
 */
bool process_input_args(const int argc, const char *argv[], bool *help)
{
	if (argc > 1) {
		// option --log
		if (strcmp(argv[1], "--log") == 0) {
			// argument --log must be followed by argument X and N
			if (argc != 4) {
				*help = true;
				return true;
			}

			//convert arguments and do validation
			bool error = false;
			double x = value_of_arg_to_double(argv[2], "--log X", &error);
			if (error) {
				return false;
			}
			unsigned int n = value_of_arg_to_unsigned_int(argv[3], "--log N", &error);
			if (error) {
				return false;
			}
			if (n == 0) {
				PRINT_ERR("Value of argument --log N must be greater than 0, given %i.", n);
				return false;
			}

			// print results to stdout
			print_log_results(x, n);

			return true;

		// option --pow
		} else if (strcmp(argv[1], "--pow") == 0) {
			// argument --pow must be followed by argument X, Y and N
			if (argc != 5) {
				*help = true;
				return true;
			}

			//convert arguments and do validation
			bool error = false;
			double x = value_of_arg_to_double(argv[2], "--pow X", &error);
			if (error) {
				return false;
			}
			double y = value_of_arg_to_double(argv[3], "--pow Y", &error);
			if (error) {
				return false;
			}
			unsigned int n = value_of_arg_to_unsigned_int(argv[4], "--pow N", &error);
			if (error) {
				return false;
			}
			if (n == 0) {
				PRINT_ERR("Value of argument --pow N must be greater than 0, given %i.", n);
				return false;
			}

			// print results to stdout
			print_pow_results(x, y, n);

			return true;
		}
	}

	// unknown option
	*help = true;
	return true;
}


int main(const int argc, const char *argv[])
{
	// start processing arguments
	bool help = false, result = process_input_args(argc, argv, &help);
	if (help) {
		// to be listed help/usage string
		printf("usage:\n%s\n", usage_string);
	}
	// return value of this program (of this main function)
	// depends on return value of function process_input_args
	return result ? EXIT_SUCCESS : EXIT_FAILURE;
}
