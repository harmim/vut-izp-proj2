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


/**
 * constant string - help
 */
const char usage_string[] =
	"	proj2 --log X N (calculating the natural logarithm of the number X in\n"
	"		N iterations)\n"
	"	proj2 --pow X Y N (calculate exponential function of the number Y with\n"
	"		the general basis X in N iterations)";


/**
 * calculate the base to the exponent power (only for exponent >= 0)
 * recursive function
 *
 * @param base
 * @param exponent (greater than or equal to 0)
 * @return the base to the exponent power
 */
double pow_positive(double base, double exponent)
{
	if (exponent == 0.0) {
		return 1.0;
	}

	return base * pow_positive(base, exponent - 1.0);
}


/**
 * calculate factorial of number greater than or equal to zero
 * recursive function
 *
 * @param n integer >= 0
 * @return factorial of number `n`
 */
unsigned long long factorial(unsigned int n)
{
	if (n == 1 || n == 0) {
		return 1;
	}

	return n * factorial(n - 1);
}


/**
 * calculating the natural logarithm of the number `x` in `n` iterations using Taylor polynomial
 *
 * @param x real number // TODO: <= 0
 * @param n number of members in polynomial // TODO  = 0
 * @return logarithm of the number x
 */
double taylor_log(double x, unsigned int n)
{
	// TODO: x <= 0
	if (x == 0.0) {
		return -INFINITY;
	} else if (x < 0) {
		return NAN;
	}

	double result = 0.0;
	if (x <= 1) {
		// algorithm for (0,1>
		x = 1 - x;
		for (unsigned int i = 1; i <= n; i++) {
			result -= pow_positive(x, i) / i;
		}
	} else {
		// algorithm for (1,INF)
		for (; n > 0; n--) {
			result += pow_positive((x - 1.0) / x, n) / n;
		}
	}

	return result;
}


/**
 * calculating the natural logarithm of the number `x` in `n` iterations using continued fractions
 *
 * @param x real number // TODO: <= 0
 * @param n number of steps in continued fraction // TODO: = 0
 * @return logarithm of the number x
 */
double cfrac_log(double x, unsigned int n)
{
	// TODO: x <= 0
	if (x == 0.0) {
		return -INFINITY;
	} else if (x < 0) {
		return NAN;
	} else if (x == 1.0) {
		return  0.0;
	}

	x = (x - 1.0) / (x + 1.0);
	double frac = 1.0, pow2x = x * x;
	for (; n > 0; n--) {
		frac = 2.0 * n - 1.0 - n * n * pow2x / frac;
	}

	return 2.0 * x / frac;
}


/**
 * calculate exponential function of the number `y` with the general basis `x` in `n` iterations
 * for the calculation of the natural logarithm uses function taylor_log
 *
 * @see taylor_log
 * @param x base // TODO: <= 0
 * @param y exponent
 * @param n number of members in polynomial // TODO: = 0
 * @return the base to the exponent power
 */
double taylor_pow(double x, double y, unsigned int n)
{
	double result = 0.0, log_x = taylor_log(x, n);
	for (; n > 0; n--) {
		result += pow_positive(y, n) * pow_positive(log_x, n) / factorial(n);
	}

	return 1.0 + result;
}


/**
 * calculate exponential function of the number `y` with the general basis `x` in `n` iterations
 * for the calculation of the natural logarithm usees function cfrac_log
 *
 * @see cfrac_log
 * @param x base // TODO: <= 0
 * @param y exponent
 * @param n number of members in polynomial // TODO: = 0
 * @return the base to the exponent power
 */
double taylorcf_pow(double x, double y, unsigned int n)
{
	double result = 0.0, log_x = cfrac_log(x, n);
	for (; n > 0; n--) {
		result += pow_positive(y, n) * pow_positive(log_x, n) / factorial(n);
	}

	return 1.0 + result;
}


/**
 * TODO
 *
 * @param x
 * @return
 */
double mylog(double x)
{
	(void) x;
	return 0.0;
}


/**
 * TODO
 *
 * @param x
 * @param y
 * @return
 */
double mypow(double x, double y)
{
	(void) x;
	(void) y;
	return 0.0;
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
 * converts value of argument in string format to int
 * if in value there is non-convertible part, or minus sign, print it to stderr
 *
 * @param value value of argument in string format (from argv)
 * @param arg name of argument that value belongs to
 * @param error set to true if in value there is non-convertible part or minus sign
 * @return converted int value
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
 * prints results of logarithm functions to stdout
 *
 * for params @see taylor_log and cfrac_log
 */
void print_log_results(double x, unsigned int n)
{
	char *formated_x, *formated_ls;
	asprintf(&formated_x, "%g", x);
	const short min_number_of_chars_ls = 12 + strlen(formated_x);

	//log
	asprintf(&formated_ls, "log(%s)", formated_x);
	printf("%*s = %.*g\n", min_number_of_chars_ls, formated_ls, RESULT_PRECISION_FOR_OUTPUT, log(x));
	//crac_log
	asprintf(&formated_ls, "cfrac_log(%s)", formated_x);
	printf("%*s = %.*g\n", min_number_of_chars_ls, formated_ls, RESULT_PRECISION_FOR_OUTPUT, cfrac_log(x, n));
	//taylor_log
	asprintf(&formated_ls, "taylor_log(%s)", formated_x);
	printf("%*s = %.*g\n", min_number_of_chars_ls, formated_ls, RESULT_PRECISION_FOR_OUTPUT, taylor_log(x, n));
}


/**
 * prints results of exponencial functions to stdout
 *
 * for params @see taylor_pow and taylorcf_pow
 */
void print_pow_results(double x, double y, unsigned int n)
{
	char *formated_x, *formated_y, *formated_ls;
	asprintf(&formated_x, "%g", x);
	asprintf(&formated_y, "%g", y);
	const short min_number_of_chars_ls = 15 + strlen(formated_x) + strlen(formated_y);

	//pow
	asprintf(&formated_ls, "pow(%s,%s)", formated_x, formated_y);
	printf("%*s = %.*g\n", min_number_of_chars_ls, formated_ls, RESULT_PRECISION_FOR_OUTPUT, pow(x, y));
	//crac_log
	asprintf(&formated_ls, "taylor_pow(%s,%s)", formated_x, formated_y);
	printf("%*s = %.*g\n", min_number_of_chars_ls, formated_ls, RESULT_PRECISION_FOR_OUTPUT, taylor_pow(x, y, n));
	//taylor_log
	asprintf(&formated_ls, "taylorcf_pow(%s,%s)", formated_x, formated_y);
	printf("%*s = %.*g\n", min_number_of_chars_ls, formated_ls, RESULT_PRECISION_FOR_OUTPUT, taylorcf_pow(x, y, n));
}


/**
 * process input arguments and print results to stdout
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
			// TODO: x <= 0
			unsigned int n = value_of_arg_to_unsigned_int(argv[3], "--log N", &error);
			if (error) {
				return false;
			}
			// TODO: n = 0
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
			// TODO: x <= 0
			double y = value_of_arg_to_double(argv[3], "--pow Y", &error);
			if (error) {
				return false;
			}
			unsigned int n = value_of_arg_to_unsigned_int(argv[4], "--pow N", &error);
			if (error) {
				return false;
			}
			// TODO: n = 0
			if (n == 0) {
				PRINT_ERR("Value of argument --log N must be greater than 0, given %i.", n);
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
