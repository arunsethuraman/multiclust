/**
 * @file cline.c
 * @author Karin S. Dorman, kdorman@iastate.edu
 * @date Tue Mar 20 18:10:33 CDT 2012
 *
 * Some functions for processing a command line.  Note, the programmer must
 * provide the definition of a function void fprintf_usage(FILE *, const char
 * *, void *) in order to use this code.  This function, as its name suggests,
 * prints help on command-line usage for the user.
 *
 * Note: The requirement used to be that the user defined print_usage(const
 * char *).  To compile older code under that obsolete system, use command-line
 * option -DOLD_USAGE with gcc.
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "cline.h"

/**
 * Print error message, usage information, and return a failure code.
 *
 * @param argv list of comman-line arguments as (unmutable) string constants
 * @param i index of argument that could not be parsed
 * @param obj void pointer to additional usage information
 * @return (simple) error code
 */
int usage_error(const char **argv, int i, void *obj)
{
	fprintf(stderr, "ERROR -- incorrect command line usage\n");
	if (i >= 0)
		fprintf(stderr, "Could not parse command line option '%s'.\n",
			argv[i]);
	fprintf(stderr, "See correct usage below.\n");
#ifdef OLD_USAGE	/* define to compile with older source code */
	print_usage(argv[0]);
#else
	fprint_usage(stderr, argv[0], obj);
#endif
	fprintf(stderr, "\nERROR -- incorrect command line usage\n");
	if (i >= 1)
		fprintf(stderr, "Could not parse command line option '%s'.\n",
			argv[i]);
	return EXIT_FAILURE;
} /* usage_error */

/**
 * Read an integer from the next argument in the command-line argument
 * list.  This code avoids some repetition in main()
 *
 * @param argc the number of arguments
 * @param argv the list of command-line arguments strings
 * @param i the argument we are to parse at int
 * @param obj void pointer to additional usage information
 * @return the parsed integer or errno is error
 */
int read_int(int argc, const char **argv, int i, void *obj)
{
	int n;
	char *ret_ptr;

	/* check for end of argument list; assignment to errno BEFORE return */
	if (i == argc)
		return (errno = usage_error(argv, i-1, obj));
	/* It turns out errno is not set to EINVAL (in C99) when the read fails
	 * as I said in class.  The only way to detect failure to read a 
	 * numeric value is to use the second argument.   The second argument
	 * is set to point to the char immediately following the read number.
	 * If it is the same as the address in the first pointer passed, then
	 * you know the read was unsuccessful.  Note the need to pass a pointer
	 * to a pointer in order to update the value (address) stored in the
	 * pointer variable. */
	/* n = strtol(argv[i], NULL, 0); */
	n = strtol(argv[i], &ret_ptr, 0);
	if (errno || ret_ptr == argv[i])
		return (errno = usage_error(argv, i, obj));
	return n;
} /* read_int */

/**
 * Read an positive integer from the next argument in the command-line argument
 * list.  This code avoids some repetition in main()
 *
 * @param argc the number of arguments
 * @param argv the list of command-line arguments strings
 * @param i the argument we are to parse at int
 * @param obj void pointer to additional usage information
 * @return the parsed unsigned integer or errno is error
 */
unsigned int read_uint(int argc, const char **argv, int i, void *obj)
{
	unsigned int n;
	char *ret_ptr;

	/* check for end of argument list; assignment to errno BEFORE return */
	if (i == argc)
		return (errno = usage_error(argv, i-1, obj));
	/* It turns out errno is not set to EINVAL (in C99) when the read fails
	 * as I said in class.  The only way to detect failure to read a 
	 * numeric value is to use the second argument.   The second argument
	 * is set to point to the char immediately following the read number.
	 * If it is the same as the address in the first pointer passed, then
	 * you know the read was unsuccessful.  Note the need to pass a pointer
	 * to a pointer in order to update the value (address) stored in the
	 * pointer variable. */
	/* n = strtol(argv[i], NULL, 0); */
	n = strtoul(argv[i], &ret_ptr, 0);
	if (errno || ret_ptr == argv[i])
		return (errno = usage_error(argv, i, obj));
	return n;
} /* read_uint */

/**
 * Read a long from the next argument on the command-line.
 *
 * @param argc the number of arguments
 * @param argv the list of command-line arguments strings
 * @param i the argument we are to parse at int
 * @param obj void pointer to additional usage information
 * @return the parsed long or errno is error
 */
long read_long(int argc, const char **argv, int i, void *obj)
{
	long n;
	char *ret_ptr;

	if (i == argc)
		return (errno = usage_error(argv, i-1, obj));

	n = strtol(argv[i], &ret_ptr, 0);
	if (errno || ret_ptr == argv[i])
		return (errno = usage_error(argv, i, obj));

	return n;
} /* read_long */

/**
 * Read a double from the next argument in the command-line argument
 * list.  This code avoids repetition in main().
 *
 * @param argc the number of arguments
 * @param argv the list of command-line arguments strings
 * @param i the argument we are to parse at int
 * @param obj void pointer to additional usage information
 * @return the parsed integer or errno is error
 */
double read_double(int argc, const char **argv, int i, void *obj)
{
	double d;
	char *ret_ptr;

	if (i==argc)
		return (errno = usage_error(argv, i-1, obj));
	
	d = strtod(argv[i], &ret_ptr);
	if (errno || ret_ptr == argv[i])
		return (errno = usage_error(argv, i, obj));
	
	return d;
} /* read_double */

/**
 * Read a char from the next argument on the command-line.
 *
 * @param argc number of arguments on command-line
 * @param argv list of command-line argument strings
 * @param i the argument we are to parse as char
 * @param obj void pointer to additional usage information
 * @return the parsed char or errno is error
 */
char read_char(int argc, const char **argv, int i, void *obj)
{
	if (i==argc)
		return (errno = usage_error(argv, i-1, obj));
	
	return argv[i][0];
} /* read_char */
