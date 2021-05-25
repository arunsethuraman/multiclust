#ifndef __H_MESSAGE__
#define __H_MESSAGE__

#include <stdio.h>
#include <stdarg.h>
#include <errno.h>
#include <sys/sysinfo.h>
#include "message.h"

/** Types of messages.
 */
enum {	NO_MSG,		/*!< no message */
	INFO_MSG,	/*!< informative message */
	DEBUG_MSG,	/*!< debugging message */
	WARNING_MSG,	/*!< warning messag */
	ERROR_MSG	/*!< error message */
};

/** Types of errors.
 */
enum {NO_ERROR,			/*!< no error */
	CUSTOM_ERROR,		/*!< customer error */
	NO_DATA,		/*!< necessary data not provided */
	MEMORY_ALLOCATION,	/*!< memory allocation error */
	FILE_NOT_FOUND,		/*!< file not found */
	FILE_OPEN_ERROR,	/*!< file open error */
	END_OF_FILE,		/*!< premature end of file error */
	FILE_FORMAT_ERROR,	/*!< invalid file format error */
	INVALID_CMDLINE,	/*!< invalid command line, specific below */
	INVALID_CMD_OPTION,	/*!< invalid command-line option */
	INVALID_CMD_ARGUMENT,	/*!< invalid argument to command-line option */
	INVALID_USER_SETUP,	/*!< invalid user setup */
	INTERNAL_MISMATCH,	/*!< data does not match in some way it must */
	INTERNAL_ERROR,		/*!< internal error */
	CLUSTER_SIZE_OVERFLOW,	/*!< overflow in no. of clusters? */
	STATE_SPACE_OVERFLOW,	/*!< memory overflow due to state space size */
	OUT_OF_TIME,		/*!< ran out of time */
	MEMORY_USAGE_LIMIT,	/*!< ran out of memory */
	MEMCPY_ERROR,		/*!< error in call to memcpy() */
	NUM_ERRORS		/*!< total number of errors */
};

/** Level of verbosity.
 */
enum {	ABSOLUTE_SILENCE,	/*!< only output through files */
	SILENT,			/*!< no verbosity; final output only */
	QUIET,			/*!< try to be quiet */
	MINIMAL,		/*!< minimal verbosity */
	RESTRAINED,		/*!< some output */
	TALKATIVE,		/*!< talkative output */
	VERBOSE,		/*!< verbose output */
	DEBUG			/*!< debugging output */
};

#define CHECK_TIME(start_time, time_lim, file_name, fxn_name) {                \
	errno = NO_ERROR;                                                      \
	if ((time_lim) && difftime(time(NULL), (start_time)) > (time_lim)) {   \
		message(stderr, (file_name), (fxn_name), __LINE__, ERROR_MSG,  \
			OUT_OF_TIME, "", (time_lim));                          \
		errno = OUT_OF_TIME;                                           \
	}                                                                      \
}

#define CHECK_MEMORY(s_info, mem_lim, mem_new, file_name, fxn_name) {         \
	errno = NO_ERROR;                                                      \
	if ((mem_lim) && !sysinfo(&(s_info))) {                               \
		unsigned long mem_b = ((s_info).totalram - (s_info).freeram)  \
			* (s_info).mem_unit + (mem_new);                      \
		if (mem_b > (mem_lim)) {                                      \
			errno = MEMORY_USAGE_LIMIT;                           \
			message(stderr, (file_name), (fxn_name), __LINE__,    \
				ERROR_MSG, errno, "request for %.1"           \
				DOUBLE_F_NFMT "Gb exceeds %.1" DOUBLE_F_NFMT  \
				"Gb limit\n", (DOUBLE) mem_b / 1e9,           \
				(DOUBLE) (mem_lim) / 1e9);                    \
		}                                                             \
	}                                                                     \
}

int message(FILE *, const char *, const char *, int, int, int, const char *, ...);

/**
 * Print a formatted message to stderr.
 */
#define mmessage(type, err, ...) message(stderr, __FILE__, __func__,  __LINE__, (type), (err),  __VA_ARGS__)


#endif
