/**
 * @file message.c
 * @author Karin Dorman, kdorman@iastate.edu
 *
 * Functions for outputting error and debugging messages.
 */

#include "message.h"

/**
 * Write uniformly formatted message to requested stream.
 * Write uniformly formatted debug or error messages to requested stream.
 * Certain common errors/problems issue default messages so they need not be
 * typed repeatedly and nonstandardly.  Returns msg_id so that the function
 * can be conveniently called as return message(...) or exit(message(...)).
 *
 * @param fp file handle where message to be written
 * @param file_name name of offending source file
 * @param fxn_name name of offending function
 * @param line line number of offending function
 * @param msg_type urgency of message
 * @param msg_id identify one of a default set of canned messages
 * @param msg formatted string to output via vprintf
 * @return msg_id
 */
int message(FILE *fp, const char *file_name, const char *fxn_name, int line,
	int msg_type, int msg_id, const char *msg, ...) {
	int nsec;
	va_list args;
	fprintf(fp, "%s [%s::%s(%d)]: ",
		msg_type == INFO_MSG ? "INFO" : msg_type == DEBUG_MSG ? "DEBUG"
			: msg_type == WARNING_MSG ? "WARNING" : "ERROR", 
		file_name, fxn_name, line);
	if (msg_id == NO_ERROR) {
		va_start(args, msg);
		vfprintf(fp, msg, args);
		va_end(args);
	} else {
		switch(msg_id) {
			case MEMORY_ALLOCATION: 
				if (msg) {
					fprintf(fp, "could not allocate ");
					va_start(args, msg);
					vfprintf(fp, msg, args);
					va_end(args);
				} else
					fprintf(fp, "memory allocation error\n");
				break;
			case INVALID_CMD_OPTION:
				fprintf(fp, "unrecognized command option");
				if (msg) {
					fprintf(fp, ": ");
					va_start(args, msg);
					vfprintf(fp, msg, args);
					va_end(args);
				} else
					fprintf(fp, "\n");
				break;
			case INVALID_CMD_ARGUMENT:
				fprintf(fp, "invalid argument to command option");
				if (msg) {
					fprintf(fp, ": ");
					va_start(args, msg);
					vfprintf(fp, msg, args);
					va_end(args);
				} else
					fprintf(fp, "\n");
				break;
			case INVALID_CMDLINE:
				fprintf(fp, "[invalid command line] %s\n", msg);
				break;
			case INVALID_USER_SETUP:
				fprintf(fp, "[invalid user choice] %s\n", msg);
				break;
			case FILE_OPEN_ERROR:
				fprintf(fp, "could not open file \"%s\"\n", msg);
				break;
			case FILE_NOT_FOUND:
				fprintf(fp, "file \"%s\" not found\n", msg);
				break;
			case FILE_FORMAT_ERROR:
				fprintf(fp, "invalid file format");
				if (msg) {
					fprintf(fp, ": ");
					va_start(args, msg);
					vfprintf(fp, msg, args);
					va_end(args);
				} else
					fprintf(fp, "\n");
				break;
			case END_OF_FILE:
				fprintf(fp, "unexpected end of file in file \"%s\"\n", msg);
				break;
			case INTERNAL_MISMATCH:
				fprintf(fp, "[internal mismatch] %s\n", msg);
				break;
			case OUT_OF_TIME:
				fprintf(fp, "out of time");
				if (msg) {
					va_start(args, msg);
					nsec = va_arg(args, int);
					va_end(args);
					fprintf(fp, " (%slimit %02d:%02dm)",
						msg, (int)(nsec/3600),
						(int)((nsec%3600)/60));
				}
				fprintf(fp, "\n");
				break;
			case MEMORY_USAGE_LIMIT:
				fprintf(fp, "exceed memory limit");
				if (msg) {
					fprintf(fp, ": ");
					va_start(args, msg);
					vfprintf(fp, msg, args);
					va_end(args);
				} else
					fprintf(fp, "\n");
				break;
			default:
				va_start(args, msg);
				vfprintf(fp, msg, args);
				va_end(args);
		}
	}
	return msg_id;
} /* message */
