

#include "write_BB.h"

static const int BUFFER_SIZE = 10000;

FILE* ROM_BB_write_fopen(
		FILE*        fp,
		const char*  filename,
		const char*  directory)
{
	char fname[BUFFER_SIZE];
	snprintf(fname, BUFFER_SIZE, "%s/%s", directory, filename);

	fp = fopen(fname, "w");
	if( fp == NULL ) {
		printf("%s ERROR: File \"%s\" cannot be opened.\n", fname);
		exit(EXIT_FAILURE);
	}
	else {
		printf("Writing file \"%s\".\n", fname);
	}

	return fp;
}


FILE* ROM_BB_write_fopen_without_error(
		FILE*        fp,
		const char*  filename,
		const char*  directory)
{
	char fname[BUFFER_SIZE];
	snprintf(fname, BUFFER_SIZE, "%s/%s", directory, filename);

	fp = fopen(fname, "w");
	printf("Writing file \"%s\".\n", fname);

	return fp;
}


FILE* ROM_BB_write_add_fopen(
		FILE*        fp,
		const char*  filename,
		const char*  directory)
{
	char fname[BUFFER_SIZE];
	snprintf(fname, BUFFER_SIZE, "%s/%s", directory, filename);

	fp = fopen(fname, "a");
	if( fp == NULL ) {
		printf("ERROR: File \"%s\" cannot be opened.\n", fname);
		exit(EXIT_FAILURE);
	}
	else {
		printf("Writing (adding) file \"%s\".\n", fname);
	}

	return fp;
}

