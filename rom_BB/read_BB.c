

#include "read_BB.h"

static const char* CODENAME = "FE_sys/read >";
static const int BUFFER_SIZE = 10000;


FILE* ROM_BB_read_fopen(
		FILE*        fp,
		const char*  filename,
		const char*  directory)
{
	char fname[BUFFER_SIZE];
	snprintf(fname, BUFFER_SIZE, "%s/%s", directory, filename);

	fp = fopen(fname, "r");
	if( fp == NULL ) {
		printf("%s ERROR: File \"%s\" cannot be opened.\n",
				CODENAME, fname);
		exit(EXIT_FAILURE);
	}
	else {
		printf("%s Reading file \"%s\".\n", CODENAME, fname);
	}

	return fp;
}


FILE* ROM_BB_read_fopen_without_error(
		FILE*        fp,
		const char*  filename,
		const char*  directory)
{
	char fname[BUFFER_SIZE];
	snprintf(fname, BUFFER_SIZE, "%s/%s", directory, filename);
	fp = fopen(fname, "r");
	printf("%s Reading file \"%s\".\n", CODENAME, fname);

	return fp;
}

