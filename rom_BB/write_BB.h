#pragma once

#include <stdio.h>
#include <stdlib.h>

FILE* ROM_BB_write_fopen(
		FILE*        fp,
		const char*  filename,
		const char*  directory);

FILE* ROM_BB_write_fopen_without_error(
		FILE*        fp,
		const char*  filename,
		const char*  directory);

FILE* ROM_BB_write_add_fopen(
		FILE*        fp,
		const char*  filename,
		const char*  directory);
