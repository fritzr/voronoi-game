#include <opencv2/core/core.hpp>

#if CV_VERSION_MAJOR < 3 && !defined(OPENCV2_4_13)
#include <opencv2/contrib/contrib.hpp>
typedef int ColormapTypes;
#endif

#if !defined(_WIN32) && !defined(_WIN64)
#include <unistd.h>
#include <getopt.h>
extern "C" {
	extern int opterr;    /* if error message should be printed */
	extern	int optind;    /* index into parent argv vector */
	extern int optopt;    /* character checked for validity */
	extern int optreset;  /* reset getopt */
	extern char    *optarg;
	extern int getopt(int nargc, char * const nargv[], const char *ostr);
};
#else
#error no getopt support
#endif
