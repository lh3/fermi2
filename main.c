#include <stdio.h>
#include <string.h>

#define FM_VERSION "r1"

int main_diff(int argc, char *argv[]);
int main_sub(int argc, char *argv[]);

void liftrlimit(void);
double cputime(void);
double realtime(void);

int main(int argc, char *argv[])
{
	int ret, i;
	double t_start;
	liftrlimit();
	if (argc == 1) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   fermi2 <command> [arguments]\n\n");
		fprintf(stderr, "Command: diff     extract reads containing special k-mers\n");
		fprintf(stderr, "\n");
		return 1;
	}
	t_start = realtime();
	if (strcmp(argv[1], "diff") == 0) ret = main_diff(argc-1, argv+1);
	else {
		fprintf(stderr, "[E::%s] unknown command\n", __func__);
		return 1;
	}
	if (ret == 0) {
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, FM_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - t_start, cputime());
	}
	return ret;
}
