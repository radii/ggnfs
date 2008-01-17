#include "if.h"
#include "poly_stage2.h"

#define START_MESSAGE \
"----------------------------------------------------\n"\
"|    pol51opt GNFS polynomial selection program    |\n"\
"| This program is copyright (c) 2005, by Thorsten  |\n"\
"| Kleinjung and Jens Franke, and is subject to the |\n"\
"| terms of the GNU GPL version 2.                  |\n"\
"| This program is part of gnfs4linux.              |\n"\
"----------------------------------------------------\n"

/*------------------------------------------------------------------------*/
static char *
get_options(int argc, char **argv, poly_stage2_t *data)
{
	char c;
	int i = 0;
	char *base_name = NULL;

	while (i < argc) {
		if (argv[i][0] != '-') {
			i++;
			continue;
		}

		c = argv[i++][1];
		switch (c) {
		case 'b':
			base_name = argv[i++];
			break;
		case 'n':
			if (sscanf(argv[i++], "%lf", &data->max_norm_1) != 1)
				complain("Bad argument to -n!\n");
			break;
		case 'N':
			if (sscanf(argv[i++], "%lf", &data->max_norm_2) != 1)
				complain("Bad argument to -n!\n");
			break;
		case 'F':
			if (sscanf(argv[i++], "%lf", &data->bound0) != 1)
				complain("Bad argument to -n!\n");
			break;
		case 'f':
			if (sscanf(argv[i++], "%lf", &data->bound1) != 1)
				complain("Bad argument to -n!\n");
			break;
		case 'P':
			if (sscanf(argv[i++], "%u", &data->p_bound) != 1)
				complain("Bad argument to -n!\n");
			break;
		case 'A':
			if (sscanf(argv[i++], "%lf", &data->area) != 1)
				complain("Bad argument to -n!\n");
			break;
		case 'e':
			if (sscanf(argv[i++], "%lf", &data->min_e) != 1)
				complain("Bad argument to -n!\n");
			break;
		case 'v':
			verbose++;
			break;
		default:
			fprintf(stderr, "Bad Option %c\n", (char) c);
			Schlendrian("");
		}
	}

	if (base_name == NULL)
		complain("argument '-b base_name' is necessary\n");
	return base_name;
}

/*------------------------------------------------------------------------*/
static void
read_data(poly_stage2_t *data, char *infile)
{
	char buf[256];
	FILE *fp;

	sprintf(buf, "%s.data", infile);
	if ((fp = fopen(buf, "r")) == NULL)
		complain("File not found: %s\n", buf);

	while (1) {
		fgets(buf, sizeof(buf), fp);
		if (buf[0] == 'N') {
			mpz_set_str(data->gmp_N, buf + 2, 10);
			break;
		}
		if (feof(fp)) {
			complain("error: no input number found");
		}
	}
	fclose(fp);
}

/*------------------------------------------------------------------------*/
static void
open_inputfile(poly_stage2_t *data, char *infile)
{
	char buf[256];

	sprintf(buf, "%s.51.m", infile);
	if ((data->infile = fopen(buf, "r")) == NULL)
		complain("File not found: %s\n", buf);
}

/*------------------------------------------------------------------------*/
static void
open_outputfile(poly_stage2_t *data, char *infile)
{
	char buf[256];

	sprintf(buf, "%s.cand", infile);
	if ((data->outfile = fopen(buf, "a")) == NULL)
		complain("Cannot open '%s'\n", buf);
}

/*------------------------------------------------------------------------*/
static void
close_inputfile(poly_stage2_t *data)
{
	fclose(data->infile);
}

/*------------------------------------------------------------------------*/
static void
close_outputfile(poly_stage2_t *data)
{
	fclose(data->outfile);
}

/*------------------------------------------------------------------------*/
int
main(int argc, char **argv)
{
	poly_stage2_t data;
	char *basename;

	printf("%s\n", START_MESSAGE);

	poly_stage2_init(&data);

	basename = get_options(argc, argv, &data);

	setbuf(stdout, NULL);
	read_data(&data, basename);
	open_inputfile(&data, basename);
	open_outputfile(&data, basename);

	if (poly_stage2_run(&data))
		printf("success\n");

	poly_stage2_free(&data);
	close_inputfile(&data);
	close_outputfile(&data);
	return 0;
}
