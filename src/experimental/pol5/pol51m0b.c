#include "if.h"
#include "stage1/poly_stage1.h"

#define START_MESSAGE \
"----------------------------------------------------\n"\
"|    pol51m0b GNFS polynomial selection program    |\n"\
"| This program is copyright (c) 2005, by Thorsten  |\n"\
"| Kleinjung and Jens Franke, and is subject to the |\n"\
"| terms of the GNU GPL version 2.                  |\n"\
"| This program is part of gnfs4linux.              |\n"\
"----------------------------------------------------\n"

#define A5_SCALE_FACTOR 1000

/*------------------------------------------------------------------------*/
static char *
get_options(int argc, char **argv, poly_stage1_t *data)
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
			if (sscanf(argv[i++], "%lf", &data->norm_max) != 1)
				complain("Bad argument to -n!\n");
			break;
		case 'p':
			if (sscanf(argv[i++], "%u", &data->npr_in_p) != 1)
				complain("Bad argument to -p!\n");
			break;
		case 'a':
			mpz_set_str(data->gmp_a5_begin, argv[i++], 10);
			mpz_mul_ui(data->gmp_a5_begin, data->gmp_a5_begin, 
					A5_SCALE_FACTOR);
			break;
		case 'A':
			mpz_set_str(data->gmp_a5_end, argv[i++], 10);
			mpz_mul_ui(data->gmp_a5_end, data->gmp_a5_end, 
					A5_SCALE_FACTOR);
			break;
		case 'l':
			if (sscanf(argv[i++], "%u", &data->p0_limit) != 1)
				complain("Bad argument to -l!\n");
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
	if (data->npr_in_p < 4)
		complain("argument '-p' must be >=4\n");
	if (mpz_cmp_ui(data->gmp_a5_end, 0) == 0)
		complain("argument '-A' not specified\n");

	mpz_add_ui(data->gmp_a5_begin, data->gmp_a5_begin, 1);
	return base_name;
}

/*------------------------------------------------------------------------*/
void
read_data(poly_stage1_t *data, char *infile)
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
open_outputfile(poly_stage1_t *data, char *infile)
{
	char buf[256];

	sprintf(buf, "%s.51.m", infile);
	if ((data->outfile = fopen(buf, "a")) == NULL)
		complain("File not found: %s\n", buf);
}

/*------------------------------------------------------------------------*/
void
close_outputfile(poly_stage1_t *data)
{
	fclose(data->outfile);
}

/*------------------------------------------------------------------------*/
int
main(int argc, char **argv)
{
	poly_stage1_t data;
	char *infile_name;

	printf("%s\n", START_MESSAGE);

	poly_stage1_init(&data);

	infile_name = get_options(argc, argv, &data);

	setbuf(stdout, NULL);
	read_data(&data, infile_name);
	open_outputfile(&data, infile_name);
	setbuf(data.outfile, NULL);

	if (poly_stage1_run(&data))
		printf("success\n");

	poly_stage1_free(&data);
	close_outputfile(&data);
	return 0;
}
