#include <unistd.h>
#include <sys/fcntl.h>

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include <string.h>

#include "loran.h"

struct lc_chain lorchain;

void loran_info(struct lc_chain *chain)
{
	for(size_t i = 0; i < chain->station_cnt; ++i) {
		printf("LORSTA %u\n"
		       "\toffset: %u\n"
		       "\tComb filter P: %d\n",
		       i, chain->sta[i].offset, chain->sta[i].comb_kp);
	}
}

int main(int argc, char **argv)
{
	const char		*filename;

	if(argc < 2) {
		fprintf(stderr, "Usage: loran_datatest <filename>\n"
			"Where filename is path to file with samples,\n"
			"\tsample values are separated by whitespace\n");
		return EXIT_FAILURE;
	}
	filename = argv[1];

	FILE			*dfile;
	dfile = fopen(filename, "r");

	if(!dfile) {
		perror("data file fopen()");
		return EXIT_FAILURE;
	}

	/* Initialize Loran-C structs */
	//lorchain.gri = 8000;
	lorchain.gri = 8000;
	lorchain.refedge = &lc_refedge_rsdn310;
	lorchain.station_cnt = 8;
	lc_init(&lorchain);

	/* Initialize buffer */
	size_t			samplebufsize	= 128;
	lc_type_sample		*samplebuf, *sampleptr, *samplebufend;
	samplebuf = malloc(samplebufsize * sizeof(lc_type_sample));
	if(!samplebuf) {
		fprintf(stderr, "malloc() failed");
		return EXIT_FAILURE;
	}
	samplebufend = samplebuf + samplebufsize;
	sampleptr = samplebuf;

	size_t skipcnt	= 160;//22000;
	printf("Skipping %u first samples\n", skipcnt);

	/* Read all the samples into memory */
	printf("Reading samples from file..\n");
	while(1) {
		int ret;
		/* Read from file, fill buffer */
		while(sampleptr < samplebufend) {
			unsigned int val;
			ret = fscanf(dfile, "%u", &val);
			if(ret > 0) {
				if(skipcnt) {
					skipcnt--;
					continue;
				}
				*sampleptr++ = val;
			} else {
				break;
			}
		}
		size_t		samplecnt;
		samplecnt = sampleptr - samplebuf;
		printf("Already read %u samples...\n", samplecnt);
		if(ret <= 0) {
			if(ret == 0)
				fprintf(stderr, "Error in input data format\n");
			if(errno != 0)
				perror("data file fscanf()");
			break;
		}

		/* Increase buffer size by factor of 2 */
		samplebufsize*= 2;
		samplebuf = realloc(samplebuf,
				    samplebufsize * sizeof(lc_type_sample));

		if(!samplebuf) {
			fprintf(stderr, "realloc() failed");
			return EXIT_FAILURE;
		}

		samplebufend = samplebuf + samplebufsize;
		sampleptr = samplebuf + samplecnt;
	}

	size_t			total_samples = sampleptr - samplebuf;
	printf("Total samples read: %u\n", total_samples);
	printf("Allocated buffer size (samples): %u\n", samplebufsize);
	fclose(dfile);

	size_t			fri_cnt;
	fri_cnt = total_samples / lorchain.frin;
	printf("%u full FRIs available for testing\n", fri_cnt);

	size_t			fri_tested = 35000;
	printf("We are going to run algorithm for %u FRIs..\n", fri_tested);
	size_t			partsize = 440;
	printf("..passing samples in chunks of %u samples each\n", partsize);

	size_t			fri_counter = 0;
	lc_type_sample		*sample_last;
	sample_last = samplebuf + fri_cnt * lorchain.frin;
	sampleptr = samplebuf;

	while(1) {
		size_t samplecnt = sample_last - sampleptr > (int) partsize ?
				    partsize : (size_t) (sample_last - sampleptr);
		lc_new_samples(&lorchain, sampleptr, samplecnt);
		sampleptr+= samplecnt;
		if(sampleptr == sample_last) {
			fri_counter+= fri_cnt;
			printf("[Passed %6u FRIs]\n", fri_counter);
			if(fri_counter >= fri_tested)
				break;
			sampleptr = samplebuf;
		}
	}

//	/* Output windows content to files... */
//	char			fname[128];
//	for(size_t sidx = 0; sidx < lorchain.station_cnt; ++sidx) {
//		snprintf(fname, 127, "data_win_%u.txt", sidx);
//		dfile = fopen(fname, "wb");
//		if(!dfile) {
//			perror("fopen() while writing windows contents");
//			continue;
//		}
//		for(size_t i = 0; i < LC_GRCNT * LC_STA_GRIWCNT; ++i) {
//			struct lc_subwin *win = &lorchain.sta[sidx].wnd[i];

//			for(size_t j = 0; j < LC_STA_WINSZ; ++j)
//				fprintf(dfile, "%hd ", win->data[j]);
//			fputs("\n", dfile);
//		}
//		fclose(dfile);
//	}
	loran_info(&lorchain);
	printf("\nEND\n");
	return EXIT_SUCCESS;
}

