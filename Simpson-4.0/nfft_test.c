/*
 * nfft_test.c
 *
 *  Created on: 10.2.2012
 *      Author: zdenek
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "fftw3.h"
#include "nfft3util.h"
#include "nfft3.h"
#include "complx.h"


int simpson_nfft_test(void) {

	FILE *fp_re, *fp_im, *fp_or;
	char line[512];
	double d;
	int i,j, num_threads;
	nfsft_plan sourceplan, targetplan;
	int Jsmax, Jsovers, N1, N2;
	//double complex *z;
	complx *z;
	double *or1, *or2;

	int dim = 4;
	const char name1[] = "D:\\Tosner\\Arhus12\\interpolation\\LEBh295_degree.dat";
	const char name2[] = "D:\\Tosner\\Arhus12\\interpolation\\ROSELEBh11617_degree.dat";

	// read source orientations
	fp_or = fopen(name1,"r");
	if (!fp_or) {
		fprintf(stderr,"Error opening file %s\n",name1);
		return 1;
	}
	fgets(line,512,fp_or);
	if (sscanf(line,"%d",&N1) != 1) {
		fprintf(stderr,"Error when reading # of elements from %s\n",name1);
		return 1;
	}
	or1 = (double*)malloc(N1*3*sizeof(double));
	for (i=0; i<N1; i++) {
		for (j=0; j<3; j++) {
			if (fscanf(fp_or,"%lg",&d) != 1) {
				fprintf(stderr,"Error when reading %s, element (%d, %d)\n",name1,i,j);
				return 1;
			}
			or1[i*3+j] = d;
		}
	}
	fclose(fp_or);
	printf("Source orientations read (%d, first = %lg, last = %lg\n",N1,or1[0],or1[(N1-1)*3+2]);
	// read target orientations
	fp_or = fopen(name2,"r");
	if (!fp_or) {
		fprintf(stderr,"Error opening file %s\n",name2);
		return 1;
	}
	if (fscanf(fp_or,"%d",&N2) != 1) {
		fprintf(stderr,"Error when reading # of elements from %s\n",name2);
		return 1;
	}
	or2 = (double*)malloc(N2*3*sizeof(double));
	for (i=0; i<N2; i++) {
		for (j=0; j<3; j++) {
			if (fscanf(fp_or,"%lg",&d) != 1) {
				fprintf(stderr,"Error when reading %s, element (%d, %d)\n",name2,i,j);
				return 1;
			}
			or2[i*3+j] = d;
		}
	}
	fclose(fp_or);
	printf("Target orientations read (%d, first = %lg, last = %lg\n",N2,or2[0],or2[(N2-1)*3+2]);

	// read data to interpolate (from Baltzar)
	fp_re = fopen("D:\\Tosner\\Arhus12\\interpolation\\source_propprod.dat","r");
	//fp_im = fopen("D:\\Tosner\\Arhus12\\interpolation\\q_lams_LEB217_im.dat","r");
	//if (!fp_re || !fp_im) {
	//	fprintf(stderr,"Error opening files...\n");
	//	return 1;
	//}

	//z = (double complex *)malloc(N1*sizeof(double complex));
	z = (complx *)malloc(N1*sizeof(complx));
	//memset(z, 0, N1*sizeof(double complex));
	memset(z, 0, N1*sizeof(complx));
	for (i=0; i<N1; i++) {
		for (j=0;j<6;j++) {
			if (fscanf(fp_re,"%lg",&d) != 1) {
				fprintf(stderr,"Error when reading real element (%d, %d)\n",i,j);
				return 1;
			}
			//if (j == 4) z[i] += d;
			if (j == 4) z[i].re = d;
			//if (fscanf(fp_im,"%lg",&d) != 1) {
			//	fprintf(stderr,"Error when reading imaginary element (%d, %d)\n",i,j);
			//	return 1;
			//}
			//if (j == 5) z[i]  += d*I;
			if (j == 5) z[i].im  = d;
		}
	}
	fclose(fp_re);
	//fclose(fp_im);

	Jsmax = ((int)(sqrt(6*(N1-1))) - 1) / 2;
	Jsovers = 2 << (int)(floor(log2((double)Jsmax)));
	printf("Jsmax = %d, Jsovers = %d\n",Jsmax, Jsovers);



	nfsft_precompute(Jsovers, 1000.0, 0U, 0U);

	/* NFSFT_NORMALIZED shall not be included*/
	nfsft_init_guru(&sourceplan, Jsovers, N1, NFSFT_NO_DIRECT_ALGORITHM | NFSFT_MALLOC_X | NFSFT_MALLOC_F |
			NFSFT_MALLOC_F_HAT | NFSFT_DESTROY_X | NFSFT_DESTROY_F | NFSFT_DESTROY_F_HAT,
			PRE_PHI_HUT | PRE_PSI | FFTW_INIT | FFT_OUT_OF_PLACE, 6);
	printf("sourceplan.M_total = %d, N_total = %d, N = %d\n",sourceplan.M_total, sourceplan.N_total, sourceplan.N);

	for (i=0; i<N1; i++) {
		d = or1[i*3];
		if (d < 180.0 ) d /= 360.0; else d = d/360.0 - 1.0;
		sourceplan.x[2*i]   = d;
		sourceplan.x[2*i+1] = or1[i*3+1]/360.0;
	}

	/* Do pre-computation for nodes. */
	nfsft_precompute_x(&sourceplan);

	/* used for initializing the odd rank terms*/
	//memset(sourceplan.f_hat, 0, sourceplan.N_total * sizeof(double complex));
	memset(sourceplan.f_hat, 0, sourceplan.N_total * sizeof(complx));



	nfsft_init_guru(&targetplan, Jsovers, N2, NFSFT_NO_DIRECT_ALGORITHM | NFSFT_MALLOC_X | NFSFT_MALLOC_F |
			NFSFT_MALLOC_F_HAT |  NFSFT_DESTROY_X | NFSFT_DESTROY_F | NFSFT_DESTROY_F_HAT,
			PRE_PHI_HUT | PRE_PSI | FFTW_INIT | FFT_OUT_OF_PLACE, 6);
	printf("targetplan.M_total = %d, N_total = %d, N = %d\n",targetplan.M_total, targetplan.N_total, targetplan.N);

	for (i=0; i<N2; i++) {
		d = or2[i*3];
		if (d < 180.0 ) d /= 360.0; else d = d/360.0 - 1.0;
		targetplan.x[2*i]   = d;
		targetplan.x[2*i+1] = or2[i*3+1]/360.0;
	}

	/* Do pre-computation for nodes. */
	nfsft_precompute_x(&targetplan);
	/* used for initializing the odd rank terms*/
	//memset(targetplan.f_hat, 0, targetplan.N_total * sizeof(double complex));
	memset(targetplan.f_hat, 0, targetplan.N_total * sizeof(complx));

	// copy data to be interpolated
	//memcpy(sourceplan.f, z, sourceplan.M_total*sizeof(double complex));
	memcpy(sourceplan.f, z, sourceplan.M_total*sizeof(complx));
	nfft_vpr_complex(sourceplan.f,10,"original source, vector f(1:10)");
	// do the adjoint transform on source
	nfsft_adjoint(&sourceplan);
	/* Copy only even rank terms, the rest shall remain zero*/

	for( i = 0; i <= Jsmax; i += 2) {
		double normalization = 2. * ((double) i) + 1.;
		//double normalization = 1.0;
		for( j = -i; j <= i; j++) {
			//targetplan.f_hat[NFSFT_INDEX(i, j, &targetplan)] = sourceplan.f_hat[NFSFT_INDEX(i, j, &sourceplan)] * normalization;
			//printf("t.f_hat(%d).re = %g * %g = %g\n",NFSFT_INDEX(i, j, &targetplan),sourceplan.f_hat[NFSFT_INDEX(i, j, &sourceplan)][0],normalization,targetplan.f_hat[NFSFT_INDEX(i, j, &targetplan)][0]);
			targetplan.f_hat[NFSFT_INDEX(i, j, &targetplan)][0] = sourceplan.f_hat[NFSFT_INDEX(i, j, &sourceplan)][0] * normalization;
			targetplan.f_hat[NFSFT_INDEX(i, j, &targetplan)][1] = sourceplan.f_hat[NFSFT_INDEX(i, j, &sourceplan)][1] * normalization;
			//printf("t.f_hat(%d).im = %g * %g = %g\n",NFSFT_INDEX(i, j, &targetplan),sourceplan.f_hat[NFSFT_INDEX(i, j, &sourceplan)][1],normalization,targetplan.f_hat[NFSFT_INDEX(i, j, &targetplan)][1]);
		}
	}
	/*
	for( j = 0; j <= sourceplan.N_total; j++) {
		targetplan.f_hat[j][0] = sourceplan.f_hat[j][0];
				//printf("t.f_hat(%d).re = %g * %g = %g\n",NFSFT_INDEX(i, j, &targetplan),sourceplan.f_hat[NFSFT_INDEX(i, j, &sourceplan)][0],normalization,targetplan.f_hat[NFSFT_INDEX(i, j, &targetplan)][0]);
		targetplan.f_hat[j][1] = sourceplan.f_hat[j][1];
				//printf("t.f_hat(%d).im = %g * %g = %g\n",NFSFT_INDEX(i, j, &targetplan),sourceplan.f_hat[NFSFT_INDEX(i, j, &sourceplan)][1],normalization,targetplan.f_hat[NFSFT_INDEX(i, j, &targetplan)][1]);
	}
	*/
	nfft_vpr_complex(sourceplan.f_hat,10,"source fHAT");
	//printf("3: %lg, %lg\n",creal(sourceplan.f_hat[3]),cimag(sourceplan.f_hat[3]));
	nfft_vpr_complex(targetplan.f_hat,10,"target fHAT");
	// do the transformation on target
	nfsft_trafo(&targetplan);
	// copy out the result from targetplan.f, size is targetplan.M_total of complx
	nfft_vpr_complex(targetplan.f,20,"interpolated target, vector f(1:20)");


	nfsft_finalize(&sourceplan);
	nfsft_finalize(&targetplan);

	nfsft_forget();




	free(or1);
	free(or2);
	free(z);

	return 0;
}


void simpson_fftw_test(void) {

    fftw_complex *in, *out;
    fftw_plan p;
    int i, N = 32;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for (i=0; i<N; i++) {
    	in[i][0] = sin((double)i);
    	in[i][1] = cos((double)i);
    }
    printf("FFTW in filled\n");

    fftw_execute(p); /* repeat as needed */

    printf("transform done\n");

    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);

}
