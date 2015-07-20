#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include "Maxfiles.h"
#include "MaxSLiCInterface.h"

#define TEST_MAIN
#include "su3.h"
#include "global.h"


/************* Calls to DFE ********************/
void init_dfe(void);
void transfer_spinors_to_dfe (spinor *in, int address);
void transfer_gauges_to_dfe (su3 *in, int address);
void apply_dirac (int in_address, int out_address, int gauge_address, int p_address, int ieo, int doSub, int isign);
void transfer_spinors_to_cpu (spinor *out, int address);
void unload_dfe(void);

/************* Verifying Functions ***************/
int verify_results (spinor* dfe_out, spinor *expected_out, int V);
int AreNotSameComplex(complex float a, float complex b);
int compare_spinor (spinor *a, spinor *b);

/************ Utility Functions *********************/
void create_random_input(spinor* s, su3* u);
void create_random_spinor(spinor * s);
void create_random_su3vector(su3_vector *v);
void create_random_su3(su3 *s);
void create_random_complex(float complex *a);
void print_spinors (spinor* s);
void print_gauges (su3* s);

/****************** Functions for reading data from files ***************/
void read_spinor(char * filename, spinor *out);
void read_gauge(char * filename, su3 *s);

/****************** Functions for reordering data in memory ***************/
void reorganize_gauge(su3 const * const in, su3 * const out, int ieo);


static max_file_t *maxfile;
static max_engine_t *engine;
static int burstsPerGaugeT, burstsPerSpinorT;
static double beta_s, beta_t_b, beta_t_f, mass;


int main(void) {

	T = S_LQCD_T;
	LX = S_LQCD_LX;
	LY = S_LQCD_LY;
	LZ  = S_LQCD_LZ;
	VOLUME = LZ * LY * LX * T;

	NUM_PIPES = S_LQCD_numPipes;
	LOOP_OFFSET = S_LQCD_loopOffset;

	beta_s   = -0.5;
	beta_t_f = 0.3;
	beta_t_b = 0.3;
	mass     = 0.1;

	spinor *in, *out, *out_dfe, *out_expected, *tmp;
	su3 *u0, *u0_re;
	su3 *u1, *u1_re;

	printf("Allocating memory for data ...");

	in = malloc(VOLUME * sizeof(spinor));
	out = &in[VOLUME / 2];
	out_dfe = malloc(VOLUME/2 * sizeof(spinor));
	out_expected = malloc(VOLUME/2 * sizeof(spinor));


	tmp = malloc(VOLUME/2 * sizeof(spinor));

	u0 = malloc (VOLUME * 8 * sizeof(su3));
	u0_re = &u0[VOLUME/2 * 8];
	u1 = malloc (VOLUME * 8 * sizeof(su3));
	u1_re = &u1[VOLUME/2 * 8];

	printf("Done!\n");

	/*printf("Creating random spinor and gauge inputs ...");

	for (int i=0 ; i<VOLUME/2 ; i++ ) {
		create_random_spinor(in + i);
	}
	for (int i=0 ; i<VOLUME*2 ; i++ ) {
		create_random_su3(ueven + i);
		create_random_su3(uodd + i);
	}

	printf("Done!\n");*/

	printf("Reading spinor and gauge inputs ...");

	read_spinor("tmp_spinor.txt", tmp);
	read_spinor("in_spinor.txt", in);
	read_spinor("out_spinor.txt", out_expected);
	read_gauge("in_gauge0.txt", u0);
	read_gauge("in_gauge1.txt", u1);

	printf("Done!\n");
	printf("Data reordering and adding necessary halos ...");

	reorganize_gauge(u0, u0_re, 0);
	reorganize_gauge(u1, u1_re, 1);

	printf("Done!\n");

	init_dfe();

	int address_g0  = 0;
	int address_g1  = T * burstsPerGaugeT;
	int address_p   = 2 * T * burstsPerGaugeT;
	int address_mp  = 2 * T * burstsPerGaugeT + T * burstsPerSpinorT;
	int address_mmp = 2 * T * burstsPerGaugeT + 2 * T * burstsPerSpinorT;
	int address_tmp = 2 * T * burstsPerGaugeT + 3 * T * burstsPerSpinorT;

	printf("Transferring Gauge to Memory  ...");
	transfer_gauges_to_dfe(u0_re, address_g0);
	transfer_gauges_to_dfe(u1_re, address_g1);
	printf("Done!\n");

	printf("Transferring Input Spinors to Memory  ...");
	transfer_spinors_to_dfe (in, address_p) ;
	printf("Done!\n");

	struct timeval start_time, end_time;

	gettimeofday(&start_time, NULL);

	for (int i = 0; i<1 ; i++) {
		printf("Calculating on DFE ...");
		apply_dirac (address_p,   address_tmp, address_g1, 0,          1, 0, 0);
		apply_dirac (address_tmp, address_mp,  address_g0, address_p,  0, 1, 0);
		apply_dirac (address_mp,  address_tmp, address_g1, 0,          1, 0, 1);
		apply_dirac (address_tmp, address_mmp, address_g0, address_mp, 0, 1, 1);
		printf("Done!\n");
	}

	gettimeofday(&end_time, NULL);

	double elapsedTime = (end_time.tv_sec - start_time.tv_sec) * 1000.0;      // sec to ms
	elapsedTime += (end_time.tv_usec - start_time.tv_usec) / 1000.0;   // us to ms

	printf("DFE time elapsed: %f ms\n", elapsedTime);
	printf("DFE Throughput: %g GFLOPS\n", (1320.0*(double)VOLUME/2*4*1.0)*1000.0/elapsedTime);

	printf("Transferring Input Spinors to Memory  ...");
	transfer_spinors_to_cpu (out_dfe, address_mmp);
	printf("Done!\n");

	unload_dfe();

	printf("Verifying LQCD output ...\n");

	return  verify_results(out_dfe, out_expected, VOLUME/2);
}
/************* Calls to DFE ********************/
void init_dfe(void) {
	maxfile = S_LQCD_init();
	engine = max_load(maxfile, "*");
	int burstSize = max_get_burst_size(maxfile, 0);
	burstsPerGaugeT  = (LX*LY*LZ/2*8*sizeof(su3)) / burstSize;
	burstsPerSpinorT = (LX*LY*LZ/2*sizeof(spinor)) / burstSize;
}

void transfer_spinors_to_dfe (spinor *in, int address) {
	max_actions_t *act = max_actions_init(maxfile, 0);
	max_queue_input(act, "spinor_in", in, VOLUME/2 * sizeof(spinor));
	max_set_ticks(act, "spWriteCmdKernel", T*burstsPerSpinorT/S_LQCD_spCmdSize);
	max_set_uint64t(act, "spWriteCmdKernel", "startAddress", address);
	max_set_uint64t(act, "spWriteCmdKernel", "halos", 0);
    max_route(act, "sptoLmemMux_fromCPU", "sptoLmemMux");
    max_lmem_set_interrupt_on(act, "sptoLmem");

	max_ignore_kernel(act, "diracKernel");
	max_ignore_kernel(act, "gReadCmdKernel");
	max_ignore_kernel(act, "gWriteCmdKernel");
	max_ignore_kernel(act, "spReadCmdKernel0");
	max_ignore_kernel(act, "spReadCmdKernel1");
    max_ignore_route(act, "spfromLmem0Demux");

	max_run(engine, act);
	max_actions_free(act);
}

void transfer_gauges_to_dfe (su3 *in, int address) {
	max_actions_t *act = max_actions_init(maxfile, 0);
	max_queue_input(act, "gauge_in", in, VOLUME/2 * 8 * sizeof(su3));
	max_set_ticks(act, "gWriteCmdKernel", T*burstsPerGaugeT/S_LQCD_gCmdSize);
	max_set_uint64t(act, "gWriteCmdKernel", "startAddress", address);
	max_lmem_set_interrupt_on(act, "gtoLmem");

	max_ignore_kernel(act, "diracKernel");
	max_ignore_kernel(act, "gReadCmdKernel");
	max_ignore_kernel(act, "spWriteCmdKernel");
	max_ignore_kernel(act, "spReadCmdKernel0");
	max_ignore_kernel(act, "spReadCmdKernel1");
	max_ignore_route(act, "sptoLmemMux");
	max_ignore_route(act, "spfromLmem0Demux");

	max_run(engine, act);
	max_actions_free(act);
}

void apply_dirac (int in_address, int out_address, int gauge_address, int p_address, int ieo, int doSub, int isign) {
	max_actions_t *act = max_actions_init(maxfile, 0);

	max_set_double(act, "diracKernel", "ieo", ieo);
	max_set_double(act, "diracKernel", "doSub", doSub);
	max_set_double(act, "diracKernel", "isign", isign);
	if (doSub == 0) {
		max_ignore_scalar(act, "diracKernel", "alpha");
		max_set_double(act, "diracKernel", "beta_s",   beta_s);
		max_set_double(act, "diracKernel", "beta_t_b", beta_t_b);
		max_set_double(act, "diracKernel", "beta_t_f", beta_t_f);
	} else {
		max_set_double(act, "diracKernel", "alpha", 4+mass);
		max_set_double(act, "diracKernel", "beta_s",   beta_s / (16+4*mass));
		max_set_double(act, "diracKernel", "beta_t_b", beta_t_b / (16+4*mass));
		max_set_double(act, "diracKernel", "beta_t_f", beta_t_f / (16+4*mass));
	}
	max_set_ticks(act, "diracKernel", 16/NUM_PIPES * ( (T+2)*LX*LY*LZ/2 +  LOOP_OFFSET) + 2 );
	max_set_ticks(act, "gReadCmdKernel", (T+2)*burstsPerGaugeT/S_LQCD_gCmdSize);
	max_set_uint64t(act, "gReadCmdKernel", "startAddress", gauge_address);
	max_set_ticks(act, "spReadCmdKernel0", (T+2)*burstsPerSpinorT/S_LQCD_spCmdSize);
	max_set_uint64t(act, "spReadCmdKernel0", "startAddress", in_address);
	max_set_uint64t(act, "spReadCmdKernel0", "halos", 1);
	if (doSub == 0) {
		max_ignore_kernel(act, "spReadCmdKernel1");
	} else  {
		max_set_ticks(act, "spReadCmdKernel1", (T)*burstsPerSpinorT/S_LQCD_spCmdSize);
		max_set_uint64t(act, "spReadCmdKernel1", "startAddress", p_address);
		max_set_uint64t(act, "spReadCmdKernel1", "halos", 0);
	}
	max_set_ticks(act, "spWriteCmdKernel", T*burstsPerSpinorT/S_LQCD_spCmdSize);
	max_set_uint64t(act, "spWriteCmdKernel", "startAddress", out_address);
	max_set_uint64t(act, "spWriteCmdKernel", "halos", 0);
	max_route(act, "sptoLmemMux_fromKernel", "sptoLmemMux");
	max_route(act, "spfromLmem0Demux", "spfromLmem0Demux_toKernel");
	max_lmem_set_interrupt_on(act, "sptoLmem");

	max_ignore_kernel(act, "gWriteCmdKernel");

	max_run(engine, act);
	max_actions_free(act);
}

void transfer_spinors_to_cpu (spinor *out, int address) {
	max_actions_t *act = max_actions_init(maxfile, 0);
	max_queue_output(act, "spinor_out", out,  VOLUME/2 * sizeof(spinor));
	max_set_ticks(act, "spReadCmdKernel0", T*burstsPerSpinorT/S_LQCD_gCmdSize);
	max_set_uint64t(act, "spReadCmdKernel0", "startAddress", address);
	max_set_uint64t(act, "spReadCmdKernel0", "halos", 0);
	max_route(act, "spfromLmem0Demux", "spfromLmem0Demux_toCPU");

	max_ignore_kernel(act, "diracKernel");
	max_ignore_kernel(act, "gReadCmdKernel");
	max_ignore_kernel(act, "gWriteCmdKernel");
	max_ignore_kernel(act, "spWriteCmdKernel");
	max_ignore_kernel(act, "spReadCmdKernel1");
	max_ignore_route(act, "sptoLmemMux");

	max_run(engine, act);
	max_actions_free(act);

}

void unload_dfe(void) {
	max_unload(engine);
}

/************* Verifying Functions ***************/

int verify_results (spinor* dfe_out, spinor *expected_out, int V) {
	for (int i=0 ; i<V ; i++ ) {
		int error = compare_spinor(expected_out+i , dfe_out+i);
		if (error) {
			printf("Wrong %d! %d\n", error,i);
			spinor *a = expected_out + i;
			spinor *b = dfe_out + i;
			printf("%f %f    %f %f ", creal(a->s0.c0), cimag(a->s0.c0),
					                  creal(b->s0.c0), cimag(b->s0.c0));
			return 1;
		}
	}
	printf("Good.\n");
	return 0;

}

int AreNotSameComplex(complex float a, float complex b)
{
    return ( fabs(creal(a) - creal(b)) > 0.001 ) ||
    	   ( cimag(creal(a) - cimag(b)) > 0.001 );
}

int compare_spinor (spinor *a, spinor *b) {
	if (AreNotSameComplex(a->s0.c0 ,b->s0.c0)) return 1;
	if (AreNotSameComplex(a->s0.c1 ,b->s0.c1)) return 2;
	if (AreNotSameComplex(a->s0.c2 ,b->s0.c2)) return 3;
	if (AreNotSameComplex(a->s1.c0 ,b->s1.c0)) return 4;
	if (AreNotSameComplex(a->s1.c1 ,b->s1.c1)) return 5;
	if (AreNotSameComplex(a->s1.c2 ,b->s1.c2)) return 6;
	if (AreNotSameComplex(a->s2.c0 ,b->s2.c0)) return 7;
	if (AreNotSameComplex(a->s2.c1 ,b->s2.c1)) return 8;
	if (AreNotSameComplex(a->s2.c2 ,b->s2.c2)) return 9;
	if (AreNotSameComplex(a->s3.c0 ,b->s3.c0)) return 10;
	if (AreNotSameComplex(a->s3.c1 ,b->s3.c1)) return 11;
	if (AreNotSameComplex(a->s3.c2 ,b->s3.c2)) return 12;
	return 0;
}

/************ Utility Functions *********************/

void create_random_input(spinor* s, su3* u) {
	for (int i = 0; i < VOLUME / 2; i++) {
		create_random_spinor(s + i);
	}
	for (int i = 0; i < VOLUME * 4; i++) {
		create_random_su3(u + i);
	}
}

void create_random_spinor(spinor * s) {
	create_random_su3vector(&s->s0);
	create_random_su3vector(&s->s1);
	create_random_su3vector(&s->s2);
	create_random_su3vector(&s->s3);
}

void create_random_su3vector(su3_vector *v) {
	create_random_complex(&v->c0);
	create_random_complex(&v->c1);
	create_random_complex(&v->c2);
}

void create_random_su3(su3 *s) {
	create_random_complex(&s->c00);
	create_random_complex(&s->c01);
	create_random_complex(&s->c02);
	create_random_complex(&s->c10);
	create_random_complex(&s->c11);
	create_random_complex(&s->c12);
	create_random_complex(&s->c20);
	create_random_complex(&s->c21);
	create_random_complex(&s->c22);
}

void create_random_complex(float complex *a) {
	float r[2];
	for (int i = 0; i < 2; i++) {
		r[i] = (float) (rand()) / RAND_MAX * 2;
	}
	*a = r[0] + I * r[1];
}

void print_spinors (spinor* s) {
	for (int i=0 ; i < VOLUME/2 ; i++ ) {
		printf("%f %f\n", creal(s->s0.c0), cimag(s->s0.c0));
		s++;
	}
}

void print_gauges (su3* s) {
	for (int i=0 ; i < VOLUME/2 * 4 ; i++ ) {
		printf("%f %f\n", creal(s->c00), cimag(s->c00));
		s++;
	}
}

/****************** Functions for reading data from files ***************/
void read_spinor(char * filename, spinor *out) {
	FILE * fptr = fopen(filename, "r" );
	char temp[64];

	for (int t=0 ; t<T ; t++ ) {
		for (int z=0 ; z<LZ ; z++ ) {
			for (int y=0 ; y<LY ; y++ ) {
				for (int x=0 ; x<LX/2 ; x++ ) {

					spinor *s = &out [ (((((t*LZ)+z)*LY)+y)*LX/2)+x ];

					fgets (temp, 64, fptr);
					float a, b;

					fscanf(fptr, "%f %f", &a, &b);
					s->s0.c0 = a + I * b;
					fscanf(fptr, "%f\n%f\n", &a, &b);
					s->s0.c1 = a + I * b;
					fscanf(fptr, "%f\n%f\n", &a, &b);
					s->s0.c2 = a + I * b;

					fscanf(fptr, "%f\n%f\n", &a, &b);
					s->s1.c0 = a + I * b;
					fscanf(fptr, "%f\n%f\n", &a, &b);
					s->s1.c1 = a + I * b;
					fscanf(fptr, "%f\n%f\n", &a, &b);
					s->s1.c2 = a + I * b;

					fscanf(fptr, "%f\n%f\n", &a, &b);
					s->s2.c0 = a + I * b;
					fscanf(fptr, "%f\n%f\n", &a, &b);
					s->s2.c1 = a + I * b;
					fscanf(fptr, "%f\n%f\n", &a, &b);
					s->s2.c2 = a + I * b;

					fscanf(fptr, "%f\n%f\n", &a, &b);
					s->s3.c0 = a + I * b;
					fscanf(fptr, "%f\n%f\n", &a, &b);
					s->s3.c1 = a + I * b;
					fscanf(fptr, "%f\n%f\n", &a, &b);
					s->s3.c2 = a + I * b;

				}
			}
		}
	}
	fclose(fptr);
}

void read_gauge(char * filename, su3 *s) {
	FILE * fptr = fopen(filename, "r" );
	char temp[64];
	for (int i = 0 ; i < VOLUME/2 * 8 ; i++ ) {
		fgets (temp, 64, fptr);
		float a, b;

		fscanf(fptr, "%f\n%f\n", &a, &b);
		s->c00 = a + I * b;
		fscanf(fptr, "%f\n%f\n", &a, &b);
		s->c01 = a + I * b;
		fscanf(fptr, "%f\n%f\n", &a, &b);
		s->c02 = a + I * b;

		fscanf(fptr, "%f\n%f\n", &a, &b);
		s->c10 = a + I * b;
		fscanf(fptr, "%f\n%f\n", &a, &b);
		s->c11 = a + I * b;
		fscanf(fptr, "%f\n%f\n", &a, &b);
		s->c12 = a + I * b;

		fscanf(fptr, "%f\n%f\n", &a, &b);
		s->c20 = a + I * b;
		fscanf(fptr, "%f\n%f\n", &a, &b);
		s->c21 = a + I * b;
		fscanf(fptr, "%f\n%f\n", &a, &b);
		s->c22 = a + I * b;

		s++;
	}
	fclose(fptr);
}

/****************** Functions for reordering data in memory ***************/
void reorganize_gauge(su3 const * const in, su3 * const out, int ieo) {
	int i = 0;
	for (int t=0 ; t<T ; t++ ) {
		for (int z=0 ; z<LZ ; z++ ) {
			for (int y=0 ; y<LY ; y++ ) {
				for (int x=0 ; x<LX/2 ; x++ ) {
					for (int mu=0; mu<4 ; mu++ ) {
						for (int f=-1; f<=1 ; f+= 2) {
							su3 tmp = in[i];

							int isOddRow = (t & 1) ^ (z & 1) ^ (y & 1) ^ ieo;

							/*int mu_ = (mu+1)%4;
							int t_ = t;               // converting from checkerboarded
							int z_ = z/2;             // coordinates of qphix along x-axis
							int y_ = y;               // to tmLQCD checkerboarding along
							int x_ = (2*x)+isOddRow;  // along y-axis*/

							if (f == 1) {

								int xx = (mu==0)?( isOddRow ? x+1 :x ):x;
								int yy = (mu==1)?y+1:y;
								int zz = (mu==2)?z+1:z;
								int tt = (mu==3)?t+1:t;

								tt = (tt+T) % T;
								zz = (zz+LZ) % LZ;
								yy = (yy+LY) % LY;
								xx = (xx+LX) % (LX/2);

								out[ ((((((tt*LZ)+zz)*LY)+yy)*LX/2)+xx)*8+mu*2+1 ] = tmp;

							} else {

								int xx = (mu==0)?( isOddRow ? x : x-1 ):x;
								int yy = (mu==1)?y-1:y;
								int zz = (mu==2)?z-1:z;
								int tt = (mu==3)?t-1:t;

								tt = (tt+T) % T;
								zz = (zz+LZ) % LZ;
								yy = (yy+LY) % LY;
								xx = (xx+LX) % (LX/2);

								out[ ((((((tt*LZ)+zz)*LY)+yy)*LX/2)+xx)*8+mu*2+0 ] = tmp;

							}

							i++;
						}
					}
				}
			}
		}
	}
}





/************* These functions are not needed anymore ************************/

void reorganize_ueven (su3 *out, su3 *in) {
	int i = 0;
	for (int t=0 ; t<T ; t++ ) {
		for (int x=0 ; x<LX ; x++ ) {
			for (int y=0 ; y<LY ; y++ ) {
				for (int z=0 ; z<LZ/2 ; z++ ) {
					int isOddRow = (t & 1) ^ (x & 1) ^ (y & 1);
					for (int mu=0; mu<4 ; mu++ ) {
						int tt = (mu==0)?t+1:t;
						int xx = (mu==1)?x+1:x;
						int yy = (mu==2)?y+1:y;
						int zz = (mu==3)?( (isOddRow)?(z):z+1 ):z;

						tt = (tt+T) % T;
						xx = (xx+LX) % LX;
						yy = (yy+LY) % LY;
						zz = (zz+LZ) % (LZ/2);

						out[i] = in [ ((((((tt*LX)+xx)*LY)+yy)*LZ/2)+zz)*4+mu ];
						i++;
					}
				}
			}
		}
	}
}

void reorganize_back_ueven (su3 *out, su3 *in) {
	int i = 0;
	for (int t=0 ; t<T ; t++ ) {
		for (int x=0 ; x<LX ; x++ ) {
			for (int y=0 ; y<LY ; y++ ) {
				for (int z=0 ; z<LZ/2 ; z++ ) {
					int isOddRow = (t & 1) ^ (x & 1) ^ (y & 1);
					for (int mu=0; mu<4 ; mu++ ) {
						int tt = (mu==0)?t+1:t;
						int xx = (mu==1)?x+1:x;
						int yy = (mu==2)?y+1:y;
						int zz = (mu==3)?( (isOddRow)?(z):z+1 ):z;

						tt = (tt+T) % T;
						xx = (xx+LX) % LX;
						yy = (yy+LY) % LY;
						zz = (zz+LZ) % (LZ/2);

						out[((((((tt*LX)+xx)*LY)+yy)*LZ/2)+zz)*4+mu ] = in [ i ];
						i++;
					}
				}
			}
		}
	}
}

/* Converting from qphix style gauge for the whole lattice, to even/odd separated
 * tmLQCD style gauge fields
 */
void devide_gauge_to_oddeven(su3 const * const in, su3 * const even, su3 * const odd, int ieo) {
	int i = 0;
	for (int t=0 ; t<T ; t++ ) {
		for (int z=0 ; z<LZ ; z++ ) {
			for (int y=0 ; y<LY ; y++ ) {
				for (int x=0 ; x<LX/2 ; x++ ) {
					for (int mu=0; mu<4 ; mu++ ) {
						for (int f=-1; f<=1 ; f+= 2) {
							su3 tmp = in[i];

							int isOddRow = (t & 1) ^ (z & 1) ^ (y & 1) ^ ieo;

							/*int mu_ = (mu+1)%4;
							int t_ = t;               // converting from checkerboarded
							int z_ = z/2;             // coordinates of qphix along x-axis
							int y_ = y;               // to tmLQCD checkerboarding along
							int x_ = (2*x)+isOddRow;  // along y-axis*/

							if (f == 1) {

								int xx = (mu==0)?( isOddRow ? x+1 :x ):x;
								int yy = (mu==1)?y+1:y;
								int zz = (mu==2)?z+1:z;
								int tt = (mu==3)?t+1:t;

								tt = (tt+T) % T;
								zz = (zz+LZ) % LZ;
								yy = (yy+LY) % LY;
								xx = (xx+LX) % (LX/2);

								even[ ((((((tt*LZ)+zz)*LY)+yy)*LX/2)+xx)*4+mu ] = tmp;

							} else {

								int xx = (mu==0)?( isOddRow ? x : x-1 ):x;
								int yy = (mu==1)?y-1:y;
								int zz = (mu==2)?z-1:z;
								int tt = (mu==3)?t-1:t;

								tt = (tt+T) % T;
								zz = (zz+LZ) % LZ;
								yy = (yy+LY) % LY;
								xx = (xx+LX) % (LX/2);

								odd[ ((((((tt*LZ)+zz)*LY)+yy)*LX/2)+xx)*4+mu ] = tmp;

							}

							i++;
						}
					}
				}
			}
		}
	}
}


void add_1d_halos_spinor(spinor* with_halos, spinor* orig, int halos) {

	for (int t = -halos ; t < T+halos ; t++ ) {
		for (int z = 0 ; z < LZ ; z++ ) {
			for (int y = 0 ; y < LY ; y++ ) {
				for (int x = 0 ; x < LX/2 ; x++ ) {
					int tt = (t+T)%T;

					with_halos[ ((( (t+halos) * LZ + z) * LY + y) * LX/2 ) + x] =
							orig[ (((tt*LZ + z) * LY + y) * LX/2 ) + x];
				}
			}
		}
	}

}

void add_1d_halos_gauge(su3* with_halos, su3* orig, int halos) {

	for (int t = -halos ; t < T+halos ; t++ ) {
		for (int z = 0 ; z < LZ ; z++ ) {
			for (int y = 0 ; y < LY ; y++ ) {
				for (int x = 0 ; x < LX/2 ; x++ ) {
					int tt = (t+T)%T;

					for (int i=0 ; i < 8 ; i++ ){
						with_halos[ (((( (t+halos) * LZ + z) * LY + y) * LX/2 ) + x) * 8 + i] =
								orig[ ((((tt*LZ + z) * LY + y) * LX/2 ) + x) * 8 + i];
					}
				}
			}
		}
	}

}


void add_4d_halos_spinor(spinor* with_halos, spinor* orig, int halos, int eo) {
	int lhx, lhy, lhz;
	lhx = LX+2*halos;
	lhy = LY+2*halos;
	lhz = LZ/2+halos;

	eo = eo ^ 1;

	for (int t = -halos ; t < T+halos ; t++ ) {
		for (int x = -halos ; x < LX+halos ; x++ ) {
			for (int y = -halos ; y < LY+halos ; y++ ) {
				for (int z = -(halos/2) ; z < LZ/2+(halos+1)/2 ; z++ ) {
					int tt = (t+T)%T;
					int xx = (x+LX)%LX;
					int yy = (y+LY)%LY;
					int isOddRow = (tt & 1) ^ (xx & 1) ^ (yy & 1) ^ eo;
					int zz;
					if (halos%2 == 0) {
						zz = (z+LZ/2)%(LZ/2);
					} else {
						zz = (z+LZ*2-isOddRow*halos)%(LZ/2);
					}

					with_halos[ (((( (t+halos) * lhx +
					                 (x+halos)) * lhy +
					                 (y+halos)) * lhz ) +
					                 (z+halos/2))] =
							orig[ (((tt*LX + xx) * LY + yy) * LZ/2 ) + zz];
				}
			}
		}
	}
}

void add_4d_halos_gauge(su3* with_halos, su3* orig, int halos, int eo) {
	int lhx, lhy, lhz;
	lhx = LX+2*halos;
	lhy = LY+2*halos;
	lhz = LZ/2+halos;

	eo = eo ^ 1;

	for (int t = -halos ; t < T+halos ; t++ ) {
		for (int x = -halos ; x < LX+halos ; x++ ) {
			for (int y = -halos ; y < LY+halos ; y++ ) {
				for (int z = -(halos/2) ; z < LZ/2+(halos+1)/2 ; z++ ) {
					int tt = (t+T)%T;
					int xx = (x+LX)%LX;
					int yy = (y+LY)%LY;
					int isOddRow = (tt & 1) ^ (xx & 1) ^ (yy & 1) ^ eo;
					int zz;
					if (halos%2 == 0) {
						zz = (z+LZ/2)%(LZ/2);
					} else {
						zz = (z+LZ*2-isOddRow*halos)%(LZ/2);
					}

					for (int i=0 ; i < 4 ; i++ ){
						with_halos[ ((((( (t+halos) * lhx +
				                          (x+halos)) * lhy +
				                          (y+halos)) * lhz ) +
				                          (z+halos/2)) * 4) + i] =
						      orig[ ((((tt*LX + xx) * LY + yy) * LZ/2 ) + zz) * 4 + i];
					}
				}
			}
		}
	}
}

