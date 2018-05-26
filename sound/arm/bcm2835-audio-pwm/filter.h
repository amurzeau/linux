#ifndef FILTER_H
#define FILTER_H

#include <asm/string.h>

typedef int16_t sample_t;
typedef int32_t sample_hires_t;

#define FILTER_TAP_NUM 128

static int16_t filter_taps[FILTER_TAP_NUM] = {
  85,
  -31,
  -30,
  -31,
  -34,
  -37,
  -40,
  -42,
  -43,
  -43,
  -41,
  -37,
  -30,
  -21,
  -9,
  5,
  21,
  38,
  56,
  74,
  91,
  106,
  119,
  127,
  132,
  131,
  124,
  110,
  90,
  64,
  31,
  -6,
  -49,
  -94,
  -142,
  -189,
  -234,
  -274,
  -308,
  -334,
  -348,
  -350,
  -337,
  -308,
  -261,
  -197,
  -116,
  -16,
  100,
  230,
  374,
  529,
  691,
  858,
  1025,
  1190,
  1349,
  1497,
  1631,
  1749,
  1847,
  1922,
  1974,
  2000,
  2000,
  1974,
  1922,
  1847,
  1749,
  1631,
  1497,
  1349,
  1190,
  1025,
  858,
  691,
  529,
  374,
  230,
  100,
  -16,
  -116,
  -197,
  -261,
  -308,
  -337,
  -350,
  -348,
  -334,
  -308,
  -274,
  -234,
  -189,
  -142,
  -94,
  -49,
  -6,
  31,
  64,
  90,
  110,
  124,
  131,
  132,
  127,
  119,
  106,
  91,
  74,
  56,
  38,
  21,
  5,
  -9,
  -21,
  -30,
  -37,
  -41,
  -43,
  -43,
  -42,
  -40,
  -37,
  -34,
  -31,
  -30,
  -31,
  85
};

/*
struct filter_iir_biquad_t {
	double b_coefs[3];
	double a_coefs[3];
	double s1;
	double s2;
};

static void filter_iir_init(struct filter_iir_biquad_t* filter, const double a_coefs[3], const double b_coefs[3]) {
	filter->s1 = filter->s2 = 0;
	memcpy(filter->b_coefs, b_coefs, sizeof(filter->b_coefs));
	memcpy(filter->a_coefs, a_coefs, sizeof(filter->a_coefs));
}

static double filter_iir_put(struct filter_iir_biquad_t* f, double input) {
	double* const a = f->a_coefs;
	double* const b = f->b_coefs;
	
	double y = b[0] * input + f->s1;
	f->s1 = f->s2 + b[1] * input - a[1] * y;
	f->s2 = b[2] * input - a[2] * y;
	
	return y;
}*/

struct filter_sample_t {
	sample_t left;
	sample_t right;
};

struct filter_sample_hires_t {
	sample_hires_t left;
	sample_hires_t right;
};

struct filter_fir_t {
	sample_t coefs[128];
	struct filter_sample_hires_t history[128];
	int last_index;
	int upsampling_ratio;
	int downsampling_ratio;
	int modpos;
};

static inline void filter_fir_init(struct filter_fir_t* f, const sample_t coefs[128], int coef_size, int upsampling_ratio, int downsampling_ratio) {
	memcpy(f->coefs, coefs, sizeof(*f->coefs) * coef_size);
	memset(f->history, 0, sizeof(f->history));
	f->last_index = 0;
	f->upsampling_ratio = upsampling_ratio;
	f->downsampling_ratio = downsampling_ratio;
	f->modpos = 0;
}

static inline void filter_fir_put(struct filter_fir_t* f, sample_hires_t input_left, sample_hires_t input_right) {
	f->history[++f->last_index & (128-1)] = (struct filter_sample_hires_t){input_left, input_right};
}

static inline void filter_fir_get(struct filter_fir_t* f, struct filter_sample_hires_t* outputs, int* output_size) {
	#define TAP 128
	int upsampling_ratio = 25;
	int downsampling_ratio = 3;
	int n = f->last_index;
	int startout = f->modpos;
	int k;
	sample_t* coefs = f->coefs;
	int tap_lefts = TAP;
	
	memset(outputs, 0, sizeof(*outputs) * (*output_size));
	*output_size = TAP > upsampling_ratio ? (upsampling_ratio - startout - 1)/downsampling_ratio + 1 : (TAP - startout - 1)/downsampling_ratio + 1;

	for(k = 0; k < TAP; k += upsampling_ratio, n--, tap_lefts -= upsampling_ratio) {

		struct filter_sample_hires_t hval = f->history[n & (128-1)];
		int l, o;
		int end = upsampling_ratio;
		if(unlikely(tap_lefts < upsampling_ratio))
			end = tap_lefts;
		for(o = 0, l = startout; l < end; l += downsampling_ratio, o++) {
			sample_t coef = coefs[k + l];
			outputs[o].left += hval.left * coef;
			outputs[o].right += hval.right * coef;
		}
	};
	
	startout += ((downsampling_ratio - (upsampling_ratio % downsampling_ratio)) % downsampling_ratio);
	if(startout >= downsampling_ratio)
		startout -= downsampling_ratio;
	f->modpos = startout;
}

#endif
