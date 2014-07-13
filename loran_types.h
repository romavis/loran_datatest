#ifndef LORAN_TYPES_H
#define LORAN_TYPES_H

#include <stdint.h>
#include <stddef.h>

/*******************************************************************************
 * Data types:
 *
 * lc_type_sample - ADC sample type, usually unsigned integer :-)
 * lc_type_data_s - data type on output from DC blocker (usually signed counter-
 *	part of lc_type_sample)
 * lc_type_comb_s - data type inside comb filter buffers, usually - fixed point
 *	with LC_BITS_COMBF fractional part length (signed, of course)
 *	We need fixed point here because just integral type introduces high
 *	quantization error when using high P values
 */

typedef uint16_t	lc_type_sample;
typedef int16_t		lc_type_data_s;
typedef int32_t		lc_type_comb_s;

#endif // LORAN_TYPES_H
