////////////////////////////////////////////////////////////////////////////////
/**
*@file      filter.h
*@brief     Various filter designs
*@author    Ziga Miklosic
*@date      02.01.2021
*@version   V0.0.1
*/
////////////////////////////////////////////////////////////////////////////////
/**
*@addtogroup FILTER_API
* @{ <!-- BEGIN GROUP -->
*
*/
////////////////////////////////////////////////////////////////////////////////

#ifndef __FILTER_H
#define __FILTER_H

////////////////////////////////////////////////////////////////////////////////
// Includes
////////////////////////////////////////////////////////////////////////////////
#include "project_config.h"
#include "stdint.h"

////////////////////////////////////////////////////////////////////////////////
// Definitions
////////////////////////////////////////////////////////////////////////////////

/**
 * 	Filter status
 */
typedef enum
{
	eFILTER_OK 			= 0x00,		/**<Normal operation */
	eFILTER_ERROR		= 0x01,		/**<General error */
} filter_status_t;

/**
 * 	RC filter instance type
 */
typedef struct filter_rc_s * p_filter_rc_t;

/**
 * 	CR filter instance type
 */
typedef struct filter_cr_s * p_filter_cr_t;

/**
 * 	FIR filter instance type
 */
typedef struct filter_fir_s * p_filter_fir_t;

/**
 * 	IIR filter instance type
 */
typedef struct filter_iir_s * p_filter_iir_t;

////////////////////////////////////////////////////////////////////////////////
// Functions
////////////////////////////////////////////////////////////////////////////////
filter_status_t filter_rc_init		(p_filter_rc_t * p_filter_inst, const float32_t fc, const float32_t dt, const uint8_t order, const float32_t init_value);
float32_t 		filter_rc_update	(p_filter_rc_t filter_inst, const float32_t x);
filter_status_t filter_cr_init		(p_filter_cr_t * p_filter_inst, const float32_t fc, const float32_t dt, const uint8_t order);
float32_t 		filter_cr_update	(p_filter_cr_t filter_inst, const float32_t x);
filter_status_t filter_fir_init		(p_filter_fir_t * p_filter_inst, const float32_t * p_a, const uint32_t order);
float32_t		filter_fir_update	(p_filter_fir_t filter_inst, const float32_t x);
filter_status_t filter_iir_init		(p_filter_iir_t * p_filter_inst, const float32_t * p_a, const float32_t * p_b, const uint32_t a_size, const uint32_t b_size);
float32_t		filter_iir_update	(p_filter_iir_t filter_inst, const float32_t x);


#endif // __FILTER_H


////////////////////////////////////////////////////////////////////////////////
/**
* @} <!-- END GROUP -->
*/
////////////////////////////////////////////////////////////////////////////////
