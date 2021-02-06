////////////////////////////////////////////////////////////////////////////////
/**
*@file      filter.c
*@brief     Various filter designs
*@author    Ziga Miklosic
*@date      02.01.2021
*@version   V0.0.1
*
*@section   Description
*   
*   This module contains different kind of digital filter
*   implementation. All of the following filter types are
*   being simulated in python for validation purposes.
*
*@section 	Dependencies
*
* 	Some implementation of filters uses ring buffers.
*
*/
////////////////////////////////////////////////////////////////////////////////
/*!
* @addtogroup FILTER
* @{ <!-- BEGIN GROUP -->
*/
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Includes
////////////////////////////////////////////////////////////////////////////////
#include "filter.h"
#include "stdint.h"
#include "stdlib.h"
#include "math.h"

#include "middleware/ring_buffer/ring_buffer.h"

////////////////////////////////////////////////////////////////////////////////
// Definitions
////////////////////////////////////////////////////////////////////////////////

/**
 * 	RC Filter data
 */
typedef struct filter_rc_s
{
	float32_t *	p_y;		/**<Output of filter + previous values */
	float32_t	alpha;		/**<Filter smoothing factor */
	uint8_t		order;		/**<Filter order - number of cascaded filter */
} filter_rc_t;

/**
 * 	CR Filter data
 */
typedef struct filter_cr_s
{
	float32_t * p_y;		/**<Output of filter + previous values */
	float32_t * p_x;		/**<Input of filter + previous values */
	float32_t  	alpha;		/**<Filter smoothing factor */
	uint8_t 	order;		/**<Filter order - number of cascaded filter */
} filter_cr_t;

/**
 * 	FIR Filter data
 */
typedef struct filter_fir_s
{
	p_ring_buffer_t  	p_x;		/**<Previous values of input filter */
	float32_t 		*	p_a;		/**<Filter coefficients */
	uint32_t			order;		/**<Number of FIR filter taps - order of filter */
} filter_fir_t;

/**
 * 	IIR Filter data
 */
typedef struct filter_iir_s
{
	p_ring_buffer_t   p_y;		/**<Previous values of filter outputs */
	p_ring_buffer_t   p_x;		/**<Previous values of filter inputs*/
	float32_t 		* p_a;		/**<Filter poles */
	float32_t		* p_b;		/**<Filter zeros */
	uint32_t		  a_size;	/**<Number of filter poles */
	uint32_t	 	  b_size;	/**<Number of filter zeros */
} filter_iir_t;


////////////////////////////////////////////////////////////////////////////////
// Variables
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Function prototypes
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Functions
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**
*   Initialize RC filter
*
*@Note: Order of RC filter is represented as number of cascaded RC
*		analog equivalent circuits!
*
* @param[in] 	p_filter_inst	- Pointer to RC filter instance
* @param[in] 	fc				- Filter cutoff frequency
* @param[in] 	dt				- Sample time
* @param[in] 	order			- Order of filter (number of cascaded filter)
* @param[in] 	init_value		- Initial value
* @return 		status			- Status of operation
*/
////////////////////////////////////////////////////////////////////////////////
filter_status_t filter_rc_init(p_filter_rc_t * p_filter_inst, const float32_t fc, const float32_t dt, const uint8_t order, const float32_t init_value)
{
	filter_status_t status = eFILTER_OK;
	uint8_t i;

	if (( NULL != p_filter_inst ) && ( order > 0UL ))
	{
		// Allocate space
		*p_filter_inst 			= malloc( sizeof( filter_rc_t ));
		(*p_filter_inst)->p_y 	= malloc( order  * sizeof( float32_t ));

		// Check if allocation succeed
		if 	(	( NULL != *p_filter_inst )
			&& 	( NULL != (*p_filter_inst)->p_y ))
		{
			// Calculate coefficient
			(*p_filter_inst)->alpha = (float32_t) ( dt / ( dt + ( 1.0f / ( M_TWOPI * fc ))));

			// Store order
			(*p_filter_inst)->order = order;

			// Initial value
			for ( i = 0; i < order; i++)
			{
				(*p_filter_inst)->p_y[i] = init_value;
			}
		}
		else
		{
			status = eFILTER_ERROR;
		}
	}
	else
	{
		status = eFILTER_ERROR;
	}

	return status;
}

////////////////////////////////////////////////////////////////////////////////
/**
*   Update RC filter
*
* @param[in] 	filter_inst	- Pointer to RC filter instance
* @param[in] 	x			- Input value
* @return 		y			- Output value
*/
////////////////////////////////////////////////////////////////////////////////
float32_t filter_rc_update(p_filter_rc_t filter_inst, const float32_t x)
{
	float32_t y = 0.0f;
	uint8_t n;

	if ( NULL != filter_inst )
	{
		for ( n = 0; n < filter_inst->order; n++)
		{
			if ( 0 == n )
			{
				filter_inst->p_y[0] = ( filter_inst->p_y[0] + ( filter_inst->alpha * ( x - filter_inst->p_y[0] )));
			}
			else
			{
				filter_inst->p_y[n] = ( filter_inst->p_y[n] + ( filter_inst->alpha * ( filter_inst->p_y[n-1] - filter_inst->p_y[n] )));
			}
		}

		y = filter_inst->p_y[ filter_inst->order - 1U ];
	}

	return y;
}

////////////////////////////////////////////////////////////////////////////////
/**
*   Initialize CR filter
*
*@Note: Order of CR filter is represented as number of cascaded CR
*		analog equivalent circuits!
*
* @param[in] 	p_filter_inst	- Pointer to CR filter instance
* @param[in] 	fc				- Filter cutoff frequency
* @param[in] 	dt				- Sample time
* @param[in] 	order			- Order of filter (number of cascaded filter)
* @return 		status			- Status of operation
*/
////////////////////////////////////////////////////////////////////////////////
filter_status_t filter_cr_init(p_filter_cr_t * p_filter_inst, const float32_t fc, const float32_t dt, const uint8_t order)
{
	filter_status_t status = eFILTER_OK;
	uint8_t i;

	if (( NULL != p_filter_inst ) && ( order > 0UL ))
	{
		// Allocate space
		*p_filter_inst 			= malloc( sizeof( filter_cr_t ));
		(*p_filter_inst)->p_y 	= malloc( order  * sizeof( float32_t ));
		(*p_filter_inst)->p_x 	= malloc( order  * sizeof( float32_t ));

		// Check if allocation succeed
		if 	(	( NULL != *p_filter_inst )
			&& 	( NULL != (*p_filter_inst)->p_y )
			&&	( NULL != (*p_filter_inst)->p_x ))
		{
			// Calculate coefficient
			(*p_filter_inst)->alpha = (float32_t) (( 1.0f / ( M_TWOPI * fc )) / ( dt + ( 1.0f / ( M_TWOPI * fc ))));

			// Store order
			(*p_filter_inst)->order = order;

			// Initial value
			for ( i = 0; i < order; i++)
			{
				(*p_filter_inst)->p_y[i] = 0.0f;
				(*p_filter_inst)->p_x[i] = 0.0f;
			}
		}
		else
		{
			status = eFILTER_ERROR;
		}
	}
	else
	{
		status = eFILTER_ERROR;
	}

	return status;
}

////////////////////////////////////////////////////////////////////////////////
/**
*   Update CR filter
*
* @param[in] 	filter_inst	- Pointer to CR filter instance
* @param[in] 	x			- Input value
* @return 		y			- Output value
*/
////////////////////////////////////////////////////////////////////////////////
float32_t filter_cr_update(p_filter_cr_t filter_inst, const float32_t x)
{
	float32_t y = 0.0f;
	uint8_t n;

	if ( NULL != filter_inst )
	{
		for ( n = 0; n < filter_inst->order; n++)
		{
			if ( 0 == n )
			{
				filter_inst->p_y[0] = (( filter_inst->alpha * filter_inst->p_y[0] ) + ( filter_inst->alpha * ( x - filter_inst->p_x[0] )));
				filter_inst->p_x[0] = x;
			}
			else
			{
				filter_inst->p_y[n] = (( filter_inst->alpha * filter_inst->p_y[n] ) + ( filter_inst->alpha * ( filter_inst->p_y[n-1] - filter_inst->p_x[n] )));
				filter_inst->p_x[n] = filter_inst->p_y[n-1];
			}
		}

		y = filter_inst->p_y[ filter_inst->order - 1U ];
	}

	return y;
}

////////////////////////////////////////////////////////////////////////////////
/**
*   Initialize FIR filter
*
*@Note: Order of FIR filter is tap number of filter.
*
* @param[in] 	p_filter_inst	- Pointer to FIR filter instance
* @param[in] 	p_a				- Pointer to FIR coefficients
* @param[in] 	order			- Number of taps
* @return 		status			- Status of operation
*/
////////////////////////////////////////////////////////////////////////////////
filter_status_t filter_fir_init(p_filter_fir_t * p_filter_inst, const float32_t * p_a, const uint32_t order)
{
	filter_status_t status = eFILTER_OK;

	if (( NULL != p_filter_inst ) && ( order > 0UL ))
	{
		// Allocate filter space
		*p_filter_inst = malloc( sizeof( filter_fir_t ));

		// Allocation succeed
		if ( NULL != *p_filter_inst )
		{
			// Create ring buffer
			if ( eRING_BUFFER_OK == ring_buffer_init( &(*p_filter_inst)->p_x, order  ))
			{
				// Get filter coefficient & order
				(*p_filter_inst)->p_a = (float32_t*) p_a;
				(*p_filter_inst)->order = order;
			}
			else
			{
				status = eFILTER_ERROR;
			}
		}
		else
		{
			status = eFILTER_ERROR;
		}
	}
	else
	{
		status = eFILTER_ERROR;
	}

	return status;
}


////////////////////////////////////////////////////////////////////////////////
/**
*   Update FIR filter
*
* @param[in] 	filter_inst	- Pointer to FIR filter instance
* @param[in] 	x			- Input value
* @return 		y			- Output value
*/
////////////////////////////////////////////////////////////////////////////////
float32_t filter_fir_update(p_filter_fir_t filter_inst, const float32_t x)
{
	float32_t y = 0.0f;
	uint32_t i;

	// Add new sample to buffer
	ring_buffer_add_f( filter_inst->p_x, x );

	// Make convolution
	for ( i = 0; i < filter_inst -> order; i++ )
	{
		y += ( filter_inst->p_a[i] * ring_buffer_get_f( filter_inst->p_x, -i ));
	}

	return y;
}


filter_status_t filter_iir_init(p_filter_iir_t * p_filter_inst, const float32_t * p_a, const float32_t * p_b, const float32_t a_size, const float32_t b_size)
{
	filter_status_t status = eFILTER_OK;


	return status;
}


float32_t filter_iir_update(p_filter_iir_t filter_inst, const float32_t x)
{
	float32_t y = 0.0f;



	return y;
}


////////////////////////////////////////////////////////////////////////////////
/**
* @} <!-- END GROUP -->
*/
////////////////////////////////////////////////////////////////////////////////
