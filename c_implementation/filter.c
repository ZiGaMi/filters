////////////////////////////////////////////////////////////////////////////////
/**
*@file      filter.c
*@brief     Various filter designs
*@author    Ziga Miklosic
*@date      01.01.2021
*@version   V0.0.1
*
*@mainpage
*@section   Description
*   
*   This module contains different kind of digital filter
*   implementation. All of the following filter types are
*   being simulated in python for validation purposes.
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

////////////////////////////////////////////////////////////////////////////////
// Definitions
////////////////////////////////////////////////////////////////////////////////

/**
 * 	RC Filter data
 */
typedef struct filter_rc_s
{
	float32_t *			p_y;		/**<Output of filter + previous value*/
	float32_t			alpha;		/**<Filter smoothing factor*/
	uint8_t				order;		/**<Filter order */
} filter_rc_t;


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
* @param[in] 	p_filter	- Pointer to RC filter instance
* @param[in] 	fc			- Filter cutoff frequency
* @param[in] 	dt			- Sample time
* @param[in] 	order		- Order of filter
* @param[in] 	init_value	- Initial value
* @return 		status		- Status of operation
*/
////////////////////////////////////////////////////////////////////////////////
filter_status_t filter_rc_init(p_filter_rc_t * p_filter_inst, const float32_t fc, const float32_t dt, const uint8_t order, const float32_t init_value)
{
	filter_status_t status = eFILTER_OK;
	uint8_t i;

	if (( NULL != p_filter_inst ) && ( order > 0 ))
	{
		// Allocate space
		*p_filter_inst 			= malloc( sizeof( filter_rc_t ));
		(*p_filter_inst)->p_y 	= malloc( order  * sizeof( float32_t ));

		// Check if allocation succeed
		if (( NULL != *p_filter_inst ) && ( NULL != (*p_filter_inst)->p_y ))
		{
			// Calculate coefficient
			(*p_filter_inst)->alpha = (float32_t) ( dt / ( dt + ( 1.0f / fc )));

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
* @param[in] 	p_filter	- Pointer to RC filter instance
* @param[in] 	x			- Input value
* @return 		y			- Output value
*/
////////////////////////////////////////////////////////////////////////////////
float32_t filter_rc_update(p_filter_rc_t p_filter_inst, const float32_t x)
{
	float32_t y = 0.0f;
	uint8_t n;

	if ( NULL != p_filter_inst )
	{
		for ( n = 0; n < p_filter_inst->order; n++)
		{
			if ( 0 == n )
			{
				p_filter_inst->p_y[0] = ( p_filter_inst->p_y[0] + ( p_filter_inst->alpha * ( x - p_filter_inst->p_y[0] )));
			}
			else
			{
				p_filter_inst->p_y[n] = ( p_filter_inst->p_y[n] + ( p_filter_inst->alpha * ( p_filter_inst->p_y[n-1] - p_filter_inst->p_y[n] )));
			}
		}

		y = p_filter_inst->p_y[ p_filter_inst->order - 1U ];
	}

	return y;
}


////////////////////////////////////////////////////////////////////////////////
/**
* @} <!-- END GROUP -->
*/
////////////////////////////////////////////////////////////////////////////////
