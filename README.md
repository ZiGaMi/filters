## Description

This repository contains simulations in Python using Matplotlib and C implementation tested on STM32 platfrom. All filters evaluated in simulations are realized in C for embeeded application use. 

---
## List of all filters
 - RC filter (LPF)
 - CR filter (HPF)
 - FIR

---
## Simulations
For signals written in CSV file use ***filter_csv.py*** script with -f argument to select file to analyse. 

```python
# Example of invocation from console
>>>py filter_csv.py -f ..\..\test_data\steward_imu_roll_pitch_full_range_1_0Hz.csv
```

In simulations signal is being shown in time and frequency domain. By inspecting frequency signature of raw filter cutoff frequency of selected filter can be easily determine. 

#### Example of accelerometer data filtering
This example shows signal acquire from accelerometer and RC low pass filter in work. Upper picture shows raw and filtered signals in time domain and bottom shows same signals in frequency domain.

![](simulations/pics/filter_analysis_example.png)

---
## C implementation
#### RC filter
There are only two functions being a part of RC filter API:
 - *filter_status_t ***filter_rc_init***(p_filter_rc_t * p_filter_inst, const float32_t fc, const float32_t dt, const uint8_t order, const float32_t init_value)*
 - *float32_t ***filter_rc_update***(p_filter_rc_t p_filter_inst, const float32_t x)*

 ##### Example of usage
```
// 1. Declare filter instance
p_filter_rc_t my_filter_inst;

/* 
*   2. Init RC filter with following parameters:
*   - fc = 10Hz
*   - order = 1
*   - inititial value = 0
*/ 
if ( eFILTER_OK != filter_rc_init( &my_filter_instance, 10.0f, SAMPLE_TIME, 1, 0 ))
{
    // Filter init failed
    // Further actions here...
}

// 3. Apply filter in period of SAMPLE_TIME
loop @SAMPLE_TIME
{
    // Update filter
    filtered_signal = filter_rc_update( my_filter_inst, raw_signal );
}
```


---
## TODO
 - [x] Evaluation of RC filter in python
 - [x] Evaluation of CR filter in python
 - [x] Implementation of RC filter in C   
 - [ ] Implementation of CR filter in C   
 - [ ] Evaluation of FIR filter in python   
 - [ ] Implementation of FIR filter in C   
 - [ ] Evaluation of washout filter in python
 - [ ] Implementation of washout filter in C

    
