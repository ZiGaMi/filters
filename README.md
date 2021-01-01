# Description

This repository contains simulations in Python using Matplotlib and C implementation tested on STM32 platfrom. All filters evaluated in simulations are realized in C for embeeded application use. 

## List of all filters
 - RC filter (LPF)
 - CR filter (HPF)
 - FIR

## Simulations
For signals written in CSV file use ***filter_csv.py*** script with -f argument to select file to analyse. 

```python
# Example of invocation from console
>>>py filter_csv.py -f ..\..\test_data\steward_imu_roll_pitch_full_range_1_0Hz.csv
```

#### Example of accelerometer data filtering
[example_pic]: simulations/pics/filter_analysis_example.png "Example pic"


## Todo
 - [x] Evaluation of RC filter in python
 - [x] Evaluation of CR filter in python
 - [ ] Implementation of RC filter in C   
 - [ ] Implementation of CR filter in C   
 - [ ] Evaluation of FIR filter in python   
 - [ ] Implementation of FIR filter in C   
 - [ ] Evaluation of washout filter in python
 - [ ] Implementation of washout filter in C

    
