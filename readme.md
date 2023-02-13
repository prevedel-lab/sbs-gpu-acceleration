# Stimulated Brillouin Spectroscopy (SBS)
This library is designed to process data acquired from Stimulated Brillouin Spectroscopy (SBS) experiments.
It provides a user-friendly interface for analyzing SBS data and includes various functions for fitting, visualization and processing of data. 
The use of Matlab ensures compatibility with a wide range of data formats and processing tools, while the incorporation of gpufit allows for fast and efficient fitting of complex models even for large datasets. 
The library is designed to be flexible and easily extendable, making it a valuable tool for researchers in the field of SBS.

# Gpufit custom functions
Follow these instructions to customize gpufit : https://gpufit.readthedocs.io/en/latest/customization.html 

Added for SBS:  
- double 1D lorentzian (absolute shifts) :
    Sum of 2 lorentzians. They each have independent shifts.
    ```    
        f(x) = A_1 * gamma_sqarred_1 / (gamma_squarred_1 + (x - x_0_1)^2) 
            + A_2 * gamma_sqarred_2 / (gamma_squarred_2 + (x - x_0_2)^2)
            + offset
    ```
- double 1D lorentzian (relative shifts) :
    Sum of 2 lorentzian. One has an absolute shift, the other's shift is in respect to that one.
    ```   
    f(x) = A_1 * gamma_sqarred_1 / (gamma_squarred_1 + (x - x_0_1)^2) 
        + A_2 * gamma_sqarred_2 / (gamma_squarred_2 + (x - (x_0_1 + delta))^2)
        + offset 
    ```

From previous Brillouin work : 
- cauchy lorentz 1D :
    ``` f(x) =  A / (1 + [(x-x0)/gamma]**2) + offset ```
- poly 2 :
    Second degree polynomial ```f(x) =  a*x*x + b*x + c```
- stokes :
    Stokes function. Used for an old microscope configuration.
- antistokes :
    Antistokes function. Used for an old microscope configuration.

See this (https://github.com/prevedel-lab/Gpufit) project for gpufit fork already integrating these functions.

# Matlab pipeline
All auxiliary functions needed to process the SBS data are inside the folder *SBS_pipeline*
This includes : 
    - pre-processing steps
    - single peak-fitting
    - double peak-fitting
    - statistics over the fit
    - double peak criteria
    - visualization tool to explore the data set

To use them, you need to have the Matlab binding from gpufit inside a folder called *gpufit-matlab*.
Because of different CUDA versions, you might not be able to directly use the folder here. 
Consider compiling again gpufit on your computer and with a compatible CUDA version to be sure it's working.

## Example
*process_data_xyscan.m* and *process_data_xzscan.m* provides 2 examples of how to use the pipeline.
They are based on xxx dataset. #TODO : Fan approval ?

## Manual correction
Due to a glitch in our setup, the first 2000 points we acquire need their X-coordinate slightly shifted.

## Constraints
For the single peak fitting, the constraints are selected by physical considerations, 
a priori knowledge of the sample and how well the fit works. You can change the constraints
outside of the fitting pipeline.

For the double peak fitting, by default the same constraints as for the single peak are used (with some slight 
changes to adapt for the existence of a second peak). If you need to modify them or add some custom one, 
look into the file *pipeline_double_peak.m*.

**Warning** : There was a bug in gpufit, where only the first set of constraints was used for all the data points. 
That's fixed in the newer versions, but be sure to check which one you're using.
At the time of the development of this tool, the bug was still there, so that's why we have the same constraints
for each point.

See https://github.com/gpufit/Gpufit/issues/105 for more informations.

## Displaying 
The GUI tool can be used to display :
- a single peak fit : ```volume_display(experiment_settings, volume_single_peak, [], []);```
- a double peak fit : ```volume_display(experiment_settings, volume_single_peak, volume_double_peak, []);```
- a double peak fit with statistics : ```volume_display(experiment_settings, volume_single_peak, volume_double_peak, fit_statistics);```

## Relative vs absolute double peak
To switch between the absolute double peak lorentzian and the relative double peak lorentzian
change the file *pipeline_double_peak.m* lines 58-59 :
``` matlab
% Either use absolute or use relative lorentzians
[double_peak(z).parameters, double_peak(z).states, double_peak(z).chi_squares, double_peak(z).number_iterations, double_peak(z).execution_time] = ...
    ... gpufit_double_lorentzian_constrained(data_X(:,:,z), data_Y(:,:,z), experiment_settings.freq_begin, experiment_settings.freq_end, double_weights, previous_peak_fit, double_constraints);
    gpufit_double_lorentzian_relative_constrained(data_X(:,:,z), data_Y(:,:,z), experiment_settings.freq_begin, experiment_settings.freq_end, double_weights, previous_peak_fit, double_constraints);
```
You'll also have to change the constraints for L2 (a bit above) : 

``` matlab
        % Contraints ranges for L_2
        %temp = sort(previous_peak_fit(z).parameters(3,:));
        %delta_x0 = temp(round(0.90 * size(previous_peak_fit(z).parameters, 2))) ;
        delta_x0 = 0.5;
        double_constraints(7, :) = constraints(1, :) ; %/2; % Min amplitude
        double_constraints(8, :) = constraints(2, :) ; %/2; % Max amplitude
        double_constraints(9, :) = -delta_x0 ; %constraints(3, :) ;% -delta_x0 ;% Min x_0 
        double_constraints(10, :) = delta_x0 ; %constraints(4, :) ;% delta_x0 ;% Max x_0
        double_constraints(11, :) = constraints(5, :).^2 ; %/2; % Min gamma
        double_constraints(12, :) = constraints(6, :).^2 ; %/2; % Max gamma
```
Use *+/- delta* for relative double peak, and use *contraints(3,:); contraints(4,:)* for the absolute double peak.

## Definitions for this tool :
- **error** :
    The chi square value for the fit, as returned by gpufit. 
    Look here for more information : https://gpufit.readthedocs.io/en/latest/gpufit_api.html (under *output_chi_squares*)

- **states** : 
    A value returned by gpufit, indicating how the fit ended. 0 is good, everything else is a problem.
    Look here for more information : https://gpufit.readthedocs.io/en/latest/gpufit_api.html (under *output_states*) 


- **discrepancy** :
    The difference between the first and the second fit of a single peak lorentzian model. There can be a difference if you're using different data points or weights. 
    ``` matlab
    discrepancy = sum((previous_peak_fit(z).parameters - single_peak(z).parameters).^2,1)
    ```

- **noise_level** :
    the standard deviation of the recorded signal outside the peak region.  
    ``` matlab
    % peak_width = 2 by default
    low_lim = single_peak(z).parameters(2,n) - peak_width *single_peak(z).parameters(3,n);
    high_lim = single_peak(z).parameters(2,n) + peak_width *single_peak(z).parameters(3,n);
    std_dev = std([data_Y(data_X(:,n) < low_lim ,n)' data_Y(data_X(:,n) > high_lim,n)']);
    ```

- **snr** : The amplitude of the fitted peak, divided by the noise_level.
    ``` matlab
    volume_single_peak.snr(:, :, :) = volume_single_peak.amplitude ./ volume_single_peak.noise_level ;
     ```


# saveastiff_4.5
Multipage TIFF stack
Version 4.5.0.0 (8.28 KB) by Yoon-Oh Tak
Much faster than 'imwrite' and supports complex number images.

Yoon-Oh Tak (2023). Multipage TIFF stack (https://www.mathworks.com/matlabcentral/fileexchange/35684-multipage-tiff-stack), MATLAB Central File Exchange. Retrieved February 13, 2023. 