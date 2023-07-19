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