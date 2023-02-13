#ifndef GPUFIT_DOUBLE_LORENTZIAN_1D_INCLUDED
#define GPUFIT_DOUBLE_LORENTZIAN_1D_INCLUDED

/* Description of the calculate_double_lorentzian_1d function
* ===================================================
* Added by Sebastian Hambura the 05/2022
* 
* This function calculates the values for a double lorentzian function :
* it's the sum of 2 Lorentzian functions. See cauchy_lorentz_1d for more information 
* about that kind of functions (or https://en.wikipedia.org/wiki/Cauchy_distribution)
*
* The expression of the fitted function is :
* f(x) = A_1 * gamma_sqarred_1 / (gamma_squarred_1 + (x - x_0_1)^2) 
*           + A_2 * gamma_sqarred_2 / (gamma_squarred_2 + (x - x_0_2)^2)
*           + offset
*
* This function makes use of the user information data to pass in the
* independent variables (X values) corresponding to the data.  The X values
* must be of type REAL.
*
* Note that if no user information is provided, the (X) coordinate of the
* first data value is assumed to be (0.0).  In this case, for a fit size of
* M data points, the (X) coordinates of the data are simply the corresponding
* array index values of the data array, starting from zero.
*
* There are three possibilities regarding the X values:
*
*   No X values provided:
*
*       If no user information is provided, the (X) coordinate of the
*       first data value is assumed to be (0.0).  In this case, for a
*       fit size of M data points, the (X) coordinates of the data are
*       simply the corresponding array index values of the data array,
*       starting from zero.
*
*   X values provided for one fit:
*
*       If the user_info array contains the X values for one fit, then
*       the same X values will be used for all fits.  In this case, the
*       size of the user_info array (in bytes) must equal
*       sizeof(REAL) * n_points.
*
*   Unique X values provided for all fits:
*
*       In this case, the user_info array must contain X values for each
*       fit in the dataset.  In this case, the size of the user_info array
*       (in bytes) must equal sizeof(REAL) * n_points * nfits.
*
* Parameters:
*
* parameters: An input vector of model parameters
*             p[0]: Amplitude A_1
*             p[1]: Shift/center x0_1
*             p[2]: gamma_squarred_1
*             p[3]: Amplitude A_2
*             p[4]: Shift/center x0_2
*             p[5]: gamma_squarred_2
*             p[6]: offset
*
* n_fits: The number of fits.
*
* n_points: The number of data points per fit.
*
* value: An output vector of model function values.
*
* derivative: An output vector of model function partial derivatives.
*
* point_index: The data point index.
*
* fit_index: The fit index.
*
* chunk_index: The chunk index. Used for indexing of user_info.
*
* user_info: An input vector containing user information.
*
* user_info_size: The size of user_info in bytes.
*
* Calling the calculate_linear1d function
* =======================================
*
* This __device__ function can be only called from a __global__ function or an other
* __device__ function.
*
*/

__device__ void calculate_double_lorentzian_1d(
    REAL const* parameters,
    int const n_fits,
    int const n_points,
    REAL* value,
    REAL* derivative,
    int const point_index,
    int const fit_index,
    int const chunk_index,
    char* user_info,
    std::size_t const user_info_size)
{
    // indices

    REAL* user_info_float = (REAL*)user_info;
    REAL x = 0;
    if (!user_info_float)
    {
        x = point_index;
    }
    else if (user_info_size / sizeof(REAL) == n_points)
    {
        x = user_info_float[point_index];
    }
    else if (user_info_size / sizeof(REAL) > n_points)
    {
        int const chunk_begin = chunk_index * n_fits * n_points;
        int const fit_begin = fit_index * n_points;
        x = user_info_float[chunk_begin + fit_begin + point_index];
    }

    // parameters
    REAL A_1                = parameters[0];
    REAL x0_1               = parameters[1];
    REAL gamma_squarred_1   = parameters[2];
    REAL A_2                = parameters[3];
    REAL x0_2               = parameters[4];
    REAL gamma_squarred_2   = parameters[5];
    REAL offset             = parameters[6];

    // 1st lorentzian
    REAL delta_1 = x - x0_1;
    REAL denominator_1 = gamma_squarred_1 + delta_1 * delta_1;
    REAL L_1 = A_1 * gamma_squarred_1 / denominator_1;

    // 2nd lorentzian
    REAL delta_2 = x - x0_2;
    REAL denominator_2 = gamma_squarred_2 + delta_2 * delta_2;
    REAL L_2 = A_2 * gamma_squarred_2 / denominator_2;

    // value
    //REAL denominator = gamma * gamma + (x - x0) * (x - x0);
    value[point_index] = L_1 + L_2 + offset;

    // derivatives
    REAL squarred_denominator_1 = denominator_1 * denominator_1;
    REAL squarred_denominator_2 = denominator_2 * denominator_2;
    REAL* current_derivatives = derivative + point_index;

    current_derivatives[0 * n_points] = gamma_squarred_1 / denominator_1;                                   // derivative A_1
    current_derivatives[1 * n_points] = 2 * A_1 * gamma_squarred_1  * delta_1 / squarred_denominator_1;  // derivative x0_1
    current_derivatives[2 * n_points] = A_1 * delta_1 * delta_1 / squarred_denominator_1;             // derivative gamma_squarred_1
    current_derivatives[3 * n_points] = gamma_squarred_2 / denominator_2;                                   // derivative A_2
    current_derivatives[4 * n_points] = 2 * A_2 * gamma_squarred_2 * delta_2 / squarred_denominator_2;   // derivative x0_2
    current_derivatives[5 * n_points] = A_2 * delta_2 * delta_2 / squarred_denominator_2;             // derivative gamma_squarred_2
    current_derivatives[6 * n_points] = 1;                                                                  // derivative offset
}

#endif
