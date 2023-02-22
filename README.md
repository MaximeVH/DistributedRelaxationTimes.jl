# DRT.jl

A Julia package to calculate the **Distribution of Relaxation Times (DRT)** of a given set of electrochemical impedance spectroscopy measurements and frequency values. The current implementation uses the Tikhonov regularisation approach with radial basis function discretization. 

# Usage 
The DRT is a method to analyse **electrochemical impedance spectroscopy (EIS)** measurements. Two of its main benefits are 

 - that it does not require a specific prior equivalent electrical circuit model assumption.
 - it is capable of resolving polarisation processes with similar time constants, which are indistinguishable when using traditional EIS data representation methods, such as Nyquist and Bode plots. 
 
The main function exported in this module is the  `compute_DRT` function. Its two essential arguments are a set of frequencies and the complex-valued impedances measured at those frequencies.

```julia

using DRT, DataFrames, CSV

# Read the measurements.

#assuming the measurements in working directory.
Z_data = CSV.read("example_measurements.csv",DataFrame) 

  

#Obtain the frequencies and measurements from the loaded data.

frequencies = Z_data[:,3]

measurements = Z_data[:,1] .+ (Z_data[:,2] .* im)

  

#Calculacte the DRT using default keyword argument values.

peak_taus, peak_amplitudes, taus_out, drt =  compute_DRT(frequencies, measurements)

```
The `compute_DRT` function has four outputs, these are the relaxation times of the DRT peaks, their amplitudes, the time values over which the DRT is calculated, and the calculated DRT, respectively. Apart from the two essential arguments, the `compute_DRT` function also accepts several keyword arguments:

 - `method` : Determines whether the DRT is calculated using the real part of the measurements ("re"), the imaginary part of the measurements ("im"), or both ("re_im"). The default value is "re_im".
 
 -  `rbf_kernel` : The radial basis function used for the discretisation of the DRT.  The default value is KernelFunctions.jl's `SqExponentialKernel()`.
- `width_coeff` : hyperparameter related to the shape factor of the radial basis functions, lower values lead to wider DRT peaks. The default value is 0.08.
- `λ` :  The hyperparameter tuning the regularisation during ridge regression. The default value is 1e-2.

While there is no universal way to automatically find the most suitable value for the regularisation hyperparameter `λ`,  the discrepancy minimisation and Re-Im crossvalidation methods proposed by  [Saccoccio et al.](https://www.sciencedirect.com/science/article/abs/pii/S0013468614018763) are implemented in this package as the `find_optimal_lambda` function.

Some basic plotting functionality is provided to visualise the outputs of the `compute_DRT` function.  The function `plot_DRT` does this as follows (using the outputs calculated earlier):

 ```julia
#With optional arguments for a label "example" and a green color of the plot.

plot_DRT(peak_taus, peak_amplitudes, taus_out, drt, lab = "example", color = "green")
```