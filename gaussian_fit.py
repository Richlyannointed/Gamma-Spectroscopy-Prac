"""
It does the fit
for the Gauss
Give it Data

"""

from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import scienceplots
import numpy as np
import csv

plt.rcParams['font.size'] = 25
plt.style.use('science')

def import_data(filename: str) -> tuple:
    """
    Returns tuple of data array and elapsed time in [s]
    """
    elapsed_time = extract_elapsed_live_time(filename) # For count rate conversion
    if elapsed_time is None:
        return None
    
    skip_rows = 22 # USX Header
    data = []
    with open(filename, 'r') as file:
        csv_reader = csv.reader(file)
        for _ in range(skip_rows):
            next(csv_reader)  # Skip rows
        for row in csv_reader:
        # Skip the middle column (assuming the array has 3 columns)
            channel = float(row[0]) if row[0] else None
            counts_per_sec = float(row[2]) / elapsed_time if row[2] else None  
            u_counts_per_sec = np.sqrt(counts_per_sec) / elapsed_time if counts_per_sec else None
            data.append((channel, counts_per_sec, u_counts_per_sec))
        
    # Convert list to NumPy array
    try:
        numpy_array = np.array(data).astype(float)
    except TypeError:
        print('Failed to cast data to float')
        return None
    return numpy_array, elapsed_time


def extract_elapsed_live_time(filename):
    """Scrapes 'live time' value from USX csv file
    
    parameters
    ----------
    filename : str
        Name of raw csv spectrum file

    returns
    -------
    elapsed_live_time : float
    """
    with open(filename, 'r') as file:
        csv_reader = csv.reader(file)

        # Skip rows until 'Calibration Coefficients' is reached
        for row in csv_reader:
            if 'Elapsed Live Time:' in row[0]:
                try:
                    elapsed_live_time = float(row[1])
                except ValueError:
                    print(f'Could not convert {row[1]} to float')
                    return 1

                return elapsed_live_time


def peak_finder(data: np.ndarray) -> np.ndarray:
    """ This is just scipy's `find_peaks` function with `width=10`

    Returns
    -------
    peaks : ndarray
        Indices of peak positions
        NB! indices are w.r.t. the passed `data` array NOT the uncropped array.
    """
    peaks, _ = find_peaks(data, width=10)
    return peaks


def find_peak_ranges(data: np.ndarray) -> list:
    """
    Returns
    -------
    ranges : tuple
        peak ranges based on `CONSTANT`: avg_peak_width
    """
    peaks = peak_finder(data[:,1])
    ranges = []
    for peak in peaks:
        ranges.append((peak - avg_peak_width, peak + avg_peak_width))
    return ranges


def gauss(x, A, mu, sigma, C) -> float:
        return A * np.exp(- (x - mu)**2 / (2 * sigma**2)) + C
    

def gaussian_fit(data: np.ndarray, peak_range: tuple) -> tuple:
    """
    Fit a Gaussian function to the data within the specified peak_range.

    Parameters
    ------------
    data : ndarray
        3-column data array in format: x | y | u(y)

    Returns
    -------
    popt : tuple
        fit parameters

    perr : tuple
        parameters' fit uncertainty with chi^2/dof correction
    """
    #Optimisation

    x_range = data[peak_range[0]:peak_range[1], 0]
    y_range = data[peak_range[0]:peak_range[1], 1]
    sigma = data[peak_range[0]:peak_range[1], 2]

    # Initial guesses based on data within the specified range
    A0 = np.max(y_range) - np.min(y_range)
    mu0 = x_range[np.argmax(y_range)]
    sigma0 = 1.0  # Initial guess for sigma
    C0 = np.min(y_range)
    p0 = [A0, mu0, sigma0, C0] # Guess List ;)

    try:
        # Adding lower bound for C parameter
        bounds = ([0, 0, 0, 0], [np.inf, np.inf, np.inf, np.inf])
        popt, pcov = curve_fit(gauss, x_range, y_range, p0=p0, sigma=sigma, absolute_sigma=True, bounds=bounds)

        # Evaluation 1
        dymin = (y_range - gauss(x_range, *popt)) / sigma
        min_chisq = sum(dymin * dymin)
        dof = len(x_range) - len(popt)
        perr = np.sqrt(np.diag(pcov)) * np.sqrt(min_chisq / dof)
        return popt, perr
    
    except RuntimeError:
        print("Fit failed for peak range:", peak_range)
        return None, None
    

def show_fits(data: np.ndarray, fit_params_dict: dict,  source_name: str) -> None:
    """
    Parameters 
    ----------
    fit_dict : dict
        dictionary of Gaussian fit parameters w.r.t. `data` axes

    Returns
    ------- 
    None : None 
        Spectrum plot with fits of all found peaks    
    
    """
    fig, ax = plt.subplots(1, figsize=(8, 5))
    ax.errorbar(x=data[:,0], y=data[:,1], yerr=data[:,2], 
                capsize=2, 
                alpha=0.2, 
                fmt='.', 
                color='k',
                ecolor='r', 
                label='Source Spectrum')
    
    for peak_name, fit_params in fit_params_dict.items():
        A = fit_params['A']
        mean = fit_params['mu']
        sigma = fit_params['sigma']
        C = fit_params['C']
        # Generate x values within the range of the data
        x_range = np.linspace(mean - sigma*4 , mean + sigma*4, 1000)
        # Generate y values for the best-fit Gaussian curve
        y_fit = gauss(x_range, A, mean, sigma, C)
        # Plot the best-fit Gaussian curve for the current peak
        ax.plot(x_range, y_fit, 
                linewidth=1.5,
                alpha=0.8)
    
    ax.set(xlabel='Channel Number',
           ylabel='Count Rate [Counts/s]',
           yscale='linear',
           xlim=(30, 400), 
           ylim=(0, 40))
    plt.legend()
    # plt.savefig(source_name + '.png', dpi=400)
    plt.show(block=True)    


def write_fit_params_to_csv(fit_params_dict, elapsed_time, source_name, output_filename):
    # Extract parameter names from the first fit_params entry
    param_names = list(fit_params_dict.values())[0].keys()

    # Write header to CSV file
    with open(output_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        header = ['Peak Name', 'Elapsed Time', 'Source Name'] + list(param_names)
        writer.writerow(header)
        # Write data rows
        for peak_name, fit_params in fit_params_dict.items():
            row = [peak_name, elapsed_time, source_name] + [fit_params[param] for param in param_names] 
            writer.writerow(row)


#CONSTANTS | Change as necessary
avg_peak_width = 30 #Average width of energy peaks [in channels] 

# best values:
'''
Mn54 : 37
Cs37 : 
Na22 : 
Co60 :
Co57 :
Unk : 
Background : NO PEAKS FOUND
'''

def main():
    # USAGE: Change `source_filename` to your source and run file.
    # OUTPUT: .csv of fit parameters and (TODO) fit uncertainties
    source_filename = 'Cs137.csv'

    source_name = source_filename.strip('.csv')
    spec, elapsed_time = import_data(source_filename)

    # GET PEAK RANGES BASED ON AVERAGE PEAK WIDTH
    ranges = find_peak_ranges(spec)
    
    # CREATE FIT PARAMETER DICTIONARY FOR EACH PEAK
    fit_params_dict = {}
    for i, peak_range in enumerate(ranges):
        peak_name = f'Peak_{i+1}' 
        fit_params, fit_error = gaussian_fit(spec, peak_range)
        if fit_params is not None:
            fit_params_dict[peak_name] = {'A' : fit_params[0], 
                                          'u(A)' : fit_error[0],
                                          'mu' : fit_params[1], 
                                          'u(mu)' : fit_error[1],
                                          'sigma' : fit_params[2],
                                          'u(sigma)' : fit_error[2],
                                          'C' : fit_params[3],
                                          'u(C)' : fit_error[3]}

    show_fits(spec, fit_params_dict, source_name)
    write_fit_params_to_csv(fit_params_dict, elapsed_time, source_name, source_name + '_fit.csv')
    
    
if __name__ =='__main__':
    main()
