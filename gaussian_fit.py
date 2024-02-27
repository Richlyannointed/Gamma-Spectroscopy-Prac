"""
It does the fit
for the Gauss
Give it Data

"""

from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import numpy as np
import csv

#CONSTANTS
avg_peak_width = 35 #Average width of energy peaks [in channels]

def import_data(filename: str) -> np.ndarray:
    skip_rows = 22 # USX Header
    data = []
    with open(filename, 'r') as file:
        csv_reader = csv.reader(file)
        for _ in range(skip_rows):
            next(csv_reader)  # Skip rows
        for row in csv_reader:
        # Skip the middle column (assuming the array has 3 columns)
            filtered_row = [row[0], row[2]]  # Keep first and third column, skip the middle one
            data.append(filtered_row)

    # Convert list to NumPy array
    try:
        numpy_array = np.array(data).astype(float)
    except TypeError:
        print('Failed to cast data to float')
        return None
    return numpy_array



def display_spectrum(data:np.ndarray, source_name: str, yscale='linear'):
    peaks = peak_finder(data[:,1])
    fig, ax = plt.subplots(1, figsize=(9, 5))
    x = data[:,0]
    y = data[:,1]
    ax.scatter(x=x, y=y, marker='.', alpha=0.4, color='black', label='Source Spectrum')
    ax.vlines(peaks, ymin=0, ymax=max(y), color='red', linestyle='--')
    ax.set(title=f'Spectrum of {source_name.upper()}',
           xlabel='Channel No.',
           yscale=yscale)
    plt.show(block=False)

    """
    # Ask user for input
    while True:
        plt.show(block=False)  # Show the plot
        start = input("Enter the start index of the data range (or 'exit' to close the plot): ")
        if start.lower() == 'exit':
            plt.close()  # Close the plot window
            break
        end = input("Enter the end index of the data range: ")
        
        # Validate input
        try:
            start = int(start)
            end = int(end)
            if start < 0 or end < 0 or end <= start or end >= len(x):
                print("Invalid input. Please try again.")
            else:
                plt.plot(x[start:end+1], y[start:end+1])  # Plot the selected range
                plt.xlabel('X-axis')
                plt.ylabel('Y-axis')
                plt.title('Selected Range')
                plt.show()  # Show the selected range plot
                break
        except ValueError:
            print("Invalid input. Please enter integers only.")
    """

def peak_finder(data: np.ndarray) -> list:
    """
    Returns indices of peak positions
    NB! indices are w.r.t. the passed `data` array.
    """
    peaks, _ = find_peaks(data, width=10)
    return peaks


def find_peak_ranges(data: np.ndarray) -> tuple:
    """

    """
    
    peaks = peak_finder(data[:,1])
    ranges = []
    for peak in peaks:
        ranges.append([peak - avg_peak_width, peak + avg_peak_width])
    return ranges


def gauss(x, A, mu, sigma, C) -> float:
        return A * np.exp(- (x - mu)**2 / (2 * sigma**2)) + C
    

def gaussian_fit(data: np.ndarray, peak_range: tuple) -> tuple:
    """
    Fit a Gaussian function to the data within the specified peak_range.
    """
    #Optimisation

    x_range = data[peak_range[0]:peak_range[1], 0]
    y_range = data[peak_range[0]:peak_range[1], 1]

    # Initial guesses based on data within the specified range
    A0 = np.max(y_range) - np.min(y_range)
    mu0 = x_range[np.argmax(y_range)]
    sigma0 = 1.0  # Initial guess for sigma
    C0 = np.min(y_range)
    p0 = [A0, mu0, sigma0, C0]

    try:
        popt, pcov = curve_fit(gauss, x_range, y_range, p0=p0)
        return popt, pcov
    except RuntimeError:
        print("Fit failed for peak range:", peak_range)
        return None, None
    

def show_fits(data: np.ndarray, fit_params_dict: dict,  source_name: str) -> None:
    """
    Returns spectrum plot with fits of all found peaks    
    
    fit_dict : dictionary of Gaussian fit parameters w.r.t. `data` axes
    """
    fig, ax = plt.subplots(1, figsize=(9, 5))
    ax.scatter(x=data[:,0], y=data[:,1], marker='^', alpha=0.4, color='cyan', label='Source Spectrum')
    for peak_name, fit_params in fit_params_dict.items():
        # Generate y values for the best-fit Gaussian curve
        y_fit = gauss(data[:, 0], *fit_params)

        # Plot the best-fit Gaussian curve for the current peak
        ax.plot(data[:, 0], y_fit, alpha=0.6, label=f'{peak_name} Fit')

    
    ax.set(title=f'Spectrum of {source_name.upper()}',
           xlabel='Channel No.',
           ylabel='Count Rate [Counts/sec]',
           yscale='linear')
    plt.legend()
    plt.show(block=True)    


def main():
    source_filename = 'Cs137.csv'
    source_name = source_filename.strip('.csv')
    spec = import_data(source_filename)
    sigma = np.ones_like(spec[:,1]) * (1 / np.sqrt(spec[:,1]))
    spec = np.column_stack((spec, sigma))

    #display_spectrum(spec, source_name, 'log')
    ranges = find_peak_ranges(spec)
    print(find_peak_ranges(spec))
    
    # VISUALLY INSPECT GAUSSIAN FITS
    fit_params_dict = {}
    for i, peak_range in enumerate(ranges):
        peak_name = f'Peak_{i+1}' 
        fit_params, _ = gaussian_fit(spec, peak_range)
        if fit_params is not None:
            fit_params_dict[peak_name] = fit_params

    show_fits(spec, fit_params_dict, source_name)
    
    

if __name__ =='__main__':
    main()
