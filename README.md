# Gamma Spectroscopy Practical Toolkit

Welcome to the Gamma Spectroscopy Practical Toolkit! This toolkit provides a set of Python scripts for analyzing energy spectra obtained from gamma spectroscopy experiments. It is designed to assist in the identification of radioactive samples based on their characteristic peaks.

## Features

- **Data Import:** Easily import energy spectra data from CSV files.
- **Peak Detection:** Automatically detect peaks in the spectra using signal processing techniques.
- **Gaussian Fitting:** Fit Gaussian curves to the identified peaks for peak characterization.
- **Interactive Plotting:** Visualize the original spectra along with the fitted Gaussian curves for easy analysis and interpretation.

## Usage

1. **Installation:** Clone this repository to your local machine.

   ```bash
   git clone https://github.com/yourusername/Gamma-Spectroscopy-Prac.git

2. **Dependencies** Install the required Python dependencies using pip.
   ```bash
   pip install -r requirements.txt

3. **Data Preparation** Prepare your energy spectra data in CSV format. Ensure that the data includes channel numbers and corresponding intensity values.
    
4. **Analysis** Use the provided Python scripts to analyze your data.
    ```bash
   python gaussian_fit.py source_spectrum.csv
      
6. **Visualisation** Visualize the analyzed spectra and Gaussian fits using the generated plots.

## Contributing
Contributions to the Gamma Spectroscopy Practical Toolkit are welcome! If you have any ideas, suggestions, or improvements, please feel free to open an issue or submit a pull request on GitHub.

## Acknowledgments
- This toolkit was developed as part of the Gamma Spectroscopy Practical at the University of Cape Town.
- Special thanks to Peanut Butter and Jam for their contributions and support.

