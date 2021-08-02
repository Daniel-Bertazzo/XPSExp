import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import os
from scipy.ndimage import shift


class XPS_Experiment:

    def __init__(self):
        
        self.header       = None
        self.data         = None
        self.cut_data     = None
        self.int_data     = None

        self.dim1_scale   = None
        self.dim2_scale   = None
        self.dim1_name    = None
        self.dim2_name    = None

        self.edges        = None
        self.int_bins     = None
        self.bin_centers  = None
        self.nbins        = None
        
        self.max_points   = None
        self.popt         = None
        self.f            = None
        self.deg          = None
        self.f_fix        = None

        self.roi_init     = None
        self.roi_end      = None

        self.__binned     = False
        self.__integrated = False
        self.__cut        = False
        self.__bin_int    = False
        self.__max_points = False
        self.__fit        = False


    def read_file(self, filename):
        """Reads text file and stores it in a XPS_Experiment object.

        Args:
            filename (string): name of the file to be read and loaded.

        """
        # ..:: Reading file header ::.. 

        # Reads file as csv, using "=" as a separator. This way, all the data in the file's header is structured
        # in a single-column dataframe, where which row stores the data piece name (represented by the row index) and its
        # corresponding value.
        self.header = pd.read_csv(filename, sep="=")

        # Reading the file using pandas causes the data rows in the [Data] section of the file to be stored as indices with no corresponding value (NaN).
        # Because of this, all NaN rows are dropped, resulting in a dataframe containing only the header data.
        self.header = self.header.dropna(axis='rows')

        # Splitting each dimension scale into a list of strings
        self.dim1_scale = self.header.loc["Dimension 1 scale"][0].split()
        self.dim2_scale = self.header.loc["Dimension 2 scale"][0].split()

        # Converting the list of strings into a numpy array of floats
        self.dim1_scale = np.asarray(self.dim1_scale).astype(float)
        self.dim2_scale = np.asarray(self.dim2_scale).astype(float)

        # Saves dimension names for easier use later
        self.dim1_name = self.header.loc["Dimension 1 name"][0]
        self.dim2_name = self.header.loc["Dimension 2 name"][0]



        # ..:: Reading file data ::..

        with open(filename) as f:

            # Skips file header (the goal here is to read only the [Data] section of the file)
            while "[Data" not in f.readline():
                continue

            # Reads first data line and transforms it into a numpy array of floats (from string values)
            self.data = f.readline()
            self.data = np.asarray(self.data.split()).astype(float)

            for data_line in f:
                if data_line != "\n":
                    # Transforms data line from string to float array
                    transformed_data = np.asarray(data_line.split()).astype(float)
                    # Stacks new data line into a matrix
                    self.data = np.vstack((self.data, transformed_data))

        # Removes first column (which contains duplicated values already present in the dimension 1 scale)
        self.data = self.data[:,1:]

        # Transforms the numpy array to a DataArray from xarray library
        self.data = xr.DataArray(data=self.data,
                                 dims=(self.dim1_name, self.dim2_name),
                                 coords={self.dim1_name: self.dim1_scale, self.dim2_name: self.dim2_scale}
                                )

        # Sets region of interests to entire data
        self.roi_init = self.data[self.dim2_name][0]
        self.roi_end  = self.data[self.dim2_name][-1]

    
    def plot(self,
             show_bins       = False,
             show_roi        = False,
             show_max_points = False,
             show_fit        = False):
        """Plots data contained in the object.

        Args:
            show_bins (bool, optional): if True, the bins in which the data has been divided into
                will be shown on the plot. Previous call to ``XPS_Experiment.divide()`` required. Default False.
            show_roi (bool, optional): if True, the region of interest will be shown on the plot. If the method
                ``XPS_Experiment.set_roi()`` has not been previously called, the region of interest will be defaulted to
                the data's original boundaries. Default False.
            show_max_points (bool, optional): if True, the points of maximum kinetic energy of each bin
                will be shown on the plot. Previous call to ``XPS_Experiment.calculate_max_points()`` required. Default False.
            show_fit (bool, optional): if True, the curve which best fits the max points will be shown on the plot.
                Previous call to ``XPS_Experiment.fit_max_points()`` required. Default False.

        """             

        if show_bins == True and self.__binned == False:
            print("Error: data not binned yet. Use method divide_data()")

        fig, ax = plt.subplots()

        self.data.plot.pcolormesh(cmap='inferno', ax=ax)
        
        ax.set_title("Particle position (mm) vs Kinetic energy (eV)")

        # Show bins in plot
        if show_bins == True:

            ax.set_xlim(left=self.data[self.dim2_name][0], right=self.data[self.dim2_name][-1])

            if self.edges[-1] == self.data.loc[:, self.roi_init : self.roi_end][self.dim2_name].data[-1]:
                # [ax.axvline(x=edge, color='red') for edge in self.edges[1:-1]]
                [ax.axvline(x=self.cut_data[self.dim2_name][edge], color='red') for edge in self.edges[1:-1]]
            else:
                # [ax.axvline(x=edge, color='red') for edge in self.edges[1:]]
                [ax.axvline(x=self.cut_data[self.dim2_name][edge], color='red') for edge in self.edges[1:-1]]
        
        if show_roi == True:
            ax.axvline(x=self.roi_init, color='cyan')
            ax.axvline(x=self.roi_end , color='cyan')

        if show_max_points == True:
            # ax.scatter(x=self.bin_centers, y=self.max_points, color='white')
            x = self.cut_data[self.dim2_name][self.bin_centers]
            y = self.cut_data[self.dim1_name][self.max_points]
            ax.scatter(x=x, y=y, color='white')
            

        if show_fit == True:
            x = np.linspace(self.roi_init, self.roi_end, 200)
            y = self.f_plot(x)
            ax.plot(x, y, color='white')

        plt.draw()
        plt.show()


    def integrate(self):
        """Sums all values for each level of kinetic energy and stores it in a 1-D array.

        """
        if self.__cut == True:
            self.int_data = self.cut_data.sum(dim=self.dim2_name)
        else:
            self.int_data = self.data.sum(dim=self.dim2_name)

        self.__integrated = True

        # For plotting later
        self.int_data.attrs["long_name"] = "Photon count"


    def plot_integrated(self):
        """Plots the integrated data (kinetic energy vs photon count).

        """
        if self.__integrated == False:
            print("Error: data not integrated yet. Use method integrate_data()")
            return

        fig, ax = plt.subplots()

        self.int_data.plot.line(ax=ax)
        
        ax.grid(linestyle='dotted')
        ax.set_title("Photon count as function of kinetic energy (eV)")

        plt.draw()
        plt.show()


    def divide(self, nbins=10):
        """Divides the data into bins (columns).

        Args:
            nbins (int, optional): number of bins to divide the data into. Default 10.

        """
        if nbins <= 1:
            print("Error: Invalid number of bins. Must be greater than 1.")
            return

        if self.__cut == False:
            print("Error: Region of interest not set yet. Use method XPS_Experiment.set_roi()")
            return

        self.nbins = nbins

        # Get bins size
        bin_size  = self.cut_data.shape[1] // nbins
        remainder = self.cut_data.shape[1] % nbins


        self.edges = []
        for i in range(0, self.cut_data[self.dim2_name].shape[0], bin_size):
            self.edges.append(i)

        self.edges[-1] = self.cut_data.shape[1] - 1
        
        self.__binned = True


    def integrate_bins(self):
        """Sums all values for each kinetic energy level, separating by bin.

        """
        if self.__binned == False:
            print("Error: data not binned yet. Use method XPS_Experiment.divide()")
            return
        
        # Creates DataArray that will hold the integrated bins
        self.int_bins = xr.DataArray(np.zeros((len(self.edges)-1, self.data.shape[0])),
                                     dims=("Bins", self.dim1_name),
                                     coords={"Bins": range(len(self.edges)-1), self.dim1_name: self.dim1_scale}
                                    )

        self.bin_centers = np.zeros((self.int_bins.shape[0]), dtype=int)

        # Integrates each column    
        for i in range(len(self.edges)-1):
            self.int_bins[i, :] = self.cut_data[:, self.edges[i]:self.edges[i+1]].sum(dim=self.dim2_name)
            self.bin_centers[i] = (self.edges[i] + self.edges[i+1]) // 2
        
        self.__bin_int = True


    def plot_columns(self, separate=True, vdist=1.0):
        """Plots each previously integrated bin.

        Args:
            separate (bool, optional): if True, the curves will be plotted in different,
                separated levels. If False, the curves will be plotted on the same scale,
                overlapping eachother. Default True.
            vdist (float, optional): Factor that changes the distance between the curves,
                if ``separate`` is True. Default 1.0.

        """
    
        # Manipulates data to separate the lines and improve visualization
        if separate == True:
            # Finds max value in data
            pad_value = np.max(self.int_bins.data)
        
            # Uses the max value, along with the vertical distance, to separate the lines of integrated data for better visualization
            for i in range(self.int_bins.shape[0]):
                self.int_bins[i, :] += i * vdist * pad_value
        
            # Sets plot size
            fig, ax = plt.subplots(figsize=(5, 8))
        
        # Sets different size for plot when not separating the columns
        else:
            fig, ax = plt.subplots(figsize=(12, 6))
                

        # Plotting
                
        for i in range(self.int_bins.shape[0]):
            label_text = f'Bin {i+1}'
            self.int_bins[i, :].plot.line(ax=ax, label=label_text)

        ax.grid(linestyle="dotted", axis='x')
        ax.set_title("Bin intensity offset")
        ax.set_yticks([])
    
        ax.legend(fontsize="large",
                bbox_to_anchor=(1.05, 1.0),
                loc="upper left")
            
        plt.draw()
        plt.show()
        
        # Reverts data back to normal if it was previously manipulated
        if separate == True:
            for i in range(self.int_bins.shape[0]):
                self.int_bins[i, :] -= i * vdist * pad_value
    

    def set_roi(self, init, end):
        """Sets region of interest to the data, wherein calculations will be performed. Data outside this region will be ignored by most methods.

        Args:
            init (float): new initial value (beginning of the slice).
            end (float): new last value (end of the slice).

        """
        self.roi_init = init
        self.roi_end  = end

        # Slices data to only consider the interior of the region of interest
        self.cut_data = self.data.loc[:, self.roi_init : self.roi_end]
        
        self.__cut = True

  
    def calculate_max_points(self):
        """Calculates points of maximum kinetic energy for each bin.

        """

        if self.__bin_int == False:
            print("Error: bins not integrated yet. Use XPS_Experiment.integrate_bins()")
            return

        # Finding max value of each integrated bin
        self.max_points = np.zeros((self.int_bins.shape[0]), dtype=int)

        for i in range(self.int_bins.shape[0]):
            self.max_points[i] = self.int_bins[i].argmax()
        
        self.__max_points = True


    def fit_max_points(self, deg):
        """Fits a curve to the coordinates of the previously calculated max points.

        Args:
            deg (int, positive): degree of polynomial function to be fitted into the points.

        """

        if self.__max_points == False:
            print("Error: max points not calculated yet. Use XPS_Experiment.calculate_max_points()")
            return

        # Fits function in the pixel scale so that it can be fixed with shifting later on
        self.popt = np.polyfit(x=self.bin_centers, y=self.max_points, deg=deg)
        self.f    = np.poly1d(self.popt)
        self.deg  = deg
        # Ignores the 0th power coefficient
        self.f_fix = np.poly1d(np.append(self.popt[:-1], 0.0))

        # Fits function in dimensions scale to enable plotting
        popt = np.polyfit(x=self.cut_data[self.dim2_name][self.bin_centers],
                          y=self.cut_data[self.dim1_name][self.max_points],
                          deg=deg
                         )
        self.f_plot = np.poly1d(popt)

        self.__fit = True


    def fix_distortion(self):
        """Fixes the distortion contained in the data. Only works for distortions of up to 2nd degree (linear and quadratic distortions).

        """

        if self.__fit == False:
            print("Error: function not fit yet. Use XPS_Experiment.fit_max_points()")
            return

        for y in range(self.data.shape[1]):
            self.data[:, y] = shift(self.data[:, y], -self.f_fix(y), prefilter=False, order=0, mode='constant')

        # If distortion is of 2nd degree, apply correction one more time to fix remaining linear distortion
        if self.deg == 2:
            
            self.divide(nbins=self.nbins)
            self.integrate_bins()
            self.calculate_max_points()
            self.fit_max_points(deg=1)
            
            for y in range(self.data.shape[1]):
                self.data[:, y] = shift(self.data[:, y], -self.f_fix(y), prefilter=False, order=0, mode='wrap')


    def add_offset(self, value):
        """Shifts data to adjust any misalignments caused by fixing the distortion.

        Args:
            value (int): value (in pixels) to be shifted.

        """
        for y in range(self.data.shape[1]):
            self.data[:, y] = shift(self.data[:, y], value, prefilter=False, order=0, mode='reflect')

    
    def save_vms(self, filename):
        """Saves data as VAMAS file.

        Args:
            filename (string): name of output vms file (format extension required, i.e. .vms).
        """

        str_not_specified = "Not Specified\n"
        

        with open(filename, 'w') as f:
            
            if self.__integrated == False:
                print("Error: data not integrated yet. Use XPS_Experiment.integrate().")
                return

            # ..:: EXPERIMENT ::..

            # Format identifier
            f.write("VAMAS Surface Chemical Analysis Standard Data Transfer Format 1988 May 4\n")
            # Institution identifier
            f.write(self.header.loc['Location'].values[0] + '\n')
            # Instrument model identifier
            f.write(self.header.loc['Instrument'].values[0] + '\n')
            # Operator identifier
            f.write(self.header.loc['User'].values[0] + '\n')
            # Experiment identifier
            f.write(str_not_specified)
            # Number of lines in comment
            f.write("1\n")
            # Comment line
            f.write(self.header.loc['Comments'].values[0] + '\n')
            # Experiment mode
            f.write("NORM\n")
            # Scan mode
            f.write("REGULAR\n")
            # Number of spectral regions
            f.write(self.header.loc['Number of Regions'].values[0] + '\n')
            # Number of experimental variables
            f.write("1\n")
            f.write("Data Set\n")
            f.write("d\n")
            # Number of entries in parameter inclusion or exclusion list
            f.write("0\n")
            # Number of manually entered items in block
            f.write("0\n")
            # Number of future upgrade experiment entries
            f.write("0\n")
            # Number of future upgrade block entries
            f.write("0\n")
            # Number of blocks
            f.write("1\n")

            # ..:: BLOCK ::..

            # Experiment identifier (block)
            f.write("Block 1\n")
            # Sample identifier
            f.write(self.header.loc['Sample'].values[0] + '\n')

            date = (self.header.loc['Date'].values[0]).split('/')
            # Year
            f.write(date[2] + '\n')
            # Month
            f.write(date[0] + '\n')
            # Day
            f.write(date[1] + '\n')

            time, period = (self.header.loc['Time'].values[0]).split()
            h, m, s      = time.split(':')
            if period == 'PM':
                if h == '12':
                    h = '00'
                else:
                    h = str(int(h) + 12)
            # Hours
            f.write(h + '\n')
            # Minutes
            f.write(m + '\n')
            # Seconds
            f.write(s + '\n')

            # Number of hours in advance of GMT
            f.write("-3\n")
            # Number of comment lines in block
            f.write("0\n")
            # Technique
            f.write("XPS\n")
            # Value of experimental variable
            f.write("0\n")
            # Analysis source label ("fonte de radiacao")
            f.write("Al\n")
            # Analysis source characteristic energy ("comprimento de onda")
            f.write("1486.61\n")
            # Analysis source strength
            f.write("100\n")
            # Analysis source beam width x
            f.write("0\n")
            # Analysis source beam width y
            f.write("0\n")
            # Analysis source polar angle of incidence
            f.write("45\n")
            # Analysis source azimuth
            f.write("0\n")
            # Analyser mode
            f.write("FAT\n")
            # Analyser pass energy or retard ratio or mass resolution
            f.write(self.header.loc['Pass Energy'].values[0] + '\n')
            # Magnification of analyser transfer lens
            f.write("1\n")
            # Analyser work function or acceptance energy of atom or ion
            f.write(str(1e37) + '\n')
            # Target bias
            f.write("0\n")
            # Analysis width x
            f.write("0\n")
            # Analysis width y
            f.write("0\n")
            # Analyser axis take off polar angle
            f.write("45\n")
            # Analyser axis ,take off azimuth
            f.write("0\n")
            # Species label
            f.write(self.header.loc['Sample'].values[0] + '\n')
            
            names = str(self.header.loc['Region Name'].values[1][0])
            names = names.split('_')
            # Transition or charge state label
            f.write(names[1] + '\n')

            # Charge of detected particle
            f.write("-1\n")
            abscissa = (self.header.loc['Dimension 1 name'].values[0]).split()
            label    = ' '.join(abscissa[:2])
            units    = abscissa[-1][1:-1]
            # Abscissa label
            f.write(label + '\n')
            # Abscissa units
            f.write(units + '\n')
            
            # Abscissa start
            f.write(self.header.loc['Low Energy'].values[0] + '\n')
            # Abscissa increment
            f.write(self.header.loc['Energy Step'].values[0] + '\n')
            # Number of corresponding variables
            f.write("1\n")
            # Corresponding variable label (contagem por canal)
            f.write("counts\n")
            # Corresponding variable units (unidade das variaveis)
            f.write('d\n')
            # Signal mode
            f.write("Pulse counting\n")
            
            # Signal collection time
            t = float(self.header.loc['Step Time']) / 100
            f.write(str(t) + '\n')

            # Number of scans to compile this block
            f.write(self.header.loc['Number of Sweeps'].values[0] + '\n')
            # Signal time correction
            f.write(str(1e37) + '\n')
            # sample normal polar angle of tilt
            f.write(str(1e37) + '\n')
            # sample normal tilt azimuth
            f.write(str(1e37) + '\n')
            # Sample rotation angle
            f.write(str(1e37) + '\n')
            # Number of additional numerical parameters
            f.write("0\n")
            """
            Aqui viria o label, units e value do parametro adicional
            """
            
            # Number of ordinate values (numero de valores obtidos)
            f.write(str(self.int_data.shape[0]) + '\n')
            # Minimum ordinate value
            f.write(str(self.int_data.min().values) + '\n')
            # Maximum ordinate value
            f.write(str(self.int_data.max().values) + '\n')

            # Ordinate values (data)
            for i in range(self.int_data.shape[0]):
                f.write(str(self.int_data[i].values) + '\n')

            f.write("end of experiment\n")


def convert_directory(path=os.getcwd()):
    """Converts all `.txt` XPS experiment files in a given directory to
    `.vms`. If this function tries to process a random `.txt` file, i.e.,
    a file not containing a XPS experiment, it will ignore it.

    Args:
        path (str): path of directory containing the `.txt` files.
            Default is the current working directory.
    """

    for file_name in os.listdir(path):
        if file_name.endswith(".txt"):

            try:
                
                print("Generating vamas file for \"", file_name, "\"", sep="")

                file_name = path + "/" + file_name

                experiment = XPS_Experiment()
                experiment.read_file(file_name)
                
                experiment.integrate()

                out_name = file_name[:-4] + ".vms"
                experiment.save_vms(out_name)
                
                print("Done!")
            
            except:
                print(f"\n\nError: Invalid file: \n{file_name}\n\n Ignoring...\n\n")


    print("All files converted!")
