import numpy
import matplotlib.pyplot as plt

import clawpack.geoclaw.units as units
import clawpack.geoclaw.surge.storm as storm

from netCDF4 import Dataset
# from ncdump import ncdump


def calculate_return_period(period, events, bins):
    r"""
    Calculate empirical return period using the following
    formulation

    T_r = T/(n * ECDF_c) = dt/(ECDF)

    http://abouthydrology.blogspot.com/2017/10/return-period_25.html

    Parameters
    ----------
    bins    : numpy array
        The bins specified
    events  : numpy array
        The events we are attempting to bin
    period  : float
        Period in which data is sampled
    n      : int
        Number of data points
    ECDF_c : numpy array
        Complement of the empirical cumulative distribution function


    Returns
    -------
    T_r    : numpy array
        Empirical return period
    a      : numpy array
        Values of the histogram
    bins   : numpy array
        Edges of the bins

    """
    hist, bin_edges = numpy.histogram(events, bins)
    n = hist.sum()

    # Complement of empirical distribution function
    #print(n)
    ECDF_c = numpy.cumsum(hist[::-1])[::-1] * 1/n
    #print(numpy.cumsum(hist[::-1])[::-1])
    #print(ECDF_c)

    return_period = period * 1/n * (1/ECDF_c)

    return hist, bin_edges, return_period


def intensity(storms, period):
    r"""
    Calculate intensity return period.

    Parameters
    ----------
    storms           : obj
        Storms over which intensity is calculated
    period           : float
        Period of the storm data
    events_intensity : numpy array
        Array containing max storm speeds
    bins_intensity   : numpy array
        Array with the binning needed for intensity


    Returns
    -------
    hist          : numpy array
        Historgram of data
    bin_edges     : numpy array
        The edges of the bins of the histogram
    return_period : numpy array
        Return period
    """


    events_intensity = numpy.ones((1,len(storms)), dtype = float)

    for (i, storm) in enumerate(storms):
        storm.max_wind_speed = units.convert(storm.max_wind_speed, 'm/s', 'knots')
        events_intensity[0,i] = numpy.max(storm.max_wind_speed)


    bins_intensity = numpy.array([0, 30, 64, 83, 96, 113, 135, 200])
    hist, bin_edges, return_period = calculate_return_period(period = period,
                                            events = events_intensity,
                                            bins = bins_intensity)
    return hist, bin_edges, return_period

def surge_return_period(storms, period):
    r"""
    Calculate storm surge return period.

    Parameters
    ----------
    storms           : obj
        Storms over which intensity is calculated
    period           : float
        Period of the storm data
    events_surge : numpy array
        Array containing max surge heights
    bins_surge   : numpy array
        Array with the binning needed for surge


    Returns
    -------
    hist          : numpy array
        Historgram of data
    bin_edges     : numpy array
        The edges of the bins of the histogram
    return_period : numpy array
        Return period
    """
    pass

def mumbai_intensity(mask_dist = None,
                     mask_coord = None,
                     mask_cat = None):
    r"""
    Plot intensities of Mumbai
    """

    # Calculate Emanuel's storm intensities
    storms_emanuel_old_mat = storm.load_emanuel_storms(path = "data/Mumbai_IO_ncep_reanal.mat",
                                   mask_distance = mask_dist,
                                   mask_coordinate= mask_coord,
                                   mask_category = mask_cat, categorization="NHC")

    period_emanuel_old_mat =1750/0.28402

    h_emanuel_old_mat, b_e_emanuel_old_mat, T_r_emanuel_old_mat = intensity(storms=storms_emanuel_old_mat, period=period_emanuel_old_mat)
    # Calculate Emanuel's storm intensities
    storms_emanuel = storm.load_emanuel_storms(path = "data/Mumbai3_io_ncep_reanalcal.mat",
                                   mask_distance = mask_dist,
                                   mask_coordinate= mask_coord,
                                   mask_category = mask_cat, categorization="NHC")

    period_emanuel = 1850/0.061545

    h_emanuel, b_e_emanuel, T_r_emanuel = intensity(storms=storms_emanuel,
                                                            period=period_emanuel)

    # Calculate Emanuel NetCDF storm intensities
    storms_emanuel_nc = storm.load_other_chaz_storms(path = "data/Emanuel_Mumbai.nc",
                                   mask_distance = mask_dist,
                                   mask_coordinate = mask_coord,
                                   mask_category = mask_cat, categorization = "NHC")


    h_emanuel_nc, b_e_emanuel_nc, T_r_emanuel_nc = intensity(storms=storms_emanuel_nc,
                                                            period=period_emanuel)


    # Calculate CHAZ's storm intensities
    storms_chaz = storm.load_other_chaz_storms(path = "data/CHAZ_knaff15_Mumbai.nc",
                                   mask_distance = mask_dist,
                                   mask_coordinate = mask_coord,
                                   mask_category = mask_cat, categorization = "NHC")


    period_chaz = 32.0*120*27

    h_chaz, b_e_chaz, T_r_chaz = intensity(storms=storms_chaz,
                                                            period=period_chaz)

    # Calculate Observations's storm intensities
    storms_obs = storm.load_other_chaz_storms(path = "data/Obs_ibtracks_knaff15_Mumbai.nc",
                                   mask_distance = mask_dist,
                                   mask_coordinate= mask_coord,
                                   mask_category = mask_cat, categorization="NHC")


    period_obs = 2016-1879+1

    h_obs, b_e_obs, T_r_obs = intensity(storms=storms_obs, period=period_obs)
    # Plot intensities on same plot
    # Initalize font sizes
    title_font = {'fontname':'Arial', 'size':'28', 'color':'black', 'weight':'normal',
              'verticalalignment':'bottom'} # Bottom vertical alignment for more space
    axis_font = {'fontname':'Arial', 'size':'20'}

    # Initialize figure
    fig = plt.figure()
    fig.set_size_inches(15.5, 7.5)
    axes = fig.add_subplot(1, 1, 1)

    axes.semilogx(T_r_emanuel_old_mat, b_e_emanuel_old_mat[:-1],
                                        label="Emanuel Old Mat")

    axes.semilogx(T_r_emanuel, b_e_emanuel[:-1],
                                        label="Emanuel Mat3")
    axes.semilogx(T_r_emanuel_nc, b_e_emanuel_nc[:-1],
                                        label="Emanuel NetCDF")
    axes.semilogx(T_r_chaz, b_e_chaz[:-1],
                                        label="CHAZ")
    axes.semilogx(T_r_obs, b_e_obs[:-1],
                                        label="OBS")

    axes.set_ylim(0)

    plot_title = "Return Period Mumbai"
    axes.set_title("%s" %plot_title,**title_font)
    axes.set_xlabel("Return Period (years)",**axis_font)
    axes.set_ylabel("Intensity (knots) ",**axis_font)
    axes.legend(prop={'size': 16})


    axes.grid(True, linestyle='-',color='0.75', linewidth=1.5)

    plt.yticks(fontsize=24)
    plt.xticks(fontsize=24)

    plt.show()

    return fig

def new_york_intensity(mask_dist = None, mask_coord = None,
                         mask_cat = None):
    r"""
    Plot intensities of New York
    """

    # Calculate CHAZ's storm intensities
    storms_chaz = storm.load_chaz_storms(path = "data/LongIsland_chaz_test.nc",
                                   mask_distance = mask_dist,
                                   mask_coordinate = mask_coord,
                                   mask_category = mask_cat,
                                   categorization="NHC")


    period_chaz = 32.0*120*27

    h_chaz, b_e_chaz, T_r_chaz = intensity(storms=storms_chaz,
                                           period=period_chaz)
    # Plot intensities on same plot
    # Initalize font sizes
    title_font = {'fontname':'Arial', 'size':'28', 'color':'black', 'weight':'normal',
              'verticalalignment':'bottom'} # Bottom vertical alignment for more space
    axis_font = {'fontname':'Arial', 'size':'20'}

    # Initialize figure
    fig = plt.figure()
    fig.set_size_inches(15.5, 7.5)
    axes = fig.add_subplot(1, 1, 1)

    axes.semilogx(T_r_chaz, b_e_chaz[:-1], label="CHAZ")

    axes.set_ylim(0)

    plot_title = "Return Period New York"
    axes.set_title("%s" %plot_title,**title_font)
    axes.set_xlabel("Return Period (years)",**axis_font)
    axes.set_ylabel("Intensity (knots) ",**axis_font)
    axes.legend(prop={'size': 16})


    axes.grid(True, linestyle='-',color='0.75', linewidth=1.5)

    plt.yticks(fontsize=24)
    plt.xticks(fontsize=24)

    plt.show()

    return fig

def ncdump(nc_fid, verb=True):
    '''
    ncdump outputs dimensions, variables and their attribute information.
    The information is similar to that of NCAR's ncdump utility.
    ncdump requires a valid instance of Dataset.

    Parameters
    ----------
    nc_fid : netCDF4.Dataset
        A netCDF4 dateset object
    verb : Boolean
        whether or not nc_attrs, nc_dims, and nc_vars are printed

    Returns
    -------
    nc_attrs : list
        A Python list of the NetCDF file global attributes
    nc_dims : list
        A Python list of the NetCDF file dimensions
    nc_vars : list
        A Python list of the NetCDF file variables
    '''
    def print_ncattr(key):
        """
        Prints the NetCDF file attributes for a given key

        Parameters
        ----------
        key : unicode
            a valid netCDF4.Dataset.variables key
        """
        try:
            print("\t\ttype:", repr(nc_fid.variables[key].dtype))
            for ncattr in nc_fid.variables[key].ncattrs():
                print('\t\t%s:' % ncattr,\
                      repr(nc_fid.variables[key].getncattr(ncattr)))
        except KeyError:
            print("\t\tWARNING: %s does not contain variable attributes" % key)

    # NetCDF global attributes
    nc_attrs = nc_fid.ncattrs()
    if verb:
        print("NetCDF Global Attributes:")
        for nc_attr in nc_attrs:
            print('\t%s:' % nc_attr, repr(nc_fid.getncattr(nc_attr)))
    nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
    # Dimension shape information.
    if verb:
        print("NetCDF dimension information:")
        for dim in nc_dims:
            print("\tName:", dim)
            print("\t\tsize:", len(nc_fid.dimensions[dim]))
            print_ncattr(dim)
    # Variable information.
    nc_vars = [var for var in nc_fid.variables]  # list of nc variables
    if verb:
        print("NetCDF variable information:")
        for var in nc_vars:
            if var not in nc_dims:
                print('\tName:', var)
                print("\t\tdimensions:", nc_fid.variables[var].dimensions)
                print("\t\tsize:", nc_fid.variables[var].size)
                print_ncattr(var)
    return nc_attrs, nc_dims, nc_vars




if __name__ == '__main__':

    from matplotlib.backends.backend_pdf import PdfPages


    #---------------Mumbai------------

    # Create one pdf file to save all intensities
    mumbai_intensities = PdfPages('output/mumbai_intensities.pdf')

    # Generating intensity plots for Mumbai
    mumbai_fig = mumbai_intensity(mask_dist = None, mask_coord = (72.8562, 19.0176),
                             mask_cat = None)
    mumbai_fig.savefig(mumbai_intensities, format='pdf')

    mumbai_fig = mumbai_intensity(mask_dist = 0.4, mask_coord = (72.8562, 19.0176),
                             mask_cat = 3)
    mumbai_fig.savefig(mumbai_intensities, format='pdf')

    # Save intensities to a pdf file
    mumbai_intensities.close()


    #---------------NYC--------------

    # Create one pdf file to save all intensities
    nyc_intensities = PdfPages('output/nyc_intensities.pdf')

    # Generating intensity plots for NYC
    nyc_fig = new_york_intensity(mask_dist = None, mask_coord = (-74.0060, 40.7128),
                             mask_cat = None)
    nyc_fig.savefig(nyc_intensities, format='pdf')

    nyc_fig = new_york_intensity(mask_dist = 0.4, mask_coord = (-74.0060, 40.7128),
                             mask_cat = 3)
    nyc_fig.savefig(nyc_intensities, format='pdf')

    # Save intensities to a pdf file
    nyc_intensities.close()

    # Show plots
    mumbai_fig.show()
    new_york_fig.show()


    # nc_fid = Dataset("data/CHAZ_knaff15_Mumbai.nc", 'r')
    # nc_fid_eman_nc = Dataset("data/Emanuel_Mumbai.nc", 'r')
    # nc_fid_chaz_ny_nc = Dataset("data/LongIsland_chaz_test.nc", 'r')
    #nc_attrs, nc_dims, nc_vars = ncdump(nc_fid)
    #nc_attrs, nc_dims, nc_vars = ncdump(nc_fid_eman_nc)
    # nc_attrs, nc_dims, nc_vars = ncdump(nc_fid_chaz_ny_nc)
