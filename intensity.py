import numpy
import matplotlib.pyplot as plt

import os

import clawpack.geoclaw.units as units
import clawpack.geoclaw.surge.storm as storm

from netCDF4 import Dataset


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

def calculate_intensity(intensity_data_descriptors = None, mask_distance = None,
                        mask_coordinate = None, mask_category = None):
    r'''

    Parameters
    ----------
    intensity_data_descriptors : list
        A list containing the path, label, and type of file


    Returns
    -------
    '''

    fig, (ax1, ax2) = plt.subplots(nrows = 2, ncols = 1, sharey=False, sharex=False)
    fig.set_size_inches(20.0, 30.0)

    bar_width = 1.0
    opacity = 1.0

    if intensity_data_descriptors is not None:
        for (index, data) in enumerate(intensity_data_descriptors):
            if data[2] == "other_chaz":
                data_path = os.path.join("data", data[0])
                storms = storm.load_other_chaz_storms(path = data_path,
                                mask_distance = mask_distance,
                                mask_coordinate = mask_coordinate,
                                mask_category = mask_category,
                                categorization = "NHC")
            elif data[2] == "mat":
                data_path = os.path.join("data", data[0])
                storms = storm.load_emanuel_storms(path = data_path,
                                mask_distance = mask_distance,
                                mask_coordinate = mask_coordinate,
                                mask_category = mask_category,
                                categorization = "NHC")
            else:
                break

            events_intensity = numpy.ones((1,len(storms)), dtype = float)
            for (i, storm_track) in enumerate(storms):
                storm_track.max_wind_speed = units.convert(storm_track.max_wind_speed, 'm/s', 'knots')
                events_intensity[0,i] = numpy.max(storm_track.max_wind_speed)


            #bins_intensity = numpy.array([1.0, 30.0, 64.0, 83.0, 96.0, 113.0, 135.0, 160.0, 180.0])
            bins_intensity = numpy.linspace(1.0, 200, 100)

            period = data[-1]

            hist, bin_edges = numpy.histogram(events_intensity, bins_intensity)
            index_edges = numpy.ones(bin_edges.shape) * (index + bar_width)
            n = hist.sum()

            # Complement of empirical distribution function
            ECDF_c = numpy.cumsum(hist[::-1])[::-1] * 1/n
            ECDF = numpy.ones(ECDF_c.shape, dtype = float) - ECDF_c

            return_period = period * 1/n * (1/ECDF_c)


            T_r = numpy.zeros(events_intensity.shape, dtype=float)

            events_intensity = numpy.sort(events_intensity)

            counter = 0

            for i in range(events_intensity.shape[1]):
                if events_intensity[0,i] < bin_edges[counter]:
                    T_r[0,i] = return_period[counter]
                else:
                    counter += 1
                    T_r[0,i] = return_period[counter]

            ax1.bar(index_edges[:-1] + bin_edges[:-1], ECDF_c, bar_width,
                            label = data[1], color = data[-2],
                            alpha = opacity)

            ax2.semilogx(T_r[0, :], events_intensity[0, :],  
                    label = data[1], color = data[-2])


    title_font = {'fontname':'Arial', 'size':'10', 'color':'black', 'weight':'normal',
              'verticalalignment':'bottom'} # Bottom vertical alignment for more space
    axis_font = {'fontname':'Arial', 'size':'8'}

    ax1.set_xlabel('Knots', **axis_font)
    ax1.set_ylabel('CDF', **axis_font)
    ax1.set_title('Mumbai 150 km', **title_font)
    
    ax1.legend()

    ax2.set_xlabel('Return Period (Years)', **axis_font)
    ax2.set_ylabel('Intensity (knots)', **axis_font)
    ax2.set_title('Mumbai 150 km', **title_font)
    ax2.set_xlim(0,5000) 
    ax2.set_ylim(20,150) 
    ax2.legend()


    plt.yticks(fontsize=10)
    plt.xticks(fontsize=10)
    plt.tight_layout()
    plt.savefig('output/Mumbai_Stats_ncdata.pdf')
    plt.show()

    return fig


if __name__ == '__main__':

    # nc_fid_chaz_ny_nc = Dataset("data/LongIsland_chaz_test.nc", 'r')
    # nc_attrs, nc_dims, nc_vars = ncdump(nc_fid_chaz_ny_nc)

    mumbai_intensity_data_descriptors = [
        ["Obs_ibtracks_knaff15_Mumbai.nc", "Obs NC", "other_chaz", "k", 2016-1879+1],
        ["Emanuel_Mumbai.nc", "Emanuel NC", "other_chaz", "b", 1850/0.061545],
        ["CHAZ_knaff15_Mumbai.nc", "CHAZ NC", "other_chaz", "r", 32.0*120*27],
        ["Mumbai_IO_ncep_reanal.mat", "Emanuel Old Mat", "mat", "c", 1850/0.061545],
        ["Mumbai3_io_ncep_reanalcal.mat", "Emanuel Mat", "mat", "m", 1850/0.061545]]

    mumbai_intensity_data_descriptors = [
        ["Obs_ibtracks_knaff15_Mumbai.nc", "Obs NC", "other_chaz", "k", 2016-1879+1]]

    calculate_intensity(intensity_data_descriptors = mumbai_intensity_data_descriptors,
                    mask_distance = None,
                    mask_coordinate = (-74.0060, 40.7128),
                    mask_category = None)
