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


def calculate_surge(surge_data_descriptors = None, mask_distance = None,
                    mask_coordinate = None, mask_category = None):
    
    r'''
    Calculate surge heights. 

    '''
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharey=False, sharex=False) 
    fig.set_size_inches(20.0, 30.0) 
    bar_width = 1.0 
    opacity = 1.0 
    
    wind_models = ["holland80", "holland10", "CLE", "SLOSH",
                   "rankine", "modified-rankine", "DeMaria"] 

    amr_levels = ["amr2", "amr5"]
 
    if surge_data_descriptors is not None:
        for (index, data) in enumerate(surge_data_descriptors):
            
            if any(model in data[0] for model in wind_models): 
                label = model 
            gauge_data = numpy.loadtxt(data[0], delimiter = ',', 
                                      skiprows = 1)
            max_surge = gauge_data[:, 2]
            mu = numpy.mean(max_surge) 
            var = numpy.var(max_surge) 
    
            events_surge = numpy.array(max_surge)

            bins_surge = numpy.linspace(0.2, 6.0, 10) 
            period = data[-1] 
            
            hist, bin_edges = numpy.histogram(events_surge, bins_surge) 
            index_edges = numpy.ones(bin_edges.shape) * (index + bar_width) 
            n = hist.sum() 
            
            # Complement of empirical distribution function
            ECDF_c = numpy.cumsum(hist[::-1])[::-1] * 1/n
            ECDF = numpy.ones(ECDF_c.shape, dtype = float) - ECDF_c

            return_period = period * 1/n * (1/ECDF_c)


            T_r = numpy.zeros(events_surge.shape, dtype=float)

            events_surge = numpy.sort(events_surge)

            counter = 0

            for i in range(events_surge.shape[1]):
                if events_surge[0,i] < bin_edges[counter]:
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


    mumbai_intensity_data_descriptors = [
        ["Obs_ibtracks_knaff15_Mumbai.nc", "Obs NC", "other_chaz", "k", 2016-1879+1],
        ["Emanuel_Mumbai.nc", "Emanuel NC", "other_chaz", "b", 1850/0.061545], 
        ["CHAZ_knaff15_Mumbai.nc", "CHAZ NC", "other_chaz", "r", 32.0*120*27]]

    calculate_intensity(intensity_data_descriptors = mumbai_intensity_data_descriptors,
                    mask_distance = None,
                    mask_coordinate = (-74.0060, 40.7128),
                    mask_category = None)
