import numpy
import matplotlib.pyplot as plt
plt.switch_backend('agg')

import os

import clawpack.geoclaw.units as units
import clawpack.geoclaw.surge.storm as storm

from load_synthetic_storms import load_other_chaz_storms 
from load_synthetic_storms import load_chaz_storms 

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

#nc_fid_chaz_ny_nc = Dataset("../data/synthetic-storm-data/LongIsland_chaz_test.nc", 'r')
#nc_attrs, nc_dims, nc_vars = ncdump(nc_fid_chaz_ny_nc)
      
#obs = [["Obs_ibtracks_knaff15_Mumbai.nc", "Obs NC", "other_chaz", "k", 2016-1879+1]]
data = [["LongIsland_chaz_test.nc", "Obs NC", "other_chaz", "k", 2016-1879+1]]



     

data_path = os.path.join("../data", "synthetic-storm-data", data[0][0])
storms = load_chaz_storms(path = data_path,
                mask_distance = None,
                mask_coordinate = None,
                mask_category = 4,
                categorization = "NHC")

test_storm = storms[0] 

fig = plt.figure()
fig.set_figwidth(fig.get_figwidth() * 2)


title_font = {'fontname':'Arial', 'size':'12', 'color':'black', 'weight':'normal',
          'verticalalignment':'bottom'} # Bottom vertical alignment for more space
axis_font = {'fontname':'Arial', 'size':'12'}

axes = fig.add_subplot(2, 2, 1)
axes.plot(test_storm.max_wind_speed, test_storm.max_wind_radius, 'ro', label="Knaff Max Wind Radius") 
axes.set_xlabel('Max Wind Speed (knots)', **axis_font)
axes.set_ylabel('Max Wind Radius (nmi)', **axis_font)
axes.legend()


axes = fig.add_subplot(2, 2, 2)
axes.plot(units.convert(test_storm.max_wind_speed, 'knots', 'm/s'),
units.convert(test_storm.max_wind_radius, 'nmi', 'km'), 'ro-', label="Knaff MWR")  
axes.set_xlabel('Max Wind Speed (m/s)', **axis_font)
axes.set_ylabel('Max Wind Radius (km)', **axis_font)
axes.legend()


axes = fig.add_subplot(2, 2, 3)
axes.plot(range(0, test_storm.max_wind_speed.shape[0]),
test_storm.max_wind_speed, 'bo-', label="Max Wind Speed")
axes.set_xlabel('Obs Count', **axis_font)
axes.set_ylabel('Max Wind Speed (knots)', **axis_font)
axes.legend()


axes = fig.add_subplot(2, 2, 4)
axes.plot(range(0, test_storm.max_wind_radius.shape[0]),
test_storm.max_wind_radius, 'ro-', label="Knaff MWR")
axes.set_xlabel('Obs Count', **axis_font)
axes.set_ylabel('Max Wind Radius (nmi)', **axis_font)
axes.legend()

plt.yticks(fontsize=10)
plt.xticks(fontsize=10)
plt.tight_layout()
plt.savefig('MaxWindRadiusCheck.pdf')
