import numpy
import matplotlib.pyplot as plt
import clawpack.geoclaw.surge.storm as storm


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
    print(numpy.cumsum(hist[::-1])[::-1])
    print(ECDF_c)

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
        events_intensity[0,i] = numpy.max(storm.max_wind_speed)


    bins_intensity = numpy.linspace(0,180, 15)
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

def mumbai_intensity():
    r"""
    Plot intensities of Mumbai
    """
    # Calculate Kerry's storm intensities
    storms_kerry = storm.load_emanuel_storms(path = "Mumbai3_io_ncep_reanalcal.mat",
                                   mask_distance = None,
                                   mask_coordinate=(72.8562, 19.0176),
                                   mask_category = None, categorization="NHC")

    period_kerry = 1850/0.061545

    hist_kerry, bin_edges_kerry, return_period_kerry = intensity(storms=storms_kerry,
                                                            period=period_kerry)


    # Calculate CHAZ's storm intensities
    storms_chaz = storm.load_other_chaz_storms(path = "CHAZ_knaff15_Mumbai.nc",
                                   mask_distance = None,
                                   mask_coordinate = (72.8562, 19.0176),
                                   mask_category = None, categorization = "NHC")


    period_chaz = 32.0*120*27

    hist_chaz, bin_edges_chaz, return_period_chaz = intensity(storms=storms_chaz,
                                                            period=period_chaz)

    # Calculate Observations's storm intensities
    storms_obs = storm.load_other_chaz_storms(path = "Obs_ibtracks_knaff15_Mumbai.nc",
                                   mask_distance = None,
                                   mask_coordinate=(72.8562, 19.0176),
                                   mask_category = None, categorization="NHC")


    period_obs = 2016-1879+1

    hist_obs, bin_edges_obs, return_period_obs = intensity(storms=storms_obs,
                                                            period=period_obs)
    # Plot intensities on same plot
    # Initalize font sizes
    title_font = {'fontname':'Arial', 'size':'28', 'color':'black', 'weight':'normal',
              'verticalalignment':'bottom'} # Bottom vertical alignment for more space
    axis_font = {'fontname':'Arial', 'size':'20'}

    # Initialize figure
    fig = plt.figure()
    fig.set_size_inches(15.5, 7.5)
    axes = fig.add_subplot(1, 1, 1)

    axes.semilogx(return_period_kerry, bin_edges_kerry[:-1], label="Kerry")
    axes.semilogx(return_period_chaz, bin_edges_chaz[:-1], label="CHAZ")
    axes.semilogx(return_period_obs, bin_edges_obs[:-1], label="OBS")

    axes.set_ylim(0)

    plot_title = "Return Period Mumbai"
    axes.set_title("%s" %plot_title,**title_font)
    axes.set_xlabel("Return Period (years)",**axis_font)
    axes.set_ylabel("Intensity (m/s) ",**axis_font)
    axes.legend(prop={'size': 16})


    axes.grid(True, linestyle='-',color='0.75', linewidth=1.5)

    plt.text(1, -0.6, r'$\sum_{i=0}^\infty x_i$', fontsize=20)
    plt.yticks(fontsize=24)
    plt.xticks(fontsize=24)

    plt.show()

    return fig

def new_york_intensity():
    r"""
    Plot intensities of New York
    """

    # Calculate CHAZ's storm intensities
    storms_chaz = storm.load_chaz_storms(path = "LongIsland_chaz_test.nc",
                                   mask_distance = None,
                                   mask_coordinate = (-74.0060, 40.7128),
                                   mask_category = None,
                                   categorization="NHC")


    period_chaz = 32.0*120*27

    hist_chaz, bin_edges_chaz, return_period_chaz = intensity(storms=storms_chaz,
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

    axes.semilogx(return_period_chaz, bin_edges_chaz[:-1], label="CHAZ")

    axes.set_ylim(0)

    plot_title = "Return Period New York"
    axes.set_title("%s" %plot_title,**title_font)
    axes.set_xlabel("Return Period (years)",**axis_font)
    axes.set_ylabel("Intensity (m/s) ",**axis_font)
    axes.legend(prop={'size': 16})


    axes.grid(True, linestyle='-',color='0.75', linewidth=1.5)

    plt.text(1, -0.6, r'$\sum_{i=0}^\infty x_i$', fontsize=20)
    plt.yticks(fontsize=24)
    plt.xticks(fontsize=24)

    plt.show()

    return fig




if __name__ == '__main__':
    mumbai_fig = mumbai_intensity()
    new_york_fig = new_york_intensity()
    mumbai_fig.show()
    new_york_fig.show()
