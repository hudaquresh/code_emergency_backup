import numpy 
import os
import sys
import matplotlib.pyplot as plt 

def calreturnP(bins=[0,1,2,3],period=0.25,events=[0,1,2,3]):
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


    Variables 
    -------
    T      : numpy array 
        Period  
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
    #events = events[events==events]
    #events = [1, 2, 1] 
    #bins = [0, 1, 2, 3]
    #print('bins:', bins)  
    a = numpy.histogram(events, bins)[0]
    print(numpy.histogram(events, bins))
    a = list(a) 
    new_a = []
    new_bins = []  
    #new_a = a
    #new_bins = bins  
    for i in range(0, len(a)): 
        if a[i] == 0:
            continue 
        else: 
            new_a.append(a[i])
            new_bins.append(bins[i]) 
    a = numpy.array(new_a) 
    #print('a:', a) 
    #returnP = period*1./a.sum()*(1./(numpy.cumsum(a[::-1])[::-1]*1.0/a.sum()))
    T = 0.25
    n = a.sum()
    ECDF_c = numpy.cumsum(a[::-1])[::-1] * 1/n
    T_r = T * 1/n * (1/ECDF_c) 
    returnP = period*1./a.sum()*(1./(numpy.cumsum(a[::-1])[::-1]*1.0/a.sum()))
    #print(T_r)
    #print(returnP) 
    #print("") 
    #print("a: ", a)
    #print("a[::-1]: ", a[::-1]) 
    #print("")
    #print("numpy.cumsum(a[::-1]):", numpy.cumsum(a[::-1]))
    #print("") 
    #print("numpy.cumsum(a[::-1])[::-1]: ", numpy.cumsum(a[::-1])[::-1]) 
    #print("")
    #print("n:", a.sum())
    #print("") 
    #print("")
    #print("numpy.cumsum(a[::-1]):", numpy.cumsum(a[::-1]))
    #print("") 
    #print("numpy.cumsum(a): ", numpy.cumsum(a)) 
    ##print(len(returnP))
    ##return a, returnP, new_bins 
    return a, T_r, new_bins 

def filtergauges(gaugedata):
    filtered_gauges = [] 
    for n in range(0, len(gauge_data[:, 2])):
        if (gauge_data[n, 2] > 0.0000 and gauge_data[n,2] < 5.0): 
            filtered_gauges.append(gauge_data[n,2])
        else: 
            #print(gauge_data[n,2]) 
            continue  

    #print('filtered_gauges:', filtered_gauges)
    filtered_gauges = numpy.array(filtered_gauges) 
    return filtered_gauges 


if __name__ == '__main__': 
    
    # Initialize the tital and axis font sizes, color, weight, and style 
    title_font = {'fontname':'Arial', 'size':'28', 'color':'black', 'weight':'normal',
              'verticalalignment':'bottom'} # Bottom vertical alignment for more space
    axis_font = {'fontname':'Arial', 'size':'20'}

    calreturnP()

    # Run through set gauges and wind models and initialize bins for return
    # period plots. The set of gauges are a simple list as well as the wind
    # models. Gauges include 1, 2, 3, 4. While the wind models are the following
    # codes: 
    # H80: Holland 1980
    # H10: Holland 2010 
    # CLE: Chavas Lin Emanuel 
    # SLOSH: Slosh   
    gauges = [1, 2, 3, 4]
    #gauges = [1, 2, 3]
    #gauges = [1]
    #wind_models = ['H80', 'H10']

    #gauges = [1]
    #wind_models = ["holland80", "holland10", "CLE", 
    #               "SLOSH", "rankine", "modified-rankine", 
    #               "DeMaria"]

    #wind_models = ["holland80"] 
    wind_models = ["rankine"] 
    #wind_models = ["holland80", "rankine"] 
    #wind_models = ["holland80", "modified-rankine"] 

    # The set of bins necessary for plotting
    bins = numpy.linspace(0, 6, 30)
    #print('a:', a) 
    #returnP = period*1./a.sum()*(1./(numpy.cumsum(a[::-1])[::-1]*1.0/a.sum()))
    print('bins:', bins)
    #bins = [0, 1, 2, 3, 4, 5, 6]
    period = 1850/0.061545
    #period = 32*120*27 
   
    # Set line styles for plots  
    line_style = ['-o','--', '-.',':', '-o', 's', 'p', '*']

    # Dictionary to contain all wind models for a particular gauge  
    gauge_comparisons = {} 
    
    # Loop to collect all wind model data for a particular gauge and group it by
    # wind model into a list of paths for that particular data. A print
    # statement is then used to print the paths so we may be sure all wind
    # models paths were collected and organized to the appropriate gauge.  
    #amr_level = 5 
    amr_level = 2 
    for g in gauges: 
        gauge_comparisons[g] = [] 
        #for wm in wind_models: 
        #    #gaugefile = "CHAZ-%s-AMR%i-Max-Surges-Gauge-%i.txt" %(wm, amr_level, g)
        #    gaugefile = "%s-AMR%i-Max-Surges-Gauge-%i.txt" %(wm, amr_level, g)
        #    gauge_path = os.path.join(os.getcwd(),
        #                                "chaz",  
        #                                gaugefile)
        #    gauge_comparisons[g].append(gauge_path)
        #    print("Gauge %i, Wind Model %s" %(g, wm), gauge_path) 
        
        for wm in wind_models: 
            #gaugefile = "CHAZ-%s-AMR%i-Max-Surges-Gauge-%i.txt" %(wm, amr_level, g)
            gaugefile = "%s-AMR%i-Max-Surges-Gauge-%i.txt" %(wm.upper(), amr_level, g)
            gauge_path = os.path.join(os.getcwd(),
                                        "kerry-surges",  
                                        gaugefile)
            gauge_comparisons[g].append(gauge_path)
            print("Gauge %i, Wind Model %s" %(g, wm), gauge_path) 
   
    # Create a return period curve for each gauge with one curve for each wind
    # model.  
    for g in gauges:
    
        # Initialize plot dimensions  
        fig = plt.figure()
        fig.set_size_inches(15.5, 7.5)
        axes = fig.add_subplot(1, 1, 1)
        
        # Collect data paths for a particular gauge 
        data_paths = gauge_comparisons[g]
        j = 0
        
        # Run through each wind model 
        for i in range(len(data_paths)):
            if int(data_paths[i][-5]) == g:
                gauge_data = numpy.loadtxt(data_paths[i], delimiter=',', skiprows=1)
                max_surges = gauge_data[:, 2]
                max_surges = filtergauges(gauge_data) 
                mu = numpy.mean(max_surges)
                var = numpy.var(max_surges)
                print(mu, var)  
                events = max_surges 
                #print('max_surges:', max_surges) 
                a, prob, new_bins = calreturnP(bins=bins,period=period,events=events)
                #print('Prob:', prob)
                #new_bins = [0, 1, 2, 3, 4, 5]
                #new_bins = bins[0:-1]
                print(new_bins)
                axes.semilogx(prob, new_bins, line_style[i],linewidth=3, label="%s Wind Model" %wind_models[i]) 
    
        axes.set_xlim(0)
        #axes.set_ylim([0, 0.1])
        axes.set_ylim(0)

 
        #title = 'CHAZ_Return_Period_Gauge-%i' %g
        title = 'KERRY_AMR%i_Return_Period_Gauge-%i' %(amr_level,g)
        plot_title = title.replace("_", " ")  
        axes.set_title("%s" %plot_title,**title_font)
        axes.set_xlabel("Return Period (years)",**axis_font)
        axes.set_ylabel("Surge Height (m) ",**axis_font)
        axes.legend(prop={'size': 16})
        
        
        axes.grid(True, linestyle='-',color='0.75', linewidth=1.5)

        plt.text(1, -0.6, r'$\sum_{i=0}^\infty x_i$', fontsize=20)
        plt.yticks(fontsize=24)
        plt.xticks(fontsize=24)
        
        plt.show() 

        return_period_plots_dir = os.path.join(os.getcwd(), "../return-period-plots")
        fig_title = os.path.join(return_period_plots_dir, "%s.png" %title) 
        fig.savefig(fig_title,format='png')
        #fig.savefig("return-period-plots/%s.png" %title,format='png')
