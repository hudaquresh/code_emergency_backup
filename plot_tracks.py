
data = [["LongIsland_chaz_test.nc", "Obs NC", "other_chaz", "k", 2016-1879+1]]



     

data_path = os.path.join("../data", "synthetic-storm-data", data[0][0])
storms = load_chaz_storms(path = data_path,
                mask_distance = None,
                mask_coordinate = None,
                mask_category = None,
                categorization = "NHC")



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
