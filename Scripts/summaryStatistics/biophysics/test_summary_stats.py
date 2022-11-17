import CellModeller
import anisotropy
import aspect_ratio
import density_calculation
import growth_rate_analysis
import pickle

# Import data
pickle_full_path = 'step-00420.pickle'
data = pickle.load(open(pickle_full_path, 'rb'))
cells = data['cellStates']   
   
# Run calculations
order_parameter = anisotropy.get_global_order_parameter(cells)
print("global_order_parameter = ", order_parameter)

aspect_ratio = aspect_ratio.get_aspect_ratio(cells)
print("aspect_ratio = ", aspect_ratio)

density_parameter = density_calculation.get_density_parameter(cells)
print("density_parameter = ", density_parameter)

dt = 0.02
growth_rate_parameter = growth_rate_analysis.get_growth_parameter(cells, dt, show_plots=False)
print("growth_rate_parameter = ", growth_rate_parameter)
