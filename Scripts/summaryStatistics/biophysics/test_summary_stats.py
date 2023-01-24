import CellModeller
import anisotropy
import aspect_ratio
import density_calculation
import growth_rate_analysis
import growth_rate_exp_deviation
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

dt = 0.025
test_directory1 = 'exp_deviation_test_data/gamma=10'
test_directory2 = 'exp_deviation_test_data/gamma=100'
exp_deviation1 = growth_rate_exp_deviation.get_exp_deviation(test_directory1, dt)
exp_deviation2 = growth_rate_exp_deviation.get_exp_deviation(test_directory2, dt)
print(f"When gamma=10, the std of the residuals from exp fit is {exp_deviation1: .2f}.\n\
When gamma=100, the std of the residuals from exp fit is {exp_deviation2: .2f}")
