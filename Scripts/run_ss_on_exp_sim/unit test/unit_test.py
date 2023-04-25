import unittest

from Scripts.run_ss_on_exp_sim.scripts.helper_functions.helperFunctions import load_cellStates

# scripts that we are testing
from Scripts.run_ss_on_exp_sim.scripts.summary_statistics.AspectRatio import calc_aspect_ratio
from Scripts.run_ss_on_exp_sim.scripts.summary_statistics.Anisotropy import get_global_order_parameter
from Scripts.run_ss_on_exp_sim.scripts.summary_statistics.density_calculation import get_density_parameter
from Scripts.run_ss_on_exp_sim.scripts.summary_statistics.growth_rate_exp_deviation import get_exp_deviation
from Scripts.run_ss_on_exp_sim.scripts.summary_statistics.CellOrientationOnBoundary import calc_cell_orientation_on_boundary
from Scripts.run_ss_on_exp_sim.scripts.summary_statistics.AgeDistanceDistribution import calcAgeDistanceDistribution
from Scripts.run_ss_on_exp_sim.scripts.summary_statistics.dyadStructure import calcDyadStructure
from Scripts.run_ss_on_exp_sim.scripts.summary_statistics.fourier_descriptor import calc_fourier_descriptor
from Scripts.run_ss_on_exp_sim.scripts.summary_statistics.convexity_smart import cal_convexity


"""
new unit test using unittest module, with the pass condition that it doesn't crash
"""

picklefile = "sep8_step-000089.pickle"  # "sep7_step-000097" "step-00200.pickle" #"circle_step-01000.pickle" 'jan3_step-000097.pickle' "sep8_step-000089.pickle"
cs_last = load_cellStates("", picklefile)
dyad_picklefile = 'sep7_step-000020.pickle'
cs_dyad = load_cellStates("", dyad_picklefile)

# constants
um_pixel_ratio = 0.144
dt = 0.05  # 3 min, 0.05=3/60
shape = 0.048
margin=10
fig_path = ""
pickle_files_directory = 'sep_scene8'
scene = 'scene8'

class Test_sum_stats(unittest.TestCase):


    def test_aspect_ratio(self):
        self.assertIsNotNone(calc_aspect_ratio(cs_last))

    def test_anisotropy(self):
        self.assertIsNotNone(get_global_order_parameter(cs_last))

    def test_density(self):
        self.assertIsNotNone(get_density_parameter(cs_last, um_pixel_ratio))

    def test_growth_rate_exp_deviation(self):
        self.assertIsNotNone(get_exp_deviation(pickle_files_directory, dt))

    def test_cell_orientaion_on_boundary(self):
        self.assertIsNotNone(calc_cell_orientation_on_boundary(cs_last, fig_path, scene, False))

    def test_AgeDistanceDistribution(self):
        self.assertIsNotNone(calcAgeDistanceDistribution(cs_last, fig_path, scene, False))

    def test_dyadStructure(self):
        self.assertIsNotNone(calcDyadStructure(cs_dyad))

    def test_fourier_descriptor(self):
        self.assertIsNotNone(calc_fourier_descriptor(cs_last, fig_name = scene, shape=shape, margin=margin))

    def test_convexity(self):
        self.assertIsNotNone(cal_convexity(cs_last, fig_path, scene, shape, margin))


if __name__ == '__main__':
    unittest.main()

# run this with: python -m unittest unit_test.Test_sum_stats

