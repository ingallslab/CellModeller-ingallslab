# this code is used to test functionality of fractal_demension.py
import pickle
import sys

sys.path.append('../../scripts')
from fractal_dimension import calc_fractal_dimension


def sim_data():
    # tun on simulation data
    sim_pickle_file = '../../input files/CellModeller_Sample/step-00200.pickle'
    # read current pickle file
    sim_bacteria_info = pickle.load(open(sim_pickle_file, 'rb'))
    cs = sim_bacteria_info['cellStates']
    # input arguments: cells (cellStates dict), fig_export_path (directory to export images), fig_name (image name)
    sim_fractal_dimension = calc_fractal_dimension(cs, 'fig/', 'sim_img')
    print("simulation data: " + str(sim_fractal_dimension))


def exp_data():
    # tun on experimental data
    exp_pickle_file = '../../input files/K12/step-000089.pickle'
    # read current pickle file
    exp_bacteria_info = pickle.load(open(exp_pickle_file, 'rb'))
    cs = exp_bacteria_info['cellStates']
    # input arguments: cells (cellStates dict), fig_export_path (directory to export images), fig_name (image name)
    exp_fractal_dimension = calc_fractal_dimension(cs, 'fig/', 'exp_img')
    print("experimental data: " + str(exp_fractal_dimension))


if __name__ == '__main__':
    sim_data()
    exp_data()
