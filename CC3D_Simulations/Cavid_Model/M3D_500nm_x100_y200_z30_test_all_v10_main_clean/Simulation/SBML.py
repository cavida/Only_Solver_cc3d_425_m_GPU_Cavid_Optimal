from cc3d.core.PySteppables import *
import sys
# sys.path.append('Y:\Jazimib1\CC3D\CompuCell3D-py3-64bit_4.2.5_GPU\lib\site-packages')

from CC3D_Simulations.Cavid_Model.M3D_500nm_x100_y200_z30_test_all_v10_main_clean import model
import roadrunner
from CC3D_Simulations.Cavid_Model.M3D_500nm_x100_y200_z30_test_all_v10_main_clean.GlobalValues import *
import random, math, time, csv
from scipy.spatial import distance
random.seed(123)


class SBML(SteppableBasePy):

    def __init__(self, frequency):

        SteppableBasePy.__init__(self, frequency)

        # Cavid parameters
        self.SBML_TimeStep = 60
        self.m_freq = 10                    # Frequency in which ECM remodeling calculations are done
        self.num_ECM_degraded = 0
        self.ECM_deg_thre = 1
        self.moving_thre = 2
        self.ECM_Fibers_thre = 0.6
        self.surface_fraction_thre = 0.9
        self.VIC_nums = 10000
        self.um_2_pixel_ratio = 0.5

        #Cavid SBML Values
        self.init_TGFb_value = 0.2

    def start(self):
        print('SBML steppable started!')

        # Sub Cellular part SBML Model
        self.model_sbml = model.antimonyToSBML(model.Mukti_model_antimony)
        # self.r_model = roadrunner.RoadRunner(self.model_sbml)

        medium_volume = self.dim.x * self.dim.y * self.dim.z * 1e-15 * (self.um_2_pixel_ratio ** 3)
        print('@@@@@@@ Meidum is: ', medium_volume)

        self.add_sbml_to_cell_types(model_string=self.model_sbml, model_name='model_M', cell_types=[self.QVEC, self.QVIC, self.AVEC, self.AVIC, self.OST],
                                    step_size=self.SBML_TimeStep, initial_conditions = {'[TGFBeta]': self.init_TGFb_value})


        # self.add_sbml_to_cell_types(model_string=self.model_sbml, model_name='model_M', cell_types=[self.QVEC, self.QVIC, self.AVEC, self.AVIC, self.OST],
        #                             step_size=self.SBML_TimeStep, initial_conditions = {'[TGFBeta]': self.init_TGFb_value, 'Medium': medium_volume})
        # self.add_antimony_to_cell_types(model_string=self.model_sbml, model_name='model_M', cell_types=[self.QVEC, self.QVIC, self.AVEC, self.AVIC, self.OST],
        #                             step_size=self.SBML_TimeStep, initial_conditions = {'[TGFBeta]': self.init_TGFb_value, 'Medium': medium_volume})

        # self.add_free_floating_sbml(model_string=self.model_sbml, model_name='model_M_free', step_size=self.SBML_TimeStep)
        # self.Medium = self.sbml.model_M_free['Medium']
        # self.Cytoplasm = self.sbml.model_M_free['Cytoplasm']
        # self.Nucleus = self.sbml.model_M_free['Nucleus']
        #
        # print('*** Model Medium Volume is: ',self.Medium)




    def step(self, mcs):

        #SBML Calculations
        tsbml0 = time.perf_counter()
        self.timestep_cell_sbml()
        tsbml1 = time.perf_counter()
        print('SBML runtime: ',tsbml1-tsbml0)

        # # Testing Medium value
        # for cell in self.cell_list_by_type(self.AVEC):
        #     print(' Medium value after applying to the cell : ', cell.sbml.model_M['Medium'] )



    # def finish(self):
    #     pass
    #
    # def on_stop(self):
    #     return
