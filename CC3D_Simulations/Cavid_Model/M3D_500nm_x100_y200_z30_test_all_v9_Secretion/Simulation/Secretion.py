from cc3d.core.PySteppables import *
import random, time, csv
from scipy.spatial import distance
random.seed(123)


class Secretion(SteppableBasePy):

    def __init__(self, frequency):

        SteppableBasePy.__init__(self, frequency)

    def start(self):
        print('CaNodule steppable started!')

        # Initializing the fields
        self.MDE_field = self.field.MDE
        # self.ECM_fibers_field = self.field.ECM_Fibers
        self.MDE_field_secretor = self.get_field_secretor("MDE")
        # self.DEAD_Cell_Field = self.field.DEAD_Cells
        # self.TGFb_field = self.field.TGFb


    def step(self, mcs):

        for cell in self.cell_list_by_type(self.AVEC):
            # self.MDE_field_secretor.secreteInsideCellConstantConcentration(cell, cell.sbml.model_M['[MMP9]'])
            # print('MMP9 ************** ', cell.sbml.model_M['[MMP9]'])
            self.MDE_field[cell.xCOM, cell.yCOM, cell.zCOM] = cell.sbml.model_M['[MMP9]'] #This implementation
            #  assigns a concentration at COM at it is faster

    # def finish(self):
    #     pass
    #
    # def on_stop(self):
    #     return
