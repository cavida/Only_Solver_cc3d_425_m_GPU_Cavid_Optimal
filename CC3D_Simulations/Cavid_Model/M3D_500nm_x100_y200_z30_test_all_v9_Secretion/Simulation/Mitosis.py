from cc3d.core.PySteppables import *
import random, time, csv
from scipy.spatial import distance
random.seed(123)


class Mitosis(MitosisSteppableBase):

    def __init__(self, frequency):

        MitosisSteppableBase.__init__(self, frequency)

    def start(self):
        print('Mitosis steppable started!')


    def step(self, mcs):

        aVIC_cells_to_divide=[]
        for cell in self.cell_list_by_type(self.AVIC):
            cell.targetVolume += 4.4  # 4.4 [voxel/min] [ This value comes from the doubling time of VICs (30 Hours) and a VIC volume of 8000 voxel
            if cell.volume > 16000:
                aVIC_cells_to_divide.append(cell)

        for cell in aVIC_cells_to_divide:
            self.divide_cell_along_minor_axis(cell)


    def update_attributes(self):
        # reducing parent target volume
        self.parent_cell.targetVolume /= 2.0

        self.clone_parent_2_child()

    # def finish(self):
    #   pass
    #
    # def on_stop(self):
    #
    #     return
