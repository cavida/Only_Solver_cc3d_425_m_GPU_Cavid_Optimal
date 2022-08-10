from cc3d.core.PySteppables import *
import random, time, csv
from scipy.spatial import distance
random.seed(123)


class CaNodule(SteppableBasePy):

    def __init__(self, frequency):

        SteppableBasePy.__init__(self, frequency)

    def start(self):
        print('CaNodule steppable started!')

        # Initializing the fields
        # self.MDE_field = self.field.MDE
        # self.ECM_fibers_field = self.field.ECM_Fibers
        # self.MDE_field_secretor = self.get_field_secretor("MDE")
        self.DEAD_Cell_Field = self.field.DEAD_Cells
        self.TGFb_field = self.field.TGFb


    def step(self, mcs):

        cell = self.fetch_cell_by_id(5)
        print('***** From CaNodule class: ',self.DEAD_Cell_Field[cell.xCOM, cell.yCOM, cell.zCOM], self.TGFb_field[cell.xCOM, cell.yCOM, cell.zCOM])


        # Nucleation for generating calcium nodules
        for cell in self.cell_list_by_type(self.OST, self.DEAD):
            if cell.dict["TIME_IN_STATE"] == 0:
                pixel_list = self.get_cell_boundary_pixel_list(cell)
                r = random.choice(list(pixel_list))  # Choosing a random pixel at the boundary of a cell
                new_nodule = self.new_cell(self.CANOD)
                new_nodule.lambdaVolume = 50
                if cell.type == self.OST:
                    new_nodule.dict['Type'] = 'OST'
                elif cell.type == self.DEAD:
                    new_nodule.dict['Type'] = 'DEAD'
                self.cell_field[r.pixel.x, r.pixel.y, r.pixel.z] = new_nodule
            else:
                cell.dict["TIME_IN_STATE"] += self.frequency

        for cell in self.cell_list_by_type(self.AVEC):
            cell.dict["LAST_COM"] = [cell.xCOM, cell.yCOM, cell.zCOM]

        # Calcium nodule growth
        for cell in self.cell_list_by_type(self.CANOD):
            CaNod_growth_rate = 10 * self.DEAD_Cell_Field[cell.xCOM, cell.yCOM, cell.zCOM]
            if cell.dict['Type'] == 'OST':
                cell.targetVolume += CaNod_growth_rate
            else:
                cell.targetVolume += CaNod_growth_rate

    # def finish(self):
    #     pass
    #
    # def on_stop(self):
    #     return
