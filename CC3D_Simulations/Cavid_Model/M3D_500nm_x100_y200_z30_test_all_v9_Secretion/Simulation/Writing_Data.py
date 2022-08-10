from cc3d.core.PySteppables import *
import random, time, csv
from scipy.spatial import distance
random.seed(123)


class Writing_Data(SteppableBasePy):

    def __init__(self, frequency):

        SteppableBasePy.__init__(self, frequency)

    def start(self):
        print('Writing_Data steppable started!')
        # Writing data to files
        ctime = time.ctime()

        # Cell count
        name = 'C:\\Users\jazimib1\CC3DWorkspace\General_Output\ ' + 'CellCount ' + ctime.replace(':', ',') + '.csv'
        self.f_CellCount = open(name, "w", newline='')
        self.f_csv_CellCount = csv.writer(self.f_CellCount)
        self.f_csv_CellCount.writerow(('MCS', 'qVEC', 'aVEC', 'ECM', 'nECM', 'qVIC', 'aVIC', 'OST', 'DEAD', 'CaNod'))

        # Cell Location
        name = 'C:\\Users\jazimib1\CC3DWorkspace\General_Output\ ' + 'qVEC location ' + ctime.replace(':', ',') + '.csv'
        self.f_qVEC = open(name, "w", newline='')
        self.f_csv_qVEC = csv.writer(self.f_qVEC)
        self.f_csv_qVEC.writerow(('MCS','Cell ID, (Location), Speed [MCS/min], [TGFBeta], [MMP9], [alphaSMA], [PECAM1] '))

        name = 'C:\\Users\jazimib1\CC3DWorkspace\General_Output\ ' + 'aVEC location ' + ctime.replace(':', ',') + '.csv'
        self.f_aVEC = open(name, "w", newline='')
        self.f_csv_aVEC = csv.writer(self.f_aVEC)
        self.f_csv_aVEC.writerow(
            ('MCS', 'Cell ID, (Location), Speed [MCS/min], [TGFBeta], [MMP9], [alphaSMA], [PECAM1] '))

        name = 'C:\\Users\jazimib1\CC3DWorkspace\General_Output\ ' + 'qVIC location ' + ctime.replace(':', ',') + '.csv'
        self.f_qVIC = open(name, "w", newline='')
        self.f_csv_qVIC = csv.writer(self.f_qVIC)
        self.f_csv_qVIC.writerow(
            ('MCS', 'Cell ID, (Location), Speed [MCS/min], [TGFBeta], [MMP9], [alphaSMA], [PECAM1] '))

        name = 'C:\\Users\jazimib1\CC3DWorkspace\General_Output\ ' + 'aVIC location ' + ctime.replace(':', ',') + '.csv'
        self.f_aVIC = open(name, "w", newline='')
        self.f_csv_aVIC = csv.writer(self.f_aVIC)
        self.f_csv_aVIC.writerow(
            ('MCS', 'Cell ID, (Location), Speed [MCS/min], [TGFBeta], [MMP9], [alphaSMA], [PECAM1] '))

        name = 'C:\\Users\jazimib1\CC3DWorkspace\General_Output\ ' + 'OST location ' + ctime.replace(':', ',') + '.csv'
        self.f_OST = open(name, "w", newline='')
        self.f_csv_OST = csv.writer(self.f_OST)
        self.f_csv_OST.writerow(
            ('MCS', 'Cell ID, (Location), Speed [MCS/min], [TGFBeta], [MMP9], [alphaSMA], [PECAM1] '))

        name = 'C:\\Users\jazimib1\CC3DWorkspace\General_Output\ ' + 'Calcium Nodule location and size ' + ctime.replace(
            ':', ',') + '.csv'
        self.f_CaNod = open(name, "w", newline='')
        self.f_csv_CaNod = csv.writer(self.f_CaNod)
        self.f_csv_CaNod.writerow(
            ('MCS', 'Cell ID, (Location), Cell Volume '))


    def step(self, mcs):
        # Writing cell data to files
        if mcs % self.frequency == 0:
            print('********* Writing data and frequency is: ', self.frequency)
            self.f_csv_CellCount.writerow(
                (mcs, len(self.cell_list_by_type(self.QVEC)), len(self.cell_list_by_type(self.AVEC)), \
                 len(self.cell_list_by_type(self.ECM)), len(self.cell_list_by_type(self.NECM)), \
                 len(self.cell_list_by_type(self.QVIC)), len(self.cell_list_by_type(self.AVIC)), \
                 len(self.cell_list_by_type(self.OST)), len(self.cell_list_by_type(self.DEAD)), \
                 len(self.cell_list_by_type(self.CANOD)) \
                 ))

            # Printing (cell ID, cell location, cell speed [pixel/MCS])
            row = [mcs]
            for cell in self.cell_list_by_type(self.QVEC):
                loc = (cell.xCOM, cell.yCOM, cell.zCOM)
                row.append((cell.id, loc, distance.euclidean(loc, cell.dict['LAST_COM']) / self.frequency,\
                            cell.sbml.model_M['[TGFBeta]'],cell.sbml.model_M['[MMP9]'],\
                            cell.sbml.model_M['[alphaSMA]'],cell.sbml.model_M['[PECAM1]']))
                cell.dict['LAST_COM'] = loc
            self.f_csv_qVEC.writerow(row)

            row = [mcs]
            for cell in self.cell_list_by_type(self.AVEC):
                loc = (cell.xCOM, cell.yCOM, cell.zCOM)
                row.append((cell.id, loc, distance.euclidean(loc, cell.dict['LAST_COM']) / self.frequency,\
                            cell.sbml.model_M['[TGFBeta]'],cell.sbml.model_M['[MMP9]'],\
                            cell.sbml.model_M['[alphaSMA]'],cell.sbml.model_M['[PECAM1]']))
                cell.dict['LAST_COM'] = loc
            self.f_csv_aVEC.writerow(row)

            row = [mcs]
            for cell in self.cell_list_by_type(self.QVIC):
                loc = (cell.xCOM, cell.yCOM, cell.zCOM)
                row.append((cell.id, loc, distance.euclidean(loc, cell.dict['LAST_COM']) / self.frequency,\
                            cell.sbml.model_M['[TGFBeta]'],cell.sbml.model_M['[MMP9]'],\
                            cell.sbml.model_M['[alphaSMA]'],cell.sbml.model_M['[PECAM1]']))
                cell.dict['LAST_COM'] = loc
            self.f_csv_qVIC.writerow(row)

            row = [mcs]
            for cell in self.cell_list_by_type(self.AVIC):
                loc = (cell.xCOM, cell.yCOM, cell.zCOM)
                row.append((cell.id, loc, distance.euclidean(loc, cell.dict['LAST_COM']) / self.frequency,\
                            cell.sbml.model_M['[TGFBeta]'],cell.sbml.model_M['[MMP9]'],\
                            cell.sbml.model_M['[alphaSMA]'],cell.sbml.model_M['[PECAM1]']))
                cell.dict['LAST_COM'] = loc
            self.f_csv_aVIC.writerow(row)

            row = [mcs]
            for cell in self.cell_list_by_type(self.OST):
                loc = (cell.xCOM, cell.yCOM, cell.zCOM)
                row.append((cell.id, loc, distance.euclidean(loc, cell.dict['LAST_COM']) / self.frequency,\
                            cell.sbml.model_M['[TGFBeta]'],cell.sbml.model_M['[MMP9]'],\
                            cell.sbml.model_M['[alphaSMA]'],cell.sbml.model_M['[PECAM1]']))
                cell.dict['LAST_COM'] = loc
            self.f_csv_OST.writerow(row)

            # Calcium nodule location and size
            row = [mcs]
            for cell in self.cell_list_by_type(self.CANOD):
                loc = (cell.xCOM, cell.yCOM, cell.zCOM)
                row.append((cell.id, loc, cell.volume))
                cell.dict['LAST_COM'] = loc
            self.f_csv_CaNod.writerow(row)


        # # This must be tha last part of code. If this part of the code executed at first, the current locations and last would be the same
        # for cell in self.cell_list_by_type(self.QVEC, self.AVEC, self.QVIC, self.AVIC, self.OST):
        #     cell.dict['LAST_COM'] = (cell.xCOM, cell.yCOM, cell.zCOM)

    def finish(self):
        self.f_CellCount.close()
        self.f_qVEC.close()
        self.f_aVEC.close()
        self.f_qVIC.close()
        self.f_aVIC.close()
        self.f_OST.close()
        self.f_CaNod.close()
        # self.f_CellInvasion.close()
        # self.f_CellDistance.close()
        # self.f_ChemicalConcentrations.close()

    def on_stop(self):
        self.f_CellCount.close()
        self.f_qVEC.close()
        self.f_aVEC.close()
        self.f_qVIC.close()
        self.f_aVIC.close()
        self.f_OST.close()
        self.f_CaNod.close()
        # self.f_CellInvasion.close()
        # self.f_CellDistance.close()
        # self.f_ChemicalConcentrations.close()
        return
