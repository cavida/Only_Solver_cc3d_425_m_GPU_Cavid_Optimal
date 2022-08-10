
from cc3d.core.PySteppables import *
from CC3D_Simulations.Cavid_Model.M3D_500nm_x100_y200_z30_test_all_v9_Secretion import model
import roadrunner
from CC3D_Simulations.Cavid_Model.M3D_500nm_x100_y200_z30_test_all_v9_Secretion.GlobalValues import *
import random, math, time, csv
from scipy.spatial import distance
random.seed(123)

class Main(SteppableBasePy):

    def __init__(self,frequency):

        SteppableBasePy.__init__(self,frequency)
        self.finalEndothelialCellCount = 0
        self.finalTumorCellCount = 0
        self.finalECMCount = 0
        self.initEndothelialCellCount = 0
        self.initTumorCellCount = 0
        self.initECMCount = 0
        self.activatedCellId = 0

        # Cavid parameters
        self.SBML_TimeStep = 60
        self.m_freq = 10                    # Frequency in which ECM remodeling calculations are done
        self.num_ECM_degraded = 0
        self.ECM_deg_thre = 1
        self.moving_thre = 2
        self.ECM_Fibers_thre = 0.6
        self.surface_fraction_thre = 0.9
        self.VIC_nums = 10000

        #Cavid SBML Values
        self.init_TGFb_value = 0.2

        # Make the global constants local to reduce run time
        self.conversionPixToMeter = conversionPixToMeter
        self.conversionPascalToMmHg = conversionPascalToMmHg
        self.universalGasConstant = universalGasConstant
        self.bodyTemperature = bodyTemperature
        self.minimumO2PressureEndothelial = minimumO2PressureEndothelial
        self.maximumCO2PressureEndothelial = maximumCO2PressureEndothelial
        self.minimumO2PressureCancer = minimumO2PressureCancer
        self.maximumCO2PressureCancer = maximumCO2PressureCancer
        self.activatedCellTimeThreshold = activatedCellTimeThreshold
        self.mediumStartCoordinate = mediumStartCoordinate
        self.ecmDegradationByMde = ecmDegradationByMde
        self.initCytokineConcentration = initCytokineConcentration
        self.maxCytokineConcentration = maxCytokineConcentration
        self.initNumberOfCadherins = initNumberOfCadherins
        self.cadherinThreshold = cadherinThreshold
        self.initNumberOfIntegrins = initNumberOfIntegrins
        self.initNumberOfLigands = initNumberOfLigands
        self.initECMConcentration = initECMConcentration

    def start(self):
        print('Main steppable started!')

        #Generating VICs
        for i in range(int (self.VIC_nums * (self.dim.x * self.dim.y * self.dim.z)/(125*10E6))):
            x = random.randint(1,self.dim.x-20)
            y = random.randint(1, self.dim.y-30) # This 150 is the place that should be subtraceted from the tissue
            z = random.randint(1, self.dim.z-20)
            QVIC_cell = self.new_cell(self.QVIC)
            self.cell_field[x:x+20,y:y+20,z:z+20] = QVIC_cell
            QVIC_cell.targetVolume = QVIC_cell.volume
            QVIC_cell.lambdaVolume = 10

        self.idList_qVIC = []
        for cell in self.cell_list_by_type(self.QVIC):
            self.idList_qVIC.append(int(cell.id))
        # print('ID list of qVICs: ',self.idList_qVIC)

        self.idList_qVEC = []
        for cell in self.cell_list_by_type(self.QVEC):
            self.idList_qVEC.append(int(cell.id))
        # print('ID list of qVECs: ',self.idList_qVEC)

        # Initializing activated VICs and activated VECs
        for i in range(10):
            # cell = self.cell_field[50,355,15]
            qVEC_cell = self.fetch_cell_by_id(random.choice(self.idList_qVEC))
            qVEC_cell.type = self.AVEC
            print('Cell lambda and target volume values: ', qVEC_cell.lambdaVolume, " ", qVEC_cell.volume)
            qVEC_cell.lambdaVolume = 100
            print('Cell id:', qVEC_cell.id)
            qVEC_cell.dict["TIME_IN_STATE"] = 0

        # Activating some VICs
        for i in range(1):
            qVIC_cell = self.fetch_cell_by_id(random.choice(self.idList_qVIC))
            qVIC_cell.type = self.AVIC

       




        for cell in self.cell_list_by_type(self.ECM):
            cell.targetVolume = cell.volume
            cell.lambdaVolume = 10000

        for cell in self.cell_list_by_type(self.QVEC, self.AVEC, self.QVIC, self.AVIC, self.OST):
            cell.lambdaVolume = 20
            cell.targetVolume = cell.volume
            # ecmField[int(round(cell.xCOM)), int(round(cell.yCOM)), int(round(cell.zCOM))] = 0
            cell.dict["INIT_COM"] = [cell.xCOM, cell.yCOM, cell.zCOM]
            cell.dict["LAST_COM"] = cell.dict["INIT_COM"]
            # The following part was in the for loop
            cell.dict["TIME_IN_STATE"] = 0
            cell.dict["TIME_EXPOSED_TO_CYTOKINES"] = 0
            cell.dict["AGE"] = 0
            cell.dict["GENERATION"] = 0






        # # Sub Cellular part SBML Model
        # self.model_sbml = model.antimonyToSBML(model.Mukti_model_antimony)
        # self.r_model = roadrunner.RoadRunner(self.model_sbml)
        #
        # self.add_sbml_to_cell_types(model_string=self.model_sbml, model_name='model_M', cell_types=[self.QVEC, self.QVIC, self.AVEC, self.AVIC],
        #                             step_size=self.SBML_TimeStep, initial_conditions = {'[TGFBeta]': self.init_TGFb_value})
        # self.add_free_floating_sbml(model_string=self.model_sbml, model_name='model_M_free', step_size=self.SBML_TimeStep)
        # self.Medium = self.sbml.model_M_free['Medium']
        # self.Cytoplasm = self.sbml.model_M_free['Cytoplasm']
        # self.Nucleus = self.sbml.model_M_free['Nucleus']
        #
        # print('*** Model Medium Volume is: ',self.Medium)




        # Setting MDE concentration
        mdeField = self.field.MDE
        mdeField[0:self.dim.x:1, 0:self.dim.y:1, 0:self.dim.z:1] = 0

        # Setting Cytokine concentration
        # dose (in nmol/mL) * 1E-09 (for convert to mole) * 0.3E-03 (mL) / (30 * 100 * 100 voxels)

        # cytokineField = self.field.Cytokines
        # cytokineField[0:self.dim.x:1, mediumStartCoordinate:self.dim.y:1, 0:self.dim.z:1] = 0.044  # 2 ng/mL

        # Setting initial ECM concentration and parameters
        # ecmField = self.field.ECMProteins



        #Applying lambda and volume
        for cell in self.cell_list_by_type(2,4):
            cell.targetVolume = cell.volume
            cell.lambdaVolume = 5



        # Initializing the fields
        self.MDE_field = self.field.MDE
        self.ECM_fibers_field = self.field.ECM_Fibers
        self.MDE_field_secretor = self.get_field_secretor("MDE")
        self.DEAD_Cell_Field = self.field.DEAD_Cells
        self.TGFb_field = self.field.TGFb



        # self.activatedCellId = self.cell_field[25, 42, 25].id  # Activating one of the quiesent cells

        # self.pW1 = self.addNewPlotWindow(_title='PECAM1 Amount', _xAxisTitle='MCS',
        #                                 _yAxisTitle='Amount', _xScaleType='linear', _yScaleType='linear')
        # self.pW1.addPlot('PECAM1', _style='Dots', _color='red', _size=5)
        #
        # self.pW2 = self.addNewPlotWindow(_title='# Integrins', _xAxisTitle='MCS',
        #                                 _yAxisTitle='Amount', _xScaleType='linear', _yScaleType='linear')
        # self.pW2.addPlot('integrins', _style='Dots', _color='green', _size=5)
        #
        # self.pW3 = self.addNewPlotWindow(_title='Alpha SMA Amount', _xAxisTitle='MCS',
        #                                 _yAxisTitle='Amount', _xScaleType='linear', _yScaleType='linear')
        # self.pW3.addPlot('alphaSMA', _style='Dots', _color='yellow', _size=5)
        #
        # self.pW4 = self.addNewPlotWindow(_title='TGFB', _xAxisTitle='MCS',
        #                                 _yAxisTitle='TGFB', _xScaleType='linear', _yScaleType='linear')
        # #self.pW4.addPlot('Tumor', _style='Dots', _color='green', _size=5)
        # self.pW4.addPlot('TGFB', _style='Dots', _color='red', _size=5)
        #
        # self.pW5 = self.addNewPlotWindow(_title='MMP9 Amount', _xAxisTitle='MCS',
        #                                 _yAxisTitle='PECAM1', _xScaleType='linear', _yScaleType='linear')
        # self.pW5.addPlot('MMP9', _style='Dots', _color='red', _size=5)

        # Writing data to files

        # ctime = time.ctime()
        #
        # name = 'C:\\Users\jazimib1\CC3DWorkspace\General_Output\ '+'EndothelialCellCount ' + ctime.replace(':', ',') + '.csv'
        # self.f_EndothelialCellCount = open(name, "w", newline='')
        # self.f_csv_EndothelialCellCount = csv.writer(self.f_EndothelialCellCount)
        # self.f_csv_EndothelialCellCount.writerow(('mcs', 'qVEC', 'aVEC', 'ECM', 'nECM','qVIC','aVIC','OST','DEAD','CaNod'))
        #
        # name = 'C:\\Users\jazimib1\CC3DWorkspace\General_Output\ '+'CellInvasion ' + ctime.replace(':', ',') + '.csv'
        # self.f_CellInvasion = open(name, "w", newline='')
        # self.f_csv_CellInvasion = csv.writer(self.f_CellInvasion)
        # self.f_csv_CellInvasion.writerow(('mcs', 'ECM','Displacement','RoG','PECAM1','alphaSMA','TGFB','MMP9'))
        #
        # name = 'C:\\Users\jazimib1\CC3DWorkspace\General_Output\ '+'CellDistance ' + ctime.replace(':', ',') + '.csv'
        # self.f_CellDistance = open(name, "w", newline='')
        # self.f_csv_CellDistance = csv.writer(self.f_CellDistance)
        # self.f_csv_CellDistance.writerow(('mcs', 'Endothelial','Tumor'))
        #
        # self.Chemicals = ('SD093', 'TGFBeta','TGFBeta_Nonspecific','TGFBeta_Cytoplasm','TGFBR_Complex','TGFBR_Complex_Cytoplasm','TGFBR1_Surface',\
        #              'TGFBR1_Cytoplasm','TGFBR1_Surface_Inhibited','TGFBR2_Surface','TGFBR2_Cytoplasm','Smad2','Smad2_Nucleus','pSmad2',\
        #              'pSmad2_Nucleus','Smad4','Smad4_Nucleus','pSmad24','pSmad24_Nucleus','pSmad22','pSmad22_Nucleus','Smad3','Smad3_Nucleus',\
        #              'pSmad3','pSmad3_Nucleus','pSmad34','pSmad34_Nucleus','pSmad33','pSmad33_Nucleus','Snail','alphaSMA','Slug','PECAM1','MMP9',\
        #              'mRNASnail','mRNAalphaSMA','mRNASlug','mRNAPECAM1','mRNAMMP9','Amino_Acids_Placeholder','Degraded_Proteins')
        # name = 'C:\\Users\jazimib1\CC3DWorkspace\General_Output\ '+'ChemicalConcentrations ' + ctime.replace(':', ',') + '.csv'
        # self.f_ChemicalConcentrations = open(name, "w", newline='')
        # self.f_csv_ChemicalConcentrations = csv.writer(self.f_ChemicalConcentrations)
        # self.f_csv_ChemicalConcentrations.writerow(['mcs']+list(self.Chemicals))



    def step(self,mcs):
        print('MCS: ', mcs)

        #Reading from the fields
        if mcs % self.m_freq == 0:
            for cell in self.cell_list_by_type(self.ECM):
                if self.MDE_field[cell.xCOM, cell.yCOM, cell.zCOM] > self.ECM_deg_thre:
                    self.delete_cell(cell)
                    self.num_ECM_degraded += 1

            # if mcs % 100 == 0:
            # Cell type change implementation
            for cell in self.cell_list_by_type(self.QVEC):
                if cell.sbml.model_M['[PECAM1]'] < 0.5 * cell.sbml.model_M["init_PECAM1"]:
                    cell.type = self.AVEC
                    cell.dict["TIME_IN_STATE"] = 0
                    cell.dict["AGE"] += self.m_freq
                else:
                    cell.dict["TIME_IN_STATE"] += self.m_freq
                    cell.dict["AGE"] += self.m_freq

            # Loop for activated endothelial cells
            for cell in self.cell_list_by_type(self.AVEC):
                neighbor_list = self.get_cell_neighbor_data_list(cell)
                surface_fraction = list(neighbor_list[0])[1] / cell.surface # neighbor_list[0] takes the first element in the list which is
                # associated with media, and list(neighbor_list[0])[1] takes second element which is common surface area
                # print('Cell area fraction: ', surface_fraction)

                if cell.sbml.model_M['[PECAM1]'] >= 0.5 * cell.sbml.model_M["init_PECAM1"] and \
                    cell.dict["TIME_IN_STATE"] < 2880: #After two days of staying in this condition, it will be superactivated mode and should not goes back to QVEC
                        # cell.type = self.QVEC
                        cell.dict["TIME_IN_STATE"] = 0
                        cell.dict["AGE"] += self.m_freq

                elif cell.dict["TIME_IN_STATE"] > 14400:      # After 10 days of being in TGFb medium, they will turn into OST type
                    cell.type = self.OST
                    cell.dict["TIME_IN_STATE"] = 0
                    cell.dict["AGE"] += self.m_freq

                elif surface_fraction < self.surface_fraction_thre or \
                        (self.ECM_fibers_field[cell.xCOM, cell.yCOM, cell.zCOM] > self.ECM_Fibers_thre and \
                        distance.euclidean(cell.dict['LAST_COM'],[cell.xCOM,cell.yCOM,cell.zCOM]) < self.moving_thre):
                    cell.type = self.DEAD
                    cell.lambdaVolume = 1000
                    cell.dict["TIME_IN_STATE"] = 0
                    cell.dict["AGE"] += self.m_freq

                else:
                    cell.dict["TIME_IN_STATE"] += self.m_freq
                    cell.dict["AGE"] += self.m_freq

            # Loop for Ostioblastic cells
            for cell in self.cell_list_by_type(self.OST):
                neighbor_list = self.get_cell_neighbor_data_list(cell)
                surface_fraction = list(neighbor_list[0])[1] / cell.surface  # neighbor_list[0] takes the first element in the list which is
                # associated with media, and list(neighbor_list[0])[1] takes second element which is common surface area
                # print('Cell area fraction: ', surface_fraction)
                if surface_fraction < self.surface_fraction_thre or \
                     (self.ECM_fibers_field[cell.xCOM, cell.yCOM, cell.zCOM] > self.ECM_Fibers_thre and \
                      distance.euclidean(cell.dict['LAST_COM'], [cell.xCOM, cell.yCOM, cell.zCOM]) < self.moving_thre):
                    cell.type = self.DEAD
                    cell.lambdaVolume = 1000
                    cell.dict["TIME_IN_STATE"] = 0
                    cell.dict["AGE"] += self.m_freq

                else:
                    cell.dict["TIME_IN_STATE"] += self.m_freq
                    cell.dict["AGE"] += self.m_freq

            #Setting chemotaxis of AVEC toward ECMs by using ECMChemo field
            for cell in self.cell_list_by_type(self.AVEC):
                cd = self.chemotaxisPlugin.addChemotaxisData(cell, "ECM_Chemo")
                if cell.yCOM < cell.dict["INIT_COM"][1] - 120:          # Cell when is deeper than 60 um inside the tissue
                    cd.setLambda(1000)
                else:
                    cd.setLambda(2000)

        # #SBML Calculations
        # tsbml0 = time.perf_counter()
        # self.timestep_cell_sbml()
        # tsbml1 = time.perf_counter()
        # print('SBML runtime: ',tsbml1-tsbml0)

        #Secretion in this part after SBML calclulation
        # for cell in self.cell_list_by_type(self.AVEC):
        #     # self.MDE_field[cell.xCOM, cell.yCOM, cell.zCOM] = 1
        #     self.MDE_field_secretor.secreteInsideCellConstantConcentration(cell,cell.sbml.model_M['[MMP9]'])
        #     # cell.sbml.model_M["MMP9"]
        # print('****** PECAM1: ', cell.sbml.model_M["[PECAM1]"])

        # # Nucleation for generating calcium nodules
        # if mcs % self.m_freq == 0:
        #     for cell in self.cell_list_by_type(self.OST, self.DEAD):
        #         if cell.dict["TIME_IN_STATE"] == 0:
        #             pixel_list = self.get_cell_boundary_pixel_list(cell)
        #             r = random.choice(list(pixel_list))  # Choosing a random pixel at the boundary of a cell
        #             new_nodule = self.new_cell(self.CANOD)
        #             new_nodule.lambdaVolume = 50
        #             if cell.type == self.OST:
        #                 new_nodule.dict['Type'] = 'OST'
        #             elif cell.type == self.DEAD:
        #                 new_nodule.dict['Type'] = 'DEAD'
        #             self.cell_field[r.pixel.x, r.pixel.y, r.pixel.z] = new_nodule
        #         else:
        #             cell.dict["TIME_IN_STATE"] += self.m_freq
        #
        #     for cell in self.cell_list_by_type(self.AVEC):
        #         cell.dict["LAST_COM"] = [cell.xCOM, cell.yCOM, cell.zCOM]
        #
        # # Calcium nodule growth
        # for cell in self.cell_list_by_type(self.CANOD):
        #     CaNod_growth_rate = 10 * self.DEAD_Cell_Field[cell.xCOM, cell.yCOM, cell.zCOM]
        #     if cell.dict['Type'] == 'OST':
        #         cell.targetVolume += CaNod_growth_rate
        #     else:
        #         cell.targetVolume += CaNod_growth_rate




        # # Quiescent VICs turn into activated VICs
        # if mcs==200:
        #     for cell in self.cell_list_by_type(self.QVIC):
        #         if cell.id == 102:
        #             cell.type = self.AVIC
        #             cell.dict["TIME_IN_STATE"] = 0
        #
        # # Activated VICs turn into quiescent VICs
        # if False:
        #     for cell in self.cell_list_by_type(self.AVIC):
        #             cell.type = self.QVIC
        #             cell.dict["TIME_IN_STATE"] = 0
        #
        #
        # # Activated VECs (aVEC) tuning into osteobalstic cells
        # idList_aVEC_aVIC = []
        # if mcs % 500 == 0:
        #     for cell in self.cell_list_by_type(self.AVEC,self.AVIC):
        #         if cell.yCOM < 300:
        #             idList_aVEC_aVIC.append(int(cell.id))
        #     print(idList_aVEC_aVIC)
        #     try:
        #         aVIC_or_aVEC_cell = self.fetch_cell_by_id(random.choice(idList_aVEC_aVIC))
        #         aVIC_or_aVEC_cell.type = self.OST
        #         aVIC_or_aVEC_cell.lambdaVolume = 100
        #         aVIC_or_aVEC_cell.dict["TIME_IN_STATE"] = 0
        #     except:
        #         print('No Active cell!')
        #
        #
        # # Activated VICs (aVIC) Tuning into osteobalstic cells
        # # if mcs==40:
        # #     for cell in self.cell_list_by_type(self.AVIC):
        # #         cell.type = self.OST
        # #         cell.dict["TIME_IN_STATE"] = 0
        #
        # # Death of Osteoblastic cells
        # if mcs==1000:
        #     for cell in self.cell_list_by_type(self.OST):
        #         cell.type = self.DEAD
        #         cell.dict["TIME_IN_STATE"] = 0

        #
        # # ECM Fibrosis *******************************************************************
        #
        # for cell in self.cell_list_by_type(self.AVIC,self.AVEC):
        #     pass
        #
        # Calcification from osteoblasts


        # cell=self.fetch_cell_by_id(16)
        # print('Cell Y Coordinate: ', cell.yCOM,'Cell Volume',cell.volume)
        #
        # #Solving
        #
        #
        #
        # # Reading values from SBML model
        #
        # #Plotting ***************************************************************************
        # if mcs % 20 == 0:
        #     self.pW1.addDataPoint("PECAM1", mcs, pecamSum)
        #     self.pW2.addDataPoint("integrins", mcs, integrinNum)
        #     self.pW3.addDataPoint("alphaSMA", mcs, alphaSMASum)
        #     self.pW4.addDataPoint("TGFB", mcs, tgfbSum)
        #     self.pW5.addDataPoint("MMP9", mcs, mmpSum)
        #
        # t4 = time.perf_counter()
        #
        # # My codes to write files
        #
        # self.f_csv_EndothelialCellCount.writerow((mcs, qEndo, pEndo, aEndo, saEndo))
        #
        # disp = (disp / (aEndo + saEndo)) if (aEndo + saEndo) != 0 else 0
        # rog = (math.pow((rog / (aEndo + saEndo)), 0.5)) if (aEndo + saEndo) != 0 else 0
        # self.f_csv_CellInvasion.writerow((mcs, numECM, disp,  rog,  pecamSum,  alphaSMASum,  tgfbSum, mmpSum))
        #
        # distE = (distE / (aEndo + saEndo + pEndo)) if (aEndo + saEndo + pEndo) != 0 else 0
        # # distT = (distT / (aTumor + saTumor + pTumor)) if (aTumor + saTumor + pTumor) != 0 else 0
        # self.f_csv_CellDistance.writerow((mcs, distE, distE))
        #
        #
        # ChemicalValues = []
        # cell = self.attemptFetchingCellById(self.activatedCellId)
        # for i in self.Chemicals:
        #     ChemicalValues.append(cell.sbml.model_M['[' + i + ']'])
        # self.f_csv_ChemicalConcentrations.writerow([mcs]+ChemicalValues)
        #
        # tend = time.perf_counter()


    def finish(self):
        pass
        # self.f_EndothelialCellCount.close()
        # self.f_CellInvasion.close()
        # self.f_CellDistance.close()
        # self.f_ChemicalConcentrations.close()

    def on_stop(self):
        # self.f_EndothelialCellCount.close()
        # self.f_CellInvasion.close()
        # self.f_CellDistance.close()
        # self.f_ChemicalConcentrations.close()
        return


        

        



