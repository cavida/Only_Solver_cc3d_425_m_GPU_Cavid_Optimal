from cc3d import CompuCellSetup

from Main import Main
from ECM_Generation import ECM_Generation
from Writing_Data import Writing_Data
from Mitosis import Mitosis
from CaNodule import CaNodule
from SBML import SBML
from Secretion import Secretion

CompuCellSetup.register_steppable(steppable=ECM_Generation(frequency=10))
CompuCellSetup.register_steppable(steppable=Main(frequency=10))
CompuCellSetup.register_steppable(steppable=Mitosis(frequency=1))
CompuCellSetup.register_steppable(steppable=CaNodule(frequency=1))

CompuCellSetup.register_steppable(steppable=SBML(frequency=1))
# Secretion Should come after the SBML solver to give new values for secreting in the fields
CompuCellSetup.register_steppable(steppable=Secretion(frequency=1))

#Writing data should come in the last part to harvest the data and assign cell.dict['LAST_COM'] the locatin
CompuCellSetup.register_steppable(steppable=Writing_Data(frequency=10))


CompuCellSetup.run()
