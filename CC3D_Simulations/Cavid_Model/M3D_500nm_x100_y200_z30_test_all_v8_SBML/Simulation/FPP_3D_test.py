from cc3d import CompuCellSetup

from Main import Main
from ECM_Generation import ECM_Generation
from Writing_Data import Writing_Data
from Mitosis import Mitosis
from CaNodule import CaNodule
from SBML import SBML

CompuCellSetup.register_steppable(steppable=ECM_Generation(frequency=30))
CompuCellSetup.register_steppable(steppable=Main(frequency=1))
CompuCellSetup.register_steppable(steppable=Mitosis(frequency=1))
CompuCellSetup.register_steppable(steppable=CaNodule(frequency=1))

CompuCellSetup.register_steppable(steppable=SBML(frequency=1))

CompuCellSetup.register_steppable(steppable=Writing_Data(frequency=10))

from Main import SecretionSteppable

CompuCellSetup.register_steppable(steppable=SecretionSteppable(frequency=1))

CompuCellSetup.run()
