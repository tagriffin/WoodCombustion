# Methanol combustion

from src.main.utils import Properties_NASA_v5 as NASA
from src.main.utils import IG_Props_Comb
from src.main.utils import A8
from src.main.utils import Combustion_Props_ed14 as Comb

Methanol = Comb.Fuel_1('Methanol',['CH3OH'],[1.0])
