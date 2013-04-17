
from uhg import *
import uhg

use("/home/qiu/Downloads/Pb_0_5_Glb_collected.db")

e("dN/dydpT(0.5)(pion_p_hydro)")

e(" v_2[2]( linspace(0,2,30) )(pion_p_hydro) ")

sqrt(mean(abs(uhg._storedEbeDBReader.get_diff_V_n(particleName="pion_p_hydro", order=2, pTs=linspace(0,2,30))**2,0)))


