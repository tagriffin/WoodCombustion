# Benzin combustion

import Combustion_Props_ed14 as Comb
import A8

Benzin = Comb.Fuel_2('fuel',['C',0.837, 'H', 0.143, 'S', 0.0, 'O', 0.02],42.9, 0.0)

# Energy Balance Example; air and methane enter at 25 °C, lambda = 1.2
# Exhaust gas at 150 °C: determine the heat lost per kg fuel (1% heat loss)

Q1 = (Benzin.H_out(298.15, 413.15,1.2) - 0.99*Benzin.H_in(288.15,1.2))/1000
print('Heat lost per kg fuel is %.2f MJ/kg' %(-Q1))

# calculating with the average heat capacities
Q2 = 0.99*Benzin.Hu_dry*1000 - (Benzin.lmin() + 1)*Benzin.Cp_ave_stExhaust(0,150)*(125)- \
    0.2*Benzin.lmin() * A8.Cp_ave_Air(150)*(125)
print('Heat lost per kg fuel is %.2f MJ/kg' %(Q2/1000))

print('Difference between both is %.2f percent' %(abs(Q2/1000-(-Q1))/abs(Q1) * 100))
