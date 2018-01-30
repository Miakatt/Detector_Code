import matplotlib.pyplot as plt
import numpy as np


Volts = range(20,30)
print Volts
I_minus11C = [1.3e-7, 1.3e-7, 1.3e-7, 1.7e-7, 1.7e-7, 3.2e-6, 10.4e-6, 20e-4, 31.4e-4, 45e-4]
print len(Volts), len(I_minus11C)

xerr = 0 
yerr = 1e-6
fig = plt.figure(facecolor='white', figsize=(8,6))
plt.errorbar(Volts, I_minus11C, yerr, xerr,'r-o',  label='Std Output')

plt.legend(loc='best')
plt.xlabel("Voltage (V)")
plt.ylabel("Current (uA)")
plt.xlim(min(Volts)-1, max(Volts)+1)
plt.draw()
plt.pause(0.001)
#if (CLLB):
#	plt.savefig('CLLB_IV_Curve.pdf', fmt='pdf')
#elif (CLLBC):
#	plt.savefig('CLLBC_IV_Curve.pdf', fmt='pdf')
plt.show()