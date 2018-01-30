import matplotlib.pyplot as plt

def PlotDecays():
	fig = plt.figure(facecolor='white', figsize=([8,6]))
	Temperature   = [-20,      -10,        0,         20,        40,         50]
	PromptNeutron = [-2.39e-2, -1.97e-2,   -1.80e-2,   -1.69e-2,  -1.69e-2,  -1.70e-2]
	DelayNeutron  = [-5.16e-3, -6.24e-3,   -7.17e-3,   -7.55e-3,  -7.83e-3,  -8.13e-3]
	PromptGamma   = [-2.31e-2, -1.76e-2,   -1.68e-2,   -1.65e-2,  -1.61e-2,  -1.48e-2]
	DelayGamma    = [-5.52e-3, -5.76e-3,   -6.31e-3,   -7.49e-3,  -7.78e-3,  -8.13e-3]

	plt.errorbar(Temperature, PromptNeutron,  1e-3, 0, 'r--o', label='Neutron Prompt $\\tau$')
	plt.errorbar(Temperature, PromptGamma,  1e-3, 0, 'g--o', label='Gamma Prompt $\\tau$')
	plt.xlabel('Temperature $\\circ$C')
	plt.ylabel('Prompt decay time constant $\\tau$')
	plt.xlim([-25, 55])
	plt.legend(loc='best')
		
	fig2 = plt.figure(facecolor='white', figsize=([8,6]))
	
	plt.errorbar(Temperature, DelayNeutron, 1e-3, 0,'r--d', label='Neutron Delay $\\tau$')
	plt.errorbar(Temperature, DelayGamma,  1e-3, 0, 'g--d', label='Gamma Delay $\\tau$')
	plt.xlabel('Temperature $\\circ$C')
	plt.ylabel('Delayed decay time constant $\\tau$')
	plt.xlim([-25, 55])
	plt.legend(loc='best')


	plt.show()



PlotDecays()
