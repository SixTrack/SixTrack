import compareSix2Mad as csm
import latticeConstructor as lc

# A simple, expandable test case.
def exampleTest():
	nbr_turns = 10
	E0_GeV = 15
	mad_init_coords = (0.001, 0.002, 0.003, 0.004, 0.005, 0.0)

	l = lc.Lattice(20)
	l.addMultipoleDef(name='qf', order=2, KN=0.11755705)
	l.addMultipoleDef(name='qd', order=2, KN=-0.11755705)
	l.addRFCavityDef(name='cav', VOLT=100, LAG=0.0, L=0.0, HARMON=100, FREQ=0)
	l.addElement(name='qd', pos=10)
	l.addElement(name='cav', pos=10.001)
	l.addElement(name='qf', pos=19.9999)
	l_s = l.getLatticeDefinition()

	diff = csm.compare('element_test', nbr_turns, E0_GeV, mad_init_coords, printToFile=1,
		lattice=l_s, norm='', verbose=0, all_files=0, trk_element=('cav', 'cav'))

	print(diff)

if __name__ == '__main__':
	exampleTest()
