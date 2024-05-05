from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

platform = Platform.getPlatformByName('CUDA')

print("Hello Dave...")

#load structure file in PDB format
pdb = PDBFile("step4_equilibration.pdb")
psf = CharmmPsfFile("step3_pbcsetup.psf")

print("loading PDB File...")

crd = CharmmCrdFile("step4_equilibration.crd")
params =CharmmParameterSet("toppar/par_all36_na.prm","toppar/top_all36_na.rtf","toppar/toppar_water_ions.str","toppar/toppar_ions_won.str")

psf.setBox(57*angstroms, 57*angstroms, 57*angstroms)


#specify forcefield CHARMM36, TIP3P water
#forcefield = ForceField('charmm36.xml', 'charmm36/water.xml')

#remove waters from PDB file, add unlisted H
#print("Removing water and adding hydrogens to PDB...")
#modeller = Modeller(psf.topology, crd.positions)	
#modeller.deleteWater()
#residues = modeller.addHydrogens(forcefield)
#

#add solvent
#modeller.addSolvent(forcefield, padding=1.0*nanometer)

#create Langevin (constant T) integrator, combine topology and forcefield into simulation object
system	= psf.createSystem(params, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer, constraints=HBonds)
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(psf.topology, system, integrator,platform)
simulation.context.setPositions(pdb.positions)

print("Minimizing Energy Now...")
simulation.minimizeEnergy()

print("Energy Minimized... Creating output files.")

positions = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, positions,'minimizedStruct.pdb')


#add reporting
simulation.reporters.append(DCDReporter('output.dcd', 1000,enforcePeriodicBox=True))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
        potentialEnergy=True, temperature=True, volume=True))
simulation.reporters.append(StateDataReporter("md_log.txt", 100, step=True,
        potentialEnergy=True, temperature=True, volume=True))

print("Running NVT... see md_log.txt for progress.")
simulation.step(100000)

system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))
simulation.context.reinitialize(preserveState=True)


print("Running NPT... see md_log.txt for progress.")
simulation.step(100000000)



