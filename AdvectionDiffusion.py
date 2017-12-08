#!/usr/bin/env python

# Intialise OpenCMISS-Iron
from opencmiss.iron import iron

# parameters.parse()

# Set problem parameters
height = 1.0
width = 2.0
length = 3.0
diff_coeff = 1.0
initial_conc= 0.5
v_strength=-0.01
start_time = 0.0 
end_time = 1.0
time_step = 0.01
screen_output_freq = 2 #how many time steps between outputs to screen

(coordinateSystemUserNumber,
    regionUserNumber,
    basisUserNumber,
    generatedMeshUserNumber,
    meshUserNumber,
    decompositionUserNumber,
    geometricFieldUserNumber,
    equationsSetFieldUserNumber,
    dependentFieldUserNumber,
    independentFieldUserNumber,
    materialFieldUserNumber,
    equationsSetUserNumber,
    problemUserNumber) = range(1,14)

numberGlobalXElements = 5
numberGlobalYElements = 5
numberGlobalZElements = 5


numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber = iron.ComputationalNodeNumberGet()

# Create a RC coordinate system
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.dimension = 3
coordinateSystem.CreateFinish()

# Create a region
region = iron.Region()
region.CreateStart(regionUserNumber,iron.WorldRegion)
region.label = "DiffusionRegion"
region.coordinateSystem = coordinateSystem
region.CreateFinish()

# Create a tri-linear lagrange basis
basis = iron.Basis()
basis.CreateStart(basisUserNumber)
basis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
basis.numberOfXi = 3
basis.interpolationXi = [
        iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE] * 3
basis.quadratureNumberOfGaussXi = [3] * 3
basis.CreateFinish()

# Create a generated mesh
generatedMesh = iron.GeneratedMesh()
generatedMesh.CreateStart(generatedMeshUserNumber,region)
generatedMesh.type = iron.GeneratedMeshTypes.REGULAR
generatedMesh.basis = [basis]
generatedMesh.extent = [width,height,length]
generatedMesh.numberOfElements = [numberGlobalXElements,numberGlobalYElements,numberGlobalZElements]

mesh = iron.Mesh()
generatedMesh.CreateFinish(meshUserNumber,mesh)

numberOfElements = mesh.numberOfElements
print("number of elements: " + str(numberOfElements))

# Create a decomposition for the mesh
decomposition = iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber,mesh)
decomposition.type = iron.DecompositionTypes.CALCULATED
decomposition.numberOfDomains = numberOfComputationalNodes
decomposition.CreateFinish()

# Create a field for the geometry
geometricField = iron.Field()
geometricField.CreateStart(geometricFieldUserNumber,region)
geometricField.meshDecomposition = decomposition
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
geometricField.CreateFinish()

# Set geometry from the generated mesh
generatedMesh.GeometricParametersCalculate(geometricField)

# Create standard Diffusion equations set
equationsSetField = iron.Field()
equationsSet = iron.EquationsSet()
equationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
        iron.EquationsSetTypes.ADVECTION_DIFFUSION_EQUATION,
        iron.EquationsSetSubtypes.NO_SOURCE_STATIC_ADVEC_DIFF]
equationsSet.CreateStart(equationsSetUserNumber,region,geometricField,
        equationsSetSpecification,equationsSetFieldUserNumber,equationsSetField)
equationsSet.CreateFinish()

# Create dependent field
dependentField = iron.Field()
equationsSet.DependentCreateStart(dependentFieldUserNumber,dependentField)
dependentField.DOFOrderTypeSet(iron.FieldVariableTypes.U,iron.FieldDOFOrderTypes.SEPARATED)
dependentField.DOFOrderTypeSet(iron.FieldVariableTypes.DELUDELN,iron.FieldDOFOrderTypes.SEPARATED)
equationsSet.DependentCreateFinish()

# Create material field
materialField = iron.Field()
equationsSet.MaterialsCreateStart(materialFieldUserNumber,materialField)
equationsSet.MaterialsCreateFinish()
## I believe this will change the diffusion coeff
materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,diff_coeff)
materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,diff_coeff)
materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3,diff_coeff)

# Initialise dependent field
dependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,initial_conc)

# Create independent field
independentField = iron.Field()
equationsSet.IndependentCreateStart(independentFieldUserNumber,independentField)
#independentField.DOFOrderTypeSet(iron.FieldVariableTypes.U,iron.FieldDOFOrderTypes.SEPARATED)
#independentField.DOFOrderTypeSet(iron.FieldVariableTypes.DELUDELN,iron.FieldDOFOrderTypes.SEPARATED)
equationsSet.IndependentCreateFinish()

# Initialise dependent field
#in this mesh y is up
independentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,0.0)
independentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,0.0)
independentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3,v_strength)

# Create equations
equations = iron.Equations()
equationsSet.EquationsCreateStart(equations)
equations.sparsityType = iron.EquationsSparsityTypes.SPARSE
equations.outputType = iron.EquationsOutputTypes.NONE
equationsSet.EquationsCreateFinish()

# Create Diffusion equation problem
problem = iron.Problem()
problemSpecification = [iron.ProblemClasses.CLASSICAL_FIELD,
        iron.ProblemTypes.ADVECTION_DIFFUSION_EQUATION,
        iron.ProblemSubtypes.NO_SOURCE_STATIC_ADVEC_DIFF]
problem.CreateStart(problemUserNumber, problemSpecification)
problem.CreateFinish()

# Create control loops
problem.ControlLoopCreateStart()
problem.ControlLoopCreateFinish()

# Create problem solver
solver = iron.Solver()
problem.SolversCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
solver.outputType = iron.SolverOutputTypes.SOLVER
solver.linearType = iron.LinearSolverTypes.ITERATIVE
solver.linearIterativeAbsoluteTolerance = 1.0E-12
solver.linearIterativeRelativeTolerance = 1.0E-12
problem.SolversCreateFinish()

## Create solver equations and add equations set to solver equations
solver = iron.Solver()
solverEquations = iron.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
solver.SolverEquationsGet(solverEquations)
solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

## Create boundary conditions and set first and last nodes to 0.0 and 1.0
boundaryConditions = iron.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)
firstNodeNumber=1
nodes = iron.Nodes()
region.NodesGet(nodes)
lastNodeNumber = nodes.numberOfNodes
#firstNodeDomain = decomposition.NodeDomainGet(firstNodeNumber,1)

#lastNodeDomain = decomposition.NodeDomainGet(lastNodeNumber,1)

#if firstNodeDomain == computationalNodeNumber:
for i in range(1,37):
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,i,1,iron.BoundaryConditionsTypes.FIXED,0.0)
for i in range(181,187):
#if lastNodeDomain == computationalNodeNumber:
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,i,1,iron.BoundaryConditionsTypes.FIXED,1.0)
for i in range(211,217):
###if lastNodeDomain == computationalNodeNumber:
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,i,1,iron.BoundaryConditionsTypes.FIXED,1.0)
boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,205,1,iron.BoundaryConditionsTypes.FIXED,1.0)
boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,199,1,iron.BoundaryConditionsTypes.FIXED,1.0)
boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,193,1,iron.BoundaryConditionsTypes.FIXED,1.0)
boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,187,1,iron.BoundaryConditionsTypes.FIXED,1.0)

boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,198,1,iron.BoundaryConditionsTypes.FIXED,1.0)
boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,192,1,iron.BoundaryConditionsTypes.FIXED,1.0)
boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,204,1,iron.BoundaryConditionsTypes.FIXED,1.0)
boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,210,1,iron.BoundaryConditionsTypes.FIXED,1.0)
solverEquations.BoundaryConditionsCreateFinish()

## Solve the problem

#iron.DiagnosticsSetOn(1,[1,2,3,4,5],"Diagnostics",["all"])
problem.Solve()
#iron.DiagnosticsSetOff()
# Export results
baseName = "Diffusion"
dataFormat = "PLAIN_TEXT"
fml = iron.FieldMLIO()
fml.OutputCreate(mesh, "", baseName, dataFormat)
fml.OutputAddFieldNoType(baseName+".geometric", dataFormat, geometricField,
    iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
fml.OutputAddFieldNoType(baseName+".phi", dataFormat, dependentField,
    iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
fml.OutputWrite("DiffusionExample.xml")
fml.Finalise()

# Export results
fields = iron.Fields()
fields.CreateRegion(region)
fields.NodesExport("DiffusionResults","FORTRAN")
fields.ElementsExport("DiffusionResults","FORTRAN")
fields.Finalise()

iron.Finalise()
