import numpy as np
from scipy.spatial import Delaunay
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import math
import random
import itertools

porosity=0.3
perm_over_vis =0.8
initial_conc= 0.0
time_step = 0.01
screen_output_freq = 2 #how many time steps between outputs to screen

#Set problem parameters:
V = float(428)
t = float(2.48)
E = float(1.68)
n = input('no. of grid nodes (682000 for optimal mesh):')

#Use parameters to find a, b and c in (x^2/a^2 + y^2/b^2 + z^2/c^2 = 1)
pi = math.pi
a = np.sqrt((6*V)/(4*pi*E*t))
b = np.sqrt((6*V*E)/(4*pi*t))
c = t/2
x, y, z = (2*a), (2*b), (2*c)
xLen, yLen, zLen = x,y,z

#Calculate nodes per unit length:
nodeSpacing = (n/(x*y*z)) **(1./3)

#Set up edge node vectors:
xVector = np.linspace(-x/2, x/2, x*nodeSpacing)
yVector = np.linspace(-y/2, y/2, y*nodeSpacing)
zVector = np.linspace(-z/2, z/2, z*nodeSpacing)

#Use these vectors to make a uniform cuboid grid
nodes = np.vstack(np.meshgrid(xVector,yVector,zVector)).reshape(3,-1).T

#Store nodes in elliptoid
Enodes = np.zeros((n*2, 3))
count = 0
for i in range(len(nodes)):
  if (((nodes[i][0]**2)/a**2) + ((nodes[i][1]**2)/b**2) + ((nodes[i][2]**2)/c**2)) <= 1:
    Enodes[count,:] = nodes[i,:]
    count = count + 1
Enodes.resize(count,3)

#Move surf nodes onto surf:
#Find xy columns, set extreme z nodes in these columns to z value on surf of ellipsoid.
xyList = Enodes[:,[0,1]]
xyListUnique = np.vstack({tuple(row) for row in xyList})
for xyColumn in xyListUnique:
  #List all nodes in each column
  xyNodes = np.where(np.all(xyList == xyColumn, axis = 1))[0]
  if len(xyNodes) > 1:
    #Set last node (highest z value) to the surface
    Enodes[xyNodes[len(xyNodes) - 1],2] = np.sqrt(c**2 * ( 1 - (Enodes[xyNodes[len(xyNodes) - 1],0]**2/a**2) - (Enodes[xyNodes[len(xyNodes) - 1],1]**2/b**2)))
    #Set first node (lowest z value) to the surface
    Enodes[xyNodes[0],2] = -1 * np.sqrt(c**2 * ( 1 - (Enodes[xyNodes[0],0]**2/a**2) - (Enodes[xyNodes[0],1]**2/b**2)))


#Perform Delaunay Triangulation on the nodes:
pyMesh = Delaunay(Enodes)

#Build arrays to pass into openCMISS conversion:
node_coordinates = pyMesh.points
element_nodes_array_raw = pyMesh.simplices


#CHECK ELEMENTS FOR 0 VOLUME:

#Initialise tolerance and loop variables:
tol = 0.00001
index = 0
indexArr = []

for element in element_nodes_array_raw:
  xa = []
  ya = []
  za = []
  for node in element:
    #find coordinates of nodes
    xa.append(node_coordinates[node][0])
    ya.append(node_coordinates[node][1])
    za.append(node_coordinates[node][2])
  #Use coordinates to calculate volume of element
  Vmat = np.vstack((xa,ya,za,[1.0,1.0,1.0,1.0]))
  Volume = (1/6.0) * abs(np.linalg.det(Vmat))
  #update index list of good elements
  if Volume > tol:
    indexArr.append(index)
  index = index+1

#update arrays without 0 volume elements, to pass into openCMISS
element_nodes_array = element_nodes_array_raw[indexArr,:]
for i in range(len(element_nodes_array)):
  element_nodes_array[i] = [x+1 for x in element_nodes_array[i]]
element_array = range(1, len(element_nodes_array)+1)
node_array = range(1, len(node_coordinates)+1)



# Initialise OpenCMISS-Iron
from opencmiss.iron import iron

# Set problem parameters
(coordinateSystemUserNumber,
    regionUserNumber,
    basisUserNumber,
    generatedMeshUserNumber,
    meshUserNumber,
    decompositionUserNumber,
    geometricFieldUserNumber,
    equationsSetFieldUserNumber,
    dependentFieldUserNumber,
    materialFieldUserNumber,
    equationsSetUserNumber,
    problemUserNumber) = range(1,13)

iron.DiagnosticsSetOn(iron.DiagnosticTypes.IN, [1, 2, 3, 4, 5], "Diagnostics",
                     ["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"])

numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber = iron.ComputationalNodeNumberGet()

number_of_dimensions = 3
number_of_mesh_components = 1
total_number_of_elements = len(element_array)
total_number_of_nodes = len(node_array)
mesh_component_number = 1
nodes_per_elem = 4  # for a tet mesh

# Create a RC coordinate system
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.dimension = 3
coordinateSystem.CreateFinish()

# Create a region
region = iron.Region()
region.CreateStart(regionUserNumber,iron.WorldRegion)
region.label = "DarcyRegion"
region.coordinateSystem = coordinateSystem
region.CreateFinish()

# Create a tri-linear simplex basis
basis = iron.Basis()
basis.CreateStart(basisUserNumber)
basis.TypeSet(iron.BasisTypes.SIMPLEX)
basis.numberOfXi = 3
basis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_SIMPLEX]*3
basis.CreateFinish()

# Start the creation of the imported mesh in the region
mesh = iron.Mesh()
mesh.CreateStart(meshUserNumber, region, number_of_dimensions)
mesh.NumberOfComponentsSet(number_of_mesh_components)
mesh.NumberOfElementsSet(total_number_of_elements)

# Define nodes for the mesh
nodes = iron.Nodes()
nodes.CreateStart(region, total_number_of_nodes)

# Refers to nodes by their user number as described in the original mesh
nodes.UserNumbersAllSet(node_array)
nodes.CreateFinish()

elements = iron.MeshElements()
elements.CreateStart(mesh, mesh_component_number, basis)

# Set the nodes pertaining to each element
for idx, elem_num in enumerate(element_array):
    elements.NodesSet(idx + 1, element_nodes_array[idx])

# Refers to elements by their user number as described in the original mesh
elements.UserNumbersAllSet(element_array)
elements.CreateFinish()

mesh.CreateFinish()

# Create a decomposition for the mesh
decomposition = iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber, mesh)
decomposition.type = iron.DecompositionTypes.CALCULATED
decomposition.numberOfDomains = numberOfComputationalNodes
decomposition.CreateFinish()

# Create a field for the geometry
geometricField = iron.Field()
geometricField.CreateStart(geometricFieldUserNumber, region)
geometricField.meshDecomposition = decomposition
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, 1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 3, 1)
geometricField.CreateFinish()

# Update the geometric field parameters
geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

for idx, node_num in enumerate(node_array):
    [x, y, z] = node_coordinates[idx]

    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1,
                                            int(node_num), 1, x)
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1,
                                            int(node_num), 2, y)
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1,
                                            int(node_num), 3, z)

geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

# Create standard Darcy equations set
equationsSetField = iron.Field()
equationsSet = iron.EquationsSet()
equationsSetSpecification = [iron.EquationsSetClasses.FLUID_MECHANICS,
        iron.EquationsSetTypes.DARCY_EQUATION,
        iron.EquationsSetSubtypes.STANDARD_DARCY]
equationsSet.CreateStart(equationsSetUserNumber,region,geometricField,
        equationsSetSpecification,equationsSetFieldUserNumber,equationsSetField)
equationsSet.CreateFinish()

# Create dependent field
dependentField = iron.Field()
equationsSet.DependentCreateStart(dependentFieldUserNumber,dependentField)
dependentField.VariableLabelSet(iron.FieldVariableTypes.U, "Dependent")
dependentField.DOFOrderTypeSet(iron.FieldVariableTypes.U,iron.FieldDOFOrderTypes.SEPARATED)
dependentField.DOFOrderTypeSet(iron.FieldVariableTypes.DELUDELN,iron.FieldDOFOrderTypes.SEPARATED)
equationsSet.DependentCreateFinish()

# Create material field
materialField = iron.Field()
equationsSet.MaterialsCreateStart(materialFieldUserNumber,materialField)
materialField.VariableLabelSet(iron.FieldVariableTypes.U, "Material")
equationsSet.MaterialsCreateFinish()

kFileRaw = open('kData.txt', 'r')
kFile = kFileRaw.readlines()
startLines = range(0,len(kFile))

for i in range(len(kFile)):
  kFile[i] = kFile[i].split()

kList = []

for i in startLines:
  nodeK = []
  nodeK.append(float(kFile[i][0]))
  nodeK.append((float(kFile[i][4])/0.003)) #perm over vis = k/0.003
  nodeK.append(float(kFile[i][5])) #porosity
  kList.append(nodeK)
  print nodeK
kFileRaw.close()

materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,porosity)
materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,perm_over_vis)

for knode in kList:
  materialField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, int(knode[0]),1, knode[2]) # set porosity
  materialField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, int(knode[0]),2, knode[1]) # set perm_over_vis

# Initialise dependent field

dependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,initial_conc)

# Create equations
equations = iron.Equations()
equationsSet.EquationsCreateStart(equations)
equations.sparsityType = iron.EquationsSparsityTypes.SPARSE
equations.outputType = iron.EquationsOutputTypes.NONE
equationsSet.EquationsCreateFinish()

# Create Darcy equation problem
problem = iron.Problem()
problemSpecification = [iron.ProblemClasses.FLUID_MECHANICS,
        iron.ProblemTypes.DARCY_EQUATION,
        iron.ProblemSubtypes.STANDARD_DARCY]
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

# Create boundary conditions
boundaryConditions = iron.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)

#list all x-y couples in ellipsoid
xyList = node_coordinates[:,[0,1]]
#Remove duplicates
xyListUnique = np.vstack({tuple(row) for row in xyList})
#Set blood vessels properties:
arteries = 40
veins = 40
BVCount = arteries + veins
#Randomly generate unique xy positions of blood vessels.
random.seed(500)
bv_y = []
bv_x = np.array(random.sample(np.linspace(-xLen/2,xLen/2, 100000), BVCount))
max_y = np.sqrt(b**2 * (1 - (bv_x**2/a**2)))
for maxY in max_y:
  bv_y.append(random.choice(np.linspace(-maxY,maxY, 100000)))
# get bv_x and bv_y in form of xy:
bv_xy = np.zeros((len(bv_x), 2))
bv_xy[:,0] = bv_x
bv_xy[:,1] = bv_y
bloodVessels = []
for i in range(0,len(bv_xy)): #for each blood vessel:
  #Cycle through xyListUnique to find closest nodes.
  validNodes = []
  for nodeX in xyListUnique:
    xfound = nodeX[0] < (bv_xy[i][0] + nodeSpacing*2) and nodeX[0] > (bv_xy[i][0] - nodeSpacing*2)
    yfound = nodeX[1] < (bv_xy[i][1] + nodeSpacing*2) and nodeX[1] > (bv_xy[i][1] - nodeSpacing*2)
    if  xfound and yfound:
      validNodes.append(nodeX)

  #Now find closest:
  distance = []
  for nodeColumn in validNodes:
    distance.append(math.sqrt((bv_xy[i][0] - nodeColumn[0])**2 + (bv_xy[i][1] - nodeColumn[1])**2))
  closestNode = validNodes[np.argmin(distance)]
  bloodVessels.append(closestNode)
  #search node list for highest z closest node and set to pressure.
  xyNodes = np.where(np.all(xyList == closestNode, axis = 1))[0]
  xyNodes = [x+1 for x in xyNodes]

  if i < arteries:
    artNode = xyNodes[len(xyNodes) - 1]
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,artNode,4,iron.BoundaryConditionsTypes.FIXED,100.0)
  else:
    veinNode = xyNodes[len(xyNodes) - 1]
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,veinNode,4,iron.BoundaryConditionsTypes.FIXED,0.0)


#Set all other surface nodes to 0 velocity:
for pair in xyListUnique:
  bvNode = False
  for bv in bloodVessels:
    if ((pair[0] - bv[0]) < tol)  & ((pair[1] - bv[1]) < tol):
      bvNode = True

  #Set v to 0
  #Find top and bottom node:
  xyNodes = np.where(np.all(xyList == pair, axis = 1))[0]
  xyNodes = [x+1 for x in xyNodes]
  topNode = xyNodes[len(xyNodes) - 1]
  bottomNode = xyNodes[0]

  #Set bottom node to 0 velocity
  boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,bottomNode,3,iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
  boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,bottomNode,2,iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
  boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,bottomNode,1,iron.BoundaryConditionsTypes.FIXED_WALL,0.0)

  #If NOT thin and NOT BV - set top node to 0 velocity
  if (topNode != bottomNode) & (bvNode == False):
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,topNode,1,iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,topNode,2,iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,topNode,3,iron.BoundaryConditionsTypes.FIXED_WALL,0.0)


solverEquations.BoundaryConditionsCreateFinish()
# Solve the problem
problem.Solve()


# Export results
#baseName = "Darcy"
#dataFormat = "PLAIN_TEXT"
#fml = iron.FieldMLIO()
#fml.OutputCreate(mesh, "", baseName, dataFormat)
#fml.OutputAddFieldNoType(baseName+".geometric", dataFormat, geometricField,
#    iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
#fml.OutputAddFieldNoType(baseName+".phi", dataFormat, dependentField,
#    iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
#fml.OutputWrite("DarcyResults.xml")
#fml.Finalise()

# Export results
fields = iron.Fields()
fields.CreateRegion(region)
fields.NodesExport("DarcyResults","FORTRAN")
#fields.ElementsExport("DarcyResults","FORTRAN")
fields.Finalise()


iron.Finalise()