#!/usr/bin/env python

import placentagen as pg
import numpy as np
import time
import matplotlib.pyplot as plt

####################
#PARAMETERS
####################
  
volume=4.7e5 #mm^3#volume of placenta
thickness=22.1 #mm#thickness of placenta (z-axis dimension)
ellipticity=1.00 #no units#ellipticity of placenta - ratio of y to x axis dimensions
volume_frac = 0.40 # measured volume fraction in histology#Volume fraction at terminals
spacing = 2.00 #mm size of sampling grid
num_int_gens = 3 #number of generations of mature intermediate villi
num_convolutes = 10 #Number of terminal convolutes per intermediate villous
len_int = 1.5 #mm #length of an intermedaie villous
rad_int = 0.03 #mm #radius of an intermediate villous
len_convolute = 3.0 #mm #length of a terminal convolute
rad_convolute = 0.025 #mm #radius of a terminal convolute
start_radius=9.32 #radius of first branch
start_elem=125;#125#this is the number of element where the stem villi start
el_type=2#Put 1 for linear element, 2 for quadrantic element

pi = np.pi
z_radius = thickness / 2.0
x_radius = np.sqrt(volume * 3.0 / (4.0 * pi * ellipticity * z_radius))
y_radius = ellipticity * x_radius


nel_x=int(np.floor((x_radius*2)/3))#number of element in x aixs in pl mesh
nel_y=int(np.floor((y_radius*2)/3))#number of element in y axis in pl mesh
nel_z=int(np.floor((z_radius*2)/3))#number of element in z axis in pl mesh

strahler_order=1.53#strahler order
el_file='full.exelem'
node_file='full.exnode'
stem_file='stem_xy.txt'
velocity=26.5#mm3 per sec inlet and outlet vel of spiral artery and decidual vein
####################################################################################################################################
total_start = time.time()

#Read in exnode and exelem files desrcibing whole placenta tree
nodes_input=pg.import_exnode_tree(node_file)
elems_input=pg.import_exelem_tree(el_file)

#Define element radii by Strahler order
elem_radius = pg.define_radius_by_order(nodes_input['nodes'],elems_input['elems'],'strahler',0,start_radius,strahler_order)
pg.export_exfield_1d_linear(elem_radius, 'arteries', 'radius', 'artery_radius')

#Calculate location of terminal branches
terminal_branches = pg.calc_terminal_branch(nodes_input['nodes'],elems_input['elems'])

#define rectangular mesh
rectangular_mesh = pg.gen_rectangular_mesh(volume, thickness, ellipticity,spacing,spacing,spacing)
pg.export_ex_coords(rectangular_mesh['nodes'],'samplingmesh','sampling_mesh','exnode')
pg.export_exelem_3d_linear(rectangular_mesh['elems'],'samplingmesh','sampling_mesh')

start = time.time()
#define the placental volume in each element of the sampling grid
placental_volume= pg.ellipse_volume_to_grid(rectangular_mesh, volume, thickness, ellipticity,20)
pg.export_exelem_3d_linear_list(rectangular_mesh['elems'], placental_volume['non_empty_rects'],'samplingmesh','sm_nonempty')
end = time.time()

#Locate terminals in sampling grid
start = time.time()
terminals_in_grid = pg.terminals_in_sampling_grid_fast(rectangular_mesh,terminal_branches,nodes_input['nodes'])
end = time.time()
 
print ('Total time for terminals in grid = '+ str((end-start)/60.0) + ' mins')
pg.export_exfield_3d_linear(terminals_in_grid['terminals_in_grid'],'samplingmesh','terminal_no','terminal_no')
start = time.time()

#Define volume of tissue appended to a terminal element

term_vol = pg.terminal_villous_volume(num_int_gens,num_convolutes,len_int,rad_int,len_convolute,rad_convolute)
term_diam = pg.terminal_villous_diameter(num_int_gens,num_convolutes,len_int,rad_int,len_convolute,rad_convolute)
 
#define volume of space occupied by a terminal element
volume_per_term = term_vol/volume_frac
term_sampling = pg.terminal_volume_to_grid(rectangular_mesh,terminal_branches,nodes_input['nodes'],volume,thickness,ellipticity,volume_per_term,term_vol,term_diam)

end = time.time()
print ('Total time for placental volume = '+ str((end-start)/60.0) + ' mins')
pg.export_exfield_3d_linear(placental_volume['pl_vol_in_grid'],'samplingmesh','volume','volonmesh')

#map branches to the sampling grid

start = time.time()
br_vol_samp_gr=pg.cal_br_vol_samp_grid(rectangular_mesh,nodes_input['nodes'],elems_input['elems'],elem_radius,volume,thickness,ellipticity,start_elem)

end = time.time()
print ('Total time for branch volume = '+ str((end-start)/60.0) + ' mins')
pg.export_exfield_3d_linear(br_vol_samp_gr['br_vol_in_grid'],'samplingmesh','volume','volo_br_mesh')
total_end = time.time()

#sum tissue volume in each element of the sampling grid including branches and terminals
tissue_vol = pg.tissue_vol_in_samp_gr(term_sampling['term_vol_in_grid'],br_vol_samp_gr['br_vol_in_grid'])
vol_frac = pg.vol_frac_in_samp_gr(tissue_vol,placental_volume)
weighted_diameter = pg.weighted_diameter_in_samp_gr(term_sampling['term_diameter_in_grid'],br_vol_samp_gr['br_diameter_in_grid'],tissue_vol)
conductivity = pg.conductivity_samp_gr(vol_frac,weighted_diameter,placental_volume['non_empty_rects'])
porosity = pg.porosity(vol_frac)
 
pg.export_exfield_3d_linear_list(vol_frac,placental_volume['non_empty_rects'],'samplingmesh','volumefrac','volumefrac')
pg.export_exfield_3d_linear_list(conductivity,placental_volume['non_empty_rects'],'samplingmesh','conductivity','conductivity')

pl_mesh=pg.gen_placental_mesh(nel_x,nel_y,nel_z,volume,thickness,ellipticity,el_type)
 
pg.export_ex_coords(pl_mesh['placental_node_coor'],'placenta','placenta_mesh','exnode')
#pg.export_exelem_3d_linear(pl_mesh['placental_el_con'],'placenta','placenta_mesh')#Use this for linear
pg.export_exelem_3d_quadrantic(pl_mesh['placental_el_con'],'placenta','placenta_mesh')#Use this for quadrantic
 
pl_node_in_grid=pg.darcynode_in_sampling_grid(rectangular_mesh, pl_mesh['placental_node_coor'])
mapped_con_por=pg.mapping_darcy_sampl_gr(pl_node_in_grid,placental_volume['non_empty_rects'],conductivity,porosity)

#FROM THIS ONWARDS, it applies for quardrantic element only
#######################################GETTING THE SURFACE NODE THAT NEED TO DEFINE IN DARCY SOLVER
#For left and right surface
nnod_x =  int((nel_x*2)+1)
nnod_y =  int((nel_y*2)+1)
nnod_z =  int((nel_z*2)+1)

sIEN=np.zeros((9,nel_y*nel_z),dtype=int)
e=0
for k in range( 1,nnod_x*nnod_y*(nnod_z-1),(nnod_x*nnod_y)*2):#go up
        for j in range(  1,nnod_x*(nnod_y-1),2*nnod_x):#go left         
          
            sIEN[0,e] = j+(k-1) #1st node
            sIEN[1,e] = sIEN[0,e]+(nnod_x)*(nnod_y)#2nd node
            sIEN[2,e] = sIEN[1,e]+(nnod_x)*(nnod_y)#3rd node
            sIEN[3,e] = sIEN[0,e]+nnod_x#4th node
            sIEN[4,e] = sIEN[1,e]+nnod_x#5th node
            sIEN[5,e] = sIEN[2,e]+nnod_x#6th node
            sIEN[6,e] = sIEN[3,e]+nnod_x#7th node
            sIEN[7,e]=sIEN[4,e]+nnod_x#8th node
            sIEN[8,e]=sIEN[5,e]+nnod_x#9th node
            e = e+1            
   
        left=np.unique(sIEN)
        right=np.unique(sIEN.T+(nnod_x-1))

#For front and back surface

sIEN=np.zeros((9,nel_x*nel_z),dtype=int)
e=0
for k in range (1,nnod_x*nnod_y*(nnod_z-2),(nnod_x*nnod_y)*2):#go up
        for i in range( 1,nnod_x-1,2):#go right            
            sIEN[0,e] = i+(k-1)           
            sIEN[1,e] = sIEN[0,e]+1
            sIEN[2,e] = sIEN[0,e]+2                    
            sIEN[3,e]=sIEN[0,e]+(nnod_x*nnod_y)
            sIEN[4,e] = sIEN[3,e]+1
            sIEN[5,e] = sIEN[3,e]+2
            sIEN[6,e] = sIEN[3,e]+(nnod_x*nnod_y)
            sIEN[7,e] = sIEN[6,e]+1    
            sIEN[8,e] = sIEN[6,e]+2            
            e = e+1      
    
        front=np.unique(sIEN)
        back=np.unique(sIEN.T+(nnod_x*(nnod_y-1)))

#For top and bottom surface
sIEN=np.zeros((9,nel_x*nel_y),dtype=int)
e=0
for j in range( 1,nnod_x*(nnod_y-1),nnod_x*2):#go up
        for i in range (  1,nnod_x-1,2):#go back                                   
            sIEN[0,e] = i+(j-1)
            sIEN[1,e] = sIEN[0,e]+1
            sIEN[2,e] = sIEN[0,e]+2
            sIEN[3,e] = sIEN[0,e]+nnod_x          
            sIEN[4,e] = sIEN[3,e]+1
            sIEN[5,e] = sIEN[3,e]+2
            sIEN[6,e] = sIEN[3,e]+nnod_x          
            sIEN[7,e] = sIEN[6,e]+1
            sIEN[8,e] = sIEN[6,e]+2
            e = e+1

        bottom=np.unique(sIEN)    
        top=np.unique(sIEN.T+(nnod_x*nnod_y)*(nnod_z-1))

surfacenode=np.hstack((front,back,left,right,bottom,top))
surfacenode=np.unique(surfacenode)#collection of surface node

######################################################GETTING VESSEL NODE THAT MAPPED WITH STEM VESSEL

ellipsoid_coor=pl_mesh['placental_node_coor']
xyList=np.zeros((len(surfacenode),2))
count=0
for i in range (0,len(surfacenode)):#taking only x and y coordinates
 if ellipsoid_coor[surfacenode[i]-1,2]>0:#take if uppersurface
  xyList[count,0]=ellipsoid_coor[surfacenode[i]-1,0]
  xyList[count,1]=ellipsoid_coor[surfacenode[i]-1,1]
  count=count+1

xyList=xyList[0:count,:]

vesselnode_temp = np.vstack({tuple(row) for row in xyList})#Remove duplicates#mignt not need this one

#reading in the stem vessel to map the spiral artery location
stem_xy = open(stem_file, 'r')
stem_coor = stem_xy.readlines()#readlines
startLines = range(0,len(stem_coor))

for i in range(len(stem_coor)):
     stem_coor[i] = stem_coor[i].split()
stem_xyList = []
for i in startLines:
     node = []
     node.append(float(stem_coor[i][0]))#x coor of stem villi
     node.append((float(stem_coor[i][1]))) #y coor of stem villi     
     stem_xyList.append(node)
stem_xy.close()

vessel_mapped_stem=stem_xyList#this is the x,y location where we want to put spiral artery
spiral_array=np.zeros((len(vessel_mapped_stem)),dtype=int)#store the node nuber of spiral artery
decidual_array=np.zeros((len(vessel_mapped_stem)),dtype=int)#store the node number of decidual vein

check=ellipsoid_coor[:,0:2]  
 
for i in range(0,len(vessel_mapped_stem)): #for each blood vessel,Cycle through to find closest nodes
     
     distance = []
     for nodeX in vesselnode_temp:
       distance.append(np.sqrt((vessel_mapped_stem[i][0] - nodeX[0])**2 + (vessel_mapped_stem[i][1] - nodeX[1])**2))#distance from the nodes
         
     A=sorted(distance)[0]#taking the nearest node
     V=np.random.choice(sorted(distance)[1:])#choosing random 
          
     arterynode= vesselnode_temp[np.where(distance==A)]#this is artery node coor        
     veinnode=vesselnode_temp[np.where(distance==V)]

     xyNodes_A = np.where(np.all(check == arterynode[0], axis = 1))[0]#collection of nodes number that have same x,y of close nodes (multiple layer)
     xyNodes_A = [x+1 for x in xyNodes_A]#adding 1 to get the right node number
     
     xyNodes_V = np.where(np.all(check == veinnode[0], axis = 1))[0]#collection of nodes number that have x,y of close nodes (multiple layer)
     xyNodes_V = [x+1 for x in xyNodes_V]#adding 1 to get the right node number

     spiral_array[i] = xyNodes_A[len(xyNodes_A) - 1]#just taking the last node (cos we want from top surface)
     decidual_array[i] = xyNodes_V[len(xyNodes_V) - 1]#just taking the last node (cos we want from top surface)
     #remove the taken nodes so that won't repeat
     vesselnode_temp=np.delete(vesselnode_temp,np.where(np.all(arterynode[0] ==vesselnode_temp, axis=1)),0)#remove taken artery node             
     vesselnode_temp=np.delete(vesselnode_temp,np.where(np.all(veinnode [0]==vesselnode_temp, axis=1)),0)#remove taken vein node

vesselnode=np.hstack((spiral_array,decidual_array))#array of vessel nodes


############################calculating the velocity vx vy and vz from velocity normal to the surface

spiral_coor= ellipsoid_coor[spiral_array-1]#coor of spiral artery node
decidual_coor= ellipsoid_coor[decidual_array-1]#coor of decidual vein node

spiral_vxyz=np.zeros((len(spiral_array),3))
for j in range (0,len(spiral_array)):
 x_coor=spiral_coor[j,0]
 y_coor=spiral_coor[j,1]
 z_coor=spiral_coor[j,2]

 #normal to the surface
 n_x=(2*x_coor)/x_radius**2#nx 
 n_y=(2*y_coor)/y_radius**2#ny
 n_z=(2*z_coor)/z_radius**2#nz
 #to avoid singular matrix
 if n_x==0:
    n_x==0.00001
 if n_x==-0:
    n_x=-0.00001
 if n_y==0:
    n_y==0.00001
 if n_y==-0:
    n_y=-0.00001
 if n_z==0:
    n_z==0.00001
 if n_z==-0:
    n_z=-0.00001
 #tangent 1
 t1_x=-x_coor
 t1_y=-y_coor
 t1_z_comp1=(x_coor**2/x_radius**2)+(y_coor**2/y_radius**2)+(z_coor**2/z_radius**2)
 t1_z_comp2=z_radius**2/z_coor
 t1_z=(t1_z_comp2*(t1_z_comp1))-z_coor
 #to avoid singular matrix
 if t1_x==0:
    t1_x==0.00001
 if t1_x==-0:
    t1_x=-0.00001
 if t1_y==0:
    t1_y==0.00001
 if t1_y==-0:
    t1_y=-0.00001
 if t1_z==0:
    t1_z==0.00001
 if t1_z==-0:
    t1_z=-0.00001
 #tangent 2
 t2_x_comp1=t1_z_comp1
 t2_x_comp2=x_radius**2/x_coor
 t2_x=(t2_x_comp2*(t2_x_comp1))-x_coor
 t2_y=-y_coor
 t2_z=-z_coor
 #3 equations for 3 unknown (vx,vy,vz)
 A=np.array([[n_x,n_y,n_z],[t1_x,t1_y,t1_z],[t2_x,t2_y,t2_z]])
 B=np.array([-velocity,0,0])
 Vxyz = np.linalg.solve(A, B)#solving equations

 spiral_vxyz[j,0]=Vxyz[0]#velocity x
 spiral_vxyz[j,1]=Vxyz[1]#velocity y
 spiral_vxyz[j,2]=Vxyz[2]#velocity z

#same for decidual vein velocity x,y,z
decidual_vxyz=np.zeros((len(decidual_array),3))
for j in range (0,len(decidual_array)):
 x_coor=decidual_coor[j,0]
 y_coor=decidual_coor[j,1]
 z_coor=decidual_coor[j,2]

 n_x=(2*x_coor)/x_radius**2
 n_y=(2*y_coor)/y_radius**2
 n_z=(2*z_coor)/z_radius**2
 if n_x==0:
    n_x==0.00001
 if n_x==-0:
    n_x=-0.00001
 if n_y==0:
    n_y==0.00001
 if n_y==-0:
    n_y=-0.00001
 if n_z==0:
    n_z==0.00001
 if n_z==-0:
    n_z=-0.00001
 t1_x=-x_coor
 t1_y=-y_coor
 t1_z_comp1=(x_coor**2/x_radius**2)+(y_coor**2/y_radius**2)+(z_coor**2/z_radius**2)
 t1_z_comp2=z_radius**2/z_coor
 t1_z=(t1_z_comp2*(t1_z_comp1))-z_coor
 if t1_x==0:
    t1_x==0.00001
 if t1_x==-0:
    t1_x=-0.00001
 if t1_y==0:
    t1_y==0.00001
 if t1_y==-0:
    t1_y=-0.00001
 if t1_z==0:
    t1_z==0.00001
 if t1_z==-0:
    t1_z=-0.00001
 t2_x_comp1=t1_z_comp1
 t2_x_comp2=x_radius**2/x_coor
 t2_x=(t2_x_comp2*(t2_x_comp1))-x_coor
 t2_y=-y_coor
 t2_z=-z_coor

 A=np.array([[n_x,n_y,n_z],[t1_x,t1_y,t1_z],[t2_x,t2_y,t2_z]])
 B=np.array([velocity,0,0])
 Vxyz = np.linalg.solve(A, B)
 
 decidual_vxyz[j,0]=Vxyz[0]
 decidual_vxyz[j,1]=Vxyz[1]
 decidual_vxyz[j,2]=Vxyz[2]


###########################################surfacenod excluding vessel node

index=[]
for i in range (0,len(vesselnode)):
  idx= np.where(surfacenode==vesselnode[i])
  index.append(idx)

surfnode_ex_vessel = np.delete(surfacenode, index)#array of surface node excluding vessel node


#########################################SOLVINGDARCY################################################################################################
import os

# Intialise OpenCMISS-Iron
from opencmiss.iron import iron

#parameters.parse()

#-----------------------------------------------------------------------------------------------------------
# SET PROBLEM PARAMETERS
#-----------------------------------------------------------------------------------------------------------


porosity = 0.3
perm_over_vis = 0.8
initial_conc = 0.0


element_array=pl_mesh['element_array']
node_array=pl_mesh['node_array']

element_nodes_array=pl_mesh['placental_el_con'][:,1:28].astype(np.int32)
node_coordinates=np.copy(pl_mesh['placental_node_coor'])
for i in range(len(element_nodes_array)):
  element_nodes_array[i] = [x+1 for x in element_nodes_array[i]]


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
nodes_per_elem = 27  # 

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
basis.TypeSet(iron.BasisTypes.LAGRANGE_HERMITE_TP)
basis.numberOfXi = 3
basis.interpolationXi = [iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE]*3
basis.quadratureLocalFaceGaussEvaluate = True
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


materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,0.8)
materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,62.5)

for knode in range (0,len(node_coordinates)):
  
  materialField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, int(knode+1),1, mapped_con_por[knode,2]) # set porosity
  materialField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, int(knode+1),2, mapped_con_por[knode,1]/0.003) # set perm_over_vis

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

#INLET BC
for i in range (0,len(spiral_array)):
 boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,spiral_array[i],1,iron.BoundaryConditionsTypes.FIXED,spiral_vxyz[i,0])
 boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,spiral_array[i],2,iron.BoundaryConditionsTypes.FIXED,spiral_vxyz[i,1])
 boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,spiral_array[i],3,iron.BoundaryConditionsTypes.FIXED,spiral_vxyz[i,2])

#OUTLET BC
for i in range (0,len(decidual_array)):
 boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,decidual_array[i],1,iron.BoundaryConditionsTypes.FIXED,decidual_vxyz[i,0])
 boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,decidual_array[i],2,iron.BoundaryConditionsTypes.FIXED,decidual_vxyz[i,1])
 boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,decidual_array[i],3,iron.BoundaryConditionsTypes.FIXED,decidual_vxyz[i,2])
 

#set the rest surface nodes with velocity boundary conditon (Vx,Vy,Vz as zero)

for i in range (0,len(surfnode_ex_vessel)):
   boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,surfnode_ex_vessel[i],3,iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
   boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,surfnode_ex_vessel[i],2,iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
   boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,surfnode_ex_vessel[i],1,iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
  
'''
#TOP
for i in range (0,len(top)):
 if not top[i] in vesselnode:
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,top[i],3,iron.BoundaryConditionsTypes.FIXED,0.0)
    #boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,top[i],1,iron.BoundaryConditionsTypes.FIXED,0.0)
    #boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,top[i],2,iron.BoundaryConditionsTypes.FIXED,0.0)
#BOTTOM
for i in range (0,len(bottom)):
 
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,bottom[i],3,iron.BoundaryConditionsTypes.FIXED,0.0)
    #boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,bottom[i],1,iron.BoundaryConditionsTypes.FIXED,0.0)
    #boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,bottom[i],2,iron.BoundaryConditionsTypes.FIXED,0.0)

#FRONT
for i in range (0,len(front)):
 if not front[i] in vesselnode:#???????????????????????
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,front[i],3,iron.BoundaryConditionsTypes.FIXED,0.0)
    #boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,front[i],1,iron.BoundaryConditionsTypes.FIXED,0.0)
    #boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,front[i],2,iron.BoundaryConditionsTypes.FIXED,0.0)
#BACK

for i in range (0,len(back)):
 if not back[i] in vesselnode:#???????????????????????
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,back[i],3,iron.BoundaryConditionsTypes.FIXED,0.0)
    #boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,back[i],1,iron.BoundaryConditionsTypes.FIXED,0.0)
    #boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,back[i],2,iron.BoundaryConditionsTypes.FIXED,0.0)
#left

for i in range (0,len(left)):
 if not left[i] in vesselnode:#???????????????????????
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,left[i],3,iron.BoundaryConditionsTypes.FIXED,0.0)
    #boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,left[i],1,iron.BoundaryConditionsTypes.FIXED,0.0)
    #boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,left[i],2,iron.BoundaryConditionsTypes.FIXED,0.0)

#RIGHT

for i in range (0,len(right)):
 if not right[i] in vesselnode:#???????????????????????
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,right[i],3,iron.BoundaryConditionsTypes.FIXED,0.0)
   # boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,right[i],1,iron.BoundaryConditionsTypes.FIXED,0.0)
    #boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,right[i],2,iron.BoundaryConditionsTypes.FIXED,0.0)
'''
solverEquations.BoundaryConditionsCreateFinish()
# Solve the problem
problem.Solve()


if not os.path.exists('./output'):
    os.makedirs('./output')

# Export results
fields = iron.Fields()
fields.CreateRegion(region)
fields.NodesExport("output/StaticDarcy","FORTRAN")
fields.ElementsExport("output/StaticDarcy","FORTRAN")
fields.Finalise()


iron.Finalise()
raise SystemExit
####################PLOTTING
plt.figure()
#number of terminals per grid element
plt.hist(terminals_in_grid['terminals_in_grid'][np.nonzero(placental_volume['pl_vol_in_grid'])],bins=20)
plt.savefig("test_histo1.eps", dpi=150)

#branch volume in sampliing grid
plt.figure()
plt.xscale('log')
plt.hist(br_vol_samp_gr['br_vol_in_grid'][np.nonzero(placental_volume['pl_vol_in_grid'])], bins = np.logspace(-6,2,50))
plt.savefig("test_histo2.eps", dpi=150)


#terminal volume in sampling grid
plt.figure()
plt.xscale('log')
plt.hist(term_sampling['term_vol_in_grid'][np.nonzero(term_sampling['term_vol_in_grid'])],bins = np.logspace(-6,2,50))
plt.savefig("test_histo3.eps", dpi=150)

#total tissue volumevolume in sampling grid
plt.figure()
plt.hist(tissue_vol[np.nonzero(tissue_vol)],bins = np.linspace(0,10,50))
plt.savefig("test_histo4.eps", dpi=150)

plt.figure()
plt.hist(vol_frac[np.nonzero(placental_volume['pl_vol_in_grid'])], bins = np.linspace(0,1,50))
plt.savefig("test_histo5.eps", dpi=150)
 

plt.figure()
plt.hist(weighted_diameter[np.nonzero(placental_volume['pl_vol_in_grid'])], bins = 50)
plt.savefig("test_histo6.eps", dpi=150)
 

plt.figure()
plt.xscale('log')
plt.hist(conductivity[np.nonzero(placental_volume['pl_vol_in_grid'])], bins = np.logspace(-8,2,50))
plt.savefig("test_histo7.eps", dpi=150)

print ('Total time for script= '+ str((total_end-total_start)/60.0) + ' mins')
