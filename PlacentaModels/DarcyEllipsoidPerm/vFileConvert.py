#Read in .exnode file
exnodeFile = open('DarcyResults.part0.exnode', 'r')
exnode = exnodeFile.readlines()
exnodeFile.close()
nodeStarts = range(23,len(exnode)-9,14)

nodes = []
nodeNum = 1

for i in nodeStarts:
  nodeV = []
  nodeV.append(int(nodeNum))
  nodeV.append(float(exnode[i]))
  nodeV.append(float(exnode[i+1]))
  nodeV.append(float(exnode[i+2]))
  nodes.append(nodeV)
  nodeNum = nodeNum+1

#Write out as file with [node, v1, v2, v3]
vOutput = open('vData.txt', 'w')

for node in nodes:
  line = str(node[0]) + ' ' + str(1) + ' ' + str(node[1]) +'\n'
  vOutput.write(line)
  line = str(node[0]) + ' ' + str(2) + ' ' + str(node[2]) +'\n'
  vOutput.write(line)
  line = str(node[0]) + ' ' + str(3) + ' ' + str(node[3]) +'\n'
  vOutput.write(line)

vOutput.close()
