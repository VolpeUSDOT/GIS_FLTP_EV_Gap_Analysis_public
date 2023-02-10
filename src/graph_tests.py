import pandas as pd
import numpy as np
import networkx as nx
import math
from rtree import index
import math
import matplotlib.pyplot as plt
from scipy import spatial
from scipy.sparse import coo_matrix
from scipy.sparse import csgraph
import geopandas as gpd

class graphFromLinesAndPoints(object):
    
    def __init__(self,indexAndPoints):
        self.nodeIDPoints = indexAndPoints
        self.G = None#Graph()#nx.Graph()
        self.K  = 10
        self.epsilon = .001
        self.allNodesMatch = {}
        self.allNodesMatchWithLines = {}
        self.rtIndex = None
        #self.totalNetworkDistance=0
        self.weights = []
        self.nodesTemporary = []
        self.edgesTemporary = []
        self.pntIndexToNodeID = {}
        self.NodeIDtopntIndex = {}
        self.pntLst = []
        self.lowerDM = None
        self.baseforsparse = []
        self.pntMat = None
        
    def distance(self,x1,y1,x2,y2):
        dx = x1-x2
        dy = y1-y2
        return math.sqrt(dx*dx+dy*dy)
    
    def distanceNode(self,node1,node2):
        dx = node1[0]-node2[0]
        dy = node1[1]-node2[1]
        return math.sqrt(dx*dx+dy*dy)
    
    #need to contend with duplicate points...
    def developGraphFromTriangulation(self):
        #self.edgesTemporary = []
        self.pntLst = [i for i in range(0,len(self.nodeIDPoints.keys()))]
        self.pntIndexToNodeID = {}
        self.NodeIDtopntIndex = {}
        print ("Remove Duplicates")
        for i,k in enumerate(self.nodeIDPoints.keys()):
            self.pntLst[i] = self.nodeIDPoints[k]
            
        
            
        XC = np.array(self.pntLst)
        self.pointMat = XC
        print ("creating triangulation")
        tri = spatial.Delaunay(XC)
        data_l = []
        col_l = []
        row_l = []
        existingPairs = []
        ids = []
        for t in tri.simplices:
            
            ds = spatial.distance.pdist(np.vstack([XC[t[0]],XC[t[1]],XC[t[2]]]),'euclidean')
            ids.extend(t)
            self.weights.append(ds[0])
            self.weights.append(ds[1])
            self.weights.append(ds[2])
            pair = (t[0],t[1])
            if not pair in existingPairs:
                row_l.append(t[0])
                col_l.append(t[1])
                data_l.append(ds[0])
                existingPairs.append(pair)
            pair = (t[0],t[2])
            if not pair in existingPairs:
                row_l.append(t[0])
                col_l.append(t[2])
                data_l.append(ds[1])
                existingPairs.append(pair)
            pair = (t[1],t[2])
            if not pair in existingPairs:
                row_l.append(t[1])
                col_l.append(t[2])
                data_l.append(ds[2])
                existingPairs.append(pair)
        self.baseforsparse = zip(row_l,col_l,data_l)
        

        print ("creating spatial index from unique points")
        triids = list(set(ids))
        self.rtIndex = spatial.KDTree(XC[triids])
        print ("match points to their index")
        for i,k in enumerate(self.nodeIDPoints.keys()):
            pnt = np.array(self.nodeIDPoints[k])
            trq = self.rtIndex.query(pnt,k=1)
            closestID = triids[trq[1]]
            self.pntIndexToNodeID[closestID]=k
            self.NodeIDtopntIndex[k]=closestID

        print ("Sparse Matrix")
        self.G = coo_matrix((np.array(data_l),(np.array(row_l),np.array(col_l))),shape=(len(XC),len(XC)))
        print ("Lower Distance Matrix")
        #self.lowerDM = np.tril(csgraph.shortest_path(self.G,directed=False))
        self.lowerDM = np.tril(csgraph.dijkstra(self.G,directed=False))                     

        ##self.G.add_nodes_from(self.nodeIDPoints.keys())
        #print "creating graph"
        #for t in tri.simplices:
            #k1= self.pntLstNodeCheck[t[0]]
            #k2=self.pntLstNodeCheck[t[1]]
            #k3=self.pntLstNodeCheck[t[2]]
            #node1 = self.nodeIDPoints[k1]
            #node2 = self.nodeIDPoints[k2]
            #node3 = self.nodeIDPoints[k3]
            #self.edgesTemporary.append([k1,k2,self.distanceNode(node1,node2)])
            #self.edgesTemporary.append([k1,k2,self.distanceNode(node1,node2)])
            #self.weights.append(self.distanceNode(node1,node2))
            #self.edgesTemporary.append([k1,k3,self.distanceNode(node1,node3)])
            #self.edgesTemporary.append([k1,k3])
            #self.weights.append(self.distanceNode(node1,node3))
            #self.edgesTemporary.append([k2,k3,self.distanceNode(node2,node3)])
            #self.edgesTemporary.append([k2,k3])
            # self.weights.append(self.distanceNode(node2,node3))
            
        #self.G.add_weighted_edges_from(self.edgesTemporary,weight='distance')
        #self.G.add_vertices(self.nodeIDPoints.keys())
        #self.G.add_edges(self.edgesTemporary)
        #self.G.es["distance"] = self.weights

        #if len(self.pntLst) <=1000:
            #nx.draw(self.G,pos=self.nodeIDPoints)
            #layout = self.G.layout("fr")
            #plot(self.G, layout = layout)
        


    
    #def buildGraphDistanceMatrix(self,K):
        #get the nearest neighbors list
        #print "find cutoff"
        #dist,nkn = self.rtIndex.query(self.pntLst,K*2)
        #print "create dm"
        #nodeIndex = {v.index:v['name']}
        #self.distanceIndexNodeID = {v['name']:v.index}
    
        #for i in self.pntLstNodeCheck.keys():
            #maxDist = max(dist[i])
            #sourceNode = self.pntLstNodeCheck[i]
            #self.sparseDistanceMatrixNodeID[sourceNode] = nx.single_source_dijkstra_path_length(self.G,sourceNode,weight='distance',cutoff=maxDist)
        #print "looking vertices"
        #for v in self.G.vs:
            #nodeID = v['name']
            #nodeDM = self.G.shortest_paths()
            #if self.numpyDistanceMatrix == None:
                #self.numpyDistanceMatrix = np.array(nodeDM[0])
            #else:
                #self.numpyDistanceMatrix=np.vstack([self.numpyDistanceMatrix,nodeDM[0]])
            #del nodeDM
    def getNodeDistances(self,nodeID1,nodeID2):
        ind1 = self.NodeIDtopntIndex[nodeID1]
        ind2 = self.NodeIDtopntIndex[nodeID2]
        if ind1 > ind2:
            dist = self.lowerDM[ind1][ind2]
        else:
            dist = self.lowerDM[ind2][ind1]
        if dist == 0 or dist<0 or dist==np.Inf:
            return None
        else:
            return dist
    
        return None
    
    
    def getMaxKNearest(self,nodeID,K):
        #try:
        #stuff = []
        #keys = set(self.p.keys())
        i = self.NodeIDtopntIndex[nodeID]
        #print "in node id %s"%nodeID
        #print "pnt index %s"%i
        #stuff = []
        #for j in self.pntIndexToNodeID.keys():
            #if j!=i:
                #if i > j:
                    #dist = self.lowerDM[i][j]
                #else:
                    #dist = self.lowerDM[j][i]
                #if dist != 0 or dist!=np.Inf:
                    #stuff.append(dist)
        stuff = [x for x in self.lowerDM[i+1:,i] if x != np.inf]+[x for x in self.lowerDM[i,0:i] if x != np.inf]
        #print stuff
        stuff = sorted(stuff)
        #stuff = sorted([x for x in self.lowerDM[i+1:,i]])#[self.p[nodeID][k] for k in self.p[nodeID].keys() if k in keys and k!=nodeID]
        #sorts = sorted(stuff)
        if max(stuff[0:K]) == 0:
            return 1.0

        return max(stuff[0:K])
       #except Exception as e:
        #    print "error"
        #    print nodeID
        #    print str(e)
        #    return 1.0
    
    def getMaxKNearDict(self,keys, K):
        maxLst = {}
        for n in keys:
            maxLst[n] = self.getMaxKNearest(n,K)
        return maxLst