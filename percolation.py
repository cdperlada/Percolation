"""
Code for Percolation (W.Kinzel/G.Reents, Physics by Computer)

Simulate percolation and label clusters using the Hoshen-Kopelman algorithm

NOTE: This code takes a really long time to finish running.
"""
import numpy as np
import matplotlib.pyplot as plt

class Percolation(object):

    #Create a lattice
    def lattice(self,L,p): 
        return (np.random.random([L,L])<p)*1
        
    #Count the number and size of percolating clusters
    def mass(self,label): 
        elements, counts = np.unique(label, return_counts = 1)      
        return np.mean(counts[elements!=0])
        
    #Hoshen-Kopelman algorithm
    def labelling(self, n, p): #Create the clusters and label them
        data = self.lattice(n,p)
        label = np.zeros_like(data, dtype = int)

        #Set the label of clustering sites that thouches the edge to -1    
        label[:,0] = data[:,0]*-1
        label[0,:] = data[0,:]*-1
        label[:,n-1] = data[:,n-1]*-1
        label[n-1,:] = data[n-1,:]*-1
        for i in xrange(1,n-1):
            for j in xrange(1,n-1):            
            # occupied
                if data[i,j]: 
                    if data[i-1, j] and not(data[i,j-1]): #up occupied
                        label[i,j] = label[i-1,j]
                    elif data[i, j-1] and not(data[i-1, j]): # left occupied
                        label[i,j] = label[i, j-1]
                        
                    elif not(data[i,j-1]) and not(data[i-1, j]): # none occupied
                        label[i,j] = np.max(label) + 1
                        
                    elif data[i-1, j] and data[i, j-1]: # both occupied
                        if label[i-1,j] < label[i, j-1]:
                            label[i,j] = label[i-1,j] # choose smaller label
                            label[label == label[i, j-1]] = label[i-1,j] 
                            
                        elif label[i-1,j] > label[i, j-1]:
                            label[i,j] = label[i,j-1] # choose smaller label
                            label[label == label[i-1, j]] = label[i,j-1]
                            
                        else:
                            label[i,j] = label[i, j-1] 
                            
                    if i == n-2: # if touching the edge
                        if label[i+1, j] == -1:
                            label[label == label[i,j]] = -1
                            
                    if j == n-2: # if touching the edge
                        if label[i, j+1] == -1:
                            label[label == label[i,j]] = -1
        label[label==-1] = 0
        clustersize = self.mass(label)
        return data, label, clustersize
    
    #percolation threshold
    def thresh(self,n = 50,trials = 10):
        pl = np.arange(0, 1.1, 0.01)
        massl = np.empty_like(pl)
        for w in xrange(0, 100):
            p = pl[w]
            massp = []
            for _ in xrange(trials):
                print 'trial', len(massp)
                print 'p= ', p
                __, __, clustersize = self.labelling(n,p)
                massp.append(clustersize)
            if len(massp) != 0:
                massl[w] = np.mean(massp)
            else:
                massl[w] = 0
        return pl, massl
    
    def labellingplot(self):
        data, label, _ = self.labelling(500,0.59275)
        plt.figure()
        plt.imshow(data, cmap = 'gray')
        plt.title('Percolation at p = %s' %(0.59275))

        plt.figure()
        plt.imshow(label)
        plt.title('Size clusters at p = %s' %(0.59275))


    def Percothreshplot(self):
        pl,massl = self.thresh(n = 500, trials = 10)
        pthreshold = pl[massl == np.max(massl)]
        print 'p_c= ', pthreshold 
        plt.figure()        
        plt.plot(pl, massl,'--ro', lw = 2, ms = 6)
        plt.xlabel('$\mathbf{p}$', fontsize = 12) #Concentration
        plt.ylabel('$\mathbf{M}$', fontsize = 12) #Cluster Mass            
 
        

if __name__ == "__main__":
    Sim = Percolation()
    Sim.labellingplot()
    Sim.Percothreshplot()
        
        



