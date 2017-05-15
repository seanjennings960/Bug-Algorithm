import numpy as np
import matplotlib.pyplot as plt
from distanceFunctions import *

x1 = 3
y1 = 4
x2 = 5
y2 = 7
q = np.array([5,6])

print('Distance from point to Segment is: ',computeDistancePointToSegment(x1,y1,x2,y2,q))
plt.figure(1)
plt.plot([x1,x2],[y1,y2])
plt.plot(q[0],q[1],'rx')


P = np.array([[0,0],[1,0],[1,1],[0,1]])
q = np.array([0,0.5])
(d,I) = computeDistancePointToPolygon(P,q)
u = computeTangentVectorToPolygon(P,q)
print('D is: ',d)
print('I is: ',I)


plt.figure(2)
plt.fill(P[:,0],P[:,1])
plt.plot(q[0],q[1],'rx')
plt.quiver(q[0],q[1],u[0],u[1],scale=10)
plt.legend(['Start Point','Obstacle','TangentVector'])
plt.axis([-1, 6, -1, 6])




plt.show()
