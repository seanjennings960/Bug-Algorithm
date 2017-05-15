import numpy as np
import matplotlib.pyplot as plt
from computeBug import *
from distanceFunctions import *


P1 = np.array([[1,2],[1,0],[3,0]])
P2 = np.array([[2,3],[4,1],[5,2]])
ObstacleList = [P1,P2]
p_goal = np.array([5,3])
p_start = [0.,0.5]
step_size = .1
(outcome,path,d_goal) = computeBug(p_start,p_goal,ObstacleList,step_size)
print(outcome)
plotObstacles(ObstacleList,p_start,p_goal,path)
plt.figure(2)
plt.plot(d_goal)
plt.title('Distance to Goal Vs. Time')
plt.xlabel('Time (steps)')
plt.ylabel('Distance')
plt.figure(3)
plt.plot(path[:,0],path[:,1])
plt.title('Path plot without obstacles')
totalPathLength = len(d_goal)*step_size
plt.show()
