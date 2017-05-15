import numpy as np
import matplotlib.pyplot as plt
from distanceFunctions import *
def plotObstacles(ObstaclesList,pos_start,pos_goal,path):
    #plots Obstacles, start and goal positions, and path to goal
    #plot start and end positions
    plt.plot(pos_start[0],pos_start[1],'bx')
    plt.plot(pos_goal[0],pos_goal[1],'ro')
    #Plot path
    plt.plot(path[:,0],path[:,1],'b')
    n = len(ObstaclesList)
    #Initilize plot bounds with first obstacle
    P = ObstaclesList[0]
    plt.fill(P[:,0],P[:,1])
    xmin = min(P[:,0])
    xmax = max(P[:,0])
    ymin = min(P[:,1])
    ymax = max(P[:,1])
    for i in range(1,n):                                                    #For each obstacle, plot and compare min and max x,y values
        P = ObstaclesList[i]
        plt.fill(P[:,0],P[:,1])
        if min(P[:,0])<xmin:
            xmin = min(P[:,0])
        if max(P[:,0])>xmax:
            xmax = max(P[:,0])
        if min(P[:,1])<ymin:
            ymin = min(P[:,1])
        if max(P[:,1])>ymax:
            ymax = max(P[:,1])
    xmin = min(xmin,pos_start[0],pos_goal[0])-1
    xmax = max(xmax,pos_start[0],pos_goal[0])+1
    ymin = min(ymin,pos_start[1],pos_goal[1])-1
    ymax = max(ymax,pos_start[1],pos_goal[1])+1
    plt.axis([xmin,xmax,ymin,ymax])

def computeBug(pos_start,pos_goal,obstacleList,step_size):
    #Given: start position and goal position as numpy arrays
    # obstacleList = [P1 P2 ... Pl] where Pi is a 2d numpy array with vertices listed in counterclockwise order
    # step_size - distance the bug moves each step of algorithm
    #Output: Outcome - string describing success or failure of algorithm
    # Path - 2d array of positions at each step
    # d_goalList - List of distance to goal at each step
    #
    #computeBug calculates the path of Bug algorithm which moves towards goal until an obstacle is hit
    #When an obstacle is hit,  it circumnavigates the edge of the obstacle and stores the position that is closest to goal
    #After circumnavigating, it returns to the position closest to the goal and continues directly towards the goal until another obstacle is hit
    #If the obstacle is still in the path at the position closest to the goal, then no path exists
    pos_current = pos_start
    path = np.array([pos_start])
    num_steps = 0
    d_goal = np.linalg.norm(pos_goal-pos_current)                          
    d_goalList = np.array([d_goal])                                         #Initialize distance to goal list to hold distance value for each step
    num_obstacles = len(obstacleList)
    hit_obstacle = 0
    while np.linalg.norm(pos_goal-pos_current)>step_size:                   #Continue until goal is reached
        if not hit_obstacle:
            P = obstacleList[0]
            d_closest = computeDistancePointToPolygon(P,pos_current)[0]
            I_closest = 0
            for i in range(1,num_obstacles):                                #Go through the list of polygons
                P = obstacleList[i]
                d = computeDistancePointToPolygon(P,pos_current)[0]         
                if d<d_closest:                                             #
                    d_closest = d
                    I_closest = i
            if d_closest<0.99*step_size:                                    #If bug has hit obstacle. 0.99*step_size makes sure it doesn't think it hits right after leaving
                tangent = computeTangentVectorToPolygon(obstacleList[I_closest],pos_current)
                u = np.array([-tangent[1],tangent[0]])                      #Compute normal vector of obstacle to move directly onto the edge
                pos_current = pos_current + d_closest*u                     #Move exactly distance to obstacle
                hit_obstacle = 1                                            #Keep track of when you're on the obstacle
                pos_hit = pos_current
                d_goal_min = np.linalg.norm(pos_goal-pos_current)
                pos_closest = pos_current                                   #Keep track of closest position to goal to decide where to leave obstacle
                circ_nav = 0                                                #Has not circumnavigated obstacle
                at_vert = 0
            else:                                                           #If you haven't hit obstacle just move towards goal
                u = (pos_goal-pos_current)/np.linalg.norm(pos_goal-pos_current)
                pos_current = pos_current + step_size*u
        # If you have hit an obstacle, circumnavigate the entire obstacle,
        # navigate to the point closest to the goal, and then attempt to continue towards goal
        else:
            d_closest = np.linalg.norm(pos_current-pos_closest)             #Distance to closest point along obstacle to goal
            if circ_nav and (d_closest<step_size):                          #If circumnavigated and at closest point
                tangent = computeTangentVectorToPolygon(obstacleList[I_closest],pos_current)
                normal = np.array([-tangent[1],tangent[0]])                 #Inner normal direction
                u = (pos_goal-pos_current)/np.linalg.norm(pos_goal-pos_current)
                if np.dot(u,normal)>0:                                      #If path to goal is in the same direction as inner normal vector
                    outcome = 'Failure: path does not exist'                #even at the closest point to the goal, then no path exists
                    return (outcome,path,d_goalList)
                else:
                    hit_obstacle = 0
                    pos_current = pos_current +step_size*u                  #Move towards goal
            else:
                #Take step around boundary counterclockwise, update closest
                #point and whether obstacle has been circumnavigated
                P_hit = obstacleList[I_closest]
                u = computeTangentVectorToPolygon(P_hit,pos_current)
                d_vert = minDistanceToVertex(P_hit,pos_current)
                if (d_vert<step_size)and(not at_vert):                      #If the distance to vertex is less than the step size
                                                                            #and bug is not already at vertex
                    pos_current = pos_current + d_vert*u                    #Move exactly the distance to the vertex
                    at_vert = 1                                             #Mark that you are at vert so you can move off
                else:
                    pos_current = pos_current + step_size*u
                    at_vert = 0
                d_goal = np.linalg.norm(pos_goal-pos_current)
                if d_goal<d_goal_min:                                       #If distance to the goal is less than the minimum,
                    d_goal_min = d_goal                                     #update minimum
                    pos_closest = pos_current                               #and the closest position
                if np.linalg.norm(pos_current-pos_hit)-step_size < -1e-6:   #If close to the point where the obstacle was originally hit,
                    circ_nav = 1                                            #then obstacle has been circumnavigated

        # each step of while loop
        # For debugging:
        #print('Current Pos: ',pos_current)
        #print('Velocity Vect: ',u)
        num_steps += 1
        path = np.insert(path,num_steps,pos_current,axis=0)
        d_goal = np.linalg.norm(pos_goal-pos_current)
        d_goalList = np.insert(d_goalList,num_steps,d_goal,axis=0)
        if np.size(path,0)>10000: #check that it is not in a continuous loop
            outcome = 'Error: Continuous Loop'
            return (outcome,path,d_goalList)
    path = np.insert(path,num_steps,pos_goal,axis=0)
    outcome = 'Success'
    return (outcome,path,d_goalList)
    
    
