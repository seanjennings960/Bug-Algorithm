import math
import numpy as np
def computeDistancePointToSegment(x1,y1,x2,y2,q):
    # x1,y1,x2,y2 = endpoints of segments
    # q = vector of point coordinates

    # Finds distance from point q to segment with endpoints (x1,y1) and (x2,y2)
    # Uses tangential coordinate system to determine whether to use distance to line or endpoints

    len_segment = math.sqrt((x2-x1)**2+(y2-y1)**2)
    t_hat = [(x2-x1)/len_segment, (y2-y1)/len_segment]                  #Define tangential unit vector
    t1 = x1*t_hat[0] + y1*t_hat[1]                                      #Find the tangential coordinates of each point
    t2 = x2*t_hat[0] + y2*t_hat[1]
    tq = q[0]*t_hat[0] + q[1]*t_hat[1]                                  #tangential coordinate of point q
    # Check if point has tangential components between endpoints
    if (t1<t2)&(tq>t1)&(tq<t2):
        d = computeDistancePointToLine(x1,y1,x2,y2,q)                   #If it is between, use Point to Line formula
    else:                                                               #If outside use distance to closest endpoint
        d1 = math.sqrt((q[0]-x1)**2+(q[1]-y1)**2)
        d2 = math.sqrt((q[0]-x2)**2+(q[1]-y2)**2)
        d = min([d1,d2])
    return d

def computeLineThroughTwoPoints(x1,y1,x2,y2):
    #Given two points (x1,y1) and (x2,y2)
    #Output: (a,b,c) that satisfies equation of the line ax+by=c between given points
    #Values are normalized so [a,b] is normal unit vector
    d = math.sqrt((y2-y1)**2+(x2-x1)**2)
    if d==0:
        a = 0
        b = 0
        c = 0
    else:
        a = (y1-y2)/d
        b = (x2-x1)/d                                                   #Rearrangement of Point-Slope formula for line
        c = ((y2-y1)*x1-(x2-x1)*y1)/d
    return [a,b,c]

def computeDistancePointToLine(x1,y1,x2,y2,q):
    # (x1,y1) and (x2,y2) are points defining line
    # q is list containing x and y coordinates of point q
    # Output: d - normal distance from the point to the line
    line_vals = computeLineThroughTwoPoints(x1,y1,x2,y2)
    n = [line_vals[0], line_vals[1]]
    v = [x1-q[0],y1-q[1]]                                               #Vector between q and one point on line
    d = abs(n[0]*v[0]+n[1]*v[1])                                        #Dot product of normal and vector to line gives distance
    return d

def computeDistancePointToPolygon(P,q):
    #P is 2d array of Polygon vertices of the form [[x1,y1],[x2,y2],...]
    #Vertices must be listed in counterclockwise order
    n = np.size(P,0)                                                    #Number of edges
    D = []                                                              #List to hold distance to each edge
    for i in range(n):
        if i<n-1:
            D.append(computeDistancePointToSegment(P[i][0],P[i][1],P[i+1][0],P[i+1][1],q))
        else:
            D.append(computeDistancePointToSegment(P[n-1][0],P[n-1][1],P[0][0],P[0][1],q))
    d = min(D)
    I = D.index(d)
    return (d,I)
    
def computeTangentVectorAtVertex(P,q):
    #P is polygon with vertices listed in counterclockwise order
    # q is a point listed as a numpy array
    # if q is at the vertex of P, then u is the unit tangent vector that points to the next counterclockwise vertex 
    n = np.size(P,0)
    u = 0
    for i in range(n):
        if (abs(P[i,0]-q[0])<1e-6) and (abs(P[i,1]-q[1])<1e-6):
            if i == n-1:
                w = P[0,:]-q
                u = w/np.linalg.norm(w)
            else:
                w = P[i+1,:]-q                                          #Tangent vector should point to the next vertex
                u = w/np.linalg.norm(w)
    return u

def computeTangentVectorToPolygon(P,q):
    #P is polygon with vertices list in counterclockwise order
    n = np.size(P,0)
    (d,I) = computeDistancePointToPolygon(P,q)
    vert = 0                                                            #0 if closest point is segment, 1 if vertex is closest
    for i in range(n):
        d_vert = np.linalg.norm(P[i,:]-q)                               #Compute the distance to each vertice
        if abs(d-d_vert)<1e-6:                                          #d==d_vert with error accounted for
            j = i                                                       #Index of closest vertex
            vert = 1

    if vert == 0:                                                       #If point is closest to segment
        if I==n-1:
            w = P[0,:]-P[n-1,:]
        else:
            w = P[I+1]-P[I,:]                                           #Vector towards next vertex
        u = w/np.linalg.norm(w)                                         

    else:                                                               #If the point is closest to vertex
        v = P[j,:]-q                                                    #Vector between point and vertex
        v_mag = np.linalg.norm(v)
        if v_mag<1e-6:                                                  #If point is exactly at vertex,
            u = computeTangentVectorAtVertex(P,q)                       #use function to find tangent vector
        else:
            u = np.array([v[1]/v_mag,-v[0]/v_mag])

    return u

def minDistanceToVertex(P,q):
    #Given: Polygon P and point q
    #Output: d_min - the minimum distance to the closest vertex of P
    n = np.size(P,0)
    d_min = np.linalg.norm(P[0,:]-q)
    for i in range(1,n):
        d = np.linalg.norm(P[i,:]-q)
        if d<d_min:
            d_min = d
    return d_min

    
