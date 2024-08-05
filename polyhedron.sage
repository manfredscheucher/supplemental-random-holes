from itertools import *
from sys import *


def test_point_inside_polyhedron(X):
    for n,F0 in bounding_planes:
        assert(n*F0 < n*M)
        if n*F0 > n*X:
            return False

    return True


def random_point_in_polyhedron():
    while True: 
        P = vector(RR(random()*(xmax-xmin)+xmin) for i in range(3))
        if test_point_inside_polyhedron(P): 
            return P


def plane_normal_vector(P,Q,R):
    return (P-Q).cross_product(R-Q)


def area_triangle(P,Q,R):
    return plane_normal_vector(P,Q,R).norm()/2


def intersect_segment_with_plane(u,v,n,point_in_face):
    l = v-u
    assert(l*u < l*v)
    assert(l*n != 0) # line not parallel with probability 1, this is negligible
    d = ((point_in_face-u)*n)/(l*n)
    q = u+l*d
    assert((q-point_in_face)*n < 10^-10) # should lie in plane    
    assert (l*q >= l*u and l*q <= l*v) # should be on segment

    return q 


def intersect_triangle_with_plane(triangle,n,F0):
    assert(n*F0 < n*M)

    P,Q,R = triangle
    if n*Q < n*R: Q,R = R,Q
    if n*P < n*Q: P,Q = Q,P
    if n*Q < n*R: Q,R = R,Q
    assert(n*P >= n*Q and n*Q >= n*R) # points sorted with decreasing distince to face

    if n*R > n*F0: # if all 3 inside, return
        return [(P,Q,R)]
    elif n*Q > n*F0: # R is single outside point, PQ inside
        P2 = intersect_segment_with_plane(P,R,n,F0)
        Q2 = intersect_segment_with_plane(Q,R,n,F0)
        return [(P,Q,Q2),(P,P2,Q2)]
    elif n*P > n*F0: # P is single inside point, QR outside
        Q2 = intersect_segment_with_plane(P,Q,n,F0)
        R2 = intersect_segment_with_plane(P,R,n,F0)
        return [(P,Q2,R2)]
    else: # no point side 
        return [] 


def area_flipped_C_triangle(A,B,C):
    C2 = A+B-C
    if test_point_inside_polyhedron(C2): 
        return area_triangle(A,B,C)
    else:
        triangles = [(A,B,C2)]
        for n,F0 in bounding_planes:
            new_triangles = []
            for T in triangles:
                new_triangles += intersect_triangle_with_plane(T,n,F0)
            triangles = new_triangles

        return sum(area_triangle(*T) for T in triangles)



sample = 0
total_sum = 0
all_samples = []

shape = argv[1]

if shape == "t": shape = "tetrahedron"
if shape == "c": shape = "cube"
if shape == "o": shape = "octahedron"
if shape == "d": shape = "dodecahedron"
if shape == "i": shape = "icosahedron"

if shape == "tetrahedron":
    vertices = [vector((0,0,0)),vector((1,0,0)),vector((0,1,0)),vector((0,0,1))]
elif shape == "cube":
    vertices = [vector((x,y,z)) for x in [0,1] for y in [0,1] for z in [0,1]]
elif shape == "octahedron":
    vertices = []
    vertices += [vector((a,0,0)) for a in [-1,1]]
    vertices += [vector((0,a,0)) for a in [-1,1]]
    vertices += [vector((0,0,a)) for a in [-1,1]]
elif shape == "dodecahedron":
    phi = RR((1+sqrt(5))/2)
    vertices = []
    vertices += [vector((a,b,c)) for a in [-1,1] for b in [-1,1] for c in [-1,1]]
    vertices += [vector((0,a*phi,b/phi)) for a in [-1,1] for b in [-1,1]]
    vertices += [vector((a*phi,b/phi,0)) for a in [-1,1] for b in [-1,1]]
    vertices += [vector((b/phi,0,a*phi)) for a in [-1,1] for b in [-1,1]]
elif shape == "icosahedron":
    phi = RR((1+sqrt(5))/2)
    vertices = []
    vertices += [vector((0,a*1,b*phi)) for a in [-1,1] for b in [-1,1]]
    vertices += [vector((a*1,b*phi,0)) for a in [-1,1] for b in [-1,1]]
    vertices += [vector((b*phi,0,a*1)) for a in [-1,1] for b in [-1,1]]
else:
    print "ERROR: shape",shape,"not implemented"
    exit()


xmin = min(x for V in vertices for x in V)
xmax = max(x for V in vertices for x in V)
assert(xmin < xmax)

M = sum(vertices)/len(vertices) # mid point

bounding_planes = []
for F in combinations(vertices,3):
    n = plane_normal_vector(*F)
    eps = 10^-6
    F_all = [V for V in vertices if abs(n*V - n*F[0])<eps]
    if F != tuple(F_all[:3]): continue # avoid that planes containing more than 3 points are overcounted

    if n*F[0] > n*M: 
        n = -n

    assert(n*F[0] <= n*M)

    is_bounding = True
    for V in vertices:
        if n*F[0] > n*V+eps:
            is_bounding = False

    if is_bounding:
        bounding_planes.append((n,F[0]))


print "shape has",len(vertices),"vertices and",len(bounding_planes),"faces/bounding planes"
print "vertices",vertices

outfile = "samples_"+shape+".txt"
print "write to outfile",outfile
#samples_file = open(outfile,"w")

next_output = 1

while True:
    A = random_point_in_polyhedron()
    B = random_point_in_polyhedron()
    C = random_point_in_polyhedron()
    area = area_triangle(A,B,C)
    area_C = area_flipped_C_triangle(A,B,C)
    area_B = area_flipped_C_triangle(C,A,B)
    area_A = area_flipped_C_triangle(B,C,A)
    total_area = area+area_A+area_B+area_C
    assert(area_A/area < 1+10^-6)
    assert(area_B/area < 1+10^-6)
    assert(area_C/area < 1+10^-6)

    ratio = total_area/area

    sample += 1
    total_sum += QQ(ratio)
    
    all_samples.append(RR(ratio))
    #samples_file.write(str(RR(ratio))+"\n")

    if sample == next_output:
        next_output += next_output
        m = RR(total_sum/sample)
        s = sqrt(sum((x-m)^2 for x in all_samples)/(sample-1))
        print "sample:",sample,"\tm:",m,"+-",RR(s/sqrt(sample)),"\ts:",s
