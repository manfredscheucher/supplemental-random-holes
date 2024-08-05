
def random_point_in_sphere():
    while True: 
        P = vector(RR(2*random()-1) for i in range(3))
        if P.norm() < 1: return P


def plane_normal_vector(P,Q,R):
    return (P-Q).cross_product(R-Q)


def area_triangle(P,Q,R):
    return plane_normal_vector(P,Q,R).norm()/2


def intersect_segment_with_sphere(U,V):
    t = var('t')
    line = ((1-t)*U+t*V)
    poly = sum(x^2 for x in line)-1
    coeff = poly.coefficients()
    for i in range(len(coeff)): assert(coeff[i][1] == i)
    a = coeff[2][0]
    b = coeff[1][0]
    c = coeff[0][0]
    assert(a > 0)
    r = (-b + sqrt(b^2-4*a*c)) / (2*a)
    assert(r > 0 and r < 1)
    return (1-r)*U+r*V


def area_flipped_C_triangle(A,B,C):
    C2 = A+B-C
    if C2.norm() < 1: 
        return area_triangle(A,B,C)
    else:
        A2 = intersect_segment_with_sphere(A,C2)
        B2 = intersect_segment_with_sphere(B,C2)
        assert(abs(A2.norm()-1)<10^-6)
        assert(abs(B2.norm()-1)<10^-6)
        part1 = area_triangle(A,A2,B)
        part2 = area_triangle(A2,B,B2)
        v = plane_normal_vector(A,B,C).normalized()
        offset = v*A2
        A3 = (A2-offset*v)
        B3 = (B2-offset*v)
        rad3 = A3.norm()
        assert((rad3 - B3.norm()) < 10^-6)
        phi = acos(A3*B3/rad3^2)
        part3u = area_triangle(A2,B2,C2)
        part3 = rad3^2*phi/2-area_triangle(origin,A3,B3)
        assert(part3 < part3u+10^-6)
        return part1+part2+part3



    
origin = vector((0,0,0))

sample = 0
total_sum = 0
all_samples = []

#samples_file = open("samples_ball.txt","w")

next_output = 1

while True:
    A = random_point_in_sphere()
    B = random_point_in_sphere()
    C = random_point_in_sphere()
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



