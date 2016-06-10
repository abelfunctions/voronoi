from sage.all import QQbar, QQ, CC, CDF, AA

from voronoi import *

# list of test curves
R = QQ['x,y']
x,y = R.gens()
f1 = (x**2 - x + 1)*y**2 - 2*x**2*y + x**4
f2 = -x**7 + 2*x**3*y + y**3
f3 = (y**2-x**2)*(x-1)*(2*x-3) - 4*(x**2+y**2-2*x)**2
f4 = y**2 + x**3 - x**2
f5 = (x**2 + y**2)**3 + 3*x**2*y - y**3
f6 = y**4 - y**2*x + x**2
f7 = y**3 - (x**3 + y)**2 + 1
f8 = x**2*y**6 + 2*x**3*y**5 - 1
f9 = 2*x**7*y + 2*x**7 + y**3 + 3*y**2 + 3*y
f10 = (x**3)*y**4 + 4*x**2*y**2 + 2*x**3*y - 1
f11 = y**7 - (x*((x-1)**2))
f12 = x*y**3 + y**3 + x
f_elliptic = y**2 - (x-2)*(x-1)*(x+1)*(x+2)

def test_branch_points():
    b = branch_points(f_elliptic)
    b_actual = map(
        lambda p: tuple(map(AA, p)),
        [(-2,0), (-1,0), (1,0), (2,0)])
    assert set(b) == set(b_actual)

def test_boundary_points():
    # tests that the Euclidean distance of each boundary point is larger than
    # the Euclidean distance of any branch point
    b = branch_points(f_elliptic)
    bdry = boundary_points(b)

    for xb,yb in bdry:
        normb = xb**2 + yb**2
        for x,y in b:
            norm = x**2 + y**2
            assert normb > norm

def test_voronoi():
    # for now, just tests that no errors are raised
    v = voronoi(f_elliptic)

def test_vertex_lifts_degree():
    # test that the correct number of roots are computed with various rings
    v = voronoi(f_elliptic)
    lifts_CC = vertex_lifts(f_elliptic, v, CC)
    for lift in lifts_CC.values():
        assert len(lift) == f_elliptic.degree(y)

    lifts_CDF = vertex_lifts(f_elliptic, v, CDF)
    for lift in lifts_CDF.values():
        assert len(lift) == f_elliptic.degree(y)

