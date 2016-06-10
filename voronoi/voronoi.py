r"""
Voronoi
=======

For now, a single script that encompasses all of the worksheet functionality.
"""

from sage.all import QQbar, QQ, CC, I
from scipy.spatial import Voronoi, voronoi_plot_2d

def branch_points(f):
    r"""Returns the branch locus of f.

    Parameters
    ----------
    f : algebraic curve

    Returns
    -------
    points : list
        A list of the branch points over QQbar as x- y-coordinates.
    """
    x,y = f.parent().gens()
    res = f.resultant(f.derivative(y),y).univariate_polynomial()
    rts = res.roots(ring=QQbar, multiplicities=False)
    points = [(z.real(), z.imag()) for z in rts]
    return points

def boundary_points(branch_pts):
    r"""Returns a list of boundary / ghost points to add to the branch locus.

    The Voronoi diagram of just the branch points will return rays going off to
    infinity. We add points "around" the branch points such that the Voronoi
    edges can be used to construct monodromy cycles.

    Parameters
    ----------
    branch_points : list, QQbar
        A list of branch points.

    Returns
    -------
    boundary : list, QQbar
        A list of boundary points
    """
    # First try on box.
    maxx = 2*max(x for x,y in branch_pts)
    minx = 2*min(x for x,y in branch_pts)
    maxy = 2*max(y for x,y in branch_pts)
    miny = 2*min(y for x,y in branch_pts)

    # degenerate case: a boundary box edge lies on the real axis, in which case
    # some branch points may lie on an edge
    if abs(maxx) < 10**(-5):
        maxx = -minx
    if abs(minx) < 10**(-5):
        minx = -maxx
    if abs(maxy) < 10**(-5):
        maxy = -miny
    if abs(miny) < 10**(-5):
        miny = -maxy

    # degenerate case: boundary box consists of an edge or point
    #
    # note: this code, at the moment, assumes that the curve has at least two
    # finite branch points. this is a very strong assumption and should be made
    # more general
    if (abs(minx - maxx) < 10**(-5)):
        x_shift = (abs(maxy) + abs(miny))/2
        minx = branch_pts[0][0] - x_shift
        maxx = branch_pts[0][0] + x_shift

    if (abs(miny - maxy) < 10**(-5)):
        y_shift = (abs(maxx) + abs(minx))/2
        miny = branch_pts[0][1] - y_shift
        maxy = branch_pts[0][1] + y_shift

    # compute the box vertices and return
    boundary_pts =  [(x,y) for x in [minx,maxx] for y in [miny,maxy]]
    return boundary_pts

def voronoi(f):
    r"""Returns the Voronoi diagram of the branch locus of the curve.

    Parameters
    ----------
    f : algebraic curve

    Returns
    -------
    v : Scipy voronoi diagram.
        See :func:`scipy.spatial.Voronoi` for details.
    """
    branch_pts = branch_points(f)
    boundary_pts = boundary_points(branch_pts)
    points = branch_pts + boundary_pts
    v = Voronoi(points)
    return v


def vertex_lifts(f, voronoi, ring=CC):
    r"""For each vertex of the Voronoi diagram, compute the lift of the vertex on
    the Riemann surface.

    Since the vertices of the Voronoi diagram are, by design, far from any
    branch points of the curve each root lying above the vertex corresponds to
    a place on the Riemann surface.

    Parameters
    ----------
    f : algebraic curve
    voronoi
        The voronoi diagram of the branch locus of `f`.
    ring : Sage Field
        (Default: `CC`) The ring or field over which to compute the lifts.

    Returns
    -------
    vertex_lifts : dictionary

        A dictionary of the lifts. The keys are the vertices, values are the
        lifts, themselves.
    """
    R = f.parent()
    x,y = R.gens()

    # make sure that the input ring contains Q[I]
    I = ring.gen()
    if not I.imag():
        raise ValueError('The ring %s must contain I.'%ring)

    # coerce the vertices to complex number in the ring
    vertices = [ring(x0) + I*ring(y0) for x0,y0 in voronoi.vertices]
    lifts = [f(z0, y).univariate_polynomial().roots(ring=ring,
                                                    multiplicities=False)
             for z0 in vertices]
    vertex_lifts = dict(zip(vertices, lifts))
    return vertex_lifts
