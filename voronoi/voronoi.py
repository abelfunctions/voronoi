r"""
Voronoi
=======

For now, a single script that encompasses all of the worksheet functionality.
"""

from sage.all import QQbar, QQ
from scipy.spatial import Voronoi, voronoi_plot_2d


def branch_points(f):
    r"""Returns the branch locus of f.

    Parameters
    ----------
    f : algebraic curve

    Returns
    -------
    rts : list
        A list of the branch points over QQbar.
    """
    x,y = f.parent().gens()
    res = f.resultant(f.derivative(y),y).univariate_polynomial()
    rts = res.roots(ring=QQbar, multiplicities=False)
    return rts

def boundary_points(branch_points):
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
    # Defining the boundary points for the Voronoi diagram. In the future, it
    # may be useful to try using a circle (or rather an approximation to one).
    # The below code will ensure that our branch points are actually contained
    # in the interior of the box, but there may be better designs for the box.

    # First try on box.
    maxx = 2*max(x for x,y in branch_points)
    minx = 2*min(x for x,y in branch_points)
    maxy = 2*max(y for x,y in branch_points)
    miny = 2*min(y for x,y in branch_points)

    # If some of the coordinates are zero, then our points will lie on the
    # boundary of the box. We don't want this. There is really only one case.
    # If the max is zero, we know the min is negative, so we can just flip the
    # sign. A similar argument works for if the min is zero, we just flip
    # signs. Then we can do it for the y as well.
    if abs(maxx) < 10^(-5):
        maxx = -minx
    if abs(minx) < 10^(-5):
        minx = -maxx
    if abs(maxy) < 10^(-5):
        maxy = -miny
    if abs(miny) < 10^(-5):
        miny = -maxy

    # Also, they are actually supposed to form a box!! So if one of the
    # intervals collapses to a single point, we have to fix that.
    #
    # In our situation, we will always have at least two points- Otherwise, our
    # algebraic curve has genus zero. This is not an interseting case for us.
    # This means that only one of the intervals could be trivial. For instance,
    # if all of the roots are real, they all lie on the x-axis.
    #
    # To fix this, we just make the interval the same as the one which is not
    # trivial. Shifted accordingly so that the points are centered at it.
    #
    # There may be a better way to do this, but I'm not sure what that is. This
    # is just a temporary fix anyways.
    if (abs(minx - maxx) < 10^(-5)):
        minx = branch_points[0][0] - (abs(maxy) + abs(miny))/2
        maxx = branch_points[0][0] + (abs(maxy) + abs(miny))/2

    if (abs(miny - maxy) < 10^(-5)):
        miny = branch_points[0][1] - (abs(maxx) + abs(minx))/2
        maxy = branch_points[0][1] + (abs(maxx) + abs(minx))/2

    boundary_points =  [(x,y) for x in [minx,maxx] for y in [miny,maxy]]
    return boundary_points

def voronoi(f):
    r"""Returns the Voronoi diagram of the branch locus of the curve.

    Parameters
    ----------
    f : algebraic curve

    Returns
    -------
    voronoi : Scipy voronoi diagram.
        See :func:`scipy.spatial.Voronoi` for details.
    """
    branch_points = branch_points(f)
    boundary_points = boundary_points(branch_points)
    points = branch_points + boundary_points
    voronoi = Voronoi(points)
    return voronoi
