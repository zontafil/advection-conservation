RIGHT = 1
LEFT = -1


def inside_convex_polygon(point, vertices):
    """Check if a point is inside a convex polygon

    Arguments:
        point {list} -- coordinates of the point
        vertices {list} -- list of the polygon's vertices coordinates

    Returns:
        bool -- True if the point is inside the polygon
    """
    previous_side = None
    n_vertices = len(vertices)
    for n in range(n_vertices):
        a, b = vertices[n], vertices[(n+1) % n_vertices]
        affine_segment = v_sub(b, a)
        affine_point = v_sub(point, a)
        current_side = get_side(affine_segment, affine_point)
        if current_side is None:
            return False  # outside or over an edge
        elif previous_side is None:  # first segment
            previous_side = current_side
        elif previous_side != current_side:
            return False
    return True


def get_side(a, b):
    """check if the point is on the left or right side from the segment

    Arguments:
        a {list} -- coordinates of the vector of the segmnent (0,a)
        b {list} -- coordinates of the point

    Returns:
        int -- return 1 if right, -1 left
    """
    x = x_product(a, b)
    if x < 0:
        return -1
    elif x > 0:
        return 1
    else:
        return None


def v_sub(a, b):
    """return the coordinates of the vector joining point a and b

    Arguments:
        a {list} -- coordinates of point a
        b {list} -- coordinates of point b

    Returns:
        list -- coordinates of the vector
    """
    return (a[0]-b[0], a[1]-b[1])


def x_product(a, b):
    """Two dimensional cross product of two vectors

    Arguments:
        a {list} -- coordinates of vector a
        b {list} -- coordinates of vector b

    Returns:
        float -- cross product
    """
    return a[0]*b[1]-a[1]*b[0]
