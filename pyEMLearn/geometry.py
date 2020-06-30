
def circle(x0,y0,r):
    return [
        lambda x,y: (x-x0)**2 + (y-y0)**2 < r
    ]

def square_bbox(b,l,t,r):
    return [
        lambda x,y: y > b,
        lambda x,y: y < b,
        lambda x,y: x > l,
        lambda x,y: y < r
    ]
