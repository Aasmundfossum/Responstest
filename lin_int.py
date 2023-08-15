# Lineær interpolering

def lin_int(x,x1,x2,y1,y2):
    y = y1+(x-x1)*(y2-y1)/(x2-x1)
    return y