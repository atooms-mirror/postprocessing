def linear_grid(min,max,delta):
    """Linear grid."""
    if type(delta) is int:
        n = delta
        if n > 1:
            delta = (max - min) / (n-1)
        else:
            delta = 0.0
    else:
        n = int((max-min)/delta)+1
    list = [min+i*delta for i in range(n)]
    return list

def logx_grid(x1, x2, n):
    """Create a list of n numbers in logx scale from x1 to x2."""
    # the shape if a*x^n. if n=0 => a=x1, if n=N => x1*x^N=x2 
    if x1 > 0:
        xx = (x2/x1)**(1.0/n)
        return [x1] + [x1 * xx**(i+1) for i in range(1,n)]
    else:
        xx = (x2)**(1.0/n)
        return [x1] + [xx**(i+1)-1 for i in range(1,n)]
