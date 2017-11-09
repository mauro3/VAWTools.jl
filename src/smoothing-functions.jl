"""
    hyperbola2(x, s1, s2, Δy, xoff=0.0, yoff=0.0)

One branch of a hyperbola as a function of x given by

- s1: asymptote slope for x<0
- s2: asymptote slope for x>0
- Δy: absolute value of abcissa at x=0
- optionally shift the origin by xoff, yoff

Ref https://maurow.bitbucket.io/docs/hyperbola.html
"""
function hyperbola(x, s1, s2, Δy, xoff=0.0, yoff=0.0)
    x = x-xoff
    # Deal with degenerate cases and select the right branch of the hyperbola
    if Δy==0 # return two lines
        i1 = x.<0
        i2 = x.>=0
        return vcat(s1*x[i1], s2*x[i2]) + yoff
    elseif s1==s2 # return one line
        branch = 0
    elseif s1>s2
        branch = -1
        s1,s2 = s2,s1
    else
        branch = +1
    end
    B = (s1+s2)/2
    C = -Δy^2
    D = - (s2-s1)^2/4
    return B*x + branch * sqrt(-D*x.^2-C) + yoff
end

"""
    hyperbola01(x, Δy, xoff=0.0, yoff=0.0)

A hyperbola with asymptotes y=0 for x<0 and y=x for x>0, and
abcissa Δy.  Its origin can be moved with xoff and yoff.
"""
function hyperbola01(x, Δy, xoff=0.0, yoff=0.0)
    x = x-xoff
    B = 1/2
    C = -Δy^2
    D = -1/4
    return B*x + sqrt(-D*x^2-C) + yoff
end

"""
    max_smooth(x0, x, delta)

to replace

    max(x0, x)

for a fixed x0 and varying x.  `delta` gives the difference:
delta = max_smooth(x0,x0,delta)-x0,
which is the maximal difference for any x.

Note, all returned values >= x0 (== for delta==0).
"""
max_smooth(x0, x, delta) =
    delta==0 ? max(x0,x) : hyperbola01(x, delta, x0, x0)
"""
    min_smooth(x0, x, delta)

see max_smooth
"""
min_smooth(x0, x, delta) = -max_smooth(-x0, -x, delta)


####
# Sigmoid smoothing

"""
    sigmoid(x, x0, w)

Smooth transition function from 0 to 1.  At x0+w its value is about 0.99,
at x0-w about 0.01.
"""
sigmoid(x, x0, w) = 1/(1+exp(-(x-x0)/w*5))


"""
    fn_cap(x0, width, fn, x, args...)

Instead of

    max(fn(x0,args...), fn(x, args...))

where at x0+/-width the value is to 1% correct.

There is also the not recommended:

    fn_cap_(f0, width, fn, args...)

where at f0+/-width the value is to 1% correct.
"""
function fn_cap(x0, width, fn, x, args...)
    ff = fn(x, args...)
    f0 = fn(x0, args...)
    si = sigmoid(x, x0, width)
    f0*(1-si) + ff*si
end
function fn_cap_(f0, width, fn, args...)
    ff = fn.(args...)
    si = sigmoid.(ff, f0, width)
    f0.*(1-si) + ff.*si
end

"""
    val_cap(x0,x,width)
instead of
    max(x0, x)

with the return value within 1% for x=x0+/-width.
"""
function val_cap(x0,x,width)
    si = sigmoid(x,x0,width)
    x0*(1-si)+x*si
end
