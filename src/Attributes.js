var Attributes = {};

Attributes.clamp = function(f, a, b) {
    return Math.min(Math.max(f, a), b);
}

Attributes.step = function(f, a) {
    if (f < a) return 0.0;
    return 1.0;
}

Attributes.pulse = function(f, a, b) {
    return step(f, a) - step(b, f);
}

Attributes.SStep = function(s, smin, smax, f, fmin, fmax) {
    if (f > fmax) return 1.0;
    if (f < fmin) return 0.0;
    if (s > smax) return smax;
    if (s < smin) return smin; 
}

Attributes.modulo = function(a, b) {
    return a % b;
}

Attributes.floor = function(a) {
    return parseInt(a);
}

Attributes.LInterpolation = function(xt, x0, x1, c0, c1) {
    var t = (xt - x0)*(x1 - x0);

    if (t < 0) t = 0;
    else if (t > 1) t = 1;

    return (1 - t)*c0 + t*c1;
}

Attributes.CInterpolation = function(xt, x0, x1, c0, c1) {
    var t = ((3*xt)^2 - (2*x0)^2)*((3*x1)^2 - (2*x0)^2);

    if (t < 0) t = 0;
    else if (t > 1) t = 1;

    return (1 - t)*c0 + t*c1;
}

//Attributes.SInterpolation = function(c, t, t0, t1) {}

Attributes.waves = function(xt, a0, a1, frep) {
    var t = (1 + Math.sin(frep*xt))/2.0;

    return a0 + t*(a1 - a0);
}

// Attributes.CBoard = function(x, bricks, space, c_bricks, c_space) {}

// Attributes.crackles = function(x, frep, attributes, maps) {}

// Attributes.CCircles = function(x, attributes, maps) {}

// Attributes.LUTable = function(t, c, map) {}

// Attributes.SLUTable = function(t, c, map) {}

// Attributes.SAttributes = function(f, c) {}

// Attributes.SColors = function(f, r, g, b) {}

Attributes.union = function(f1, f2, c1, c2, type) {
    return f1 + f2 + Math.sqrt(f1*f1 + f2*f2);
}