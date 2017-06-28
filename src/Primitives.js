var Primitives = {};


Primitives.sphere = function(x, y, z, center, radius) {
    var xt = x - center[0];
    var yt = y - center[1];
    var zt = z - center[2];

    return radius*radius - xt*xt - yt*yt - zt*zt;
}


Primitives.ellipsoid = function(x, y, z, center, a, b, c) {
    var xt = (x - center[0]) / a;
    var yt = (y - center[1]) / b;
    var zt = (z - center[2]) / c;

    return 1.0 - xt*xt - yt*yt - zt*zt;
}


Primitives.ellCylX = function(x, y, z, center, a, b) {
    var yt = (y - center[1]) / a;
    var zt = (z - center[2]) / b;

    return 1.0 - yt*yt - zt*zt;
}


Primitives.ellCylY = function(x, y, z, center, a, b) {
    var xt = (x - center[0]) / a;
    var zt = (z - center[1]) / b;

    return 1.0 - xt*xt - zt*zt;
}


Primitives.ellCylZ = function(x, y, z, center, a, b) {
    var xt = (x - center[0]) / a;
    var yt = (y - center[1]) / b;

    return 1.0 - xt*xt - yt*yt;
}

// add Primitives
Primitives.cylX = function(x, y, z, center, radius) {
    var yt = y - center[1];
    var zt = z - center[2];

    return radius*radius - yt*yt - zt*zt;
}

Primitives.cylY = function(x, y, z, center, radius) {
    var xt = x - center[0];
    var zt = z - center[2];

    return radius*radius - xt*xt - zt*zt;
}

Primitives.coneX = function(x, y, z, center, radius) {
    var xt = x - center[0];
    var yt = y - center[1];
    var zt = z - center[2];

    return xt*xt - (yt*yt / radius*radius) - (zt*zt / radius*radius);
}

Primitives.ellConeX = function(x, y, z, center, a, b) {
    var xt = x - center[0];
    var yt = y - center[1];
    var zt = z - center[2];

    return xt*xt - (yt*yt / a*a) - (zt*zt / b*b);
}

Primitives.torusX = function(x, y, z, center, radius, radius0) {
    var xt = x - center[0];
    var yt = y - center[1];
    var zt = z - center[2];

    return radius0*radius0 - xt*xt - yt*yt - zt*zt - radius*radius
            + 2*radius*Math.sqrt(zt*zt + yt*yt);
}

Primitives.torusY = function(x, y, z, center, radius, radius0) {
    var xt = x - center[0];
    var yt = y - center[1];
    var zt = z - center[2];

    return radius0*radius0 - xt*xt - yt*yt - zt*zt - radius*radius
            + 2*radius*Math.sqrt(zt*zt + xt*xt);
}

Primitives.torusZ = function(x, y, z, center, radius, radius0) {
    var xt = x - center[0];
    var yt = y - center[1];
    var zt = z - center[2];

    return radius0*radius0 - xt*xt - yt*yt - zt*zt - radius*radius
            + 2*radius*Math.sqrt(xt*xt + yt*yt);
}