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

Primitives.cylZ = function(x, y, z, center, radius) {
    var xt = x - center[0];
    var yt = y - center[1];

    return radius*radius - xt*xt - yt*yt;
}

Primitives.coneX = function(x, y, z, center, radius) {
    var xt = x - center[0];
    var yt = y - center[1];
    var zt = z - center[2];

    return xt*xt - (yt*yt / radius*radius) - (zt*zt / radius*radius);
}

Primitives.coneY = function(x, y, z, center, radius) {
    var xt = x - center[0];
    var yt = y - center[1];
    var zt = z - center[2];

    return yt*yt - (xt*xt / radius*radius) - (zt*zt / radius*radius);
}

Primitives.coneZ = function(x, y, z, center, radius) {
    var xt = x - center[0];
    var yt = y - center[1];
    var zt = z - center[2];

    return zt*zt - (xt*xt / radius*radius) - (yt*yt / radius*radius);
}

Primitives.ellConeX = function(x, y, z, center, a, b) {
    var xt = x - center[0];
    var yt = y - center[1];
    var zt = z - center[2];

    return xt*xt - (yt*yt / a*a) - (zt*zt / b*b);
}

Primitives.ellConeY = function(x, y, z, center, a, b) {
    var xt = x - center[0];
    var yt = y - center[1];
    var zt = z - center[2];

    return yt*yt - (xt*xt / a*a) - (zt*zt / b*b);
}

Primitives.ellConeZ = function(x, y, z, center, a, b) {
    var xt = x - center[0];
    var yt = y - center[1];
    var zt = z - center[2];

    return zt*zt - (xt*xt / a*a) - (yt*yt / b*b);
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

Primitives.superEllipsoid = function(x, y, z, center, a, b, c, s1, s2) {
    var xt = x - center[0];
    var yt = y - center[1];
    var zt = z - center[2];

    return 1 - Math.pow((Math.pow((xt/a), 2)/s2 + Math.pow((yt/b), 2)/s2), (s2/s1))
            - Math.pow((zt/c), (2/s1));
}

Primitives.block = function(x, y, z, vertex, dx, dy, dz) {
    var xt = (x - vertex[0]) * (vertex[0] + dx - x);
    var yt = (y - vertex[1]) * (vertex[1] + dy - y);
    var zt = (z - vertex[2]) * (vertex[2] + dz - z);

    var f = (xt < yt) ? xt : yt;

    return (f < zt) ? f : zt;
}

primitives.blobby = function(x, y, z, center, a, b, t) {
    var xt = x - center[0];
    var yt = y - center[1];
    var zt = z - center[2];

    var blobby = b * Math.exp(-a*(xt*xt + yt*yt + zt*zt));

    return (blobby - t);
}

primitives.metaball = function(x, y, z, center, b, d, t) {
    var xt = x - center[0];
    var yt = y - center[1];
    var zt = z - center[2];

    var r = Math.sqrt(xt*xt + yt*yt + zt*zt);
    var rd = r / d;

    if (r <= d/3) {
        var mb = b * (1 - 3*rd*rd);
    } else var mb = (r < d) ? 1.5*b*(1 - rd)*(1 - rd) : mb;

    return (mb - t);
}

primitives.soft = function(x, y, z, center, d, t) {
    var xt = x - center[0];
    var yt = y - center[1];
    var zt = z - center[2];

    var r2 = xt*xt + yt*yt + zt*zt;
    var r4 = r2*r2;
    var d2 = d*d;
    var d4 = d2*d2;

    var soft = (Math.sqrt(r2) <= d) ? soft + 1 - 22*r2/(9*d2) + 17*r4/(9*d4) 
                - 4*r4*r2/(9*d4*d2) : soft;

    return (soft - t);
}

primitives.convPoint = function(x, y, z, vx, vy, vz, s, t) {
    var r2 = SQ(vx - x) + SQ(vy - y) + SQ(vz - z);
    var f = 1/SQ(1 + SQ(s)*r2);

    return f - t;
}

primitives.convLine = function(x, y, z, bx, by, bz, ex, ey, ez, s, t) {
    var I = sqrt(SQ(ex - bx) + SQ(ey - by) + SQ(ez - bz));

    if (I == 0.0) {
        printf("ERROR:Tips of the segment take same coordinate!\n");
        exit(EXIT_FAILURE);
    }

    var ax = (ex - bx)/I;
    var ay = (ey - by)/I;
    var az = (ez - bz)/I;

    var dx = x - bx;
    var dy = x - by;
    var dz = x - bz;

    var xx = dx*ax + dy*ay + dz*az;
    var p = Math.sqrt(1 + s*s*(dx*dx + dy*dy + dz*dz - xx*xx));
    var q = Math.sqrt(1 + s*s*(dx*dx + dy*dy + dz*dz - 2*xx));

    var f = xx/(2*p*p*(p*p + s*s*xx*xx)) + (I - xx)/(2*p*p*q*q)
            + (Math.atan(s*xx/p) + Math.atan(s*(I - xx)/p))/(2*s*p*p*p);

    return f - t;
}

primitives.convArc = function(x, y, z, center, r, theta, axis, angle, s, t) {
    var rd = Math.PI/180;
    var over_i = 0;
    var over_j = 0;
    var over_k = 1;

    var cx = center[0];
    var cy = center[1];
    var cz = center[2];

    angle += EPS;

    var i = axis[0] + EPS;
    var j = axis[1] + EPS;
    var k = axis[2] + EPS;

    var length = Math.sqrt(i*i + j*j + k*k);
    if (length < EPS) {
        length = EPS;
    }

    i /= length;
    j /= length;
    k /= length;

    var c = Math.cos(rd*(- angle));
    var s = Math.sin(rd*(- angle));

    var one_c = 1 - c;

    var ii = i*i;
    var jj = j*j;
    var kk = k*k;
    var ij = i*j;
    var jk = j*k;
    var ki = k*i;
    var is = i*s;
    var js = j*s;
    var ks = k*s;

    if (theta > 360) theta = 360;
    if (theta > 180) {
        var over_th = (theta - 180)*rd;
        theta = 180;
    
        /* rotate by -angle */
        var tempx = (c + ii*one_c)*(x - cx) + (-ks + ij*one_c)*(y - cy)
                    + (js + ki*one_c)*(z - cz);
        var tempy = (ks + ij*one_c)*(x - cx) + (c + jj*one_c)*(y - cy)
                    + (-is + jk*one_c)*(z - cz);
        var tempz = (-js + ki*one_c)*(x - cx) + (is + jk*one_c)*(y - cy)
                    + (c + kk*one_c)*(z - cz);
    
        var over_c = Math.cos(rd*(-180));
        var over_s = Math.sin(rd*(-180));
        var over_one_c = 1 - over_c;
    
        var over_ii = SQ(over_i);
        var over_jj = SQ(over_j);
        var over_kk = SQ(over_k);
        var over_ij = over_i*over_j;
        var over_jk = over_j*over_k;
        var over_ki = over_k*over_i;
        var over_is = over_i*over_s;
        var over_js = over_j*over_s;
        var over_ks = over_k*over_s;
    
        var over_x = (over_c + over_ii*over_one_c)*tempx
                     + (-over_ks + over_ij*over_one_c)*tempy
                     + (over_js + over_ki*over_one_c)*tempz;
        var over_y = (over_ks + over_ij*over_one_c)*tempx
                     + (over_c + over_jj*over_one_c)*tempy
                     + (-over_is + over_jk*over_one_c)*tempz;
        var over_z = (-over_js + over_ki*over_one_c)*tempx
                     + (over_is + over_jk*over_one_c)*tempy
                     + (over_c + over_kk*over_one_c)*tempz;
        
        var a = 2*r*s*s;
        var d2 = SQ(over_x) + SQ(over_y) + SQ(over_z);
        var b = 1 + SQ(r)*SQ(s) + SQ(s)*d2;
        var p2 = -QU(r)*QU(s) + 2*SQ(r)*SQ(s)*(SQ(s)*(d2 - 2*SQ(over_z)) - 1)
                 - SQ(1 + SQ(s)*d2);
        var p1 = (p2 < 0) ? Math.sqrt(-p2) : Math.sqrt(p2);
        var p3 = p1*p2;
    
        var f1 = (b*over_y)/(over_x*p2*(a*over_x - b))
                    + (a*(SQ(over_x) + SQ(over_y))*Math.sin(over_th) - b*over_y) 
                    /(over_x*p2*(a*(over_x*Math.cos(over_th) + over_y*Math.sin(over_th)) - b));
        
        if (p2 < 0) {
            var f2 = 2.0*b*(Math.atan(-a*over_y/p1)
                             + Math.atan((a*over_y
                                 - (a*over_x + b)*Math.tan(over_th/2.0))/p1))/p3;                 
        } else {
            var f2 = 2.0*b*(Math.atanh(a*over_y/p1)
                             + Math.atanh(((a*over_x + b)*Math.tan(over_th/2.0)
                                - a*over_y)/p1))/p3;
        }
        
        f = f1 + f2;
    }

    return f - t;
}