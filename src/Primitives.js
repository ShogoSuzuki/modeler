var Primitives = {};

SQ = function(x) {
    return x*x;
}

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

Primitives.blobby = function(x, y, z, center, a, b, t) {
    var xt = x - center[0];
    var yt = y - center[1];
    var zt = z - center[2];

    var blobby = b * Math.exp(-a*(xt*xt + yt*yt + zt*zt));

    return (blobby - t);
}

Primitives.metaball = function(x, y, z, center, b, d, t) {
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

Primitives.soft = function(x, y, z, center, d, t) {
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

Primitives.convPoint = function(x, y, z, vect, s, t) {
    var r2 = SQ(vect[0] - x) + SQ(vect[1] - y) + SQ(vect[2] - z);
    var f = 1/SQ(1 + SQ(s)*r2);

    return f - t;
}

Primitives.convLine = function(x, y, z, begin, end, s, t) {
    var l = sqrt(SQ(end[0] - begin[0])
             + SQ(end[1] - begin[1]) + SQ(end[2] - begin[2]));

    if (l == 0.0) {
        printf("ERROR:Tips of the segment take same coordinate!\n");
        exit(EXIT_FAILURE);
    }

    var ax = (end[0] - begin[0])/l;
    var ay = (end[1] - begin[1])/l;
    var az = (end[2] - begin[2])/l;

    var dx = x - begin[0];
    var dy = x - begin[1];
    var dz = x - begin[2];

    var xx = dx*ax + dy*ay + dz*az;
    var p = Math.sqrt(1 + s*s*(dx*dx + dy*dy + dz*dz - xx*xx));
    var q = Math.sqrt(1 + s*s*(dx*dx + dy*dy + dz*dz - 2*xx));

    var f = xx/(2*p*p*(p*p + s*s*xx*xx)) + (l - xx)/(2*p*p*q*q)
            + (Math.atan(s*xx/p) + Math.atan(s*(l - xx)/p))/(2*s*p*p*p);

    return f - t;
}

Primitives.convArc = function(x, y, z, center, r, theta, axis, angle, s, t) {
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

Primitives.convTriangle = function(x, y, z, vect, s, T) {
    var a1x = vect[0];
    var a1y = vect[1];
    var a1z = vect[2];
    var a2x = vect[3];
    var a2y = vect[4];
    var a2z = vect[5];
    var a3x = vect[6];
    var a3y = vect[7];
    var a3z = vect[8];

    var len1 = Math.sqrt(SQ(a2x - a1x) + SQ(a2y - a1y) + SQ(a2z - a1z));
    var len2 = Math.sqrt(SQ(a3x - a2x) + SQ(a3y - a2y) + SQ(a3z - a2z));
    var len3 = Math.sqrt(SQ(a1x - a3x) + SQ(a1y - a3y) + SQ(a1z - a3z));

    if ((len2 >= len3) && (len2 > len1)) {
        var tempx = a1x;
        var tempy = a1y;
        var tempz = a1z;
        a1x = a2x;
        a1y = a2y;
        a1z = a2z;
        a2x = a3x;
        a2y = a3y;
        a2z = a3z;
        a3x = tempx;
        a3y = tempy;
        a3z = tempz;
    }  else if ((len3 >= len2) && (len3 > len1)) {
        var tempx = a1x;
        var tempy = a1y;
        var tempz = a1z;
        a1x = a3x;
        a1y = a3y;
        a1z = a3z;
        a3x = a2x;
        a3y = a2y;
        a3z = a2z;
        a2x = tempx;
        a2y = tempy;
        a2z = tempz;
    }
    len1 = Math.sqrt(SQ(a2x - a1x) + SQ(a2y - a1y) + SQ(a2z - a1z));
    len2 = Math.sqrt(SQ(a3x - a2x) + SQ(a3y - a2y) + SQ(a3z - a2z));
    len3 = Math.sqrt(SQ(a1x - a3x) + SQ(a1y - a3y) + SQ(a1z - a3z));

    var a21x = a2x - a1x;
    var a21y = a2y - a1y;
    var a21z = a2z - a1z;
    var a13x = a1x - a3x;
    var a13y = a1y - a3y;
    var a13z = a1z - a3z;

    var t = -(a21x*a13x + a21y*a13y + a21z*a13z)/SQ(len1);
    var bx = a1x + t*a21x;
    var by = a1y + t*a21y;
    var bz = a1z + t*a21z;

    var dx = x - bx;
    var dy = y - by;
    var dz = z - bz;

    var ux = a2x - bx;
    var uy = a2y - by;
    var uz = a2z - bz;
    var ul = Math.sqrt(SQ(ux) + SQ(uy) + SQ(uz));
    ux /= ul;
    uy /= ul;
    uz /= ul;

    var vx = a3x - bx;
    var vy = a3y - by;
    var vz = a3z - bz;
    var vl = Math.sqrt(SQ(vx) + SQ(vy) + SQ(vz));
    vx /= vl;
    vy /= vl;
    vz /= vl;

    var d2 = SQ(dx) + SQ(dy) + SQ(dz);
    var u = dx*ux + dy*uy + dz*uz;
    var v = dx*vx + dy*vy + dz*vz;
    var h = Math.sqrt(SQ(a3x - bx) + SQ(a3y - by) + SQ(a3z - bz));
    var a1 = Math.sqrt(SQ(a1x - bx) + SQ(a1y - by) + SQ(a1z - bz));
    var a2 = Math.sqrt(SQ(a2x - bx) + SQ(a2y - by) + SQ(a2z - bz));

    var g = v - h;
    var m = a2*g + u*h;
    var k = u*h + a1*g;
    var C2 = 1/SQ(s) + d2 - SQ(u);
    var C = Math.sqrt(C2);
    var q = C2 - SQ(v);
    var w = C2 - 2*v*h + SQ(h);
    var A2 = SQ(a1)*w + SQ(h)*(q + SQ(u)) - 2*a1*h*u*g;
    var A = Math.sqrt(A2);
    var B2 = SQ(a2)*w + SQ(h)*(q + SQ(u)) - 2*a2*h*u*g;
    var B = Math.sqrt(B2);

    var n1 = a1 + u;
    var n2 = a2 - u;
    var n3 = a1*n1 + v*h;
    var n4 = -a1*u - g*h;
    var n5 = -a2*n2 - v*h;
    var n6 = -a2*u + g*h;

    var arc1 = k*(Math.atan(n3/A) + Math.atan(n4/A))/A;
    var arc2 = m*(Math.atan(n5/B) + Math.atan(n6/B))/B;
    var arc3 = v*(Math.atan(n1/C) + Math.atan(n2/C))/C;
    var f = (arc1 + arc2 + arc3)/(2*q*s);

    return f - T;
}

Primitives.convCurve = function(x, y, z, vect, s, t) {
    var l = Math.sqrt(SQ(vect[0] + 1 - vect[0])
             + SQ(vect[1] + 1 - vect[1]) + SQ(vect[2] + 1 - vect[2]));

    if (l == 0) {
        printf("ERROR:Tips of the segment take same coordinate!\n");
        exit(EXIT_FAILURE);
    }

    var ax = (vect[0] + 1 - vect[0])/l;
    var ay = (vect[1] + 1 - vect[1])/l;
    var az = (vect[2] + 1 - vect[2])/l;

    var dx = x - vect[0];
    var dy = y - vect[1];
    var dz = z - vect[2];

    var xx = dx*ax + dy*ay + dz*az;
    var p = Math.sqrt(1 + s*s*(dx*dx + dy*dy + dz*dz - xx*xx));
    var q = Math.sqrt(1 + s*s*(dx*dx + dy*dy + dz*dz - 2*xx));

    var f = xx/(2*p*p*(p*p + s*s*xx*xx)) + (l - xx)/(2*p*p*q*q)
            + (Math.atan(s*xx/p) + Math.atan(s*(l - xx)/p))/(2*s*p*p*p);

    return f - t;
}

Primitives.convMesh = function(x, y, z, vect, tri, s, t) {
    var a1x = vect[3*(tri[3] - 1)];
    var a1y = vect[3*(tri[3] - 1) + 1];
    var a1z = vect[3*(tri[3] - 1) + 2];
    var a2x = vect[3*(tri[3 + 1] - 1)];
    var a2y = vect[3*(tri[3 + 1] - 1) + 1];
    var a2z = vect[3*(tri[3 + 1] - 1) + 2];
    var a3x = vect[3*(tri[3 + 2] - 1)];
    var a3y = vect[3*(tri[3 + 2] - 1) + 1];
    var a3z = vect[3*(tri[3 + 2] - 1) + 2];

    var len1 = Math.sqrt(SQ(a2x - a1x) + SQ(a2y - a1y) + SQ(a2z - a1z));
    var len2 = Math.sqrt(SQ(a3x - a2x) + SQ(a3y - a2y) + SQ(a3z - a2z));
    var len3 = Math.sqrt(SQ(a1x - a3x) + SQ(a1y - a3y) + SQ(a1z - a3z));

    if ((len2 >= len3) && (len2 > len1)) {
        var tempx = a1x;
        var tempy = a1y;
        var tempz = a1z;
        a1x = a2x;
        a1y = a2y;
        a1z = a2z;
        a2x = a3x;
        a2y = a3y;
        a2z = a3z;
        a3x = tempx;
        a3y = tempy;
        a3z = tempz;
    }  else if ((len3 >= len2) && (len3 > len1)) {
        var tempx = a1x;
        var tempy = a1y;
        var tempz = a1z;
        a1x = a3x;
        a1y = a3y;
        a1z = a3z;
        a3x = a2x;
        a3y = a2y;
        a3z = a2z;
        a2x = tempx;
        a2y = tempy;
        a2z = tempz;
    }
    len1 = Math.sqrt(SQ(a2x - a1x) + SQ(a2y - a1y) + SQ(a2z - a1z));
    len2 = Math.sqrt(SQ(a3x - a2x) + SQ(a3y - a2y) + SQ(a3z - a2z));
    len3 = Math.sqrt(SQ(a1x - a3x) + SQ(a1y - a3y) + SQ(a1z - a3z));

    var a21x = a2x - a1x;
    var a21y = a2y - a1y;
    var a21z = a2z - a1z;
    var a13x = a1x - a3x;
    var a13y = a1y - a3y;
    var a13z = a1z - a3z;

    var t = -(a21x*a13x + a21y*a13y + a21z*a13z)/SQ(len1);
    var bx = a1x + t*a21x;
    var by = a1y + t*a21y;
    var bz = a1z + t*a21z;

    var dx = x - bx;
    var dy = y - by;
    var dz = z - bz;

    var ux = a2x - bx;
    var uy = a2y - by;
    var uz = a2z - bz;
    var ul = Math.sqrt(SQ(ux) + SQ(uy) + SQ(uz));
    ux /= ul;
    uy /= ul;
    uz /= ul;

    var vx = a3x - bx;
    var vy = a3y - by;
    var vz = a3z - bz;
    var vl = Math.sqrt(SQ(vx) + SQ(vy) + SQ(vz));
    vx /= vl;
    vy /= vl;
    vz /= vl;

    var d2 = SQ(dx) + SQ(dy) + SQ(dz);
    var u = dx*ux + dy*uy + dz*uz;
    var v = dx*vx + dy*vy + dz*vz;
    var h = Math.sqrt(SQ(a3x - bx) + SQ(a3y - by) + SQ(a3z - bz));
    var a1 = Math.sqrt(SQ(a1x - bx) + SQ(a1y - by) + SQ(a1z - bz));
    var a2 = Math.sqrt(SQ(a2x - bx) + SQ(a2y - by) + SQ(a2z - bz));

    var g = v - h;
    var m = a2*g + u*h;
    var k = u*h + a1*g;
    var C2 = 1/SQ(s) + d2 - SQ(u);
    var C = Math.sqrt(C2);
    var q = C2 - SQ(v);
    var w = C2 - 2*v*h + SQ(h);
    var A2 = SQ(a1)*w + SQ(h)*(q + SQ(u)) - 2*a1*h*u*g;
    var A = Math.sqrt(A2);
    var B2 = SQ(a2)*w + SQ(h)*(q + SQ(u)) - 2*a2*h*u*g;
    var B = Math.sqrt(B2);

    var n1 = a1 + u;
    var n2 = a2 - u;
    var n3 = a1*n1 + v*h;
    var n4 = -a1*u - g*h;
    var n5 = -a2*n2 - v*h;
    var n6 = -a2*u + g*h;

    var arc1 = k*(Math.atan(n3/A) + Math.atan(n4/A))/A;
    var arc2 = m*(Math.atan(n5/B) + Math.atan(n6/B))/B;
    var arc3 = v*(Math.atan(n1/C) + Math.atan(n2/C))/C;
    var f = (arc1 + arc2 + arc3)/(2*q*s);

    return f - T;
}

Primitives.noiseG = function(x, y, z, amp, frep, phase) {
    var a2 = amp;
    var a1 = frep;

    var xt = x;
    var yt = y;
    var zt = z;

    var a2x = a2*xt;
    var a2y = a2*yt;
    var a2z = a2*zt;
    var sx = Math.sin(a2x);
    var sy = Math.sin(a2y);
    var sz = Math.sin(a2z);
    var a1d = a1/1.17;
    var sx2 = a1d*Math.sin(a2x/1.35 + phase*sz);
    var sy2 = a1d*Math.sin(a2y/1.35 + phase*sx);
    var sz2 = a1d*Math.sin(a2z/1.35 + phase*sy);
    var Serx = a1*sx + sx2;
    var Sery = a1*sy + sy2;
    var Serz = a1*sz + sz2;
    var SS = Serx*Sery*Serz;

    return SS;
}