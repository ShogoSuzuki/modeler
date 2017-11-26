var Operations = {};


Operations.rUnion = function(f1, f2) {
    return f1 + f2 + Math.sqrt(f1*f1 + f2*f2);
}


Operations.rInter = function(f1, f2) {
    return f1 + f2 - Math.sqrt(f1*f1 + f2*f2);
}


Operations.rSub = function(f1, f2) {
    return f1 - f2 - Math.sqrt(f1*f1 + f2*f2);
}

// add Operations
Operations.rBlendUni = function(f1, f2, a, b, c) {
    var bt = f1 / b;
    var ct = f2 / c;

    return Operations.rUnion(f1, f2)
     + (a / (1 + bt*bt + ct*ct));
}

Operations.rScale3d = function(x, y, z, sx, sy, sz) {
    var x = x/sx;
    var y = y/sy;
    var z = z/sz;

    return 1.0;
} 

Operations.rShift3d = function(x, y, z, dx, dy, dz) {
    var x = x - dx;
    var y = y - dy;
    var z = z - dz;

    return 1.0;
}

Operations.rRotate3dX = function(x, y, z, theta) {
    var ct = Math.cos(theta);
    var st = Math.sin(theta);

    var yr = y*ct + z*st;
    var zr = y*st + z*ct;

    y = yr;
    z = zr;

    return 1.0;
}

Operations.rRotate3dY = function(x, y, z, theta) {
    var ct = Math.cos(theta);
    var st = Math.sin(theta);

    var zr = z*ct + x*st;
    var xr = z*st + x*ct;

    x = xr;
    z = zr;

    return 1.0;
}

Operations.rRotate3dZ = function(x, y, z, theta) {
    var ct = Math.cos(theta);
    var st = Math.sin(theta);

    var xr = x*ct + y*st;
    var yr = x*st + y*ct;

    x = xr;
    y = yr;

    return 1.0;
}

Operations.rTwistX = function(x, y, z, z1, z2, theta1, theta2) {
    var x2 = z2;
    var x1 = z1;

    var t = (x - x1)/(x2 - x1);
    var theta = (1 - t)*theta1 + t*theta2;
    var ct = Math.cos(theta);
    var st = Math.sin(theta);

    var yr = y*ct + z*st;
    var zr = -y*st + z*ct;

    y = yr;
    z = zr;

    return 1.0;
}

Operations.rTwistY = function(x, y, z, z1, z2, theta1, theta2) {
    var y2 = z2;
    var y1 = z1;

    var t = (y - y1)/(y2 - y1);
    var theta = (1 - t)*theta1 + t*theta2;
    var ct = Math.cos(theta);
    var st = Math.sin(theta);

    var zr = z*ct + x*st;
    var xr = -z*st + x*ct;

    x = xr;
    z = zr;

    return 1.0;
}

Operations.rTwistZ = function(x, y, z, z1, z2, theta1, theta2) {
    var z2 = z2;
    var z1 = z1;

    var t = (z - z1)/(z2 - z1);
    var theta = (1 - t)*theta1 + t*theta2;
    var ct = Math.cos(theta);
    var st = Math.sin(theta);

    var xr = x*ct + y*st;
    var yr = -x*st + y*ct;

    x = xr;
    y = yr;

    return 1.0;
}

Operations.rSTretch3d = function(x, y, z, x0, sx, sy, sz) {
    x = x0[0] + (x - x0[0])/sx;
    y = x0[1] + (y - x0[1])/sy;
    z = x0[2] + (z - x0[2])/sz;

    return 1.0;
}

Operations.rTaaperX = function(x, y, z, x1, x2, s1, s2) {
    var scale;

    if (x < x1) {
        scale = s1;
    } else {
        if (x > x2) {
            scale = s2;
        } else {
            var t = (x - x1)/(x2 - x1);
            scale = (1 - t)*s1 + t*s2;
        }
    }
    if (Math.abs(scale) < EPS) {
        scale = 1.0;
    }
    y = y/scale;
    z = z/scale;

    return 1.0;
}

Operations.rTaaperY = function(x, y, z, x1, x2, s1, s2) {
    var scale;

    if (y < x1) {
        scale = s1;
    } else {
        if (y > x2) {
            scale = s2;
        } else {
            var t = (y - x1)/(x2 - x1);
            scale = (1 - t)*s1 + t*s2;
        }
    }
    if (Math.abs(scale) < EPS) {
        scale = 1.0;
    }
    z = z/scale;
    x = x/scale;

    return 1.0;
}

Operations.rTaaperZ = function(x, y, z, x1, x2, s1, s2) {
    var scale;

    if (z < x1) {
        scale = s1;
    } else {
        if (z > x2) {
            scale = s2;
        } else {
            var t = (z - x1)/(x2 - x1);
            scale = (1 - t)*s1 + t*s2;
        }
    }
    if (Math.abs(scale) < EPS) {
        scale = 1.0;
    }
    x = x/scale;
    y = y/scale;

    return 1.0;
}

Operations.rSpaceMapCubic = function(x, y, z, original_points,
     delta_points, b) {
        var fratio = b;
        var to = delta_points;
        var from = original_points;
    
        var dtx = (to[0] - from[0]);
        var dty = (to[1] - from[1]);
        var dtz = (to[2] - from[2]);
    
        var xd = (x - to[0])/(fratio[0]*(1.0 + Math.sqrt(dtx*dtx)));
        var yd = (y - to[1])/(fratio[1]*(1.0 + Math.sqrt(dty*dty)));
        var zd = (z - to[2])/(fratio[2]*(1.0 + Math.sqrt(dtz*dtz)));
    
        var xt = xd*xd;
        var yt = yd*yd;
        var zt = zd*zd;
    
        var tmp = Math.sqrt(xt + yt + zt);
        if (tmp <= 1.0) {
            tmp = ((1 - tmp)*(1 - tmp))*((1 + tmp)*(1 + tmp));
        
            x = x - tmp*dtx;
            y = y - tmp*dty;
            z = z - tmp*dtz;
        }
    
        return 1.0;
}

Operations.rBlendInt = function(f1, f2, a, b, c) {
    var bt = f1/b;
    var ct = f2/c;

    return Operations.rInter(f1, f2)
     + (a/1 + bt*bt + ct*ct);
}

Operations.rSpaceMapExp = function(x, y, z, original_points,
    delta_points, b) {
       var fratio = b;
       var to = delta_points;
       var from = original_points;
   
       var dtx = (to[0] - from[0]);
       var dty = (to[1] - from[1]);
       var dtz = (to[2] - from[2]);
   
       var xd = (x - to[0])/(fratio[0]*(1.0 + Math.sqrt(dtx*dtx)));
       var yd = (y - to[1])/(fratio[1]*(1.0 + Math.sqrt(dty*dty)));
       var zd = (z - to[2])/(fratio[2]*(1.0 + Math.sqrt(dtz*dtz)));
   
       var xt = xd*xd;
       var yt = yd*yd;
       var zt = zd*zd;
   
       var tmp = xt + yt + zt;
       tmp = Math.exp(-tmp);
       
        x = x - tmp*dtx;
        y = y - tmp*dty;
        z = z - tmp*dtz;
   
        return 1.0;
}

Operations.rMapBlob = function(x, y, z, fobj, x0, y0, z0,
     fobj0, sigma) {
        var xt = x - x0;
        var yt = y - y0;
        var zt = z - z0;
        if (sigma > EPS) {
            var d2 = (xt*xt + yt*yt + zt*zt)/(sigma*sigma);
            fobj = fobj - fobj0/(1 + d2);
        }
    
        return fobj;
}