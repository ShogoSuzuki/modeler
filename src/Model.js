// The model to be rendered. It contains:
// - its defining function
// - its bounding-box
var model = {};

/*
model.eval = function(x, y, z) {
    var value = x*x + y*y + z*z - 25;
    return value;
}

model.boundingBox = function() {
    var bbox = {x_min: -10, x_max: 10, 
		y_min: -10, y_max: 10,
		z_min: -10, z_max: 10}

    return bbox;
}
*/


model.eval = function(x, y, z) {
    var px = x / 8.0;
    var py = y / 8.0;
    var pz = z / 8.0;

    var center = [0.0, 0.0, 0.0];

    // Hull
    var hull = Primitives.sphere(px, py, pz, center, 0.4);

    // Tail
    center = [0.0, 0.18, -0.3];
    var tail = Primitives.ellipsoid(px, py, pz, center, 0.18, 0.15, 0.50);
    tail = Operations.rSub(tail, pz);

    // Floats
    center = [0.2, -0.3, 0.0];
    var float1 = Primitives.ellipsoid(px, py, pz, center, 0.15, 0.15, 0.5);
    center = [-0.2, -0.3, 0.0];
    var float2 = Primitives.ellipsoid(px, py, pz, center, 0.15, 0.15, 0.5);

    // Stabilizer and wings
    center = [0.0, 0.2, -0.7];
    var stab = Primitives.ellCylX(px, py, pz, center, 0.05, 0.12);
    stab = Operations.rInter(stab, (px+0.4));
    stab = Operations.rInter(stab, (-px+0.4));
    center = [0.0, 0.16, 0.0];
    var wings = Primitives.ellCylX(px, py, pz, center, 0.05, 0.14);
    wings = Operations.rInter(wings, (px+0.7));
    wings = Operations.rInter(wings, (-px+0.7));

    // Propeller
    center = [0.0, 0.4, 0.0];
    var prop1 = Primitives.ellCylZ(px, py, pz, center, 0.08, 0.04);
    center = [0.0, 0.4, 0.0];
    var prop2 = Primitives.ellCylX(px, py, pz, center, 0.04, 0.08);
    center = [0.0, 0.0, 0.0];
    var prop = Operations.rUnion(prop1, prop2);
    var temp = Primitives.ellCylY(px, py, pz, center, 0.6, 0.7);
    prop = Operations.rInter(prop, temp);
    var axle = Primitives.ellCylY(px, py, pz, center, 0.05, 0.05);
    axle = Operations.rInter(axle, (-py+0.5));
    axle = Operations.rInter(axle, py);

    // Eyes
    center = [-0.16, 0.0, 0.36];
    var eye1 = Primitives.sphere(px, py, pz, center, 0.08);
    center = [0.16, 0.0, 0.36];
    var eye2 = Primitives.sphere(px, py, pz, center, 0.08);

    // Final model
    var my_model = Operations.rUnion(hull, tail);
    my_model = Operations.rUnion(my_model, float1);
    my_model = Operations.rUnion(my_model, float2);
    my_model = Operations.rUnion(my_model, stab);
    my_model = Operations.rUnion(my_model, wings);
    my_model = Operations.rUnion(my_model, prop);
    my_model = Operations.rUnion(my_model, axle);
    my_model = Operations.rUnion(my_model, eye1);
    my_model = Operations.rUnion(my_model, eye2);

    return my_model;
}


model.boundingBox = function() {
    var bbox = {x_min: -7, x_max: 7, 
		y_min: -7, y_max: 7,
		z_min: -7, z_max: 7}

    return bbox;
}
