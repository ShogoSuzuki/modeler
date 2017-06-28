var model = {};

model.eval = function(x, y, z) {
// cap
    // brim(futi)
    var center = [0.0, 0.0, 0.0];
    var px = x / 0.9;
    var py = y;
    var pz = z / 1.05;

    qy = py * Math.cos(-0.9) + pz * Math.sin(-0.9);
    qz = -py * Math.sin(-0.9) + pz * Math.cos(-0.9);

    py = qy + 3.3;
    pz = qz + 3.5;

    var brim = Primitives.torusY(px, py, pz, center, 2.9, 0.5);

    // steeple(tongari)
    px = x;
    py = y * Math.cos(-0.9) + z * Math.sin(-0.9);
    pz = -y * Math.sin(-0.9) + z * Math.cos(-0.9);

    py = py + 2.8;
    pz = pz + 3.5;

    var cap1 = Primitives.ellipsoid(px, py, pz, center, 3.1, 4.0, 3.35);
    cap1 = Operations.rSub(cap1, (py - 0.5));

    center = [0.0, 2.2, -3.3];
    var cap2 = Primitives.ellipsoid(x, 0, 0, center, 1.4, 3.2, 0.6);
    cap2 = Operations.rSub(cap2, (y - 5.5));
    var cap3 = Operations.rBUnion(cap1, cap2, 0.1, 0.1, 0.1);

    // tassel
    center = [0.0, 2.2, -3.9];
    var tassel = Primitives.sphere(x, 0, 0, center, 0.8);
    var cap = Operations.rUnion(cap3, brim);
    cap = Operations.rUnion(cap, tassel);

// ears
py = y * Math.cos(-0.3) + z * Math.sin(-0.3);
pz = -y * Math.sin(-0.3) + z * Math.cos(-0.3);

py = py + 1.7;
pz = pz + 0.27;

center = [2.5, 5.6, 1.0];
var ear1 = Primitives.ellipsoid(px, py, pz, center, 1.5, 1.5, 0.9);
center = [2.8, 5.6, 1.3];
var ear2 = Primitives.ellipsoid(px, py, pz, center, 1.9, 1.9, 0.8);

center = [-2.5, 5.6, 1.0];
var ear3 = Primitives.ellipsoid(px, py, pz, center, 1.5, 1.5, 0.9);
center = [-2.8, 5.6, 1.3];
var ear4 = Primitives.ellipsoid(px, py, pz, center, 1.9, 1.9, 0.8);

var earLeft = Operations.rSub(ear1, ear2);
var earRight = Operations.rSub(ear3, ear4);
var ear = Operations.rUnion(earLeft, earRight);

// hear
center = [0.0, 4.0, 0.5];
var head1 = Primitives.ellipsoid(x, 0, 0, center, 2.7, 3.0, 2.8);

center = [0.0, 5.0, 0.9];
var head2 = Primitives.ellipsoid(x, 0, 0, center, 2.0, 2.0, 2.1);
var head3 = 4.0 - y;
var head4 = z - 2.5;
var hear = Operations.rSub(head1, head2);
var hear2 = Operations.rInter(head3, head4);
hear =  Operations.rSub(hear, hear2);
hear = Operations.rSub(hear, (2.3 - y));

// face
center = [0.0, 4.0, 0.5];
var face1 = Primitives.ellipsoid(x, 0, 0, center, 2.3, 2.9, 2.4);
var face = Operations.rSub(face1, hear2);
face = Operations.rSub(face, (2.3 - y));

// eyes
center = [0.0, 4.0, 0.6];
var eye1 = Primitives.ellipsoid(x, 0, 0, center, 2.3, 2.9, 2.4);

center = [0.5, 5.5, 0.0];
var eye2 = Primitives.ellCylZ(x, 0, 0, center, 0.4, 0.7);
var eye3 = Operations.rInter(eye1, eye2);

center = [-0.5, 5.5, 0.0];
var eye4 = Primitives.ellCylZ(x, 0, 0, center, 0.4, 0.7);
var eye5 = Operations.rInter(eye1, eye4);

    // iris
    center = [0.5, 5.2, 2.65];
    var eye6 = Primitives.sphere(x, 0, 0, center, 0.28);
    var eye7 = Operations.rInter(eye3, -eye6);

    center = [-0.5, 5.2, 2.65];
    var eye8 = Primitives.sphere(x, 0, 0, center, 0.28);
    var eye9 = Operations.rInter(eye5, -eye8);

// cheek and mouth
px = x;
py = y / 0.7;
pz = z / 2.0;

center = [1.55, 5.7, 0.6];
var cheek1 = Primitives.torusZ(px, py, pz, center, 0.5, 0.5);
cheek1 = Operations.rSub(cheek1, (5.6 - pz));
center = [-1.55, 5.7, 0.6];
var cheek2 = Primitives.torusZ(px, py, pz, center, 0.5, 0.5);
cheek2 = Operations.rSub(cheek1, (5.6 - pz));
var cheek3 = Operations.rUnion(cheek1, cheek2);

px = x;
py = y / 0.7;
pz = z / 2.1;

center = [0.0, 5.6, 0.7];
var cheek4 = Primitives.torusZ(px, py, pz, center, 2.0, 0.5);
cheek4 = Operations.rSub(cheek4, (y - 4.4));

var cheek = Operations.rBUnion(cheek3, cheek4, 0.2, 1.0, 1.0);

    // mouth
    center = [0.0, 2.0, 3.3];
    py = y * Math.cos(-0.3) + z * Math.sin(-0.3);
    pz = -y * Math.sin(-0.3) + z * Math.cos(-0.3);
    var mouth1 = Primitives.sphere(px, py, pz, center, 0.7);
    mouth1 = Operations.rSub(mouth1, (pz - 3.4));
    var mouth2 = Primitives.sphere(px, py, pz, center, 0.5);

    var mouth3 = Operations.rSub(mouth1, mouth2);
    mouth3 = Operations.rUnion(mouth3, cheek4);
    var mouth = Operations.rBUnion(cheek3, mouth3, 0.3, 0.8, 0.8);

// nose
center = [0.0, 2.7, 0.0];
py = y * Math.cos(-0.4) + z * Math.sin(-0.4);
pz = -y * Math.sin(-0.4) + z * Math.cos(-0.4);
var nose1 = Primitives.ellipsoid(px, py, pz, center, 1.97, 1.5, 5.5);
nose1 = Operations.rSub(nose1, (1.0 - pz));

center = [0.0, 5.0, 3.9];
var nose2 = Primitives.ellipsoid(x, 0, 0, center, 0.5, 0.7, 0.5);

var nose = Operations.rUnion(nose1, nose2);

// neck
center = [0.0, 1.9, 0.0];
var neck = Primitives.torusY(x, 0, 0, center, 1.0, 0.8);

var head = Operations.rUnion(hear, face);
head = Operations.rUnion(head, eye7);
head = Operations.rUnion(head, eye9);
head = Operations.rUnion(head, neck);
head = Operations.rUnion(head, mouth);
head = Operations.rUnion(head, nose);

// breast
center = [0.0, 0.1, -1.7];
var body1 = Primitives.ellipsoid(x, 0, 0, center, 2.0, 1.6, 3.4);
body1 = Operations.rSub(body1, ((-1.7) - z));

// arm
center = [0.0, -1.8, -1.6];

px = x;
py = y / 2;
pz = z;

qz = pz * Math.cos(0.9) + px * Math.sin(0.9);
qx = -pz * Math.sin(0.9) + px * Math.cos(0.9);

px = qx;
pz = qz;
var arm1 = Primitives.torusX(px, py, pz, center, 2.0, 0.7);
arm1 = Operations.rSub(arm1, ((-2.0) - py));
arm1 = Operations.rSub(arm1, (pz + 1.6));

px = x;
py = y / 2;
pz = z;

qz = pz * Math.cos(-0.9) + px * Math.sin(-0.9);
qx = -pz * Math.sin(-0.9) + px * Math.cos(-0.9);

px = qx;
pz = qz;
var arm2 = Primitives.torusX(px, py, pz, center, 2.0, 0.7);
arm2 = Operations.rSub(arm2, ((-2.0) - py));
arm2 = Operations.rSub(arm2, (pz + 1.6));

    // cuffs
    center = [2.6, -3.0, -2.2];
    var cuff1 = Primitives.torusY(x, 0, 0, center, 0.6, 0.6);
    center = [-2.6, -3.0, -2.2];
    var cuff2 = Primitives.torusY(x, 0, 0, center, 0.6, 0.6);

    // hands
    center = [2.6, -4.3, -2.0];
    var hand1 = Primitives.ellipsoid(x, 0, 0, center, 0.8, 1.7, 2.5);
    hand1 = Operations.rSub(hand1, ((-4.6) - y));

    center = [-2.6, -4.3, -2.0];
    var hand2 = Primitives.ellipsoid(x, 0, 0, center, 0.8, 1.7, 2.5);
    hand2 = Operations.rSub(hand2, ((-4.6) - y));

    var armLeft = Operations.rUnion(arm1, cuff1);
    armLeft = Operations.rUnion(armLeft, hand1);
    var armRight = Operations.rUnion(arm2, cuff2);
    armRight = Operations.rUnion(armRight, hand2);
    var arm = Operations.rUnion(armLeft, armRight);

// lower half of body
    //　belt
    center = [0.0, -4.2, -0.5];
    var belt1 = Primitives.torusY(x, 0, 0, center, 1.9, 0.8);
    belt1 = Operations.rSub(belt1, ((-1.6) - z));
    var belt2 = Primitives.cylY(x, 0, 0, center, 1.8);

    var belt = Operations.rInter(belt1, belt2);

    // belly
    center = [0.0, -4.2, -2.0];
    var belly1 = Primitives.ellipsoid(x, 0, 0, center, 4.0, 2.9, 5.0);
    belly1 = Operations.rSub(belly1, ((-4.6) - y));
    belly1 = Operations.rSub(belly1, ((-2.0) - z));
    var belly2 = Primitives.ellipsoid(x, 0, 0, center, 4.0, 2.9, 1.0);
    belly2 = Operations.rSub(belly2, ((-4.6) - y));

    var belly = Operations.rUnion(belly1, belly2);

    var body2 = Operations.rUnion(belt, belly);
    var body = Operations.rBUnion(body1, body2, 1.0, 0.1, 0.1);

// leg
center = [0.0, -4.0, 19.8];
var leg1 = Primitives.torusY(x, 0, 0, center, 20.0, 1.5);
leg1 = Operations.rInter(leg1, (x + 4.5));
leg1 = Operations.rInter(leg1, (4.5 - x));

center = [4.5, -4.0, 0.0];
var hem1 = Primitives.torusX(x, 0, 0, center, 1.6, 0.4);
center = [-4.5, -4.0, 0.0];
var hem2 = Primitives.torusX(x, 0, 0, center, 1.6, 0.4);

center = [5.9, -3.0, 0.0];
var boots1 = Primitives.ellipsoid(x, 0, 0, center, 1.5, 3.5, 1.9);
boots1 = Operations.rSub(boots1, (x - 5.9));
boots1 = Operations.rSub(boots1, ((-6.0) - y));
center = [-5.9, -3.0, 0.0];
var boots2 = Primitives.ellipsoid(x, 0, 0, center, 1.5, 3.5, 1.9);
boots2 = Operations.rSub(boots2, ((-5.9) - x));
boots2 = Operations.rSub(boots2, ((-6.0) - y));

var leg = Operations.rUnion(leg1, hem1);
leg = Operations.rUnion(leg, hem2);
leg = Operations.rUnion(leg, boots1);
leg = Operations.rUnion(leg, boots2);

var my_model = Operations.rUnion(cap, ear);
my_model = Operations.rUnion(my_model, head);
my_model = Operations.rUnion(my_model, body);
my_model = Operations.rUnion(my_model, arm);
my_model = Operations.rUnion(my_model, leg);

return my_model;
}

model.boundingBox = function() {
    var bbox = {x_min: -7, x_max: 7, 
		y_min: -7, y_max: 7,
		z_min: -7, z_max: 7}

    return bbox;
}