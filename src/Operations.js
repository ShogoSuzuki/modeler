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
Operations.rBUnion = function(f1, f2, a, b, c) {
    var bt = f1 / b;
    var ct = f2 / c;

    return Operations.rUnion(f1, f2) + (a / (1 + bt*bt + ct*ct));
}
