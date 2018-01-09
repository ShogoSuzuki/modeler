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

Attributes.Sstep = function(s, smin, smax, f, fmin, fmax) {
    if (f > fmax) return 1.0;
    if (f < fmin) return 0.0;
    if (s > smax) return smax;
    if (s < smin) return smin; 
}

