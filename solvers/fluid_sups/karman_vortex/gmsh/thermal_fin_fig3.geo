// ============================================================
// Gmsh model corresponding to Fig. 3 in Gelsomino & Rozza (2011)
// x = x1, y = x2, z = x3
//
// Paper geometry (non-dimensional):
//   Omega01 (spreader): 0 <= x <= 1/2,
//                       0 <= y <= 3/5,
//                      -1 <= z <= 0
//   Omega02 (fin):      0 <= x <= 3/20,
//                       3/5 <= y <= 3/5 + mu2,
//                      -1 <= z <= 0
//
// To obtain dimensional coordinates in meters, set dper to the
// dimensional fin pitch/period used for non-dimensionalization.
// ============================================================

SetFactory("Built-in");

// -------------------------
// Parameters
// -------------------------
dper = 1.0;       // [m] scale factor; dper=1 reproduces Fig. 3 numerically
mu2  = 2.0;       // non-dimensional fin height L/dper; paper reference value is 2

x_base = (1.0/2.0)  * dper;
y_base = (3.0/5.0)  * dper;
x_fin  = (3.0/20.0) * dper;
h_fin  = mu2 * dper;
depth  = 1.0 * dper;

// Structured mesh resolution (number of points on each curve)
nx_fin  = 8;
nx_rest = 19;
ny_base = 21;
ny_fin  = 41;
nz      = 21;     // number of points through the depth

// -------------------------
// Points on z = 0 plane
// -------------------------
Point(1) = {0,      0,              0};
Point(2) = {x_fin,  0,              0};
Point(3) = {x_base, 0,              0};

Point(4) = {0,      y_base,         0};
Point(5) = {x_fin,  y_base,         0};
Point(6) = {x_base, y_base,         0};

Point(7) = {0,      y_base + h_fin, 0};
Point(8) = {x_fin,  y_base + h_fin, 0};

// -------------------------
// Curves
// -------------------------
// Omega01: left block (under the fin)
Line(1) = {1, 2};
Line(2) = {2, 5};
Line(3) = {5, 4};
Line(4) = {4, 1};

// Omega01: right block
Line(5) = {2, 3};
Line(6) = {3, 6};
Line(7) = {6, 5};

// Omega02: fin
Line(8)  = {4, 7};
Line(9)  = {7, 8};
Line(10) = {8, 5};

// -------------------------
// 2D block surfaces
// -------------------------
// Split Omega01 into two rectangles only for a fully structured mesh.
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Curve Loop(2) = {5, 6, 7, -2};
Plane Surface(2) = {2};

Curve Loop(3) = {8, 9, 10, 3};
Plane Surface(3) = {3};

// -------------------------
// Transfinite spacings
// -------------------------
// x direction under fin / fin thickness
Transfinite Curve {1, 3, 9} = nx_fin Using Progression 1;

// x direction on the remaining spreader width
Transfinite Curve {5, 7} = nx_rest Using Progression 1;

// y direction in spreader
Transfinite Curve {2, 4, 6} = ny_base Using Progression 1;

// y direction in fin
Transfinite Curve {8, 10} = ny_fin Using Progression 1;

Transfinite Surface {1, 2, 3};
Recombine Surface {1, 2, 3};

// -------------------------
// Extrude in x3 = z direction
// -------------------------
// Fig. 3 uses -1 <= x3 <= 0.
// Layers takes the number of elements, hence nz-1.
out[] = Extrude {0, 0, -depth} {
  Surface{1, 2, 3};
  Layers{nz-1};
  Recombine;
};

// -------------------------
// Physical groups
// -------------------------
// Use geometric queries instead of hard-coded extrusion entity IDs.
eps = 1.e-8 * dper;

// Volumes
Omega01() = Volume In BoundingBox {
  -eps, -eps, -depth-eps,
  x_base+eps, y_base+eps, eps
};

Omega02() = Volume In BoundingBox {
  -eps, y_base-eps, -depth-eps,
  x_fin+eps, y_base+h_fin+eps, eps
};

Physical Volume("Omega01_spreader") = {Omega01()};
Physical Volume("Omega02_fin")      = {Omega02()};
Physical Volume("Fluid")            = {Omega01(), Omega02()};

// Boundary naming follows Fig. 3(b).
// Gamma01: spreader face z = 0
Gamma01() = Surface In BoundingBox {
  -eps, -eps, -eps,
  x_base+eps, y_base+eps, eps
};

// Gamma02: spreader base y = 0 (uniform heat flux in the paper)
Gamma02() = Surface In BoundingBox {
  -eps, -eps, -depth-eps,
  x_base+eps, eps, eps
};

// Gamma03: spreader side x = 1/2
Gamma03() = Surface In BoundingBox {
  x_base-eps, -eps, -depth-eps,
  x_base+eps, y_base+eps, eps
};

// Gamma04: exposed top of spreader, y = 3/5 and 3/20 <= x <= 1/2
Gamma04() = Surface In BoundingBox {
  x_fin-eps, y_base-eps, -depth-eps,
  x_base+eps, y_base+eps, eps
};

// Gamma05: spreader face z = -1
Gamma05() = Surface In BoundingBox {
  -eps, -eps, -depth-eps,
  x_base+eps, y_base+eps, -depth+eps
};

// Gamma06: spreader-fin interface, y = 3/5 and 0 <= x <= 3/20
Gamma06() = Surface In BoundingBox {
  -eps, y_base-eps, -depth-eps,
  x_fin+eps, y_base+eps, eps
};

// Gamma07: fin face z = 0
Gamma07() = Surface In BoundingBox {
  -eps, y_base-eps, -eps,
  x_fin+eps, y_base+h_fin+eps, eps
};

// Gamma08: fin tip y = 3/5 + mu2
Gamma08() = Surface In BoundingBox {
  -eps, y_base+h_fin-eps, -depth-eps,
  x_fin+eps, y_base+h_fin+eps, eps
};

// Gamma09: fin face z = -1
Gamma09() = Surface In BoundingBox {
  -eps, y_base-eps, -depth-eps,
  x_fin+eps, y_base+h_fin+eps, -depth+eps
};

// Gamma10: exposed fin side x = 3/20 (Robin convection in the paper)
Gamma10() = Surface In BoundingBox {
  x_fin-eps, y_base-eps, -depth-eps,
  x_fin+eps, y_base+h_fin+eps, eps
};

// Gamma11: fin symmetry/adiabatic side x = 0
Gamma11() = Surface In BoundingBox {
  -eps, y_base-eps, -depth-eps,
  eps, y_base+h_fin+eps, eps
};

// Gamma12: spreader symmetry/adiabatic side x = 0
Gamma12() = Surface In BoundingBox {
  -eps, -eps, -depth-eps,
  eps, y_base+eps, eps
};

Physical Surface("Gamma01") = {Gamma01()};
Physical Surface("Gamma02") = {Gamma02()};
Physical Surface("Gamma03") = {Gamma03()};
Physical Surface("Gamma04") = {Gamma04()};
Physical Surface("Gamma05") = {Gamma05()};
Physical Surface("Gamma06_interface") = {Gamma06()};
Physical Surface("Gamma07") = {Gamma07()};
Physical Surface("Gamma08") = {Gamma08()};
Physical Surface("Gamma09") = {Gamma09()};
Physical Surface("Gamma10") = {Gamma10()};
Physical Surface("Gamma11") = {Gamma11()};
Physical Surface("Gamma12") = {Gamma12()};

// Convenient groups corresponding to the thermal boundary conditions in the paper
Physical Surface("Heat_flux_Gamma02") = {Gamma02()};
Physical Surface("Convection_Gamma10") = {Gamma10()};
Physical Surface("Adiabatic") = {
  Gamma01(), Gamma03(), Gamma04(), Gamma05(),
  Gamma07(), Gamma08(), Gamma09(), Gamma11(), Gamma12()
};

// External boundary only (Gamma06 is internal and is intentionally excluded)
Physical Surface("All_external_boundaries") = {
  Gamma01(), Gamma02(), Gamma03(), Gamma04(), Gamma05(),
  Gamma07(), Gamma08(), Gamma09(), Gamma10(), Gamma11(), Gamma12()
};

// -------------------------
// Mesh and output
// -------------------------
Mesh 3;
Save "thermal_fin_fig3.msh";
