// NASA TMR 2DBFS / Driver-Seegmiller backward-facing step mesh
// Gmsh .geo for stabilized FEM / VMS / Nitsche / wall-function tests
// Units: nondimensional, step height h = 1.
// Geometry follows the common NASA TMR 2DBFS layout:
//   upstream channel:   x/h = [-20, 0], y/h = [1, 9]
//   downstream channel: x/h = [0, 50],  y/h = [0, 9]
//   step at x/h = 0, vertical wall y/h = [0, 1]
// The usual BFS expansion ratio is downstream height / upstream height = 9/8 = 1.125.
// Re_H for the validation case is approximately 36000, but this file only defines the mesh.
// TMR profile comparison stations x/h = 1, 4, 6, 10 are explicit block boundaries.

SetFactory("Built-in");

// -----------------------
// Geometry
// -----------------------
h  = 1.0;       // step height
H  = 9.0*h;     // downstream channel height; upstream channel height is H-h = 8h
Lu = 110.0*h;    // inlet location x/h = -20
Ld = 50.0*h;    // outlet location x/h = 50

// Profile stations used by NASA TMR data
x1 = 1.0*h;
x4 = 4.0*h;
x6 = 6.0*h;
x10 = 10.0*h;

// Quasi-2D extrusion. Set extrude_width = O(h), nz > 1 for spanwise-resolved 3D.
extrude_width = 0.08*h;
nz = 1;

// -----------------------
// Mesh resolution controls
// -----------------------
// These are numbers of points, not elements. Elements = points - 1.
// Medium baseline: about O(1e5) hex elements for nz=1.
// For grid convergence, scale these counts consistently.
nx_up   = 1101;    // x/h in [-20,0]
nx_01   = 41;     // x/h in [0,1]
nx_14   = 121;    // x/h in [1,4]
nx_46   = 81;     // x/h in [4,6]
nx_610  = 121;    // x/h in [6,10]
nx_far  = 321;    // x/h in [10,50]

ny_upper = 129;   // y/h in [1,9]
ny_lower = 97;    // y/h in [0,1]

// Bump clusters at both ends of a wall-normal block:
// upper block clusters at y/h=1 shear-layer line and top wall y/h=9;
// lower block clusters at bottom wall y/h=0 and shear-layer line y/h=1.
bump_upper = 0.05;
bump_lower = 0.05;

// Far downstream growth from x/h=10 to outlet, clustered at x/h=10.
prog_far = 1.010;

lc = h/20;
eps = 1.0e-8*h;

// -----------------------
// Points
// -----------------------
// y = h line: upstream lower wall and downstream internal shear-layer line
Point(1)  = {-Lu, h, 0, lc};   // M0
Point(9)  = {0,   h, 0, lc};   // M1
Point(10) = {x1,  h, 0, lc};   // M2
Point(11) = {x4,  h, 0, lc};   // M3
Point(12) = {x6,  h, 0, lc};   // M4
Point(13) = {x10, h, 0, lc};   // M5
Point(14) = {Ld,  h, 0, lc};   // M6

// top wall y = H
Point(2)  = {-Lu, H, 0, lc};   // T0
Point(3)  = {0,   H, 0, lc};   // T1
Point(4)  = {x1,  H, 0, lc};   // T2
Point(5)  = {x4,  H, 0, lc};   // T3
Point(6)  = {x6,  H, 0, lc};   // T4
Point(7)  = {x10, H, 0, lc};   // T5
Point(8)  = {Ld,  H, 0, lc};   // T6

// bottom wall y = 0, downstream only
Point(15) = {0,   0, 0, lc};   // B1
Point(16) = {x1,  0, 0, lc};   // B2
Point(17) = {x4,  0, 0, lc};   // B3
Point(18) = {x6,  0, 0, lc};   // B4
Point(19) = {x10, 0, 0, lc};   // B5
Point(20) = {Ld,  0, 0, lc};   // B6

// -----------------------
// Curves
// -----------------------
// inlet and top wall
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};

// y = h: upstream wall plus downstream internal block boundaries
Line(9)  = {1, 9};
Line(10) = {9, 10};
Line(11) = {10, 11};
Line(12) = {11, 12};
Line(13) = {12, 13};
Line(14) = {13, 14};

// bottom wall downstream
Line(15) = {15, 16};
Line(16) = {16, 17};
Line(17) = {17, 18};
Line(18) = {18, 19};
Line(19) = {19, 20};

// verticals for upper block; x=0,1,4,6,10,50
Line(20) = {9, 3};
Line(21) = {10, 4};
Line(22) = {11, 5};
Line(23) = {12, 6};
Line(24) = {13, 7};
Line(25) = {14, 8};

// verticals for lower block; x=0 is the step wall, x=50 lower outlet
Line(26) = {15, 9};
Line(27) = {16, 10};
Line(28) = {17, 11};
Line(29) = {18, 12};
Line(30) = {19, 13};
Line(31) = {20, 14};

// -----------------------
// Structured block surfaces
// -----------------------
Curve Loop(1) = {1, 2, -20, -9};       Plane Surface(1) = {1};   // upstream channel [-20,0] x [1,9]

Curve Loop(2) = {20, 3, -21, -10};     Plane Surface(2) = {2};   // upper [0,1]
Curve Loop(3) = {21, 4, -22, -11};     Plane Surface(3) = {3};   // upper [1,4]
Curve Loop(4) = {22, 5, -23, -12};     Plane Surface(4) = {4};   // upper [4,6]
Curve Loop(5) = {23, 6, -24, -13};     Plane Surface(5) = {5};   // upper [6,10]
Curve Loop(6) = {24, 7, -25, -14};     Plane Surface(6) = {6};   // upper [10,50]

Curve Loop(7)  = {26, 10, -27, -15};   Plane Surface(7)  = {7};  // lower [0,1]
Curve Loop(8)  = {27, 11, -28, -16};   Plane Surface(8)  = {8};  // lower [1,4]
Curve Loop(9)  = {28, 12, -29, -17};   Plane Surface(9)  = {9};  // lower [4,6]
Curve Loop(10) = {29, 13, -30, -18};   Plane Surface(10) = {10}; // lower [6,10]
Curve Loop(11) = {30, 14, -31, -19};   Plane Surface(11) = {11}; // lower [10,50]

// -----------------------
// Transfinite spacing
// -----------------------
// Streamwise counts. Profile stations are exact block boundaries.
Transfinite Curve {2, 9}       = nx_up  Using Progression 1.0;
Transfinite Curve {3, 10, 15}  = nx_01  Using Progression 1.0;
Transfinite Curve {4, 11, 16}  = nx_14  Using Progression 1.0;
Transfinite Curve {5, 12, 17}  = nx_46  Using Progression 1.0;
Transfinite Curve {6, 13, 18}  = nx_610 Using Progression 1.0;
Transfinite Curve {7, 14, 19}  = nx_far Using Progression prog_far;

// Wall-normal counts.
Transfinite Curve {1,20,21,22,23,24,25} = ny_upper Using Bump bump_upper;
Transfinite Curve {26,27,28,29,30,31}   = ny_lower Using Bump bump_lower;

// Surface corner ordering.
Transfinite Surface {1} = {1, 2, 3, 9};
Transfinite Surface {2} = {9, 3, 4, 10};
Transfinite Surface {3} = {10, 4, 5, 11};
Transfinite Surface {4} = {11, 5, 6, 12};
Transfinite Surface {5} = {12, 6, 7, 13};
Transfinite Surface {6} = {13, 7, 8, 14};

Transfinite Surface {7}  = {15, 9, 10, 16};
Transfinite Surface {8}  = {16, 10, 11, 17};
Transfinite Surface {9}  = {17, 11, 12, 18};
Transfinite Surface {10} = {18, 12, 13, 19};
Transfinite Surface {11} = {19, 13, 14, 20};

Recombine Surface {1,2,3,4,5,6,7,8,9,10,11};
Mesh.RecombineAll = 1;

// -----------------------
// Extrude to quasi-2D 3D mesh
// -----------------------
Extrude {0, 0, extrude_width} {
  Surface{1,2,3,4,5,6,7,8,9,10,11};
  Layers{nz};
  Recombine;
}

Coherence;

// -----------------------
// 3D physical groups, selected geometrically to avoid hard-coded extrusion IDs
// -----------------------
inlet3d[]  = Surface In BoundingBox{-Lu-eps, h-eps, -eps, -Lu+eps, H+eps, extrude_width+eps};
outlet3d[] = Surface In BoundingBox{ Ld-eps, -eps, -eps,  Ld+eps, H+eps, extrude_width+eps};

topwall3d[] = Surface In BoundingBox{-Lu-eps, H-eps, -eps, Ld+eps, H+eps, extrude_width+eps};
botwall3d[] = Surface In BoundingBox{-eps, -eps, -eps, Ld+eps, eps, extrude_width+eps};
stepwall3d[] = Surface In BoundingBox{-eps, -eps, -eps, eps, h+eps, extrude_width+eps};
upstreamlower3d[] = Surface In BoundingBox{-Lu-eps, h-eps, -eps, eps, h+eps, extrude_width+eps};

front3d[] = Surface In BoundingBox{-Lu-eps, -eps, -eps, Ld+eps, H+eps, eps};
back3d[]  = Surface In BoundingBox{-Lu-eps, -eps, extrude_width-eps, Ld+eps, H+eps, extrude_width+eps};
fluid3d[] = Volume In BoundingBox{-Lu-eps, -eps, -eps, Ld+eps, H+eps, extrude_width+eps};

Physical Surface("Inlet")  = {inlet3d[]};
Physical Surface("Outlet") = {outlet3d[]};

Physical Surface("Top_wall") = {topwall3d[]};
Physical Surface("Bottom_wall") = {botwall3d[]};
Physical Surface("Step_wall") = {stepwall3d[]};
Physical Surface("Upstream_lower_wall") = {upstreamlower3d[]};

// Convenient groups for weak wall imposition or wall-function tests.
Physical Surface("Viscous_walls") = {topwall3d[], botwall3d[], stepwall3d[], upstreamlower3d[]};
Physical Surface("Nitsche_walls") = {topwall3d[], botwall3d[], stepwall3d[], upstreamlower3d[]};
Physical Surface("Wall_function_walls") = {topwall3d[], botwall3d[], stepwall3d[], upstreamlower3d[]};

// Quasi-2D planes: set as symmetry/slip/periodic depending on the solver.
Physical Surface("Front") = {front3d[]};
Physical Surface("Back")  = {back3d[]};
Physical Volume("Fluid") = {fluid3d[]};

// Use MSH2 for compatibility with many FEM solvers. Change to 4.1 if preferred.
Mesh.MshFileVersion = 2.2;
Mesh 3;
Save "nasa_tmr_2dbfs_quasi2d.msh";
