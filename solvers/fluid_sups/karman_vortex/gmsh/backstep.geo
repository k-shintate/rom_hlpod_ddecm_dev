// Backward-facing step benchmark mesh for stabilized FEM / VMS / Nitsche wall-function tests
// Units: nondimensional meters. Step height h = 1 by default.
// Geometry: inlet channel y=[h,H], downstream channel y=[0,H], step at x=0.
// Mesh: structured quadrilateral in 2D, extruded by one layer to obtain quasi-2D hexa/prism-style 3D mesh.

SetFactory("Built-in");

// -----------------------
// Benchmark-like geometry
// -----------------------
h  = 1.0;              // step height
ER = 2.0;              // expansion ratio H/h. Use 1.9423 if matching Biswas et al.-type BFS data.
H  = ER*h;             // downstream channel height
Lu = 5.0*h;            // upstream length before the step
Lr = 8.0*h;            // refined downstream block, includes expected reattachment region for many laminar/turbulent cases
Ld = 30.0*h;           // total downstream length

// Quasi-2D extrusion width. Set to e.g. 4*h and nz>1 for a spanwise-resolved 3D BFS.
extrude_width = 0.08*h;
nz = 1;

// -----------------------
// Mesh resolution controls
// -----------------------
// Recommended starting point for verification studies: run at least Coarse/Medium/Fine by scaling these counts.
nx_up   = 61;          // x points in upstream channel [-Lu,0]
nx_near = 161;         // x points in near downstream block [0,Lr]
nx_far  = 201;         // x points in far downstream block [Lr,Ld]
ny_upper = 65;         // y points in upper channel block [h,H]
ny_lower = 65;         // y points in lower/step-height block [0,h]

// Grading. Bump clusters nodes near both sides of a wall-normal block.
// Increase bump_wall for stronger clustering near walls and the separated shear layer y=h.
bump_wall = 0.08;
prog_downstream_near = 1.015;  // cluster near step, grow to x=Lr
prog_downstream_far  = 1.010;  // gentle growth toward outlet

lc = h/20;
eps = 1.0e-8*h;

// -----------------------
// Points
// -----------------------
// Upstream channel corners and partition points
Point(1)  = {-Lu, h, 0, lc};
Point(2)  = {-Lu, H, 0, lc};
Point(3)  = {0,   H, 0, lc};
Point(4)  = {Lr,  H, 0, lc};
Point(5)  = {Ld,  H, 0, lc};
Point(6)  = {Ld,  h, 0, lc};
Point(7)  = {Lr,  h, 0, lc};
Point(8)  = {0,   h, 0, lc};
Point(9)  = {0,   0, 0, lc};
Point(10) = {Lr,  0, 0, lc};
Point(11) = {Ld,  0, 0, lc};

// -----------------------
// Lines
// -----------------------
Line(1)  = {1, 2};    // inlet
Line(2)  = {2, 3};    // top wall, upstream
Line(3)  = {3, 4};    // top wall, near downstream
Line(4)  = {4, 5};    // top wall, far downstream
Line(5)  = {5, 6};    // outlet, upper segment
Line(6)  = {6, 11};   // outlet, lower segment
Line(7)  = {9, 10};   // bottom wall, near downstream
Line(8)  = {10, 11};  // bottom wall, far downstream
Line(9)  = {9, 8};    // vertical step wall
Line(10) = {8, 1};    // lower wall upstream of the step
Line(11) = {8, 7};    // internal line y=h, near downstream
Line(12) = {7, 6};    // internal line y=h, far downstream
Line(13) = {8, 3};    // internal line x=0, upper block
Line(14) = {7, 4};    // internal line x=Lr, upper block
Line(15) = {10, 7};   // internal line x=Lr, lower block

// -----------------------
// Structured block surfaces
// -----------------------
Curve Loop(1) = {1, 2, -13, 10};         // upstream channel
Plane Surface(1) = {1};

Curve Loop(2) = {13, 3, -14, -11};       // upper near downstream
Plane Surface(2) = {2};

Curve Loop(3) = {14, 4, 5, -12};         // upper far downstream
Plane Surface(3) = {3};

Curve Loop(4) = {9, 11, -15, -7};        // lower near downstream / recirculation block
Plane Surface(4) = {4};

Curve Loop(5) = {15, 12, 6, -8};         // lower far downstream
Plane Surface(5) = {5};

// -----------------------
// Transfinite spacing
// -----------------------
// Streamwise counts
Transfinite Curve {2, 10}       = nx_up   Using Progression 1.0;
Transfinite Curve {3, 11, 7}    = nx_near Using Progression prog_downstream_near;
Transfinite Curve {4, 12, 8}    = nx_far  Using Progression prog_downstream_far;

// Wall-normal counts: Bump clusters near both ends of each block.
// This gives near-wall points at solid walls and also resolves the shear-layer line y=h.
Transfinite Curve {1, 13, 14, 5} = ny_upper Using Bump bump_wall;
Transfinite Curve {9, 15, 6}     = ny_lower Using Bump bump_wall;

Transfinite Surface {1} = {1, 2, 3, 8};
Transfinite Surface {2} = {8, 3, 4, 7};
Transfinite Surface {3} = {7, 4, 5, 6};
Transfinite Surface {4} = {9, 8, 7, 10};
Transfinite Surface {5} = {10, 7, 6, 11};

Recombine Surface {1,2,3,4,5};
Mesh.RecombineAll = 1;

// -----------------------
// 2D physical groups, useful if you generate Mesh 2 instead of Mesh 3
// -----------------------
Physical Curve("Inlet_2D")  = {1};
Physical Curve("Outlet_2D") = {5,6};
Physical Curve("Top_wall_2D") = {2,3,4};
Physical Curve("Bottom_wall_2D") = {7,8};
Physical Curve("Step_wall_2D") = {9};
Physical Curve("Upstream_lower_wall_2D") = {10};
Physical Curve("Solid_walls_2D") = {2,3,4,7,8,9,10};
Physical Surface("Fluid_2D") = {1,2,3,4,5};

// -----------------------
// Extrude to quasi-2D 3D mesh
// -----------------------
Extrude {0, 0, extrude_width} {
  Surface{1,2,3,4,5};
  Layers{nz};
  Recombine;
}

Coherence;

// -----------------------
// 3D physical groups selected geometrically, avoiding fragile hard-coded extrusion IDs
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

// Convenient groups for weak imposition / wall functions in the FEM code.
Physical Surface("Nitsche_walls") = {topwall3d[], botwall3d[], stepwall3d[], upstreamlower3d[]};
Physical Surface("Wall_function_walls") = {topwall3d[], botwall3d[], stepwall3d[], upstreamlower3d[]};

// For quasi-2D simulations: use Front/Back as periodic, symmetry, or slip planes depending on the solver.
Physical Surface("Front") = {front3d[]};
Physical Surface("Back")  = {back3d[]};
Physical Volume("Fluid") = {fluid3d[]};

Mesh.MshFileVersion = 2.2;
Mesh 3;
Save "backward_facing_step_benchmark.msh";