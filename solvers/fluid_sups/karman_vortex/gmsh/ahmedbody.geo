// ============================================================
// ahmed_walllaw_nitsche_pilot.geo
//
// Target:
//   Ahmed body CFD
//   Stabilized FEM / VMS
//   Wall-model / wall-law + weak no-penetration by Nitsche
//
// Mesh concept:
//   - TET4 volume mesh
//   - TRI3 boundary mesh
//   - Ahmed body surface tangential size: 6--10 mm
//   - Wake refinement: 40--80 mm
//   - Near-body refinement: 30--50 mm
//
// Important:
//   This mesh file does NOT impose wall-law.
//   Wall-law and Nitsche are solver-side boundary treatments.
//   This file only provides boundary groups and mesh resolution.
// ============================================================

SetFactory("Built-in");

// ------------------------------------------------------------
// Mesh format
// ------------------------------------------------------------
Mesh.MshFileVersion = 2.2;

// ------------------------------------------------------------
// Stability settings
// ------------------------------------------------------------
// HXT(10) can be fast but less robust for dirty STL.
// Delaunay(1) is slower but usually more stable.
Mesh.Algorithm3D = 1;

// Small random perturbation to avoid geometric degeneracy
Mesh.RandomFactor   = 1.0e-6;
Mesh.RandomFactor3D = 1.0e-6;

// More aggressive STL cleanup
Mesh.StlRemoveBadTriangles = 2;

// Mesh optimization
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 0;

// Parallel meshing
General.NumThreads = 24;

// Avoid too much automatic extension from boundary sizes.
// Background field controls the main mesh size.
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 1;
Mesh.MeshSizeFromCurvature = 12;

// ------------------------------------------------------------
// Geometry offset
// ------------------------------------------------------------
// Small offset to avoid exact contact between body and ground.
// z0 = -gap makes ground slightly below original STL z=0.
gap = 1.0e-3;

// ============================================================
// Mesh size parameters
// ============================================================
//
// Units are meters.
//
// For wall-law + Nitsche:
//   Tangential surface size should be around 5--10 mm.
//   First wall sampling distance y_wall is handled in solver
//   from wall-face/cell geometry. If you need y+ = 30--150,
//   the near-wall element height must be checked after meshing.
//
// This is a pilot mesh, not a final production prism-layer mesh.
// ============================================================

lcFar      = 0.60;    // far field: 600 mm
lcNearBox  = 0.040;   // near body: 40 mm
lcWake     = 0.060;   // wake: 60 mm
lcSurf     = 0.008;   // Ahmed body surface: 8 mm
lcEdge     = 0.004;   // sharp/slant/rear edge target: 4 mm

// Distance transition from body surface
dMinSurf = 0.010;     // within 10 mm keep lcSurf
dMaxSurf = 0.200;     // transition to lcNearBox over 200 mm

// Distance transition from feature curves
dMinEdge = 0.000;
dMaxEdge = 0.040;     // refine within 40 mm from feature curves

// ============================================================
// 1) Import Ahmed body STL
// ============================================================

Merge "ahmed_25deg_m.stl";

// Build discrete topology from STL.
// If your STL has poor feature detection, replace this by
// ClassifySurfaces + CreateGeometry workflow.
CreateTopology;

stlPts[]   = Point{:};
stlCur[]   = Curve{:};
stlSurfs[] = Surface{:};

// Move body upward slightly to avoid ground intersection.
Translate {0, 0, +gap} {
  Point{stlPts[]};
  Curve{stlCur[]};
  Surface{stlSurfs[]};
}

// Remove duplicate mesh entities where possible.
Coherence Mesh;

// ============================================================
// 2) Outer wind-tunnel box
// ============================================================
//
// Coordinate convention assumed:
//   x: streamwise
//   y: lateral
//   z: vertical
//
// Inlet  : x = -5
// Outlet : x = 15
// Ground : z = -gap
// Top    : z = 8
// ============================================================

z0 = -gap;
H  = 8.0 + gap;

Point(1) = {-5.0, -4.0, z0, lcFar};
Point(2) = {15.0, -4.0, z0, lcFar};
Point(3) = {15.0,  4.0, z0, lcFar};
Point(4) = {-5.0,  4.0, z0, lcFar};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(10) = {1, 2, 3, 4};
Plane Surface(100) = {10};

// Extrude ground surface to create wind-tunnel boundary surfaces.
out[] = Extrude {0, 0, H} {
  Surface{100};
};

topS   = out[0];
tmpVol = out[1];

side1  = out[2];  // y = -4
side2  = out[3];  // x = 15, outlet
side3  = out[4];  // y =  4
side4  = out[5];  // x = -5, inlet

// Delete temporary box volume.
// We will rebuild fluid volume with Ahmed body as hole.
Delete {
  Volume{tmpVol};
}

// ============================================================
// 3) Fluid volume
// ============================================================

Surface Loop(200) = {
  100,
  topS,
  side1,
  side2,
  side3,
  side4
};

Surface Loop(201) = {
  stlSurfs[]
};

Volume(300) = {
  200,
  201
};

// ============================================================
// 4) Mesh size fields
// ============================================================

// ------------------------------------------------------------
// Field 1: distance from Ahmed body surface
// ------------------------------------------------------------
Field[1] = Distance;
Field[1].SurfacesList = {stlSurfs[]};
Field[1].Sampling = 40;

// ------------------------------------------------------------
// Field 2: near-wall surface refinement
// ------------------------------------------------------------
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lcSurf;
Field[2].SizeMax = lcNearBox;
Field[2].DistMin = dMinSurf;
Field[2].DistMax = dMaxSurf;

// ------------------------------------------------------------
// Field 3: feature-edge refinement
//
// This helps around slant edge, rear edge, and other STL
// topology curves, if CreateTopology detects them.
// ------------------------------------------------------------
Field[3] = Distance;
Field[3].CurvesList = {stlCur[]};
Field[3].Sampling = 80;

Field[4] = Threshold;
Field[4].InField = 3;
Field[4].SizeMin = lcEdge;
Field[4].SizeMax = lcSurf;
Field[4].DistMin = dMinEdge;
Field[4].DistMax = dMaxEdge;

// ------------------------------------------------------------
// Field 5: near-body box
// ------------------------------------------------------------
Field[5] = Box;
Field[5].VIn  = lcNearBox;
Field[5].VOut = lcFar;

Field[5].XMin = -1.5;
Field[5].XMax =  2.5;

Field[5].YMin = -1.0;
Field[5].YMax =  1.0;

Field[5].ZMin = z0;
Field[5].ZMax = 1.5;

Field[5].Thickness = 0.5;

// ------------------------------------------------------------
// Field 6: wake refinement box
// ------------------------------------------------------------
Field[6] = Box;
Field[6].VIn  = lcWake;
Field[6].VOut = lcFar;

Field[6].XMin = 2.0;
Field[6].XMax = 8.0;

Field[6].YMin = -0.9;
Field[6].YMax =  0.9;

Field[6].ZMin = z0;
Field[6].ZMax = 1.4;

Field[6].Thickness = 1.0;

// ------------------------------------------------------------
// Field 7: stronger rear wake core
// ------------------------------------------------------------
Field[7] = Box;
Field[7].VIn  = 0.040;
Field[7].VOut = lcWake;

Field[7].XMin = 0.8;
Field[7].XMax = 4.0;

Field[7].YMin = -0.55;
Field[7].YMax =  0.55;

Field[7].ZMin = z0;
Field[7].ZMax = 0.9;

Field[7].Thickness = 0.5;

// ------------------------------------------------------------
// Field 8: ground near-body refinement
//
// Useful when using wall-law also on the ground.
// ------------------------------------------------------------
Field[8] = Box;
Field[8].VIn  = 0.030;
Field[8].VOut = lcFar;

Field[8].XMin = -1.5;
Field[8].XMax =  5.0;

Field[8].YMin = -1.2;
Field[8].YMax =  1.2;

Field[8].ZMin = z0;
Field[8].ZMax = z0 + 0.20;

Field[8].Thickness = 0.3;

// ------------------------------------------------------------
// Final background field
// ------------------------------------------------------------
Field[20] = Min;
Field[20].FieldsList = {
  2,
  4,
  5,
  6,
  7,
  8
};

Background Field = 20;

// ============================================================
// 5) Physical groups
// ============================================================
//
// Keep names simple and solver-friendly.
// Wall-law/Nitsche treatment should be assigned in the solver
// to "ahmed_body" and optionally "lowerWall".
// ============================================================

Physical Volume("Fluid") = {300};

Physical Surface("lowerWall") = {100};
Physical Surface("upperWall") = {topS};

Physical Surface("Inlet")  = {side4};
Physical Surface("outlet") = {side2};

Physical Surface("frontAndBack") = {
  side1,
  side3
};

Physical Surface("ahmed_body") = {
  stlSurfs[]
};

// ============================================================
// 6) Notes for solver-side wall-law + Nitsche
// ============================================================
//
// For ahmed_body:
//
//   u_rel = u - U_wall
//   u_n   = u_rel . n
//   u_t   = u_rel - u_n n
//
//   normal:
//      impose u_n = 0 weakly by symmetric Nitsche
//
//   tangential:
//      do NOT impose no-slip penalty
//      use wall-law traction:
//
//        tau_t = beta_wall u_t
//
//      with Spalding/Reichardt/log-law wall model.
//
// Recommended initial solver parameters:
//
//   gamma_n = 20 -- 50
//   start with gamma_n = 30
//
// Monitor:
//
//   max |u . n| / U_inf
//   integrated leakage / inlet flux
//   y+ distribution
//   u_tau distribution
//   beta_wall distribution
//   C_D, C_L
//
// For wall-law validity:
//
//   target: 30 < y+ < 150 on most wall faces
//
// If y+ is too large, this isotropic tet mesh is not enough.
// Use prism layers or reduce near-wall size.
// ============================================================