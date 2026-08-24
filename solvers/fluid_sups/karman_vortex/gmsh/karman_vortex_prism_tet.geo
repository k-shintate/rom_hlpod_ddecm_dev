// ============================================================
// Karman vortex mixed mesh v9
// Re = 100 quasi-2D
//
// Mesh:
//   - TET4 : outer/core region
//   - PRI6 : annular boundary-layer region
//   - TRI3 : cylinder wall for Nitsche
//
// Main improvements from v8:
//
//   1) The outer TET cross-section is ONE single Plane Surface
//      with a circular hole.
//
//      -> No internal radial/wake partition lines.
//      -> No stripe-like mesh-density discontinuities caused
//         by Surface(5)...Surface(12).
//
//   2) Smooth mesh-size transitions are imposed by fields.
//
//   3) Quasi-2D:
//        span Lz = 0.10 D
//        4 spanwise elements
//
//   4) Cylinder/interface:
//        160 circumferential elements
//
//   5) PRI boundary layer:
//        14 radial layers
//        first wall-normal spacing ~= 0.00370 D
//
// IMPORTANT:
//   Span_top and Span_bottom are NOT no-slip walls.
//   Use periodic BC in the flow solver.
// ============================================================

SetFactory("Built-in");

// ============================================================
// Parameters
// ============================================================

d = 1.0;

// ------------------------------------------------------------
// Cylinder / prism region
// ------------------------------------------------------------

r_cyl_nominal = 0.50*d;

r_prism_outer = 0.65*d;

bl_thickness =
    r_prism_outer - r_cyl_nominal;


// ------------------------------------------------------------
// Computational domain
//
// Re=100 benchmark-like domain:
//   inlet  : -10 D
//   outlet : +25 D
//   upper  : +10 D
//   lower  : -10 D
// ------------------------------------------------------------

xmin = -10.0*d;
xmax =  25.0*d;

ymin = -10.0*d;
ymax =  10.0*d;


// ------------------------------------------------------------
// Quasi-2D span
// ------------------------------------------------------------
//
// 0.10 D / 4 = 0.025 D
//

extrude_width = -0.10*d;

zmax = 0.0;
zmin = extrude_width;


// ------------------------------------------------------------
// Interface resolution
// ------------------------------------------------------------
//
// 40 segments / quadrant
// -> 160 circumferential segments
//
// Cylinder wall:
//   pi D / 160 ~= 0.01963 D
//
// PRI/TET interface:
//   2 pi 0.65D / 160 ~= 0.02553 D
//

n_theta_quarter = 41;


// ------------------------------------------------------------
// Spanwise resolution
// ------------------------------------------------------------
//
// 5 nodes -> 4 elements
//
// dz = 0.10D / 4
//    = 0.025D
//

n_z_interface = 5;


// ------------------------------------------------------------
// PRI boundary layer
// ------------------------------------------------------------

n_bl = 14;


// ------------------------------------------------------------
// General TET sizes
// ------------------------------------------------------------

mesh_size_min = 0.020*d;

mesh_size_max = 0.25*d;


// ------------------------------------------------------------
// Geometry tolerance
// ------------------------------------------------------------

eps = 1.e-6*d;


// ------------------------------------------------------------
// 45 degree position used for circular quadrant endpoints
// ------------------------------------------------------------

theta = 45.0*Pi/180.0;

alpha = Cos(theta);


// ============================================================
// Mesh options
// ============================================================

Mesh.MeshSizeMin = mesh_size_min;

Mesh.MeshSizeMax = mesh_size_max;


// Frontal-Delaunay for TRI3 surface mesh
Mesh.Algorithm = 6;


// Delaunay TET4
Mesh.Algorithm3D = 1;


// Never globally recombine.
// Only the explicit PRI extrusion is recombined.
Mesh.RecombineAll = 0;


// Explicit background fields control TET sizes.
Mesh.MeshSizeFromPoints = 0;

Mesh.MeshSizeFromCurvature = 0;

Mesh.MeshSizeExtendFromBoundary = 0;


// ============================================================
// STEP 1
// Circular PRI/TET interface at z = 0
// ============================================================

// Center
Point(1) = {
    0,
    0,
    0
};


// ------------------------------------------------------------
// Four 45-degree interface points
//
// The arcs below run CLOCKWISE.
// This orientation is useful because the circular boundary
// represents a hole in the outer TET cross-section.
// ------------------------------------------------------------

// NE
Point(6) = {

     r_prism_outer*alpha,

     r_prism_outer*alpha,

     0
};


// NW
Point(7) = {

    -r_prism_outer*alpha,

     r_prism_outer*alpha,

     0
};


// SW
Point(8) = {

    -r_prism_outer*alpha,

    -r_prism_outer*alpha,

     0
};


// SE
Point(9) = {

     r_prism_outer*alpha,

    -r_prism_outer*alpha,

     0
};


// ------------------------------------------------------------
// Exact circular interface
// ------------------------------------------------------------

// SW -> NW
Circle(5) = {
    8,
    1,
    7
};


// NW -> NE
Circle(6) = {
    7,
    1,
    6
};


// NE -> SE
Circle(7) = {
    6,
    1,
    9
};


// SE -> SW
Circle(8) = {
    9,
    1,
    8
};


// Same angular resolution on all four quadrants.
Transfinite Curve {
    5,
    6,
    7,
    8
}
=
n_theta_quarter
Using Progression 1;


// ============================================================
// STEP 2
// Single outer rectangular boundary
// ============================================================
//
// IMPORTANT:
//
// There are NO lines connecting the cylinder/interface
// to the outer boundary.
//
// Hence there are no artificial horizontal/radial surface
// partitions that can create stripe-like TET mesh patterns.
// ============================================================


// lower-left
Point(10) = {
    xmin,
    ymin,
    0
};


// lower-right
Point(11) = {
    xmax,
    ymin,
    0
};


// upper-right
Point(12) = {
    xmax,
    ymax,
    0
};


// upper-left
Point(13) = {
    xmin,
    ymax,
    0
};


// ------------------------------------------------------------
// Outer boundary, CCW viewed from +z
// ------------------------------------------------------------

// lower boundary
Line(10) = {
    10,
    11
};


// outlet
Line(11) = {
    11,
    12
};


// upper boundary
Line(12) = {
    12,
    13
};


// inlet
Line(13) = {
    13,
    10
};


// ============================================================
// STEP 3
// ONE single 2D TET-core cross-section
// ============================================================


// Outer rectangle loop
Curve Loop(100) = {

    10,

    11,

    12,

    13
};


// Inner circular hole
//
// Curves 5-8 are clockwise.
//
Curve Loop(101) = {

    5,

    6,

    7,

    8
};


// ------------------------------------------------------------
// Single unstructured surface with one circular hole
// ------------------------------------------------------------

Plane Surface(100) = {

    100,

    101

};


// ============================================================
// STEP 4
// Geometry-only extrusion for TET4 core
// ============================================================
//
// Do NOT specify Layers.
// Do NOT Recombine.
//
// The resulting volume is later filled with unstructured TET4.
// ============================================================

oldReturn =
    Geometry.ExtrudeReturnLateralEntities;

Geometry.ExtrudeReturnLateralEntities = 0;


core[] = Extrude {

    0,

    0,

    extrude_width

} {

    Surface{100};

};


// With lateral-return disabled:
//
// core[0] = translated bottom surface
// core[1] = core volume

If (#core[] != 2)

    Error(
        "Expected one source surface -> [bottom surface, volume]."
    );

EndIf


sCoreBottom = core[0];

vCore       = core[1];


// ============================================================
// STEP 5
// Identical span-top/span-bottom TET surface meshes
// ============================================================
//
// Bottom is copied from top.
//
// This is useful for solver-side periodic BC.
//
// NOTE:
// Gmsh periodic mesh correspondence does NOT by itself impose
// the CFD periodic boundary condition.
// ============================================================

Periodic Surface {

    sCoreBottom

}
=
{

    100

}
Translate {

    0,

    0,

    extrude_width

};


// ============================================================
// STEP 6
// Identify exact cylindrical PRI/TET interface
// ============================================================

interface[] =
Surface In BoundingBox {

    -r_prism_outer-eps,

    -r_prism_outer-eps,

     zmin-eps,


     r_prism_outer+eps,

     r_prism_outer+eps,

     zmax+eps

};


If (#interface[] != 4)

    Error(
        "Expected exactly 4 cylindrical PRI/TET interface surfaces."
    );

EndIf


// ============================================================
// STEP 7
// Regularize cylindrical TRI3 source mesh
// ============================================================
//
// Circumference:
//   160 segments
//
// Span:
//   4 segments
//
// At PRI/TET interface:
//
//   ds_theta ~= 0.02553D
//   dz       = 0.02500D
//
// Hence the triangular interface mesh is nearly isotropic.
// ============================================================


// ------------------------------------------------------------
// Vertical seam curves
// ------------------------------------------------------------

// NE
vNE[] =
Curve In BoundingBox {

     r_prism_outer*alpha-eps,

     r_prism_outer*alpha-eps,

     zmin-eps,


     r_prism_outer*alpha+eps,

     r_prism_outer*alpha+eps,

     zmax+eps

};


// NW
vNW[] =
Curve In BoundingBox {

    -r_prism_outer*alpha-eps,

     r_prism_outer*alpha-eps,

     zmin-eps,


    -r_prism_outer*alpha+eps,

     r_prism_outer*alpha+eps,

     zmax+eps

};


// SW
vSW[] =
Curve In BoundingBox {

    -r_prism_outer*alpha-eps,

    -r_prism_outer*alpha-eps,

     zmin-eps,


    -r_prism_outer*alpha+eps,

    -r_prism_outer*alpha+eps,

     zmax+eps

};


// SE
vSE[] =
Curve In BoundingBox {

     r_prism_outer*alpha-eps,

    -r_prism_outer*alpha-eps,

     zmin-eps,


     r_prism_outer*alpha+eps,

    -r_prism_outer*alpha+eps,

     zmax+eps

};


Transfinite Curve {

    vNE[],

    vNW[],

    vSW[],

    vSE[]

}
=
n_z_interface
Using Progression 1;


// ------------------------------------------------------------
// Explicit bottom circular resolution
// ------------------------------------------------------------
//
// The top circular curves 5-8 are already constrained.
//
// Apply the same count to the translated bottom arcs.
// ------------------------------------------------------------

circBottom[] =
Curve In BoundingBox {

    -r_prism_outer-eps,

    -r_prism_outer-eps,

     zmin-eps,


     r_prism_outer+eps,

     r_prism_outer+eps,

     zmin+eps

};


Transfinite Curve {

    circBottom[]

}
=
n_theta_quarter
Using Progression 1;


// ------------------------------------------------------------
// TRI3 interface mesh
// ------------------------------------------------------------
//
// Do NOT recombine.
//
// TRI3 x normal extrusion -> PRI6.
//
Transfinite Surface {

    interface[]

}
Alternate;


// ============================================================
// STEP 8
// Smooth TET mesh-size fields
// ============================================================
//
// No geometric partition is introduced.
//
// All refinement is handled through background mesh-size fields.
// ============================================================


// ------------------------------------------------------------
// Near-cylinder/interface refinement
// ------------------------------------------------------------
//
// At r=0.65D:
// h ~= 0.025D
//
// Smoothly grows to h ~= 0.20D over 2D distance.
// ------------------------------------------------------------

Field[1] = Distance;

Field[1].SurfacesList = {

    interface[]

};


Field[2] = Threshold;

Field[2].InField = 1;

Field[2].SizeMin = 0.025*d;

Field[2].SizeMax = 0.20*d;

Field[2].DistMin = 0.0*d;

Field[2].DistMax = 2.0*d;


// ------------------------------------------------------------
// Main vortex-formation / near-wake region
// ------------------------------------------------------------
//
// Keep h ~= 0.030D.
//
// Transition thickness is deliberately large to avoid
// visible rectangular mesh-density bands.
// ------------------------------------------------------------

Field[3] = Box;

Field[3].VIn = 0.030*d;

Field[3].VOut = mesh_size_max;


// Slightly upstream of the cylinder center
Field[3].XMin = -0.5*d;

Field[3].XMax = 10.0*d;


// Main wake width
Field[3].YMin = -1.5*d;

Field[3].YMax =  1.5*d;


// Entire quasi-2D thickness
Field[3].ZMin = zmin-eps;

Field[3].ZMax = zmax+eps;


// Smooth transition around Box boundary
Field[3].Thickness = 1.5*d;


// ------------------------------------------------------------
// Far wake
// ------------------------------------------------------------
//
// Keep h ~= 0.060D downstream.
// ------------------------------------------------------------

Field[4] = Box;

Field[4].VIn = 0.060*d;

Field[4].VOut = mesh_size_max;


Field[4].XMin = 8.0*d;

Field[4].XMax = xmax;


Field[4].YMin = -2.0*d;

Field[4].YMax =  2.0*d;


Field[4].ZMin = zmin-eps;

Field[4].ZMax = zmax+eps;


// Even smoother far-wake transition.
Field[4].Thickness = 2.5*d;


// ------------------------------------------------------------
// Minimum of all requested sizes
// ------------------------------------------------------------

Field[5] = Min;

Field[5].FieldsList = {

    2,

    3,

    4

};


Background Field = 5;


// ============================================================
// STEP 9
// PRI6 radial boundary layer
// ============================================================
//
// Source:
//   r = 0.65D PRI/TET interface
//
// Extrusion:
//   toward the cylinder
//
// Final wall:
//   r ~= 0.50D
//
// Total thickness:
//   0.15D
//
// Number of radial layers:
//   14
//
// Growth:
//   about 1.15 away from wall
//
// Because extrusion proceeds from the OUTER interface toward
// the wall, layer widths DECREASE along extrusion direction.
//
// Final wall-normal spacing:
//   ~= 0.00370D
// ============================================================


bl[] = Extrude {

    Surface{

        interface[]

    };


    Layers {

        {

            1,1,1,1,1,1,1,

            1,1,1,1,1,1,1

        },


        {

            0.1519030346381232*bl_thickness,

            0.2839926299756217*bl_thickness,

            0.3988531476604030*bl_thickness,

            0.4987318586906476*bl_thickness,

            0.5855829117604254*bl_thickness,

            0.6611055666037106*bl_thickness,

            0.7267774403804803*bl_thickness,

            0.7838834175776713*bl_thickness,

            0.8335407890534896*bl_thickness,

            0.8767211120759404*bl_thickness,

            0.9142692190519843*bl_thickness,

            0.9469197468572400*bl_thickness,

            0.9753115101661580*bl_thickness,

            1.0000000000000000*bl_thickness

        }

    };


    Recombine;

};


Geometry.ExtrudeReturnLateralEntities =
    oldReturn;


// ------------------------------------------------------------
// Four cylindrical source surfaces:
//
// each returns
//
//   [wall surface, prism volume]
//
// -> 8 entities
// ------------------------------------------------------------

If (#bl[] != 8)

    Error(
        "Expected 4 interface surfaces -> 8 boundary-layer entities."
    );

EndIf


// ============================================================
// STEP 10
// Physical cylinder wall
// ============================================================
//
// These are TRI3 surfaces at the final boundary-layer location.
// ============================================================

Physical Surface(
    "Cylinder_wall"
)
=
{

    bl[0],

    bl[2],

    bl[4],

    bl[6]

};


// Optional debugging group:
//
// Physical Surface(
//     "Prism_Tet_interface"
// ) = {
//
//     interface[]
//
// };


// ============================================================
// STEP 11
// Physical fluid volumes
// ============================================================

Physical Volume(
    "Fluid"
)
=
{

    vCore,


    bl[1],

    bl[3],

    bl[5],

    bl[7]

};


// ============================================================
// STEP 12
// External boundary identification
// ============================================================


// ------------------------------------------------------------
// Inlet
// x = xmin
// ------------------------------------------------------------

inlet[] =
Surface In BoundingBox {

    xmin-eps,

    ymin-eps,

    zmin-eps,


    xmin+eps,

    ymax+eps,

    zmax+eps

};


// ------------------------------------------------------------
// Outlet
// x = xmax
// ------------------------------------------------------------

outlet[] =
Surface In BoundingBox {

    xmax-eps,

    ymin-eps,

    zmin-eps,


    xmax+eps,

    ymax+eps,

    zmax+eps

};


// ------------------------------------------------------------
// Upper far-field
// y = ymax
// ------------------------------------------------------------

upper_far[] =
Surface In BoundingBox {

    xmin-eps,

    ymax-eps,

    zmin-eps,


    xmax+eps,

    ymax+eps,

    zmax+eps

};


// ------------------------------------------------------------
// Lower far-field
// y = ymin
// ------------------------------------------------------------

lower_far[] =
Surface In BoundingBox {

    xmin-eps,

    ymin-eps,

    zmin-eps,


    xmax+eps,

    ymin+eps,

    zmax+eps

};


// ------------------------------------------------------------
// Span top
// z = 0
//
// Includes:
//   - TET core top
//   - PRI top cap surfaces
// ------------------------------------------------------------

span_top[] =
Surface In BoundingBox {

    xmin-eps,

    ymin-eps,

    zmax-eps,


    xmax+eps,

    ymax+eps,

    zmax+eps

};


// ------------------------------------------------------------
// Span bottom
// z = zmin
//
// Includes:
//   - TET core bottom
//   - PRI bottom cap surfaces
// ------------------------------------------------------------

span_bottom[] =
Surface In BoundingBox {

    xmin-eps,

    ymin-eps,

    zmin-eps,


    xmax+eps,

    ymax+eps,

    zmin+eps

};


// ============================================================
// STEP 13
// Physical external boundaries
// ============================================================

Physical Surface(
    "Inlet"
)
=
{

    inlet[]

};


Physical Surface(
    "Outlet"
)
=
{

    outlet[]

};


// Keep the previous names for compatibility with
// existing solver/pre-processing scripts.

Physical Surface(
    "Left_wall"
)
=
{

    upper_far[]

};


Physical Surface(
    "Right_wall"
)
=
{

    lower_far[]

};


// ------------------------------------------------------------
// Span boundaries
//
// IMPORTANT:
//
// These are NOT physical solid walls.
//
// Preferred CFD BC:
//
//   u(x,y,z=0)
//       =
//   u(x,y,z=-0.1D)
//
//   p(x,y,z=0)
//       =
//   p(x,y,z=-0.1D)
//
// i.e. periodic.
//
// If periodic BC is unavailable, use symmetry/slip,
// not no-slip.
// ------------------------------------------------------------

Physical Surface(
    "Span_top"
)
=
{

    span_top[]

};


Physical Surface(
    "Span_bottom"
)
=
{

    span_bottom[]

};


// ============================================================
// STEP 14
// Diagnostic output
// ============================================================

Printf(
    "-------------------------------------------"
);

Printf(
    "Re=100 quasi-2D mixed TET4/PRI6 mesh"
);

Printf(
    "Domain x = [%gD, %gD]",
    xmin/d,
    xmax/d
);

Printf(
    "Domain y = [%gD, %gD]",
    ymin/d,
    ymax/d
);

Printf(
    "Span Lz = %gD",
    -extrude_width/d
);

Printf(
    "Cylinder circumferential segments = %g",
    4*(n_theta_quarter-1)
);

Printf(
    "Spanwise segments = %g",
    n_z_interface-1
);

Printf(
    "Cylinder ds ~= %gD",
    Pi/(4*(n_theta_quarter-1))
);

Printf(
    "Interface ds ~= %gD",
    2*Pi*(r_prism_outer/d) /
    (4*(n_theta_quarter-1))
);

Printf(
    "dz = %gD",
    (-extrude_width/d) /
    (n_z_interface-1)
);

Printf(
    "PRI radial layers = %g",
    n_bl
);

Printf(
    "-------------------------------------------"
);


// ============================================================
// STEP 15
// Final mesh
// ============================================================

Mesh.RecombineAll = 0;


// Improve TET quality.
// Does not alter explicit PRI6 connectivity.
Mesh.Optimize = 1;


Mesh 3;


// ============================================================
// Save
// ============================================================

Save
"karman_vortex_prism_tet_v9_re100_quasi2d_smooth.msh";