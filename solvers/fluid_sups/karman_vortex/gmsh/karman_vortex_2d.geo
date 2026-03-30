// This grid is drawn in units of meters

d = 1; // Diameter of the particle
offset = 1*d;  // Distance away from particle where structured circular grid extends to
upstream_offset = 7*d ;
downstream_offset = 25*d;
vertical_offset = 7*d;
wake_factor = 1.1;  // How much higher to place the downstream points to create a growing wake region in the mesh

theta = 45 * Pi / 180.0 ;
alpha = Cos(theta);

// Center of particle point
Point(1) = {0, 0, 0};

// Define points on the surface of the particle
Point(2) = {d/2*alpha, d/2*alpha, 0};
Point(3) = {-d/2*alpha, d/2*alpha, 0};
Point(4) = {-d/2*alpha, -d/2*alpha, 0};
Point(5) = {d/2*alpha, -d/2*alpha, 0};

// Define points for circular grid around particle
Point(6) = {offset*alpha, offset*alpha, 0};
Point(7) = {-offset*alpha, offset*alpha, 0};
Point(8) = {-offset*alpha, -offset*alpha, 0};
Point(9) = {offset*alpha, -offset*alpha, 0};

// Define corner points for outer boundary
Point(10) = {-upstream_offset*alpha, -vertical_offset*alpha, 0};
Point(11) = {-upstream_offset*alpha, vertical_offset*alpha, 0};
Point(12) = {downstream_offset*alpha, vertical_offset*alpha, 0};
Point(13) = {downstream_offset*alpha, -vertical_offset*alpha, 0};

// Define points on the horizontal outer boundaries
Point(14) = {offset*alpha, vertical_offset*alpha, 0};
Point(15) = {-offset*alpha, -vertical_offset*alpha, 0};
Point(16) = {-offset*alpha, vertical_offset*alpha, 0};
Point(17) = {offset*alpha, -vertical_offset*alpha, 0};

// Define points on the vertical outer boundaries
Point(18) = {-upstream_offset*alpha, -offset*alpha, 0};
Point(19) = {-upstream_offset*alpha, offset*alpha, 0};
Point(20) = {downstream_offset*alpha, wake_factor*offset*alpha, 0}; // Make this point a little higher to capture wake region
Point(21) = {downstream_offset*alpha, -wake_factor*offset*alpha, 0}; // Make point a little lower to capture wake region

// Define circle of the particle
Circle(1) = {4, 1, 3};
Circle(2) = {3, 1, 2};
Circle(3) = {2, 1, 5};
Circle(4) = {5, 1, 4};

// Define circle around circular offset around particle
Circle(5) = {8, 1, 7};
Circle(6) = {7, 1, 6};
Circle(7) = {6, 1, 9};
Circle(8) = {9, 1, 8};

// Define outer boundary
Line(9) = {10, 18};
Line(10) = {18,19};
Line(11) = {19, 11};
Line(12) = {11,16};
Line(13) = {16,14};
Line(14) = {14,12};
Line(15) = {12,20};
Line(16) = {20,21};
Line(17) = {21,13};
Line(18) = {13,17};
Line(19) = {17,15};
Line(20) = {15,10};

// Define Lines around particle
Line(21) = {3,7};
Line(22) = {2,6};
Line(23) = {4,8};
Line(24) = {5,9};

// Define lines moving from particle to boundary
Line(25) = {8,18};
Line(26) = {7,19};
Line(27) = {7,16};
Line(28) = {6,14};
Line(29) = {8,15};
Line(30) = {9,17};
Line(31) = {6,20};
Line(32) = {9,21};

// Define Surfaces
Curve Loop(1) = {23, 5, -21, -1};
Plane Surface(1) = {1};

Curve Loop(2) = {21, 6, -22, -2};
Plane Surface(2) = {2};

Curve Loop(3) = {22, 7, -24, -3};
Plane Surface(3) = {3};

Curve Loop(4) = {4, 23, -8, -24};
Plane Surface(4) = {4};

Curve Loop(5) = {25, 10, -26, -5};
Plane Surface(5) = {5};

Curve Loop(6) = {26, 11, 12, -27};
Plane Surface(6) = {6};

Curve Loop(7) = {27, 13, -28, -6};
Plane Surface(7) = {7};

Curve Loop(8) = {28, 14, 15, -31};
Plane Surface(8) = {8};

Curve Loop(9) = {7, 32, -16, -31};
Plane Surface(9) = {9};

Curve Loop(10) = {32, 17, 18, -30};
Plane Surface(10) = {10};

Curve Loop(11) = {30, 19, -29, -8};
Plane Surface(11) = {11};

Curve Loop(12) = {29, 20, 9, -25};
Plane Surface(12) = {12};

// Define Transfinite edge spacings
Transfinite Curve {11, 27, 28, 15, 9, 29, 30, 17} = 91 Using Progression 1;
Transfinite Curve {20, 25, 26, 12} = 91 Using Progression 1;
Transfinite Curve {14, 31, 32, -18} = 151 Using Progression 1;
Transfinite Curve {10, 5, 1} = 31 Using Progression 1;
Transfinite Curve {3, 7, 16} = 31 Using Progression 1;
Transfinite Curve {19, 8, 4, 2, 6, 13} = 61 Using Progression 1;
Transfinite Curve {-21, -22, -24, -23} = 51 Using Progression 0.97;

Transfinite Surface "*";
Recombine Surface "*";

extrude_width = -0.08*d;
Extrude {0, 0, extrude_width} {
   Surface{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
   Layers{1};
   Recombine;
}

Physical Surface("Inlet") = {155, 133, 291};
Physical Surface("Cylinder_wall") = {53, 75, 97, 107};
Physical Surface("Left_wall") = {159, 177, 199};
Physical Surface("Right_wall") = {287, 265, 247};
Physical Surface("Top_wall") = {1,2,3,4,5,6,7,8,9,10,11,12};
Physical Surface("Bottom_wall") = {296,142,164,274,120,54,76,98,186,252,230,208};
Physical Surface("Outlet") = {243,225,203};

Physical Volume("Fluid") = {1,2,3,4,5,6,7,8,9,10,11,12};
Physical Surface("All_Boundaries") = {1, 2, 3, 4, 5, 6, 7};

Coherence;
Mesh 3;
Save "karman_vortex_2d.msh";
