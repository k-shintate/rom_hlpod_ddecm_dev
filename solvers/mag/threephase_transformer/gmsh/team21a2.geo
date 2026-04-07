SetFactory("OpenCASCADE");

// =====================================================
// TEAM Problem 21a-2 mesh
//
// P21a-2:
//   - no shield
//   - one non-magnetic steel plate
//   - plate height = 820 mm
//   - plate thickness = 10 mm
//   - coil-to-plate gap = 12 mm
//   - two slits in the plate
//
// Notes:
//   - material constants should be assigned in solver input:
//       NONMAG_PLATE: mu_r = 1, sigma = 1.3889e6 S/m
// =====================================================

Mesh.Algorithm3D = 1;
Mesh.Algorithm   = 6;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;

Scale = 0.001;

// ---------------- Parameters ----------------
lcCoil  = 12 * Scale;
lcPlate = 10 * Scale;
lcAir   = 30 * Scale;
lcFar   = 60 * Scale;
lcSlit  = 3  * Scale;

// Coil pack cross-section (x-y plane)
coilOuterX = 270 * Scale;
coilOuterY = 270 * Scale;
coilInnerX = 200 * Scale;
coilInnerY = 200 * Scale;
rOuter = 45 * Scale;
rInner = 10 * Scale;

// z direction stacking
coilHeight = 217 * Scale;
coilGapZ   = 24  * Scale;
coilSpanZ  = 458 * Scale;
plateHeightZ = 820 * Scale;

// Plate only
plateThk = 10 * Scale;
coilToPlateGap = 12 * Scale;
plateWidthY = 360 * Scale;

// Slits for P21a-2
slitLenZ = 660 * Scale;
slitWidthY = 10 * Scale;
slitDepthX = 10 * Scale;   // equals plate thickness: through slit

// Interpretation of "120 120" in Fig. A1-2(c):
// two slit center positions at y = -60 mm and +60 mm
slitCenter1Y = -60 * Scale;
slitCenter2Y =  60 * Scale;

// centered in z within the plate
slitZ0 = -0.5 * slitLenZ;
slitZ1 =  0.5 * slitLenZ;

// ---------------- Coil cross-section placement ----------------
ox0 = -0.5 * coilOuterX;
ox1 =  0.5 * coilOuterX;
oy0 = -0.5 * coilOuterY;
oy1 =  0.5 * coilOuterY;

ix0 = -0.5 * coilInnerX;
ix1 =  0.5 * coilInnerX;
iy0 = -0.5 * coilInnerY;
iy1 =  0.5 * coilInnerY;

// ---------------- Placement of plate ----------------
xPlate0  = ox1 + coilToPlateGap;
xPlate1  = xPlate0 + plateThk;

yPlate0 = -0.5 * plateWidthY;
yPlate1 =  0.5 * plateWidthY;

zPlate0 = -0.5 * plateHeightZ;
zPlate1 =  0.5 * plateHeightZ;

zCoil1_0 = -(coilGapZ/2 + coilHeight);
zCoil1_1 = - coilGapZ/2;
zCoil2_0 =   coilGapZ/2;
zCoil2_1 =   coilGapZ/2 + coilHeight;

// ---------------- Coil 1 with rounded corners ----------------
z0  = zCoil1_0;

// Outer rounded rectangle
p1 = newp; Point(p1) = {ox0+rOuter, oy0, z0, lcCoil};
p2 = newp; Point(p2) = {ox1-rOuter, oy0, z0, lcCoil};
p3 = newp; Point(p3) = {ox1, oy0+rOuter, z0, lcCoil};
p4 = newp; Point(p4) = {ox1, oy1-rOuter, z0, lcCoil};
p5 = newp; Point(p5) = {ox1-rOuter, oy1, z0, lcCoil};
p6 = newp; Point(p6) = {ox0+rOuter, oy1, z0, lcCoil};
p7 = newp; Point(p7) = {ox0, oy1-rOuter, z0, lcCoil};
p8 = newp; Point(p8) = {ox0, oy0+rOuter, z0, lcCoil};

cTR = newp; Point(cTR) = {ox1-rOuter, oy1-rOuter, z0, lcCoil};
cBR = newp; Point(cBR) = {ox1-rOuter, oy0+rOuter, z0, lcCoil};
cBL = newp; Point(cBL) = {ox0+rOuter, oy0+rOuter, z0, lcCoil};
cTL = newp; Point(cTL) = {ox0+rOuter, oy1-rOuter, z0, lcCoil};

l1 = newl; Line(l1) = {p1,p2};
a1 = newl; Circle(a1) = {p2,cBR,p3};
l2 = newl; Line(l2) = {p3,p4};
a2 = newl; Circle(a2) = {p4,cTR,p5};
l3 = newl; Line(l3) = {p5,p6};
a3 = newl; Circle(a3) = {p6,cTL,p7};
l4 = newl; Line(l4) = {p7,p8};
a4 = newl; Circle(a4) = {p8,cBL,p1};
clOut1 = newll; Curve Loop(clOut1) = {l1,a1,l2,a2,l3,a3,l4,a4};

// Inner rounded rectangle (hole)
q1 = newp; Point(q1) = {ix0+rInner, iy0, z0, lcCoil};
q2 = newp; Point(q2) = {ix1-rInner, iy0, z0, lcCoil};
q3 = newp; Point(q3) = {ix1, iy0+rInner, z0, lcCoil};
q4 = newp; Point(q4) = {ix1, iy1-rInner, z0, lcCoil};
q5 = newp; Point(q5) = {ix1-rInner, iy1, z0, lcCoil};
q6 = newp; Point(q6) = {ix0+rInner, iy1, z0, lcCoil};
q7 = newp; Point(q7) = {ix0, iy1-rInner, z0, lcCoil};
q8 = newp; Point(q8) = {ix0, iy0+rInner, z0, lcCoil};

dTR = newp; Point(dTR) = {ix1-rInner, iy1-rInner, z0, lcCoil};
dBR = newp; Point(dBR) = {ix1-rInner, iy0+rInner, z0, lcCoil};
dBL = newp; Point(dBL) = {ix0+rInner, iy0+rInner, z0, lcCoil};
dTL = newp; Point(dTL) = {ix0+rInner, iy1-rInner, z0, lcCoil};

m1 = newl; Line(m1) = {q1,q2};
b1 = newl; Circle(b1) = {q2,dBR,q3};
m2 = newl; Line(m2) = {q3,q4};
b2 = newl; Circle(b2) = {q4,dTR,q5};
m3 = newl; Line(m3) = {q5,q6};
b3 = newl; Circle(b3) = {q6,dTL,q7};
m4 = newl; Line(m4) = {q7,q8};
b4 = newl; Circle(b4) = {q8,dBL,q1};
clIn1 = newll; Curve Loop(clIn1) = {m1,b1,m2,b2,m3,b3,m4,b4};

s1 = news; Plane Surface(s1) = {clOut1, clIn1};
out1[] = Extrude {0,0,coilHeight} { Surface{s1}; };
coil1[] = {out1[1]};

// ---------------- Coil 2 with rounded corners ----------------
z0  = zCoil2_0;

// Outer rounded rectangle
p11 = newp; Point(p11) = {ox0+rOuter, oy0, z0, lcCoil};
p12 = newp; Point(p12) = {ox1-rOuter, oy0, z0, lcCoil};
p13 = newp; Point(p13) = {ox1, oy0+rOuter, z0, lcCoil};
p14 = newp; Point(p14) = {ox1, oy1-rOuter, z0, lcCoil};
p15 = newp; Point(p15) = {ox1-rOuter, oy1, z0, lcCoil};
p16 = newp; Point(p16) = {ox0+rOuter, oy1, z0, lcCoil};
p17 = newp; Point(p17) = {ox0, oy1-rOuter, z0, lcCoil};
p18 = newp; Point(p18) = {ox0, oy0+rOuter, z0, lcCoil};

cTR2 = newp; Point(cTR2) = {ox1-rOuter, oy1-rOuter, z0, lcCoil};
cBR2 = newp; Point(cBR2) = {ox1-rOuter, oy0+rOuter, z0, lcCoil};
cBL2 = newp; Point(cBL2) = {ox0+rOuter, oy0+rOuter, z0, lcCoil};
cTL2 = newp; Point(cTL2) = {ox0+rOuter, oy1-rOuter, z0, lcCoil};

l11 = newl; Line(l11) = {p11,p12};
a11 = newl; Circle(a11) = {p12,cBR2,p13};
l12 = newl; Line(l12) = {p13,p14};
a12 = newl; Circle(a12) = {p14,cTR2,p15};
l13 = newl; Line(l13) = {p15,p16};
a13 = newl; Circle(a13) = {p16,cTL2,p17};
l14 = newl; Line(l14) = {p17,p18};
a14 = newl; Circle(a14) = {p18,cBL2,p11};
clOut2 = newll; Curve Loop(clOut2) = {l11,a11,l12,a12,l13,a13,l14,a14};

// Inner rounded rectangle (hole)
q11 = newp; Point(q11) = {ix0+rInner, iy0, z0, lcCoil};
q12 = newp; Point(q12) = {ix1-rInner, iy0, z0, lcCoil};
q13 = newp; Point(q13) = {ix1, iy0+rInner, z0, lcCoil};
q14 = newp; Point(q14) = {ix1, iy1-rInner, z0, lcCoil};
q15 = newp; Point(q15) = {ix1-rInner, iy1, z0, lcCoil};
q16 = newp; Point(q16) = {ix0+rInner, iy1, z0, lcCoil};
q17 = newp; Point(q17) = {ix0, iy1-rInner, z0, lcCoil};
q18 = newp; Point(q18) = {ix0, iy0+rInner, z0, lcCoil};

dTR2 = newp; Point(dTR2) = {ix1-rInner, iy1-rInner, z0, lcCoil};
dBR2 = newp; Point(dBR2) = {ix1-rInner, iy0+rInner, z0, lcCoil};
dBL2 = newp; Point(dBL2) = {ix0+rInner, iy0+rInner, z0, lcCoil};
dTL2 = newp; Point(dTL2) = {ix0+rInner, iy1-rInner, z0, lcCoil};

m11 = newl; Line(m11) = {q11,q12};
b11 = newl; Circle(b11) = {q12,dBR2,q13};
m12 = newl; Line(m12) = {q13,q14};
b12 = newl; Circle(b12) = {q14,dTR2,q15};
m13 = newl; Line(m13) = {q15,q16};
b13 = newl; Circle(b13) = {q16,dTL2,q17};
m14 = newl; Line(m14) = {q17,q18};
b14 = newl; Circle(b14) = {q18,dBL2,q11};
clIn2 = newll; Curve Loop(clIn2) = {m11,b11,m12,b12,m13,b13,m14,b14};

s2 = news; Plane Surface(s2) = {clOut2, clIn2};
out2[] = Extrude {0,0,coilHeight} { Surface{s2}; };
coil2[] = {out2[1]};

// ---------------- Plate and slits ----------------
Box(10) = {xPlate0, yPlate0, zPlate0, plateThk, plateWidthY, plateHeightZ};

slit1Y0 = slitCenter1Y - 0.5 * slitWidthY;
slit1Y1 = slitCenter1Y + 0.5 * slitWidthY;
slit2Y0 = slitCenter2Y - 0.5 * slitWidthY;
slit2Y1 = slitCenter2Y + 0.5 * slitWidthY;

// two through-slits
Box(11) = {xPlate0, slit1Y0, slitZ0, slitDepthX, slitWidthY, slitLenZ};
Box(12) = {xPlate0, slit2Y0, slitZ0, slitDepthX, slitWidthY, slitLenZ};

// subtract slits from plate
plateCut[] = BooleanDifference{ Volume{10}; Delete; }{ Volume{11,12}; Delete; };
plateVol[] = {plateCut[]};

// ---------------- Air domain ----------------
marginXneg = 220 * Scale;
marginXpos = 220 * Scale;
marginY    = 220 * Scale;
marginZ    = 220 * Scale;

xAir0 = ox0 - marginXneg;
xAir1 = xPlate1 + marginXpos;
yAir0 = yPlate0 - marginY;
yAir1 = yPlate1 + marginY;
zAir0 = zPlate0 - marginZ;
zAir1 = zPlate1 + marginZ;

Box(100) = {xAir0, yAir0, zAir0, xAir1-xAir0, yAir1-yAir0, zAir1-zAir0};

allSolids[] = {coil1[], coil2[], plateVol[]};
frag[] = BooleanFragments{ Volume{100}; Delete; }{ Volume{allSolids[]}; Delete; };
Coherence;

// ---------------- Identify fragmented volumes ----------------
eps = 1e-6;

volCoil1[] = Volume In BoundingBox{ox0-eps, oy0-eps, zCoil1_0-eps,
                                   ox1+eps, oy1+eps, zCoil1_1+eps};
volCoil2[] = Volume In BoundingBox{ox0-eps, oy0-eps, zCoil2_0-eps,
                                   ox1+eps, oy1+eps, zCoil2_1+eps};
volPlate[] = Volume In BoundingBox{xPlate0-eps, yPlate0-eps, zPlate0-eps,
                                   xPlate1+eps, yPlate1+eps, zPlate1+eps};

volAllInside[] = {volCoil1[], volCoil2[], volPlate[]};
volAir[] = Volume{:};
volAir[] -= {volAllInside[]};

// ---------------- Mesh refinement ----------------
solidBnd[] = Boundary{ Volume{volAllInside[]}; };

Field[1] = Distance;
Field[1].SurfacesList = {solidBnd[]};
Field[1].Sampling = 100;

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lcPlate;
Field[2].LcMax = lcFar;
Field[2].DistMin = 12  * Scale;
Field[2].DistMax = 100 * Scale;

// local refinement in the coil-plate gap
lcGap = 5 * Scale;
gapPadY = 20 * Scale;
gapPadZ = 20 * Scale;
gapPadX = 2  * Scale;

Field[3] = Box;
Field[3].VIn  = lcGap;
Field[3].VOut = lcFar;
Field[3].XMin = ox1 - gapPadX;
Field[3].XMax = xPlate1 + gapPadX;
Field[3].YMin = oy0 - gapPadY;
Field[3].YMax = oy1 + gapPadY;
Field[3].ZMin = zCoil1_0 - gapPadZ;
Field[3].ZMax = zCoil2_1 + gapPadZ;

// extra local refinement around slit zone
Field[5] = Box;
Field[5].VIn  = lcSlit;
Field[5].VOut = lcFar;
Field[5].XMin = xPlate0 - 1 * Scale;
Field[5].XMax = xPlate1 + 1 * Scale;
Field[5].YMin = slitCenter1Y - 20 * Scale;
Field[5].YMax = slitCenter2Y + 20 * Scale;
Field[5].ZMin = slitZ0 - 10 * Scale;
Field[5].ZMax = slitZ1 + 10 * Scale;

MeshSize{ PointsOf{ Volume{volCoil1[], volCoil2[]}; } } = lcCoil;
MeshSize{ PointsOf{ Volume{volPlate[]}; } } = lcPlate;
MeshSize{ PointsOf{ Volume{volAir[]}; } } = lcAir;

Field[4] = Min;
Field[4].FieldsList = {2, 3, 5};

Background Field = 4;

// ---------------- Physical groups ----------------
Physical Volume("DOMAIN") = {volAir[]};
Physical Volume("COIL_1") = {volCoil1[]};
Physical Volume("COIL_2") = {volCoil2[]};
Physical Volume("NONMAG_PLATE") = {volPlate[]};
Physical Volume("All") = {volAir[], volCoil1[], volCoil2[], volPlate[]};

sXmin[] = Surface In BoundingBox{xAir0-eps, yAir0-eps, zAir0-eps, xAir0+eps, yAir1+eps, zAir1+eps};
sXmax[] = Surface In BoundingBox{xAir1-eps, yAir0-eps, zAir0-eps, xAir1+eps, yAir1+eps, zAir1+eps};
sYmin[] = Surface In BoundingBox{xAir0-eps, yAir0-eps, zAir0-eps, xAir1+eps, yAir0+eps, zAir1+eps};
sYmax[] = Surface In BoundingBox{xAir0-eps, yAir1-eps, zAir0-eps, xAir1+eps, yAir1+eps, zAir1+eps};
sZmin[] = Surface In BoundingBox{xAir0-eps, yAir0-eps, zAir0-eps, xAir1+eps, yAir1+eps, zAir0+eps};
sZmax[] = Surface In BoundingBox{xAir0-eps, yAir0-eps, zAir1-eps, xAir1+eps, yAir1+eps, zAir1+eps};
Physical Surface("AIR_OUTER_WALLS") = {sXmin[], sXmax[], sYmin[], sYmax[], sZmin[], sZmax[]};

allAirBnd[] = Boundary{ Volume{volAir[]}; };
interfaceBnd[] = allAirBnd[];
interfaceBnd[] -= {sXmin[], sXmax[], sYmin[], sYmax[], sZmin[], sZmax[]};
Physical Surface("INTERFACE_AIR_SOLID") = {interfaceBnd[]};

Mesh 3;