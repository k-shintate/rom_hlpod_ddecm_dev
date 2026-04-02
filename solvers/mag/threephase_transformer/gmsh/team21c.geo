SetFactory("OpenCASCADE");

// =====================================================
// TEAM Problem 21c mesh (geometry from Fig. A1-4)
// Rounded coil corners added from the drawing:
//   outer corner radius = R45
//   inner window radius = R10
// x-direction alignment corrected: coil cross-section centered at x = 0
// =====================================================

Mesh.Algorithm3D = 1;
Mesh.Algorithm   = 6;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;

Scale = 0.001;

// ---------------- Parameters ----------------
ShieldType = 1;  // 1: single shield plate (EM1/M1), 2: three-piece shield (EM2/M2)

lcCoil  = 12 * Scale;
lcPlate = 10 * Scale;
lcAir   = 30 * Scale;
lcFar   = 60 * Scale;

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
plateHeightZ = 520 * Scale;

// Plate and shield
plateThk = 10 * Scale;
shieldThk = 6 * Scale;
//coilToShieldGap = 6 * Scale;
coilToShieldGap = 12 * Scale;
// shieldWidthType1 = 360 * Scale;
shieldWidthType1 = 270 * Scale;
shieldStripW = 80 * Scale;
shieldStripGap = 15 * Scale;
plateWidthY = 360 * Scale;
shieldHeightZ = 458 * Scale;


// ---------------- Coil cross-section placement ----------------
// Center the coil at x = 0 and y = 0
ox0 = -0.5 * coilOuterX;
ox1 =  0.5 * coilOuterX;
oy0 = -0.5 * coilOuterY;
oy1 =  0.5 * coilOuterY;

ix0 = -0.5 * coilInnerX;
ix1 =  0.5 * coilInnerX;
iy0 = -0.5 * coilInnerY;
iy1 =  0.5 * coilInnerY;

// ---------------- Placement of shield and plate ----------------
xShield0 = ox1 + coilToShieldGap;
xShield1 = xShield0 + shieldThk;
xPlate0  = xShield1;
xPlate1  = xPlate0 + plateThk;

yPlate0 = -0.5 * plateWidthY;
yPlate1 =  0.5 * plateWidthY;
yShield1a0 = -0.5 * shieldWidthType1;
yShield1a1 =  0.5 * shieldWidthType1;

zPlate0 = -0.5 * plateHeightZ;
zPlate1 =  0.5 * plateHeightZ;
zShield0 = -0.5 * shieldHeightZ;
zShield1 =  0.5 * shieldHeightZ;

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

// ---------------- Magnetic steel plate ----------------
Box(10) = {xPlate0, yPlate0, zPlate0, plateThk, plateWidthY, plateHeightZ};
plateVol[] = {10};

// ---------------- Shields ----------------
shieldVols[] = {};
If (ShieldType == 1)
  Box(20) = {xShield0, yShield1a0, zShield0, shieldThk, shieldWidthType1, shieldHeightZ};
  shieldVols[] = {20};
Else
  yS0 = -0.5*(3*shieldStripW + 2*shieldStripGap);
  Box(21) = {xShield0, yS0, zShield0, shieldThk, shieldStripW, shieldHeightZ};
  Box(22) = {xShield0, yS0 + shieldStripW + shieldStripGap, zShield0, shieldThk, shieldStripW, shieldHeightZ};
  Box(23) = {xShield0, yS0 + 2*(shieldStripW + shieldStripGap), zShield0, shieldThk, shieldStripW, shieldHeightZ};
  shieldVols[] = {21,22,23};
EndIf

// ---------------- Air domain ----------------
marginXneg = 220 * Scale;
marginXpos = 220 * Scale;
marginY    = 220 * Scale;
marginZ    = 220 * Scale;

xAir0 = ox0 - marginXneg;
xAir1 = xPlate1 + marginXpos;
yAir0 = Min(yShield1a0, yPlate0) - marginY;
yAir1 = Max(yShield1a1, yPlate1) + marginY;
zAir0 = zPlate0 - marginZ;
zAir1 = zPlate1 + marginZ;

Box(100) = {xAir0, yAir0, zAir0, xAir1-xAir0, yAir1-yAir0, zAir1-zAir0};

allSolids[] = {coil1[], coil2[], plateVol[], shieldVols[]};
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

If (ShieldType == 1)
  volShield[] = Volume In BoundingBox{xShield0-eps, yShield1a0-eps, zShield0-eps,
                                      xShield1+eps, yShield1a1+eps, zShield1+eps};
Else
  volShield1[] = Volume In BoundingBox{xShield0-eps, yS0-eps, zShield0-eps,
                                       xShield1+eps, yS0+shieldStripW+eps, zShield1+eps};
  volShield2[] = Volume In BoundingBox{xShield0-eps, yS0 + shieldStripW + shieldStripGap - eps, zShield0-eps,
                                       xShield1+eps, yS0 + 2*shieldStripW + shieldStripGap + eps, zShield1+eps};
  volShield3[] = Volume In BoundingBox{xShield0-eps, yS0 + 2*(shieldStripW + shieldStripGap) - eps, zShield0-eps,
                                       xShield1+eps, yS0 + 3*shieldStripW + 2*shieldStripGap + eps, zShield1+eps};
  volShield[] = {volShield1[], volShield2[], volShield3[]};
EndIf

volAllInside[] = {volCoil1[], volCoil2[], volPlate[], volShield[]};
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

// --- local refinement only in the coil-shield gap ---
lcGap = 2 * Scale;   // まずは 2 mm 程度から試す

gapPadY = 20 * Scale;
gapPadZ = 20 * Scale;
gapPadX = 2  * Scale;

// coil outer face ~ shield face の間を囲う Box
Field[3] = Box;
Field[3].VIn  = lcGap;
Field[3].VOut = lcFar;

// x方向: コイル右面からシールド厚みまで少し広めに囲う
Field[3].XMin = ox1 - gapPadX;
Field[3].XMax = xShield1 + gapPadX;

// y方向: コイル断面付近を少し余裕を持って囲う
Field[3].YMin = oy0 - gapPadY;
Field[3].YMax = oy1 + gapPadY;

// z方向: 2つのコイル高さ全体を少し余裕を持って囲う
Field[3].ZMin = zCoil1_0 - gapPadZ;
Field[3].ZMax = zCoil2_1 + gapPadZ;

// 最小メッシュサイズを採用
Field[4] = Min;
Field[4].FieldsList = {2, 3};

Background Field = 4;

MeshSize{ PointsOf{ Volume{volCoil1[], volCoil2[]}; } } = lcCoil;
MeshSize{ PointsOf{ Volume{volPlate[], volShield[]}; } } = lcPlate;
MeshSize{ PointsOf{ Volume{volAir[]}; } } = lcAir;

// ---------------- Physical groups ----------------
Physical Volume("DOMAIN") = {volAir[]};
Physical Volume("COIL_1") = {volCoil1[]};
Physical Volume("COIL_2") = {volCoil2[]};
Physical Volume("MAGNETIC_STEEL") = {volPlate[]};
Physical Volume("SHIELD") = {volShield[]};
Physical Volume("All") = {volAir[], volCoil1[], volCoil2[], volPlate[], volShield[]};

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