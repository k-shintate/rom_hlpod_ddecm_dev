SetFactory("OpenCASCADE");

// =====================================================
// TEAM Problem 7 mesh (structured-refined version)
// Modified: split the racetrack coil into two volumes
//           COIL_1 and COIL_2 (left / right split)
// Units are meters.
// =====================================================

Mesh.Algorithm3D = 1;
Mesh.Algorithm   = 6;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;

Scale = 0.001;

// ---------------- Base mesh sizes ----------------
lcCoil = 6  * Scale;
lcHole = 6  * Scale;
lcCond = 8  * Scale;
lcNear = 12 * Scale;
lcMid  = 24 * Scale;
lcFar  = 80 * Scale;
lcCoil = 12 * Scale;
lcHole = 12 * Scale;
lcCond = 16 * Scale;
lcNear = 24 * Scale;
lcMid  = 48 * Scale;
lcFar  = 160 * Scale;

// ---------------- Geometry ----------------
aluX = 294 * Scale;
aluY = 294 * Scale;
aluZ = 19  * Scale;

holeX0 = 18  * Scale;
holeY0 = 18  * Scale;
holeX1 = 126 * Scale;
holeY1 = 126 * Scale;

coilGapZ    = 30  * Scale;
coilHeightZ = 100 * Scale;
zCoil0 = aluZ + coilGapZ;
zCoil1 = zCoil0 + coilHeightZ;

// Coil race-track from Fig.1
rOuter = 50 * Scale;
rInner = 25 * Scale;

// Exact racetrack dimensions from Problem 7 Fig.1:
// inner corner radius = 25 mm, outer corner radius = 50 mm,
// inner straight length = 150 mm in both x and y directions.
// Therefore:
//   outer overall size = 150 + 2*50 = 250 mm
//   inner overall size = 150 + 2*25 = 200 mm
// Placement follows the benchmark sketch: outer right edge at x = 294 mm
// and outer bottom edge at y = 0 mm.
ox0 = 94  * Scale;
ox1 = 294 * Scale;
oy0 = 0   * Scale;
oy1 = 200 * Scale;

ix0 = 119 * Scale;
ix1 = 269 * Scale;
iy0 = 25  * Scale;
iy1 = 175 * Scale;

// ---------------- Conductor with hole ----------------
Box(1) = {0, 0, 0, aluX, aluY, aluZ};
Box(2) = {holeX0, holeY0, 0, holeX1-holeX0, holeY1-holeY0, aluZ};
condCut[] = BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; };
condVol[] = {condCut[]};

// ---------------- Coil (split into left/right) ----------------
z0 = zCoil0;
xSplit = 0.5 * (ox0 + ox1);

// ----- outer racetrack split points -----
p1  = newp; Point(p1)  = {ox0+rOuter, oy0,         z0, lcCoil};
p2  = newp; Point(p2)  = {xSplit,     oy0,         z0, lcCoil};
p3  = newp; Point(p3)  = {ox1-rOuter, oy0,         z0, lcCoil};
p4  = newp; Point(p4)  = {ox1,        oy0+rOuter,  z0, lcCoil};
p5  = newp; Point(p5)  = {ox1,        oy1-rOuter,  z0, lcCoil};
p6  = newp; Point(p6)  = {ox1-rOuter, oy1,         z0, lcCoil};
p7  = newp; Point(p7)  = {xSplit,     oy1,         z0, lcCoil};
p8  = newp; Point(p8)  = {ox0+rOuter, oy1,         z0, lcCoil};
p9  = newp; Point(p9)  = {ox0,        oy1-rOuter,  z0, lcCoil};
p10 = newp; Point(p10) = {ox0,        oy0+rOuter,  z0, lcCoil};

cTR = newp; Point(cTR) = {ox1-rOuter, oy1-rOuter, z0, lcCoil};
cBR = newp; Point(cBR) = {ox1-rOuter, oy0+rOuter, z0, lcCoil};
cBL = newp; Point(cBL) = {ox0+rOuter, oy0+rOuter, z0, lcCoil};
cTL = newp; Point(cTL) = {ox0+rOuter, oy1-rOuter, z0, lcCoil};

// outer-right boundary
l_or_1 = newl; Line(l_or_1)   = {p2,p3};
a_or_1 = newl; Circle(a_or_1) = {p3,cBR,p4};
l_or_2 = newl; Line(l_or_2)   = {p4,p5};
a_or_2 = newl; Circle(a_or_2) = {p5,cTR,p6};
l_or_3 = newl; Line(l_or_3)   = {p6,p7};

// outer-left boundary
l_ol_1 = newl; Line(l_ol_1)   = {p7,p8};
a_ol_1 = newl; Circle(a_ol_1) = {p8,cTL,p9};
l_ol_2 = newl; Line(l_ol_2)   = {p9,p10};
a_ol_2 = newl; Circle(a_ol_2) = {p10,cBL,p1};
l_ol_3 = newl; Line(l_ol_3)   = {p1,p2};

// ----- inner racetrack split points -----
q1  = newp; Point(q1)  = {ix0+rInner, iy0,         z0, lcCoil};
q2  = newp; Point(q2)  = {xSplit,     iy0,         z0, lcCoil};
q3  = newp; Point(q3)  = {ix1-rInner, iy0,         z0, lcCoil};
q4  = newp; Point(q4)  = {ix1,        iy0+rInner,  z0, lcCoil};
q5  = newp; Point(q5)  = {ix1,        iy1-rInner,  z0, lcCoil};
q6  = newp; Point(q6)  = {ix1-rInner, iy1,         z0, lcCoil};
q7  = newp; Point(q7)  = {xSplit,     iy1,         z0, lcCoil};
q8  = newp; Point(q8)  = {ix0+rInner, iy1,         z0, lcCoil};
q9  = newp; Point(q9)  = {ix0,        iy1-rInner,  z0, lcCoil};
q10 = newp; Point(q10) = {ix0,        iy0+rInner,  z0, lcCoil};

dTR = newp; Point(dTR) = {ix1-rInner, iy1-rInner, z0, lcCoil};
dBR = newp; Point(dBR) = {ix1-rInner, iy0+rInner, z0, lcCoil};
dBL = newp; Point(dBL) = {ix0+rInner, iy0+rInner, z0, lcCoil};
dTL = newp; Point(dTL) = {ix0+rInner, iy1-rInner, z0, lcCoil};

// inner-right boundary
l_ir_1 = newl; Line(l_ir_1)   = {q2,q3};
a_ir_1 = newl; Circle(a_ir_1) = {q3,dBR,q4};
l_ir_2 = newl; Line(l_ir_2)   = {q4,q5};
a_ir_2 = newl; Circle(a_ir_2) = {q5,dTR,q6};
l_ir_3 = newl; Line(l_ir_3)   = {q6,q7};

// inner-left boundary
l_il_1 = newl; Line(l_il_1)   = {q7,q8};
a_il_1 = newl; Circle(a_il_1) = {q8,dTL,q9};
l_il_2 = newl; Line(l_il_2)   = {q9,q10};
a_il_2 = newl; Circle(a_il_2) = {q10,dBL,q1};
l_il_3 = newl; Line(l_il_3)   = {q1,q2};

// split connectors
cut_bot = newl; Line(cut_bot) = {p2,q2};
cut_top = newl; Line(cut_top) = {q7,p7};

// left coil face
clLeft = newll;
Curve Loop(clLeft) = {
  l_ol_1, a_ol_1, l_ol_2, a_ol_2, l_ol_3,
  cut_bot,
  -l_il_3, -a_il_2, -l_il_2, -a_il_1, -l_il_1,
  cut_top
};
sLeft = news; Plane Surface(sLeft) = {clLeft};

// right coil face
clRight = newll;
Curve Loop(clRight) = {
  l_or_1, a_or_1, l_or_2, a_or_2, l_or_3,
  -cut_top,
  -l_ir_3, -a_ir_2, -l_ir_2, -a_ir_1, -l_ir_1,
  -cut_bot
};
sRight = news; Plane Surface(sRight) = {clRight};

// extrude separately
outL[] = Extrude {0,0,coilHeightZ} { Surface{sLeft}; };
outR[] = Extrude {0,0,coilHeightZ} { Surface{sRight}; };

coilVol1[] = {outL[1]};
coilVol2[] = {outR[1]};
coilVol[]  = {coilVol1[], coilVol2[]};

// ---------------- Air box ----------------
xAir0 = -1353 * Scale;
xAir1 =  1647 * Scale;
yAir0 = -1353 * Scale;
yAir1 =  1647 * Scale;
zAir0 =  -300 * Scale;
zAir1 =   449 * Scale;

Box(100) = {xAir0, yAir0, zAir0, xAir1-xAir0, yAir1-yAir0, zAir1-zAir0};
frag[] = BooleanFragments{ Volume{100}; Delete; }{ Volume{condVol[], coilVol[]}; Delete; };
Coherence;

eps = 1e-6;
volCond[] = Volume In BoundingBox{-eps, -eps, -eps, aluX+eps, aluY+eps, aluZ+eps};
volCoil1[] = Volume In BoundingBox{ox0-eps, oy0-eps, zCoil0-eps, xSplit+eps, oy1+eps, zCoil1+eps};
volCoil2[] = Volume In BoundingBox{xSplit-eps, oy0-eps, zCoil0-eps, ox1+eps, oy1+eps, zCoil1+eps};
volCoil[]  = {volCoil1[], volCoil2[]};

volSolid[] = {volCond[], volCoil[]};
volAir[] = Volume{:};
volAir[] -= {volSolid[]};

// ---------------- Structured-like refinement ----------------
solidBnd[] = Boundary{ Volume{volSolid[]}; };
Field[1] = Distance;
Field[1].SurfacesList = {solidBnd[]};
Field[1].Sampling = 100;

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lcNear;
Field[2].LcMax = lcFar;
Field[2].DistMin = 10 * Scale;
Field[2].DistMax = 100 * Scale;

Field[3] = Box;
Field[3].VIn  = lcCond;
Field[3].VOut = lcFar;
Field[3].XMin = -10 * Scale;
Field[3].XMax = aluX + 10 * Scale;
Field[3].YMin = -10 * Scale;
Field[3].YMax = aluY + 10 * Scale;
Field[3].ZMin = -5  * Scale;
Field[3].ZMax = zCoil1 + 5 * Scale;

Field[4] = Box;
Field[4].VIn  = lcHole;
Field[4].VOut = lcFar;
Field[4].XMin = holeX0 - 8 * Scale;
Field[4].XMax = holeX1 + 8 * Scale;
Field[4].YMin = holeY0 - 8 * Scale;
Field[4].YMax = holeY1 + 8 * Scale;
Field[4].ZMin = -2 * Scale;
Field[4].ZMax = aluZ + 2 * Scale;

Field[5] = Box;
Field[5].VIn  = lcCoil;
Field[5].VOut = lcFar;
Field[5].XMin = ox0 - 12 * Scale;
Field[5].XMax = ox1 + 12 * Scale;
Field[5].YMin = oy0 - 12 * Scale;
Field[5].YMax = oy1 + 12 * Scale;
Field[5].ZMin = zCoil0 - 10 * Scale;
Field[5].ZMax = zCoil1 + 10 * Scale;

// Characteristic planes read from Fig.2(a),(b)
xList[] = {0,18,36,54,72,90,108,126,144,162,180,198,216,234,252,270,288,294,322,647,-353,1647};
yList[] = {0,18,36,54,72,90,108,126,144,162,180,198,216,234,252,270,288,294,322,647,-353,1647};
zList[] = {-300,0,7,13,19,27,34,41,49,149,249,449};

// Thin refinement slabs around x = const planes
For i In {0:#xList[]-1}
  Field[10+i] = Box;
  Field[10+i].VIn  = lcMid;
  Field[10+i].VOut = lcFar;
  Field[10+i].XMin = xList[i] * Scale - 2 * Scale;
  Field[10+i].XMax = xList[i] * Scale + 2 * Scale;
  Field[10+i].YMin = yAir0;
  Field[10+i].YMax = yAir1;
  Field[10+i].ZMin = zAir0;
  Field[10+i].ZMax = zAir1;
EndFor

// Thin refinement slabs around y = const planes
baseY = 10 + #xList[];
For i In {0:#yList[]-1}
  Field[baseY+i] = Box;
  Field[baseY+i].VIn  = lcMid;
  Field[baseY+i].VOut = lcFar;
  Field[baseY+i].XMin = xAir0;
  Field[baseY+i].XMax = xAir1;
  Field[baseY+i].YMin = yList[i] * Scale - 2 * Scale;
  Field[baseY+i].YMax = yList[i] * Scale + 2 * Scale;
  Field[baseY+i].ZMin = zAir0;
  Field[baseY+i].ZMax = zAir1;
EndFor

// Thin refinement slabs around z = const planes
baseZ = baseY + #yList[];
For i In {0:#zList[]-1}
  Field[baseZ+i] = Box;
  Field[baseZ+i].VIn  = lcMid;
  Field[baseZ+i].VOut = lcFar;
  Field[baseZ+i].XMin = xAir0;
  Field[baseZ+i].XMax = xAir1;
  Field[baseZ+i].YMin = yAir0;
  Field[baseZ+i].YMax = yAir1;
  Field[baseZ+i].ZMin = zList[i] * Scale - 1.5 * Scale;
  Field[baseZ+i].ZMax = zList[i] * Scale + 1.5 * Scale;
EndFor

MeshSize{ PointsOf{ Volume{volCoil[]}; } } = lcCoil;
MeshSize{ PointsOf{ Volume{volCond[]}; } } = lcCond;
MeshSize{ PointsOf{ Volume{volAir[]}; } }  = lcFar;

Field[200] = Min;
Field[200].FieldsList = {
  2, 3, 4, 5,
  10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
  21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
  32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42,
  43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53,
  54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65
};
Background Field = 200;

// ---------------- Physical groups ----------------
Physical Volume("AIR")       = {volAir[]};
Physical Volume("COIL_1")    = {volCoil1[]};
Physical Volume("COIL_2")    = {volCoil2[]};
Physical Volume("COIL")      = {volCoil1[], volCoil2[]};
Physical Volume("ALUMINUM")  = {volCond[]};
Physical Volume("All")       = {volAir[], volCoil1[], volCoil2[], volCond[]};

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
Physical Surface("AIR_SOLID_INTERFACE") = {interfaceBnd[]};

Mesh 3;