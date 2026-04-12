SetFactory("OpenCASCADE");

// =====================================================
// TEAM Problem 7 coil split into:
//   - COIL_INNER
//   - COIL_LIMB
//   - COIL_CORNER
// Units: meters
// =====================================================

Mesh.Algorithm3D = 1;
Mesh.Algorithm   = 6;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;

Scale = 0.001;

// ---------------- Base mesh sizes ----------------
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

// TEAM-7 exact dimensions
rOuter = 50 * Scale;
rInner = 25 * Scale;

ox0 = 94  * Scale;
ox1 = 294 * Scale;
oy0 = 0   * Scale;
oy1 = 200 * Scale;

ix0 = 119 * Scale;
ix1 = 269 * Scale;
iy0 = 25  * Scale;
iy1 = 175 * Scale;

// centers for corner quarter-annuli
cxBL = ox0 + rOuter;  cyBL = oy0 + rOuter;
cxBR = ox1 - rOuter;  cyBR = oy0 + rOuter;
cxTR = ox1 - rOuter;  cyTR = oy1 - rOuter;
cxTL = ox0 + rOuter;  cyTL = oy1 - rOuter;

// ---------------- Conductor with hole ----------------
Box(1) = {0, 0, 0, aluX, aluY, aluZ};
Box(2) = {holeX0, holeY0, 0, holeX1-holeX0, holeY1-holeY0, aluZ};
condCut[] = BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; };
condVol[] = {condCut[]};

// =====================================================
// 2D coil cross-section split into 9 faces
//   1 inner rectangle
//   4 limbs
//   4 corners
// =====================================================

z0 = zCoil0;

// ---- Inner rectangle ----
Rectangle(1001) = {ix0, iy0, z0, ix1-ix0, iy1-iy0};
innerFace = 1001;

// ---- Straight limbs ----
// bottom
Rectangle(1101) = {ix0 + rInner, oy0, z0, (ix1-ix0) - 2*rInner, rInner};
// right
Rectangle(1102) = {ix1, iy0 + rInner, z0, rOuter-rInner, (iy1-iy0) - 2*rInner};
// top
Rectangle(1103) = {ix0 + rInner, iy1, z0, (ix1-ix0) - 2*rInner, rOuter-rInner};
// left
Rectangle(1104) = {ox0, iy0 + rInner, z0, rOuter-rInner, (iy1-iy0) - 2*rInner};

limbFaces[] = {1101,1102,1103,1104};

// ---- Corner quarter-annuli ----
// Build each as (outer disk sector support) - (inner disk sector support) - trimming boxes
// BL corner
Disk(1201) = {cxBL, cyBL, z0, rOuter, rOuter};
Disk(1202) = {cxBL, cyBL, z0, rInner, rInner};
tmpBL[] = BooleanDifference{ Surface{1201}; Delete; }{ Surface{1202}; Delete; };
Rectangle(1211) = {ox0-rOuter, oy0-rOuter, z0, 2*rOuter, rOuter};
Rectangle(1212) = {ox0-rOuter, oy0-rOuter, z0, rOuter, 2*rOuter};
blKeep1[] = BooleanIntersection{ Surface{tmpBL[]}; Delete; }{ Surface{1211}; Delete; };
blCorner[] = BooleanIntersection{ Surface{blKeep1[]}; Delete; }{ Surface{1212}; Delete; };

// BR corner
Disk(1301) = {cxBR, cyBR, z0, rOuter, rOuter};
Disk(1302) = {cxBR, cyBR, z0, rInner, rInner};
tmpBR[] = BooleanDifference{ Surface{1301}; Delete; }{ Surface{1302}; Delete; };
Rectangle(1311) = {ox1-rOuter, oy0-rOuter, z0, 2*rOuter, rOuter};
Rectangle(1312) = {ox1,        oy0-rOuter, z0, rOuter, 2*rOuter};
brKeep1[] = BooleanIntersection{ Surface{tmpBR[]}; Delete; }{ Surface{1311}; Delete; };
brCorner[] = BooleanIntersection{ Surface{brKeep1[]}; Delete; }{ Surface{1312}; Delete; };

// TR corner
Disk(1401) = {cxTR, cyTR, z0, rOuter, rOuter};
Disk(1402) = {cxTR, cyTR, z0, rInner, rInner};
tmpTR[] = BooleanDifference{ Surface{1401}; Delete; }{ Surface{1402}; Delete; };
Rectangle(1411) = {ox1-rOuter, oy1, z0, 2*rOuter, rOuter};
Rectangle(1412) = {ox1,        oy1-rOuter, z0, rOuter, 2*rOuter};
trKeep1[] = BooleanIntersection{ Surface{tmpTR[]}; Delete; }{ Surface{1411}; Delete; };
trCorner[] = BooleanIntersection{ Surface{trKeep1[]}; Delete; }{ Surface{1412}; Delete; };

// TL corner
Disk(1501) = {cxTL, cyTL, z0, rOuter, rOuter};
Disk(1502) = {cxTL, cyTL, z0, rInner, rInner};
tmpTL[] = BooleanDifference{ Surface{1501}; Delete; }{ Surface{1502}; Delete; };
Rectangle(1511) = {ox0-rOuter, oy1, z0, 2*rOuter, rOuter};
Rectangle(1512) = {ox0-rOuter, oy1-rOuter, z0, rOuter, 2*rOuter};
tlKeep1[] = BooleanIntersection{ Surface{tmpTL[]}; Delete; }{ Surface{1511}; Delete; };
tlCorner[] = BooleanIntersection{ Surface{tlKeep1[]}; Delete; }{ Surface{1512}; Delete; };

cornerFaces[] = {blCorner[0], brCorner[0], trCorner[0], tlCorner[0]};

// optional cleanup
Coherence;

// =====================================================
// Extrude each face set separately
// =====================================================

outInner[] = Extrude {0,0,coilHeightZ} { Surface{innerFace}; };
volInner[] = {outInner[1]};

volLimb[] = {};
For i In {0:#limbFaces[]-1}
  out[] = Extrude {0,0,coilHeightZ} { Surface{limbFaces[i]}; };
  volLimb[] += {out[1]};
EndFor

volCorner[] = {};
For i In {0:#cornerFaces[]-1}
  out[] = Extrude {0,0,coilHeightZ} { Surface{cornerFaces[i]}; };
  volCorner[] += {out[1]};
EndFor

coilVol[] = {volInner[], volLimb[], volCorner[]};

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

volInnerSel[]  = Volume In BoundingBox{ix0-eps, iy0-eps, zCoil0-eps, ix1+eps, iy1+eps, zCoil1+eps};
volCoilAll[]   = Volume In BoundingBox{ox0-eps, oy0-eps, zCoil0-eps, ox1+eps, oy1+eps, zCoil1+eps};

// remove inner from all coil vols to detect outer pieces
volOuterPieces[] = {volCoilAll[]};
volOuterPieces[] -= {volInnerSel[]};

// for robustness, corners occupy the four square neighborhoods near coil corners
volCornerSel[] = {};
volCornerSel[] += Volume In BoundingBox{ox0-eps, oy0-eps, zCoil0-eps, ox0+rOuter+eps, oy0+rOuter+eps, zCoil1+eps};
volCornerSel[] += Volume In BoundingBox{ox1-rOuter-eps, oy0-eps, zCoil0-eps, ox1+eps, oy0+rOuter+eps, zCoil1+eps};
volCornerSel[] += Volume In BoundingBox{ox1-rOuter-eps, oy1-rOuter-eps, zCoil0-eps, ox1+eps, oy1+eps, zCoil1+eps};
volCornerSel[] += Volume In BoundingBox{ox0-eps, oy1-rOuter-eps, zCoil0-eps, ox0+rOuter+eps, oy1+eps, zCoil1+eps};

volLimbSel[] = {volOuterPieces[]};
volLimbSel[] -= {volCornerSel[]};

volSolid[] = {volCond[], volInnerSel[], volLimbSel[], volCornerSel[]};
volAir[] = Volume{:};
volAir[] -= {volSolid[]};

// ---------------- Mesh fields ----------------
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

Field[200] = Min;
Field[200].FieldsList = {2,3,4,5};
Background Field = 200;

// ---------------- Physical groups ----------------
Physical Volume("AIR")          = {volAir[]};
Physical Volume("ALUMINUM")     = {volCond[]};
Physical Volume("COIL_INNER")   = {volInnerSel[]};
Physical Volume("COIL_LIMB")    = {volLimbSel[]};
Physical Volume("COIL_CORNER")  = {volCornerSel[]};
Physical Volume("COIL")         = {volInnerSel[], volLimbSel[], volCornerSel[]};
Physical Volume("All")          = {volAir[], volCond[], volInnerSel[], volLimbSel[], volCornerSel[]};

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