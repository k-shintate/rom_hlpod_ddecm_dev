SetFactory("OpenCASCADE");

// =====================================================
// TEAM Problem 7 coil split into:
//   - COIL_INNER
//   - COIL_LIMB
//   - COIL_CORNER
// Hole is kept as a separate volume
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

// =====================================================
// Conductor block and hole block
// Keep both as independent volumes using variable tags
// =====================================================

vAl = newv;
Box(vAl) = {0, 0, 0, aluX, aluY, aluZ};

vHole = newv;
Box(vHole) = {holeX0, holeY0, 0, holeX1-holeX0, holeY1-holeY0, aluZ};

// =====================================================
// 2D coil cross-section split into 9 faces
//   1 inner rectangle
//   4 limbs
//   4 corners
// =====================================================

z0 = zCoil0;

// ---- Inner rectangle ----
sInner = news;
Rectangle(sInner) = {ix0, iy0, z0, ix1-ix0, iy1-iy0};
innerFace = sInner;

// ---- Straight limbs ----
sBtm = news;
Rectangle(sBtm) = {ix0 + rInner, oy0, z0, (ix1-ix0) - 2*rInner, rInner};

sRgt = news;
Rectangle(sRgt) = {ix1, iy0 + rInner, z0, rOuter-rInner, (iy1-iy0) - 2*rInner};

sTop = news;
Rectangle(sTop) = {ix0 + rInner, iy1, z0, (ix1-ix0) - 2*rInner, rOuter-rInner};

sLft = news;
Rectangle(sLft) = {ox0, iy0 + rInner, z0, rOuter-rInner, (iy1-iy0) - 2*rInner};

limbFaces[] = {sBtm, sRgt, sTop, sLft};

// ---- Corner quarter-annuli ----
// Use helper rectangles with tags created by news to avoid OCC tag conflicts

// BL corner
sBLout = news; Disk(sBLout) = {cxBL, cyBL, z0, rOuter, rOuter};
sBLin  = news; Disk(sBLin ) = {cxBL, cyBL, z0, rInner, rInner};
tmpBL[] = BooleanDifference{ Surface{sBLout}; Delete; }{ Surface{sBLin}; Delete; };
sBLh = news; Rectangle(sBLh) = {ox0, oy0, z0, 2*rOuter, rOuter};
sBLv = news; Rectangle(sBLv) = {ox0, oy0, z0, rOuter, 2*rOuter};
blKeep1[]  = BooleanIntersection{ Surface{tmpBL[]}; Delete; }{ Surface{sBLh}; Delete; };
blCorner[] = BooleanIntersection{ Surface{blKeep1[]}; Delete; }{ Surface{sBLv}; Delete; };

// BR corner
sBRout = news; Disk(sBRout) = {cxBR, cyBR, z0, rOuter, rOuter};
sBRin  = news; Disk(sBRin ) = {cxBR, cyBR, z0, rInner, rInner};
tmpBR[] = BooleanDifference{ Surface{sBRout}; Delete; }{ Surface{sBRin}; Delete; };
sBRh = news; Rectangle(sBRh) = {ox1 - 2*rOuter, oy0, z0, 2*rOuter, rOuter};
sBRv = news; Rectangle(sBRv) = {ox1 - rOuter,   oy0, z0, rOuter, 2*rOuter};
brKeep1[]  = BooleanIntersection{ Surface{tmpBR[]}; Delete; }{ Surface{sBRh}; Delete; };
brCorner[] = BooleanIntersection{ Surface{brKeep1[]}; Delete; }{ Surface{sBRv}; Delete; };

// TR corner
sTRout = news; Disk(sTRout) = {cxTR, cyTR, z0, rOuter, rOuter};
sTRin  = news; Disk(sTRin ) = {cxTR, cyTR, z0, rInner, rInner};
tmpTR[] = BooleanDifference{ Surface{sTRout}; Delete; }{ Surface{sTRin}; Delete; };
sTRh = news; Rectangle(sTRh) = {ox1 - 2*rOuter, oy1 - rOuter,   z0, 2*rOuter, rOuter};
sTRv = news; Rectangle(sTRv) = {ox1 - rOuter,   oy1 - 2*rOuter, z0, rOuter, 2*rOuter};
trKeep1[]  = BooleanIntersection{ Surface{tmpTR[]}; Delete; }{ Surface{sTRh}; Delete; };
trCorner[] = BooleanIntersection{ Surface{trKeep1[]}; Delete; }{ Surface{sTRv}; Delete; };

// TL corner
sTLout = news; Disk(sTLout) = {cxTL, cyTL, z0, rOuter, rOuter};
sTLin  = news; Disk(sTLin ) = {cxTL, cyTL, z0, rInner, rInner};
tmpTL[] = BooleanDifference{ Surface{sTLout}; Delete; }{ Surface{sTLin}; Delete; };
sTLh = news; Rectangle(sTLh) = {ox0, oy1 - rOuter,   z0, 2*rOuter, rOuter};
sTLv = news; Rectangle(sTLv) = {ox0, oy1 - 2*rOuter, z0, rOuter, 2*rOuter};
tlKeep1[]  = BooleanIntersection{ Surface{tmpTL[]}; Delete; }{ Surface{sTLh}; Delete; };
tlCorner[] = BooleanIntersection{ Surface{tlKeep1[]}; Delete; }{ Surface{sTLv}; Delete; };

// guard
If (#blCorner[] == 0 || #brCorner[] == 0 || #trCorner[] == 0 || #tlCorner[] == 0)
  Error("Corner construction failed");
EndIf

cornerFaces[] = {blCorner[0], brCorner[0], trCorner[0], tlCorner[0]};

// =====================================================
// Extrude each face set separately
// =====================================================

outInner[] = Extrude {0, 0, coilHeightZ} { Surface{innerFace}; };
volInner[] = {outInner[1]};

volLimb[] = {};
For i In {0:#limbFaces[]-1}
  out[] = Extrude {0, 0, coilHeightZ} { Surface{limbFaces[i]}; };
  volLimb[] += {out[1]};
EndFor

volCorner[] = {};
For i In {0:#cornerFaces[]-1}
  out[] = Extrude {0, 0, coilHeightZ} { Surface{cornerFaces[i]}; };
  volCorner[] += {out[1]};
EndFor

coilVol[] = {volInner[], volLimb[], volCorner[]};

// =====================================================
// Air box and fragmentation
// =====================================================

xAir0 = -1353 * Scale;
xAir1 =  1647 * Scale;
yAir0 = -1353 * Scale;
yAir1 =  1647 * Scale;
zAir0 =  -300 * Scale;
zAir1 =   449 * Scale;

vAirBox = newv;
Box(vAirBox) = {xAir0, yAir0, zAir0, xAir1-xAir0, yAir1-yAir0, zAir1-zAir0};

// Use variable tags, not fixed volume IDs
frag[] = BooleanFragments{ Volume{vAirBox}; Delete; }{ Volume{vAl, vHole, coilVol[]}; Delete; };

// Do NOT call Coherence before this point.
// One cleanup here is enough.
Coherence;

eps = 1e-6;

// =====================================================
// Select fragmented volumes
// =====================================================

// hole
volHole[] = Volume In BoundingBox{holeX0-eps, holeY0-eps, -eps, holeX1+eps, holeY1+eps, aluZ+eps};

// conductor-space fragments
volAluminumAll[] = Volume In BoundingBox{-eps, -eps, -eps, aluX+eps, aluY+eps, aluZ+eps};

// subtract hole from conductor-space set to get metal only
volCond[] = {volAluminumAll[]};
volCond[] -= {volHole[]};

// coil selections
volInnerSel[] = Volume In BoundingBox{ix0-eps, iy0-eps, zCoil0-eps, ix1+eps, iy1+eps, zCoil1+eps};
volCoilAll[]  = Volume In BoundingBox{ox0-eps, oy0-eps, zCoil0-eps, ox1+eps, oy1+eps, zCoil1+eps};

// remove inner from all coil volumes to get outer pieces
volOuterPieces[] = {volCoilAll[]};
volOuterPieces[] -= {volInnerSel[]};

// corners occupy the four square neighborhoods near coil corners
volCornerSel[] = {};
volCornerSel[] += Volume In BoundingBox{ox0-eps, oy0-eps, zCoil0-eps, ox0+rOuter+eps, oy0+rOuter+eps, zCoil1+eps};
volCornerSel[] += Volume In BoundingBox{ox1-rOuter-eps, oy0-eps, zCoil0-eps, ox1+eps, oy0+rOuter+eps, zCoil1+eps};
volCornerSel[] += Volume In BoundingBox{ox1-rOuter-eps, oy1-rOuter-eps, zCoil0-eps, ox1+eps, oy1+eps, zCoil1+eps};
volCornerSel[] += Volume In BoundingBox{ox0-eps, oy1-rOuter-eps, zCoil0-eps, ox0+rOuter+eps, oy1+eps, zCoil1+eps};

volLimbSel[] = {volOuterPieces[]};
volLimbSel[] -= {volCornerSel[]};

// Keep hole as air-like region, not solid
volSolid[] = {volCond[], volInnerSel[], volLimbSel[], volCornerSel[]};

volAir[] = Volume{:};
volAir[] -= {volSolid[]};

// =====================================================
// Mesh fields
// =====================================================

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
Field[200].FieldsList = {2, 3, 4, 5};
Background Field = 200;

// =====================================================
// Physical groups
// IDs aligned with solver-side assumptions
//
// 1 : COIL_INNER
// 2 : COIL_LIMB
// 3 : COIL_CORNER
// 4 : ALUMINUM
// 5 : AIR
// 6 : HOLE
// =====================================================

Physical Volume("COIL_INNER",  1) = {volInnerSel[]};
Physical Volume("COIL_LIMB",   2) = {volLimbSel[]};
Physical Volume("COIL_CORNER", 3) = {volCornerSel[]};
Physical Volume("ALUMINUM",    4) = {volCond[]};
Physical Volume("AIR",         5) = {volAir[]};
Physical Volume("HOLE",        6) = {volHole[]};

Physical Volume("COIL", 7) = {volInnerSel[], volLimbSel[], volCornerSel[]};
Physical Volume("All",  8) = {volAir[], volHole[], volCond[], volInnerSel[], volLimbSel[], volCornerSel[]};

// outer air-box walls
sXmin[] = Surface In BoundingBox{xAir0-eps, yAir0-eps, zAir0-eps, xAir0+eps, yAir1+eps, zAir1+eps};
sXmax[] = Surface In BoundingBox{xAir1-eps, yAir0-eps, zAir0-eps, xAir1+eps, yAir1+eps, zAir1+eps};
sYmin[] = Surface In BoundingBox{xAir0-eps, yAir0-eps, zAir0-eps, xAir1+eps, yAir0+eps, zAir1+eps};
sYmax[] = Surface In BoundingBox{xAir0-eps, yAir1-eps, zAir0-eps, xAir1+eps, yAir1+eps, zAir1+eps};
sZmin[] = Surface In BoundingBox{xAir0-eps, yAir0-eps, zAir0-eps, xAir1+eps, yAir1+eps, zAir0+eps};
sZmax[] = Surface In BoundingBox{xAir0-eps, yAir0-eps, zAir1-eps, xAir1+eps, yAir1+eps, zAir1+eps};

aluSurf[] = Boundary{ Volume{volCond[]}; };
//Physical Surface("AIR_OUTER_WALLS") = {aluSurf[]};
Physical Surface("AIR_OUTER_WALLS") = {sXmin[], sXmax[], sYmin[], sYmax[], sZmin[], sZmax[]};

// Air-solid interfaces
allAirBnd[] = Boundary{ Volume{volAir[]}; };
interfaceBnd[] = allAirBnd[];
interfaceBnd[] -= {sXmin[], sXmax[], sYmin[], sYmax[], sZmin[], sZmax[]};
Physical Surface("AIR_SOLID_INTERFACE") = {interfaceBnd[]};

Mesh 3;
