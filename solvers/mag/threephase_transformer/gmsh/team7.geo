SetFactory("OpenCASCADE");

// =====================================================
// TEAM Problem 7 coil split into:
//   - COIL_INNER  : rounded rectangle with R25
//   - COIL_LIMB   : straight parts
//   - COIL_CORNER : four quarter-annuli, R50/R25
//
// COIL_INNER is included as coil conductor.
// Units: meters
// =====================================================

Mesh.Algorithm3D = 1;
Mesh.Algorithm   = 6;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;

Scale = 0.001;

// ---------------- Base mesh sizes ----------------
// TEAM 7 精度検証用
// 必要部位を Mesh Field で局所的に細かくする
lcCoil = 6  * Scale;
lcHole = 4  * Scale;
lcCond = 4  * Scale;
lcNear = 8  * Scale;
lcMid  = 16 * Scale;
lcFar  = 80 * Scale;

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

// TEAM-7 dimensions
rOuter = 50 * Scale;
rInner = 25 * Scale;

// Outer rounded rectangle bounding box
ox0 = 94  * Scale;
ox1 = 294 * Scale;
oy0 = 0   * Scale;
oy1 = 200 * Scale;

// Inner rounded rectangle bounding box
ix0 = 119 * Scale;
ix1 = 269 * Scale;
iy0 = 25  * Scale;
iy1 = 175 * Scale;

// Centers for corner quarter-annuli
cxBL = ox0 + rOuter;  cyBL = oy0 + rOuter;
cxBR = ox1 - rOuter;  cyBR = oy0 + rOuter;
cxTR = ox1 - rOuter;  cyTR = oy1 - rOuter;
cxTL = ox0 + rOuter;  cyTL = oy1 - rOuter;

// =====================================================
// Conductor block and hole block
// =====================================================

vAl = newv;
Box(vAl) = {0, 0, 0, aluX, aluY, aluZ};

vHole = newv;
Box(vHole) = {holeX0, holeY0, 0, holeX1-holeX0, holeY1-holeY0, aluZ};

// =====================================================
// 2D coil cross-section
// =====================================================

z0 = zCoil0;

// =====================================================
// COIL_INNER
//
// IMPORTANT:
// COIL_INNER must NOT be a simple rectangle.
// It must be a rounded rectangle with R25 corners.
// Otherwise it overlaps COIL_CORNER and BooleanFragments
// splits it into 5 volumes.
// =====================================================

// Start from inner bounding rectangle
sInnerBase = news;
Rectangle(sInnerBase) = {ix0, iy0, z0, ix1-ix0, iy1-iy0};

// ---- BL cut: square outside R25 arc ----
sIBLsq = news;
Rectangle(sIBLsq) = {ix0, iy0, z0, rInner, rInner};

sIBLdisk = news;
Disk(sIBLdisk) = {cxBL, cyBL, z0, rInner, rInner};

sIBLcut[] = BooleanDifference{ Surface{sIBLsq}; Delete; }
                              { Surface{sIBLdisk}; Delete; };

// ---- BR cut ----
sIBRsq = news;
Rectangle(sIBRsq) = {ix1-rInner, iy0, z0, rInner, rInner};

sIBRdisk = news;
Disk(sIBRdisk) = {cxBR, cyBR, z0, rInner, rInner};

sIBRcut[] = BooleanDifference{ Surface{sIBRsq}; Delete; }
                              { Surface{sIBRdisk}; Delete; };

// ---- TR cut ----
sITRsq = news;
Rectangle(sITRsq) = {ix1-rInner, iy1-rInner, z0, rInner, rInner};

sITRdisk = news;
Disk(sITRdisk) = {cxTR, cyTR, z0, rInner, rInner};

sITRcut[] = BooleanDifference{ Surface{sITRsq}; Delete; }
                              { Surface{sITRdisk}; Delete; };

// ---- TL cut ----
sITLsq = news;
Rectangle(sITLsq) = {ix0, iy1-rInner, z0, rInner, rInner};

sITLdisk = news;
Disk(sITLdisk) = {cxTL, cyTL, z0, rInner, rInner};

sITLcut[] = BooleanDifference{ Surface{sITLsq}; Delete; }
                              { Surface{sITLdisk}; Delete; };

// Subtract the four outside-corner cuts from the inner rectangle
innerTmp1[] = BooleanDifference{ Surface{sInnerBase}; Delete; }
                                { Surface{sIBLcut[]}; Delete; };

innerTmp2[] = BooleanDifference{ Surface{innerTmp1[]}; Delete; }
                                { Surface{sIBRcut[]}; Delete; };

innerTmp3[] = BooleanDifference{ Surface{innerTmp2[]}; Delete; }
                                { Surface{sITRcut[]}; Delete; };

innerRounded[] = BooleanDifference{ Surface{innerTmp3[]}; Delete; }
                                   { Surface{sITLcut[]}; Delete; };

If (#innerRounded[] == 0)
  Error("COIL_INNER rounded rectangle construction failed");
EndIf

sInner = innerRounded[0];

// =====================================================
// COIL_LIMB straight parts
// =====================================================

sBtm = news;
Rectangle(sBtm) = {
  ix0 + rInner,
  oy0,
  z0,
  (ix1-ix0) - 2*rInner,
  rInner
};

sRgt = news;
Rectangle(sRgt) = {
  ix1,
  iy0 + rInner,
  z0,
  rOuter-rInner,
  (iy1-iy0) - 2*rInner
};

sTop = news;
Rectangle(sTop) = {
  ix0 + rInner,
  iy1,
  z0,
  (ix1-ix0) - 2*rInner,
  rOuter-rInner
};

sLft = news;
Rectangle(sLft) = {
  ox0,
  iy0 + rInner,
  z0,
  rOuter-rInner,
  (iy1-iy0) - 2*rInner
};

limbFaces[] = {sBtm, sRgt, sTop, sLft};

// =====================================================
// COIL_CORNER quarter-annuli
// R50 outer, R25 inner
// =====================================================

// ---- BL corner ----
sBLout = news;
Disk(sBLout) = {cxBL, cyBL, z0, rOuter, rOuter};

sBLin = news;
Disk(sBLin) = {cxBL, cyBL, z0, rInner, rInner};

tmpBL[] = BooleanDifference{ Surface{sBLout}; Delete; }
                            { Surface{sBLin};  Delete; };

sBLh = news;
Rectangle(sBLh) = {ox0, oy0, z0, 2*rOuter, rOuter};

sBLv = news;
Rectangle(sBLv) = {ox0, oy0, z0, rOuter, 2*rOuter};

blKeep1[] = BooleanIntersection{ Surface{tmpBL[]}; Delete; }
                               { Surface{sBLh};   Delete; };

blCorner[] = BooleanIntersection{ Surface{blKeep1[]}; Delete; }
                                { Surface{sBLv};     Delete; };

// ---- BR corner ----
sBRout = news;
Disk(sBRout) = {cxBR, cyBR, z0, rOuter, rOuter};

sBRin = news;
Disk(sBRin) = {cxBR, cyBR, z0, rInner, rInner};

tmpBR[] = BooleanDifference{ Surface{sBRout}; Delete; }
                            { Surface{sBRin};  Delete; };

sBRh = news;
Rectangle(sBRh) = {ox1 - 2*rOuter, oy0, z0, 2*rOuter, rOuter};

sBRv = news;
Rectangle(sBRv) = {ox1 - rOuter, oy0, z0, rOuter, 2*rOuter};

brKeep1[] = BooleanIntersection{ Surface{tmpBR[]}; Delete; }
                               { Surface{sBRh};   Delete; };

brCorner[] = BooleanIntersection{ Surface{brKeep1[]}; Delete; }
                                { Surface{sBRv};     Delete; };

// ---- TR corner ----
sTRout = news;
Disk(sTRout) = {cxTR, cyTR, z0, rOuter, rOuter};

sTRin = news;
Disk(sTRin) = {cxTR, cyTR, z0, rInner, rInner};

tmpTR[] = BooleanDifference{ Surface{sTRout}; Delete; }
                            { Surface{sTRin};  Delete; };

sTRh = news;
Rectangle(sTRh) = {ox1 - 2*rOuter, oy1 - rOuter, z0, 2*rOuter, rOuter};

sTRv = news;
Rectangle(sTRv) = {ox1 - rOuter, oy1 - 2*rOuter, z0, rOuter, 2*rOuter};

trKeep1[] = BooleanIntersection{ Surface{tmpTR[]}; Delete; }
                               { Surface{sTRh};   Delete; };

trCorner[] = BooleanIntersection{ Surface{trKeep1[]}; Delete; }
                                { Surface{sTRv};     Delete; };

// ---- TL corner ----
sTLout = news;
Disk(sTLout) = {cxTL, cyTL, z0, rOuter, rOuter};

sTLin = news;
Disk(sTLin) = {cxTL, cyTL, z0, rInner, rInner};

tmpTL[] = BooleanDifference{ Surface{sTLout}; Delete; }
                            { Surface{sTLin};  Delete; };

sTLh = news;
Rectangle(sTLh) = {ox0, oy1 - rOuter, z0, 2*rOuter, rOuter};

sTLv = news;
Rectangle(sTLv) = {ox0, oy1 - 2*rOuter, z0, rOuter, 2*rOuter};

tlKeep1[] = BooleanIntersection{ Surface{tmpTL[]}; Delete; }
                               { Surface{sTLh};   Delete; };

tlCorner[] = BooleanIntersection{ Surface{tlKeep1[]}; Delete; }
                                { Surface{sTLv};     Delete; };

// ---- Guard ----
If (#blCorner[] == 0 || #brCorner[] == 0 || #trCorner[] == 0 || #tlCorner[] == 0)
  Error("Corner construction failed");
EndIf

cornerFaces[] = {blCorner[0], brCorner[0], trCorner[0], tlCorner[0]};

// =====================================================
// Extrude each face set separately
// =====================================================

// COIL_INNER
outInner[] = Extrude {0, 0, coilHeightZ} {
  Surface{sInner};
};
volInner[] = {outInner[1]};

// COIL_LIMB
volLimb[] = {};
For i In {0:#limbFaces[]-1}
  out[] = Extrude {0, 0, coilHeightZ} {
    Surface{limbFaces[i]};
  };
  volLimb[] += {out[1]};
EndFor

// COIL_CORNER
volCorner[] = {};
For i In {0:#cornerFaces[]-1}
  out[] = Extrude {0, 0, coilHeightZ} {
    Surface{cornerFaces[i]};
  };
  volCorner[] += {out[1]};
EndFor

coilVol[] = {volInner[], volLimb[], volCorner[]};

// =====================================================
// Air box and fragmentation
// =====================================================

xAir0 = -1453 * Scale;
xAir1 =  1453 * Scale;
yAir0 = -1453 * Scale;
yAir1 =  1453 * Scale;
zAir0 = -1453 * Scale;
zAir1 =  1453 * Scale;

vAirBox = newv;
Box(vAirBox) = {
  xAir0, yAir0, zAir0,
  xAir1-xAir0,
  yAir1-yAir0,
  zAir1-zAir0
};

frag[] = BooleanFragments{ Volume{vAirBox}; Delete; }
                         { Volume{vAl, vHole, coilVol[]}; Delete; };

Coherence;

eps = 0.1 * Scale;

// =====================================================
// Select fragmented volumes
// =====================================================

// ---- Hole ----
volHole[] = Volume In BoundingBox{
  holeX0-eps, holeY0-eps, -eps,
  holeX1+eps, holeY1+eps, aluZ+eps
};

// ---- Aluminum ----
volAluminumAll[] = Volume In BoundingBox{
  -eps, -eps, -eps,
  aluX+eps, aluY+eps, aluZ+eps
};

volCond[] = {volAluminumAll[]};
volCond[] -= {volHole[]};

// ---- COIL_INNER ----
// Because COIL_INNER is now a rounded rectangle,
// this should select exactly 1 volume.
volInnerSel[] = Volume In BoundingBox{
  ix0-eps, iy0-eps, zCoil0-eps,
  ix1+eps, iy1+eps, zCoil1+eps
};

// ---- All coil-related fragments in outer box ----
volCoilAll[] = Volume In BoundingBox{
  ox0-eps, oy0-eps, zCoil0-eps,
  ox1+eps, oy1+eps, zCoil1+eps
};

// Outer pieces = all coil-region volumes minus inner
volOuterPieces[] = {volCoilAll[]};
volOuterPieces[] -= {volInnerSel[]};

// ---- COIL_CORNER selection ----
volCornerSel[] = {};

volCornerSel[] += Volume In BoundingBox{
  ox0-eps, oy0-eps, zCoil0-eps,
  ox0+rOuter+eps, oy0+rOuter+eps, zCoil1+eps
};

volCornerSel[] += Volume In BoundingBox{
  ox1-rOuter-eps, oy0-eps, zCoil0-eps,
  ox1+eps, oy0+rOuter+eps, zCoil1+eps
};

volCornerSel[] += Volume In BoundingBox{
  ox1-rOuter-eps, oy1-rOuter-eps, zCoil0-eps,
  ox1+eps, oy1+eps, zCoil1+eps
};

volCornerSel[] += Volume In BoundingBox{
  ox0-eps, oy1-rOuter-eps, zCoil0-eps,
  ox0+rOuter+eps, oy1+eps, zCoil1+eps
};

// ---- COIL_LIMB selection ----
volLimbSel[] = {volOuterPieces[]};
volLimbSel[] -= {volCornerSel[]};

// ---- Diagnostics ----
Printf("n volHole       = %g", #volHole[]);
Printf("n volCond       = %g", #volCond[]);
Printf("n volInnerSel   = %g", #volInnerSel[]);
Printf("n volLimbSel    = %g", #volLimbSel[]);
Printf("n volCornerSel  = %g", #volCornerSel[]);

If (#volInnerSel[] != 1)
  Error("volInnerSel selection failed: expected 1 volume");
EndIf

If (#volLimbSel[] != 4)
  Error("volLimbSel selection failed: expected 4 volumes");
EndIf

If (#volCornerSel[] != 4)
  Error("volCornerSel selection failed: expected 4 volumes");
EndIf

// =====================================================
// Solid and air volumes
// =====================================================

volSolid[] = {
  volCond[],
  volInnerSel[],
  volLimbSel[],
  volCornerSel[]
};

volAir[] = Volume{:};
volAir[] -= {volSolid[], volHole[]};

// =====================================================
// Mesh fields for TEAM 7 accuracy verification
// =====================================================

// -----------------------------------------------------
// Accuracy verification mesh sizes
// -----------------------------------------------------
// coarse / medium / fine の精度検証を行う場合は、
// 主に以下の値だけを変更する。

lcCoilFine   = 5.0 * Scale;   // コイル導体まわり
lcAlFine     = 3.5 * Scale;   // アルミ板まわり
lcHoleFine   = 3.0 * Scale;   // 穴まわり
lcGapFine    = 4.0 * Scale;   // コイル-アルミ間の空気
lcLineFine   = 3.0 * Scale;   // 評価線近傍
lcNearFine   = 9.0 * Scale;   // 全固体境界近傍
lcFarFine    = lcFar;

// -----------------------------------------------------
// Distance field from all solid boundaries
// -----------------------------------------------------
// アルミ、穴、コイル、空気との界面近傍を全体的に細かくする。

solidBnd[] = Boundary{ Volume{volSolid[]}; };

Field[1] = Distance;
Field[1].SurfacesList = {solidBnd[]};
Field[1].Sampling = 200;

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lcNearFine;
Field[2].LcMax = lcFarFine;
Field[2].DistMin = 5   * Scale;
Field[2].DistMax = 120 * Scale;

// -----------------------------------------------------
// Fine mesh around aluminum plate
// -----------------------------------------------------
// アルミ板本体およびその近傍。
// 渦電流密度 Jey の精度に直接効く。

Field[3] = Box;
Field[3].VIn  = lcAlFine;
Field[3].VOut = lcFarFine;
Field[3].XMin = -15 * Scale;
Field[3].XMax = aluX + 15 * Scale;
Field[3].YMin = -15 * Scale;
Field[3].YMax = aluY + 15 * Scale;
Field[3].ZMin = -8  * Scale;
Field[3].ZMax = aluZ + 8 * Scale;

// -----------------------------------------------------
// Very fine mesh around the hole in aluminum
// -----------------------------------------------------
// 穴の角・周辺は渦電流が変化しやすいため細かくする。

Field[4] = Box;
Field[4].VIn  = lcHoleFine;
Field[4].VOut = lcFarFine;
Field[4].XMin = holeX0 - 10 * Scale;
Field[4].XMax = holeX1 + 10 * Scale;
Field[4].YMin = holeY0 - 10 * Scale;
Field[4].YMax = holeY1 + 10 * Scale;
Field[4].ZMin = -4 * Scale;
Field[4].ZMax = aluZ + 4 * Scale;

// -----------------------------------------------------
// Fine mesh around coil conductor
// -----------------------------------------------------
// 励磁電流による磁界分布の精度に効く。

Field[5] = Box;
Field[5].VIn  = lcCoilFine;
Field[5].VOut = lcFarFine;
Field[5].XMin = ox0 - 15 * Scale;
Field[5].XMax = ox1 + 15 * Scale;
Field[5].YMin = oy0 - 15 * Scale;
Field[5].YMax = oy1 + 15 * Scale;
Field[5].ZMin = zCoil0 - 15 * Scale;
Field[5].ZMax = zCoil1 + 15 * Scale;

// -----------------------------------------------------
// Fine mesh in air gap between aluminum and coil
// -----------------------------------------------------
// z = aluZ から z = zCoil0 の空気領域。
// コイルとアルミ板の磁気結合に効く。

Field[6] = Box;
Field[6].VIn  = lcGapFine;
Field[6].VOut = lcFarFine;
Field[6].XMin = -20 * Scale;
Field[6].XMax = aluX + 20 * Scale;
Field[6].YMin = -20 * Scale;
Field[6].YMax = aluY + 20 * Scale;
Field[6].ZMin = aluZ - 3 * Scale;
Field[6].ZMax = zCoil0 + 3 * Scale;

// -----------------------------------------------------
// Fine mesh around Bz evaluation line A1-B1
// y = 72 mm, z = 34 mm
// -----------------------------------------------------

Field[7] = Box;
Field[7].VIn  = lcLineFine;
Field[7].VOut = lcFarFine;
Field[7].XMin = -5  * Scale;
Field[7].XMax = 295 * Scale;
Field[7].YMin = 72 * Scale - 4 * Scale;
Field[7].YMax = 72 * Scale + 4 * Scale;
Field[7].ZMin = 34 * Scale - 4 * Scale;
Field[7].ZMax = 34 * Scale + 4 * Scale;

// -----------------------------------------------------
// Fine mesh around Bz evaluation line A2-B2
// y = 144 mm, z = 34 mm
// -----------------------------------------------------

Field[8] = Box;
Field[8].VIn  = lcLineFine;
Field[8].VOut = lcFarFine;
Field[8].XMin = -5  * Scale;
Field[8].XMax = 295 * Scale;
Field[8].YMin = 144 * Scale - 4 * Scale;
Field[8].YMax = 144 * Scale + 4 * Scale;
Field[8].ZMin = 34 * Scale - 4 * Scale;
Field[8].ZMax = 34 * Scale + 4 * Scale;

// -----------------------------------------------------
// Fine mesh around Jey evaluation line A3-B3
// y = 72 mm, z = 19 mm
// upper aluminum surface
// -----------------------------------------------------

Field[9] = Box;
Field[9].VIn  = lcLineFine;
Field[9].VOut = lcFarFine;
Field[9].XMin = -5  * Scale;
Field[9].XMax = 295 * Scale;
Field[9].YMin = 72 * Scale - 4 * Scale;
Field[9].YMax = 72 * Scale + 4 * Scale;
Field[9].ZMin = aluZ - 2 * Scale;
Field[9].ZMax = aluZ + 2 * Scale;

// -----------------------------------------------------
// Fine mesh around Jey evaluation line A4-B4
// y = 72 mm, z = 0 mm
// lower aluminum surface
// -----------------------------------------------------

Field[10] = Box;
Field[10].VIn  = lcLineFine;
Field[10].VOut = lcFarFine;
Field[10].XMin = -5  * Scale;
Field[10].XMax = 295 * Scale;
Field[10].YMin = 72 * Scale - 4 * Scale;
Field[10].YMax = 72 * Scale + 4 * Scale;
Field[10].ZMin = -2 * Scale;
Field[10].ZMax =  2 * Scale;

// -----------------------------------------------------
// Combine all mesh fields
// -----------------------------------------------------
// Min を使うことで、各領域で最も細かいメッシュサイズを採用する。

Field[200] = Min;
Field[200].FieldsList = {
  2,
  3,
  4,
  5,
  6,
  7,
  8,
  9,
  10
};

Background Field = 200;

// =====================================================
// Physical volumes
// =====================================================

Physical Volume("COIL_INNER",  1) = {volInnerSel[]};
Physical Volume("COIL_LIMB",   2) = {volLimbSel[]};
Physical Volume("COIL_CORNER", 3) = {volCornerSel[]};
Physical Volume("ALUMINUM",    4) = {volCond[]};
Physical Volume("HOLE",        5) = {volHole[]};
Physical Volume("AIR",         6) = {volAir[]};

Physical Volume("COIL", 7) = {
  volInnerSel[],
  volLimbSel[],
  volCornerSel[]
};

Physical Volume("All", 8) = {
  volAir[],
  volHole[],
  volCond[],
  volInnerSel[],
  volLimbSel[],
  volCornerSel[]
};

// =====================================================
// Physical surfaces
// =====================================================

coilInnerSurf[]  = Boundary{ Volume{volInnerSel[]}; };
coilLimbSurf[]   = Boundary{ Volume{volLimbSel[]}; };
coilCornerSurf[] = Boundary{ Volume{volCornerSel[]}; };

coilSurf[] = Boundary{
  Volume{
    volInnerSel[],
    volLimbSel[],
    volCornerSel[]
  };
};

Physical Surface("COIL_INNER_SURF",  101) = {coilInnerSurf[]};
Physical Surface("COIL_LIMB_SURF",   102) = {coilLimbSurf[]};
Physical Surface("COIL_CORNER_SURF", 103) = {coilCornerSurf[]};
Physical Surface("COIL_SURF",        107) = {coilSurf[]};

// =====================================================
// Outer air-box walls
// =====================================================

sXmin[] = Surface In BoundingBox{
  xAir0-eps, yAir0-eps, zAir0-eps,
  xAir0+eps, yAir1+eps, zAir1+eps
};

sXmax[] = Surface In BoundingBox{
  xAir1-eps, yAir0-eps, zAir0-eps,
  xAir1+eps, yAir1+eps, zAir1+eps
};

sYmin[] = Surface In BoundingBox{
  xAir0-eps, yAir0-eps, zAir0-eps,
  xAir1+eps, yAir0+eps, zAir1+eps
};

sYmax[] = Surface In BoundingBox{
  xAir0-eps, yAir1-eps, zAir0-eps,
  xAir1+eps, yAir1+eps, zAir1+eps
};

sZmin[] = Surface In BoundingBox{
  xAir0-eps, yAir0-eps, zAir0-eps,
  xAir1+eps, yAir1+eps, zAir0+eps
};

sZmax[] = Surface In BoundingBox{
  xAir0-eps, yAir0-eps, zAir1-eps,
  xAir1+eps, yAir1+eps, zAir1+eps
};

Physical Surface("AIR_OUTER_WALLS", 201) = {
  sXmin[], sXmax[],
  sYmin[], sYmax[],
  sZmin[], sZmax[]
};

// =====================================================
// Air-solid interfaces
// =====================================================

allAirBnd[] = Boundary{ Volume{volAir[]}; };

interfaceBnd[] = allAirBnd[];
interfaceBnd[] -= {
  sXmin[], sXmax[],
  sYmin[], sYmax[],
  sZmin[], sZmax[]
};

Physical Surface("AIR_SOLID_INTERFACE", 202) = {interfaceBnd[]};

// =====================================================
// Mesh
// =====================================================

Mesh 3;
