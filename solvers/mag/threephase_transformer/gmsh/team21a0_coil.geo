
SetFactory("OpenCASCADE");

// =====================================================
// TEAM Problem 21a-2 mesh - turn-resolved coil model
//
// Rewritten from the bulk-coil version so that each exciting coil
// is represented by individual conductor turns.
//
// Assumptions used here to obtain exactly 300 turns / coil:
//   - bare conductor cross-section = 2.0 mm (radial) x 6.7 mm (axial)
//   - 17 radial columns x 18 axial rows = 306 slots
//   - 6 slots are intentionally left empty:
//         top row    : radial columns 0,1,2
//         bottom row : radial columns 0,1,2
//     => 306 - 6 = 300 active turns
//
// IMPORTANT:
//   1) If your actual winding layout is different, edit the skip rule
//      in the "Generate turns" loop.
//   2) This .geo resolves each turn as an independent solid volume,
//      but it does NOT impose electrical series connection by itself.
//      That must be done in the FEM solver.
//
// Material constants should be assigned in solver input:
//   NONMAG_PLATE: mu_r = 1, sigma = 1.3889e6 S/m
//   COPPER      : sigma ~ 5.8e7 S/m
// =====================================================

Mesh.Algorithm3D = 1;
Mesh.Algorithm   = 6;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;

Scale = 0.001;

// ---------------- Global mesh sizes ----------------
lcTurn  = 1.0 * Scale;   // individual turn conductor mesh size
lcPlate = 8.0 * Scale;
lcAir   = 30  * Scale;
lcFar   = 60  * Scale;
lcSlit  = 3.0 * Scale;
lcGap   = 5.0 * Scale;

// ---------------- Coil pack outer/inner envelope ----------------
// Same envelope as the original bulk model
coilOuterX = 270 * Scale;
coilOuterY = 270 * Scale;
coilInnerX = 200 * Scale;
coilInnerY = 200 * Scale;
rOuter = 45 * Scale;
rInner = 10 * Scale;

// z direction stacking of the two exciting coils
coilHeight = 217 * Scale;
coilGapZ   = 24  * Scale;

// ---------------- Physical winding specification ----------------
nTurnsTarget = 300;

turnRadial = 2.0 * Scale;   // 2.0 mm
turnAxial  = 6.7 * Scale;   // 6.7 mm

// We use 17 x 18 candidate slots (=306) and remove 6 slots symmetrically
nRad = 17;
nAx  = 18;

// Available radial build from outer to inner envelope
radialAvail = 0.5 * (coilOuterX - coilInnerX); // = 35 mm
axialAvail  = coilHeight;

// distribute unused space as equal margins/gaps
gapRad = (radialAvail - nRad * turnRadial) / (nRad + 1);
gapAx  = (axialAvail  - nAx  * turnAxial ) / (nAx  + 1);

// sanity comments (for user)
// radialAvail = 35 mm -> gapRad ~= 0.0556 mm
// axialAvail  = 217 mm -> gapAx ~= 5.0737 mm

// ---------------- Plate only ----------------
plateHeightZ = 820 * Scale;
plateThk = 10 * Scale;
coilToPlateGap = 12 * Scale;
plateWidthY = 360 * Scale;

// Slits for P21a-2
slitLenZ = 660 * Scale;
slitWidthY = 10 * Scale;
slitDepthX = 10 * Scale;   // equals plate thickness: through slit
slitCenter1Y = -60 * Scale;
slitCenter2Y =  60 * Scale;

// ---------------- Derived global coordinates ----------------
ox0 = -0.5 * coilOuterX;
ox1 =  0.5 * coilOuterX;
oy0 = -0.5 * coilOuterY;
oy1 =  0.5 * coilOuterY;

// Plate placement
xPlate0  = ox1 + coilToPlateGap;
xPlate1  = xPlate0 + plateThk;
yPlate0 = -0.5 * plateWidthY;
yPlate1 =  0.5 * plateWidthY;
zPlate0 = -0.5 * plateHeightZ;
zPlate1 =  0.5 * plateHeightZ;

// Coil placement
zCoil1_0 = -(coilGapZ/2 + coilHeight);
zCoil1_1 = - coilGapZ/2;
zCoil2_0 =   coilGapZ/2;
zCoil2_1 =   coilGapZ/2 + coilHeight;

// Slit z range
slitZ0 = -0.5 * slitLenZ;
slitZ1 =  0.5 * slitLenZ;

// ---------------- Macro: rounded rectangle curve loop at inward offset d ----------------
// Uses globals:
//   dOff, zSec, lcSec
// Returns global:
//   madeLoop
Macro MakeRoundedRectLoop
  x0 = ox0 + dOff;
  x1 = ox1 - dOff;
  y0 = oy0 + dOff;
  y1 = oy1 - dOff;
  rr = rOuter - dOff;

  // Basic validity check: rr must remain positive
  If (rr <= 0)
    Printf("ERROR: invalid inward offset %g for rounded rectangle", dOff/Scale);
  EndIf

  p1 = newp; Point(p1) = {x0 + rr, y0,      zSec, lcSec};
  p2 = newp; Point(p2) = {x1 - rr, y0,      zSec, lcSec};
  p3 = newp; Point(p3) = {x1,      y0 + rr, zSec, lcSec};
  p4 = newp; Point(p4) = {x1,      y1 - rr, zSec, lcSec};
  p5 = newp; Point(p5) = {x1 - rr, y1,      zSec, lcSec};
  p6 = newp; Point(p6) = {x0 + rr, y1,      zSec, lcSec};
  p7 = newp; Point(p7) = {x0,      y1 - rr, zSec, lcSec};
  p8 = newp; Point(p8) = {x0,      y0 + rr, zSec, lcSec};

  cTR = newp; Point(cTR) = {x1 - rr, y1 - rr, zSec, lcSec};
  cBR = newp; Point(cBR) = {x1 - rr, y0 + rr, zSec, lcSec};
  cBL = newp; Point(cBL) = {x0 + rr, y0 + rr, zSec, lcSec};
  cTL = newp; Point(cTL) = {x0 + rr, y1 - rr, zSec, lcSec};

  l1 = newl; Line(l1)   = {p1, p2};
  a1 = newl; Circle(a1) = {p2, cBR, p3};
  l2 = newl; Line(l2)   = {p3, p4};
  a2 = newl; Circle(a2) = {p4, cTR, p5};
  l3 = newl; Line(l3)   = {p5, p6};
  a3 = newl; Circle(a3) = {p6, cTL, p7};
  l4 = newl; Line(l4)   = {p7, p8};
  a4 = newl; Circle(a4) = {p8, cBL, p1};

  madeLoop = newll;
  Curve Loop(madeLoop) = {l1, a1, l2, a2, l3, a3, l4, a4};
Return

// ---------------- Macro: create one turn volume ----------------
// Uses globals:
//   d0Turn, z0Turn, hTurn, lcSec
// Returns global:
//   madeVol
Macro MakeTurnVolume
  dOff = d0Turn;
  zSec = z0Turn;
  Call MakeRoundedRectLoop;
  loopOut = madeLoop;

  dOff = d0Turn + turnRadial;
  zSec = z0Turn;
  Call MakeRoundedRectLoop;
  loopIn = madeLoop;

  sTurn = news;
  Plane Surface(sTurn) = {loopOut, loopIn};

  exTurn[] = Extrude {0, 0, hTurn} { Surface{sTurn}; };
  madeVol = exTurn[1];
Return

// ---------------- Generate individual turns ----------------
coil1Turns[] = {};
coil2Turns[] = {};
allTurnVolumes[] = {};

turnCount1 = 0;
turnCount2 = 0;

// The 6 omitted slots are:
//   j = 0 and j = nAx-1, with i = 0,1,2
// Edit this if your actual winding termination space is different.
For j In {0:nAx-1}
  zBase1 = zCoil1_0 + gapAx + j * (turnAxial + gapAx);
  zBase2 = zCoil2_0 + gapAx + j * (turnAxial + gapAx);

  For i In {0:nRad-1}
    isVoid = ((j == 0 || j == nAx-1) && (i <= 2));

    If (!isVoid)
      d0Turn = gapRad + i * (turnRadial + gapRad);
      hTurn  = turnAxial;
      lcSec  = lcTurn;

      // coil 1 turn
      z0Turn = zBase1;
      Call MakeTurnVolume;
      v1 = madeVol;
      coil1Turns[] += {v1};
      allTurnVolumes[] += {v1};
      turnCount1 += 1;

      // coil 2 turn
      z0Turn = zBase2;
      Call MakeTurnVolume;
      v2 = madeVol;
      coil2Turns[] += {v2};
      allTurnVolumes[] += {v2};
      turnCount2 += 1;
    EndIf
  EndFor
EndFor

Printf("coil1 turns generated = %g", turnCount1);
Printf("coil2 turns generated = %g", turnCount2);

// ---------------- Plate and slits ----------------
Box(10) = {xPlate0, yPlate0, zPlate0, plateThk, plateWidthY, plateHeightZ};

slit1Y0 = slitCenter1Y - 0.5 * slitWidthY;
slit2Y0 = slitCenter2Y - 0.5 * slitWidthY;

Box(11) = {xPlate0, slit1Y0, slitZ0, slitDepthX, slitWidthY, slitLenZ};
Box(12) = {xPlate0, slit2Y0, slitZ0, slitDepthX, slitWidthY, slitLenZ};

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

// ---------------- Fragment air with all solids ----------------
allSolids[] = {allTurnVolumes[], plateVol[]};

frag[] = BooleanFragments{ Volume{100}; Delete; }{ Volume{allSolids[]}; Delete; };
Coherence;

// ---------------- Identify fragmented volumes ----------------
eps = 1e-6;

// Turn volumes are already individual solids, but after fragment we
// recover them via bounding boxes.
volCoil1[] = {};
volCoil2[] = {};

For j In {0:nAx-1}
  zBase1 = zCoil1_0 + gapAx + j * (turnAxial + gapAx);
  zBase2 = zCoil2_0 + gapAx + j * (turnAxial + gapAx);

  For i In {0:nRad-1}
    isVoid = ((j == 0 || j == nAx-1) && (i <= 2));

    If (!isVoid)
      x0bb = ox0 + gapRad + i * (turnRadial + gapRad) - eps;
      x1bb = ox1 - gapRad - i * (turnRadial + gapRad) + eps;
      y0bb = oy0 + gapRad + i * (turnRadial + gapRad) - eps;
      y1bb = oy1 - gapRad - i * (turnRadial + gapRad) + eps;

      // Because the turn is a ring band, use the whole outer box of that band
      // and the specific z range to recover the volume(s).
      vtmp1[] = Volume In BoundingBox{x0bb, y0bb, zBase1-eps, x1bb, y1bb, zBase1+turnAxial+eps};
      vtmp2[] = Volume In BoundingBox{x0bb, y0bb, zBase2-eps, x1bb, y1bb, zBase2+turnAxial+eps};

      volCoil1[] += {vtmp1[]};
      volCoil2[] += {vtmp2[]};
    EndIf
  EndFor
EndFor

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
Field[2].LcMin = lcTurn;
Field[2].LcMax = lcFar;
Field[2].DistMin = 6   * Scale;
Field[2].DistMax = 80  * Scale;

// local refinement in the coil-plate gap
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

Field[4] = Min;
Field[4].FieldsList = {2, 3, 5};
Background Field = 4;

// Explicit point mesh sizes
MeshSize{ PointsOf{ Volume{volAir[]}; } } = lcAir;
MeshSize{ PointsOf{ Volume{volPlate[]}; } } = lcPlate;
MeshSize{ PointsOf{ Volume{volCoil1[], volCoil2[]}; } } = lcTurn;

// ---------------- Physical groups ----------------
Physical Volume("DOMAIN") = {volAir[]};
Physical Volume("NONMAG_PLATE") = {volPlate[]};
Physical Volume("COIL_1_ALL_TURNS") = {volCoil1[]};
Physical Volume("COIL_2_ALL_TURNS") = {volCoil2[]};
Physical Volume("ALL_SOLIDS") = {volCoil1[], volCoil2[], volPlate[]};
Physical Volume("ALL") = {volAir[], volCoil1[], volCoil2[], volPlate[]};

// Outer air walls
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