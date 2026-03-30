SetFactory("OpenCASCADE");

// =====================================================
// Unit Conversion (mm -> m)
// =====================================================
Scale = 0.001; 

// =====================================================
// Parameters (Scaled to Meters)
// =====================================================
// メッシュサイズ基準
lc    = 3.0 * Scale;  // デバイス近傍 (1mm)
lcFar = 20.0 * Scale; // 遠方境界 (20mm) - 要素数爆発防止のため粗くする

// --- Open Boundary Margins ---
// デバイスサイズが約100mmなので、片側200mm程度のマージンを取る（計500mm規模）
Margin  = 100 * Scale; 
zMargin = 100 * Scale;

// --- Device Dimensions ---
x1 = -20 * Scale;  x2 = 25 * Scale;  x3 = 45 * Scale;  x4 = 90 * Scale;
y0 = -10 * Scale;  y3 = 80 * Scale;

// Z length (Device Height)
Lz = 40 * Scale; 

leadLen = 5 * Scale;
NyLead  = 4;

// --------------------
// Coil parameters
// --------------------
lenC = 30 * Scale;
yC0  = 30 * Scale;

yStart   = yC0 - 0.5*lenC + (5 * Scale);
yCoilEnd = yStart + lenC;

y1 = yStart   - leadLen;
y2 = yCoilEnd + leadLen;

// --- Expanded Grid for Air Domain ---
// 既存の座標(x1..x4)を内包しつつ、外側にMarginを追加
xMin = x1 - Margin;
xMax = x4 + Margin;
yMin = y0 - Margin;
yMax = y3 + Margin;

// 配列を拡張
x[] = {xMin, x1, x2, x3, x4, xMax};
y[] = {yMin, y0, y1, y2, y3, yMax};

nx = #x[];
ny = #y[];

// =====================================================
// Point grid
// =====================================================
p[] = {};
For j In {0:ny-1}
  For i In {0:nx-1}
    id = newp;
    
    // メッシュサイズ制御: グリッドの外縁部は粗く(lcFar)、内部は細かく(lc)
    local_lc = lc;
    If (i == 0 || i == nx-1 || j == 0 || j == ny-1)
      local_lc = lcFar;
    EndIf
    
    Point(id) = {x[i], y[j], 0, local_lc};
    p[] += id;
  EndFor
EndFor

// =====================================================
// Lines on the grid
// =====================================================
nxSeg = nx - 1;
lh[] = {};
For j In {0:ny-1}
  For i In {0:nx-2}
    l = newl;
    Line(l) = { p[i + nx*j], p[i+1 + nx*j] };
    lh[] += l;
  EndFor
EndFor

lv[] = {};
For j In {0:ny-2}
  For i In {0:nx-1}
    l = newl;
    Line(l) = { p[i + nx*j], p[i + nx*(j+1)] };
    lv[] += l;
  EndFor
EndFor

// =====================================================
// Build block surfaces (Base Layer at Z=0)
// =====================================================
surfs[] = {};
For j In {0:ny-2}
  For i In {0:nx-2}
    bottom = lh[i + nxSeg*j];
    top    = lh[i + nxSeg*(j+1)];
    left   = lv[i + nx*j];
    right  = lv[i+1 + nx*j];

    cl = newreg;
    Curve Loop(cl) = { bottom, right, -top, -left };
    s = newreg;
    Plane Surface(s) = {cl};
    surfs[] += s;
  EndFor
EndFor

// =====================================================
// Extrude DOMAIN surfaces (3 Layers: Bottom, Mid, Top)
// =====================================================

// 1. Bottom Air (Z < 0) - Pure Air
outDomBot[] = Extrude {0,0, -zMargin} {
  Surface{surfs[]};
};

volAirBot[] = {};
For k In {0:#surfs[]-1}
  volAirBot[] += outDomBot[1 + 6*k];
EndFor

// 2. Middle Air (0 <= Z <= Lz) - Device Zone
// ここにコイルやコアが埋め込まれるため、後でBooleanDifferenceを行う
outDomMid[] = Extrude {0,0, Lz} {
  Surface{surfs[]};
};

volAirMid[] = {};
topSurfsMid[] = {}; // 上層への押し出し用に中層の上面を取得
For k In {0:#surfs[]-1}
  topSurfsMid[] += outDomMid[0 + 6*k]; // Index 0 is the top surface
  volAirMid[]   += outDomMid[1 + 6*k]; // Index 1 is the volume
EndFor

// 3. Top Air (Z > Lz) - Pure Air
// 中層の上面からさらに上へ押し出す
outDomTop[] = Extrude {0,0, zMargin} {
  Surface{topSurfsMid[]};
};

volAirTop[] = {};
For k In {0:#surfs[]-1}
  volAirTop[] += outDomTop[1 + 6*k];
EndFor


// =====================================================
// 3 coils + CORES
// =====================================================
xC[] = { 0.5*(x1+x2), 0.5*(x2+x3), 0.5*(x3+x4) };

zMid = 0.5*Lz;
rC[] = { 9 * Scale, 9 * Scale, 9 * Scale }; 

Nc      = 16;
Nr      = 6;
NrCore  = 4;
NyC     = 20;

rInFactor   = 0.5;
rCoreFactor = 0.5;

ang = Pi/4;
ct  = Cos(ang);
st  = Sin(ang);

volCore[] = {};
volCoil[] = {};

coreDnFace[] = {};
coreUpFace[] = {};

For c In {0:2}

  rOut = rC[c];
  rIn  = rInFactor * rOut;

  pc = newp; Point(pc)  = {xC[c], yStart, zMid, lc};

  pEo = newp; Point(pEo) = {xC[c] + rOut*ct, yStart, zMid + rOut*st, lc};
  pNo = newp; Point(pNo) = {xC[c] - rOut*st, yStart, zMid + rOut*ct, lc};
  pWo = newp; Point(pWo) = {xC[c] - rOut*ct, yStart, zMid - rOut*st, lc};
  pSo = newp; Point(pSo) = {xC[c] + rOut*st, yStart, zMid - rOut*ct, lc};

  pEi = newp; Point(pEi) = {xC[c] + rIn*ct, yStart, zMid + rIn*st, lc};
  pNi = newp; Point(pNi) = {xC[c] - rIn*st, yStart, zMid + rIn*ct, lc};
  pWi = newp; Point(pWi) = {xC[c] - rIn*ct, yStart, zMid - rIn*st, lc};
  pSi = newp; Point(pSi) = {xC[c] + rIn*st, yStart, zMid - rIn*ct, lc};

  aO0=newl; Circle(aO0)={pEo,pc,pNo};
  aO1=newl; Circle(aO1)={pNo,pc,pWo};
  aO2=newl; Circle(aO2)={pWo,pc,pSo};
  aO3=newl; Circle(aO3)={pSo,pc,pEo};

  aI0=newl; Circle(aI0)={pEi,pc,pNi};
  aI1=newl; Circle(aI1)={pNi,pc,pWi};
  aI2=newl; Circle(aI2)={pWi,pc,pSi};
  aI3=newl; Circle(aI3)={pSi,pc,pEi};

  rE=newl; Line(rE)={pEi,pEo};
  rN=newl; Line(rN)={pNi,pNo};
  rW=newl; Line(rW)={pWi,pWo};
  rS=newl; Line(rS)={pSi,pSo};

  rD = rCoreFactor * rIn;
  pEd = newp; Point(pEd) = {xC[c] + rD*ct, yStart, zMid + rD*st, lc};
  pNd = newp; Point(pNd) = {xC[c] - rD*st, yStart, zMid + rD*ct, lc};
  pWd = newp; Point(pWd) = {xC[c] - rD*ct, yStart, zMid - rD*st, lc};
  pSd = newp; Point(pSd) = {xC[c] + rD*st, yStart, zMid - rD*ct, lc};

  d0=newl; Line(d0)={pEd,pNd};
  d1=newl; Line(d1)={pNd,pWd};
  d2=newl; Line(d2)={pWd,pSd};
  d3=newl; Line(d3)={pSd,pEd};

  cE=newl; Line(cE)={pEd,pEi};
  cN=newl; Line(cN)={pNd,pNi};
  cW=newl; Line(cW)={pWd,pWi};
  cS=newl; Line(cS)={pSd,pSi};

  clCore=newreg; Curve Loop(clCore) = {d0,d1,d2,d3};
  sCore = newreg; Plane Surface(sCore) = {clCore};

  clIn0=newreg; Curve Loop(clIn0) = { cE,  aI0, -cN, -d0}; sIn0=newreg; Plane Surface(sIn0)={clIn0};
  clIn1=newreg; Curve Loop(clIn1) = { cN,  aI1, -cW, -d1}; sIn1=newreg; Plane Surface(sIn1)={clIn1};
  clIn2=newreg; Curve Loop(clIn2) = { cW,  aI2, -cS, -d2}; sIn2=newreg; Plane Surface(sIn2)={clIn2};
  clIn3=newreg; Curve Loop(clIn3) = { cS,  aI3, -cE, -d3}; sIn3=newreg; Plane Surface(sIn3)={clIn3};

  cl0=newreg; Curve Loop(cl0)={aI0, rN, -aO0, -rE}; s0=newreg; Plane Surface(s0)={cl0};
  cl1=newreg; Curve Loop(cl1)={aI1, rW, -aO1, -rN}; s1=newreg; Plane Surface(s1)={cl1};
  cl2=newreg; Curve Loop(cl2)={aI2, rS, -aO2, -rW}; s2=newreg; Plane Surface(s2)={cl2};
  cl3=newreg; Curve Loop(cl3)={aI3, rE, -aO3, -rS}; s3=newreg; Plane Surface(s3)={cl3};

  // Transfinite Line は保持
  Transfinite Line{aO0,aO1,aO2,aO3,aI0,aI1,aI2,aI3} = Nc + 1;
  Transfinite Line{d0,d1,d2,d3} = Nc + 1;
  Transfinite Line{rE,rN,rW,rS} = Nr + 1;
  Transfinite Line{cE,cN,cW,cS} = NrCore + 1;

  Transfinite Surface{sCore} = {pEd,pNd,pWd,pSd};
  Transfinite Surface{sIn0}  = {pEd,pEi,pNi,pNd};
  Transfinite Surface{sIn1}  = {pNd,pNi,pWi,pWd};
  Transfinite Surface{sIn2}  = {pWd,pWi,pSi,pSd};
  Transfinite Surface{sIn3}  = {pSd,pSi,pEi,pEd};

  Transfinite Surface{s0} = {pEi,pNi,pNo,pEo};
  Transfinite Surface{s1} = {pNi,pWi,pWo,pNo};
  Transfinite Surface{s2} = {pWi,pSi,pSo,pWo};
  Transfinite Surface{s3} = {pSi,pEi,pEo,pSo};

  coilOut[] = Extrude {0, lenC, 0} {
    Surface{sCore,sIn0,sIn1,sIn2,sIn3, s0,s1,s2,s3};
    Layers{NyC};
  };
    
  volCore[] += coilOut[1];
  For k In {1:8}
      volCoil[] += coilOut[1 + 6*k];
  EndFor

  leadDn[] = Extrude {0, -leadLen, 0} { Surface{sCore}; Layers{NyLead}; };
  leadUp[] = Extrude {0,  leadLen, 0} { Surface{coilOut[0]}; Layers{NyLead}; };

  volCore[] += leadDn[1];
  volCore[] += leadUp[1];

  coreDnFace[] += leadDn[0];
  coreUpFace[] += leadUp[0];

EndFor

allCoilCore[] = {volCore[], volCoil[]};

// =====================================================
// YOKES
// =====================================================
yokeLenDn = leadLen;
yokeLenUp = leadLen;
NyYoke    = NyLead;

ctY   = Cos(Pi/4);
aCore = (rCoreFactor * rInFactor * rC[0]) * ctY;

z1c = zMid - aCore;
z2c = zMid + aCore;

x0R = xC[0] + aCore;
x1L = xC[1] - aCore;
x1R = xC[1] + aCore;
x2L = xC[2] - aCore;

hCore   = (2*aCore)/Nc;
gap01L  = (x1L - x0R);
gap12L  = (x2L - x1R);

NxGap01 = Ceil(gap01L / hCore); If(NxGap01 < 1) NxGap01 = 1; EndIf
NxGap12 = Ceil(gap12L / hCore); If(NxGap12 < 1) NxGap12 = 1; EndIf

epsBB = 1e-4; 

// --- Lower Gaps ---
p0b[] = Point In BoundingBox {x0R-epsBB, y1-epsBB, z1c-epsBB, x0R+epsBB, y1+epsBB, z1c+epsBB};
p1b[] = Point In BoundingBox {x1L-epsBB, y1-epsBB, z1c-epsBB, x1L+epsBB, y1+epsBB, z1c+epsBB};
p1t[] = Point In BoundingBox {x1L-epsBB, y1-epsBB, z2c-epsBB, x1L+epsBB, y1+epsBB, z2c+epsBB};
p0t[] = Point In BoundingBox {x0R-epsBB, y1-epsBB, z2c-epsBB, x0R+epsBB, y1+epsBB, z2c+epsBB};

l0Rdn[] = Line In BoundingBox {x0R-epsBB, y1-epsBB, z1c-epsBB, x0R+epsBB, y1+epsBB, z2c+epsBB};
l1Ldn[] = Line In BoundingBox {x1L-epsBB, y1-epsBB, z1c-epsBB, x1L+epsBB, y1+epsBB, z2c+epsBB};

lBot01 = newl; Line(lBot01) = {p0b[0], p1b[0]};
lTop01 = newl; Line(lTop01) = {p1t[0], p0t[0]};

cl01 = newreg;
Curve Loop(cl01) = { lBot01, l1Ldn[0], lTop01, -l0Rdn[0] };
sGapDn01 = newreg; Plane Surface(sGapDn01) = {cl01};

Transfinite Line{lBot01,lTop01} = (NxGap01 + 1);
Transfinite Line{l0Rdn[0],l1Ldn[0]} = (Nc + 1);
Transfinite Surface{sGapDn01} = {p0b[0], p1b[0], p1t[0], p0t[0]};

q1b[] = Point In BoundingBox {x1R-epsBB, y1-epsBB, z1c-epsBB, x1R+epsBB, y1+epsBB, z1c+epsBB};
q2b[] = Point In BoundingBox {x2L-epsBB, y1-epsBB, z1c-epsBB, x2L+epsBB, y1+epsBB, z1c+epsBB};
q2t[] = Point In BoundingBox {x2L-epsBB, y1-epsBB, z2c-epsBB, x2L+epsBB, y1+epsBB, z2c+epsBB};
q1t[] = Point In BoundingBox {x1R-epsBB, y1-epsBB, z2c-epsBB, x1R+epsBB, y1+epsBB, z2c+epsBB};

l1Rdn[] = Line In BoundingBox {x1R-epsBB, y1-epsBB, z1c-epsBB, x1R+epsBB, y1+epsBB, z2c+epsBB};
l2Ldn[] = Line In BoundingBox {x2L-epsBB, y1-epsBB, z1c-epsBB, x2L+epsBB, y1+epsBB, z2c+epsBB};

lBot12 = newl; Line(lBot12) = {q1b[0], q2b[0]};
lTop12 = newl; Line(lTop12) = {q2t[0], q1t[0]};

cl12 = newreg;
Curve Loop(cl12) = { lBot12, l2Ldn[0], lTop12, -l1Rdn[0] };
sGapDn12 = newreg; Plane Surface(sGapDn12) = {cl12};

Transfinite Line{lBot12,lTop12} = (NxGap12 + 1);
Transfinite Line{l1Rdn[0],l2Ldn[0]} = (Nc + 1);
Transfinite Surface{sGapDn12} = {q1b[0], q2b[0], q2t[0], q1t[0]};

// --- Upper Gaps ---
u0b[] = Point In BoundingBox {x0R-epsBB, y2-epsBB, z1c-epsBB, x0R+epsBB, y2+epsBB, z1c+epsBB};
u1b[] = Point In BoundingBox {x1L-epsBB, y2-epsBB, z1c-epsBB, x1L+epsBB, y2+epsBB, z1c+epsBB};
u1t[] = Point In BoundingBox {x1L-epsBB, y2-epsBB, z2c-epsBB, x1L+epsBB, y2+epsBB, z2c+epsBB};
u0t[] = Point In BoundingBox {x0R-epsBB, y2-epsBB, z2c-epsBB, x0R+epsBB, y2+epsBB, z2c+epsBB};

l0Rup[] = Line In BoundingBox {x0R-epsBB, y2-epsBB, z1c-epsBB, x0R+epsBB, y2+epsBB, z2c+epsBB};
l1Lup[] = Line In BoundingBox {x1L-epsBB, y2-epsBB, z1c-epsBB, x1L+epsBB, y2+epsBB, z2c+epsBB};

lBotU01 = newl; Line(lBotU01) = {u0b[0], u1b[0]};
lTopU01 = newl; Line(lTopU01) = {u1t[0], u0t[0]};

clu01 = newreg;
Curve Loop(clu01) = { lBotU01, l1Lup[0], lTopU01, -l0Rup[0] };
sGapUp01 = newreg; Plane Surface(sGapUp01) = {clu01};

Transfinite Line{lBotU01,lTopU01} = (NxGap01 + 1);
Transfinite Line{l0Rup[0],l1Lup[0]} = (Nc + 1);
Transfinite Surface{sGapUp01} = {u0b[0], u1b[0], u1t[0], u0t[0]};

v1b[] = Point In BoundingBox {x1R-epsBB, y2-epsBB, z1c-epsBB, x1R+epsBB, y2+epsBB, z1c+epsBB};
v2b[] = Point In BoundingBox {x2L-epsBB, y2-epsBB, z1c-epsBB, x2L+epsBB, y2+epsBB, z1c+epsBB};
v2t[] = Point In BoundingBox {x2L-epsBB, y2-epsBB, z2c-epsBB, x2L+epsBB, y2+epsBB, z2c+epsBB};
v1t[] = Point In BoundingBox {x1R-epsBB, y2-epsBB, z2c-epsBB, x1R+epsBB, y2+epsBB, z2c+epsBB};

l1Rup[] = Line In BoundingBox {x1R-epsBB, y2-epsBB, z1c-epsBB, x1R+epsBB, y2+epsBB, z2c+epsBB};
l2Lup[] = Line In BoundingBox {x2L-epsBB, y2-epsBB, z1c-epsBB, x2L+epsBB, y2+epsBB, z2c+epsBB};

lBotU12 = newl; Line(lBotU12) = {v1b[0], v2b[0]};
lTopU12 = newl; Line(lTopU12) = {v2t[0], v1t[0]};

clu12 = newreg;
Curve Loop(clu12) = { lBotU12, l2Lup[0], lTopU12, -l1Rup[0] };
sGapUp12 = newreg; Plane Surface(sGapUp12) = {clu12};

Transfinite Line{lBotU12,lTopU12} = (NxGap12 + 1);
Transfinite Line{l1Rup[0],l2Lup[0]} = (Nc + 1);
Transfinite Surface{sGapUp12} = {v1b[0], v2b[0], v2t[0], v1t[0]};

// Extrude yokes
yokeDnOut[] = Extrude {0, -yokeLenDn, 0} {
  Surface{ coreDnFace[0], sGapDn01, coreDnFace[1], sGapDn12, coreDnFace[2] };
  Layers{NyYoke};
};

yokeUpOut[] = Extrude {0,  yokeLenUp, 0} {
  Surface{ coreUpFace[0], sGapUp01, coreUpFace[1], sGapUp12, coreUpFace[2] };
  Layers{NyYoke};
};

yokeVolAll[] = {};
For i In {0:4}
  yokeVolAll[] += yokeDnOut[1 + 6*i];
  yokeVolAll[] += yokeUpOut[1 + 6*i];
EndFor

volCore[] += yokeVolAll[];

// =====================================================
// DOMAIN - COILS Difference
// =====================================================
toolVols[] = {volCoil[], volCore[]};

// ブーリアン差分は、コイル等が存在する「中層 (volAirMid)」に対してのみ実行する
airVolsMid[] = BooleanDifference{ Volume{volAirMid[]}; Delete; }{ Volume{toolVols[]}; };

// 最終的な空気領域は、下層 + 穴あき中層 + 上層
finalAirVols[] = {volAirBot[], airVolsMid[], volAirTop[]};

// =====================================================
// Mesh Control (Field) - Scaled
// =====================================================
coilBnd[] = Boundary{ Volume{toolVols[]}; };

Field[1] = Distance;
Field[1].SurfacesList = { coilBnd[] };
Field[1].Sampling = 100;

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lc * 1.0; 
Field[2].DistMin = 0.5 * Scale;
Field[2].LcMax = lc * 40.0;
Field[2].DistMax = 30 * Scale;

Background Field = 2;

// =====================================================
// Physical groups definition
// =====================================================

Physical Volume("DOMAIN") = {finalAirVols[]};

volCoil1[] = {};
volCoil2[] = {};
volCoil3[] = {};

For k In {0:7}
  volCoil1[] += volCoil[k];
  volCoil2[] += volCoil[k + 8];
  volCoil3[] += volCoil[k + 16];
EndFor

Physical Volume("COIL_1") = {volCoil1[]};
Physical Volume("COIL_2") = {volCoil2[]};
Physical Volume("COIL_3") = {volCoil3[]};

Physical Volume("CORE_IRON") = {volCore[]};
Physical Volume("All") = {finalAirVols[], volCoil[], volCore[]};

// 2. Air Outer Boundaries (Scaled eps)
// 新しい拡張境界(xMin, xMax等)に基づいてBoundingBoxを設定
eps = 1e-4; 

sXmin[] = Surface In BoundingBox {xMin-eps, yMin-eps, -zMargin-eps,     xMin+eps, yMax+eps, Lz+zMargin+eps};
sXmax[] = Surface In BoundingBox {xMax-eps, yMin-eps, -zMargin-eps,     xMax+eps, yMax+eps, Lz+zMargin+eps};

sYmin[] = Surface In BoundingBox {xMin-eps, yMin-eps, -zMargin-eps,     xMax+eps, yMin+eps, Lz+zMargin+eps};
sYmax[] = Surface In BoundingBox {xMin-eps, yMax-eps, -zMargin-eps,     xMax+eps, yMax+eps, Lz+zMargin+eps};

sZmin[] = Surface In BoundingBox {xMin-eps, yMin-eps, -zMargin-eps,     xMax+eps, yMax+eps, -zMargin+eps};
sZmax[] = Surface In BoundingBox {xMin-eps, yMin-eps, Lz+zMargin-eps,   xMax+eps, yMax+eps, Lz+zMargin+eps};

allOuterBnd[] = {sXmin[], sXmax[], sYmin[], sYmax[], sZmin[], sZmax[]};
Physical Surface("AIR_OUTER_WALLS") = {allOuterBnd[]};

// 3. Air-Solid Interface
// 全ての空気ボリュームの境界を取得
allAirBnd[] = Boundary{ Volume{finalAirVols[]}; };

// 外壁を除外して、界面のみを残す（重複除去のためBoolean的な処理が必要だが、Gmshリスト操作で簡易に行う）
// 注意: リスト引き算 -= は IDベースで動作する
interfaceBnd[] = allAirBnd[];
interfaceBnd[] -= allOuterBnd[];

Physical Surface("INTERFACE_AIR_SOLID") = {interfaceBnd[]};