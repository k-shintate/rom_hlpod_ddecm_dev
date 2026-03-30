SetFactory("OpenCASCADE");

// =====================================================
// 1. Mesh Options (Stability & Quality)
// =====================================================
// HXT(10)はエラーが出やすいため、堅牢なDelaunay(1)を使用
Mesh.Algorithm3D = 1;    // Delaunay
Mesh.Algorithm   = 6;    // Frontal-Delaunay (2D)

// アスペクト比最適化 (四面体の品質向上)
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;

// =====================================================
// 2. Parameters (Scaled to Meters)
// =====================================================
Scale = 0.001; 

// --- Mesh Size Settings ---
lcRef   = 4.0 * Scale;  // 基本サイズ
lcFine  = 2.0 * Scale;  // 表面付近の目標サイズ (1mm)
lcFar   = 10.0 * Scale; // 遠方境界

// --- Dimensions ---
Margin  = 50 * Scale; 
zMargin = 50 * Scale;

x1 = -20 * Scale; x2 = 25 * Scale; x3 = 45 * Scale; x4 = 90 * Scale;
y0 = -10 * Scale; y3 = 80 * Scale;
Lz = 40 * Scale; 

leadLen = 5 * Scale;
lenC = 30 * Scale;
yC0  = 30 * Scale;

yStart   = yC0 - 0.5*lenC + (5 * Scale);
yCoilEnd = yStart + lenC;
y1 = yStart - leadLen;
y2 = yCoilEnd + leadLen;

// --- Expanded Grid for Air Domain ---
xMin = x1 - Margin; xMax = x4 + Margin;
yMin = y0 - Margin; yMax = y3 + Margin;

x[] = {xMin, x1, x2, x3, x4, xMax};
y[] = {yMin, y0, y1, y2, y3, yMax};

nx = #x[]; ny = #y[];

// =====================================================
// 3. Base Grid Generation
// =====================================================
p[] = {};
For j In {0:ny-1}
  For i In {0:nx-1}
    id = newp;
    local_lc = lcRef;
    // 境界付近は粗く
    If (i == 0 || i == nx-1 || j == 0 || j == ny-1)
      local_lc = lcFar;
    EndIf
    Point(id) = {x[i], y[j], 0, local_lc};
    p[] += id;
  EndFor
EndFor

lh[] = {};
For j In {0:ny-1}
  For i In {0:nx-2}
    l = newl; Line(l) = { p[i + nx*j], p[i+1 + nx*j] }; lh[] += l;
  EndFor
EndFor

lv[] = {};
For j In {0:ny-2}
  For i In {0:nx-1}
    l = newl; Line(l) = { p[i + nx*j], p[i + nx*(j+1)] }; lv[] += l;
  EndFor
EndFor

surfs[] = {};
For j In {0:ny-2}
  For i In {0:nx-2}
    b = lh[i + nx*(j+1)- (j+1)]; // bottom (check index logic) - simplified below
    // Re-indexing for clarity based on loops
    // lh index: i + (nx-1)*j
    // lv index: i + nx*j
    idx_b = i + (nx-1)*j;
    idx_t = i + (nx-1)*(j+1);
    idx_l = i + nx*j;
    idx_r = i+1 + nx*j;
    
    cl = newreg; 
    Curve Loop(cl) = { lh[idx_b], lv[idx_r], -lh[idx_t], -lv[idx_l] };
    s = newreg; 
    Plane Surface(s) = {cl};
    surfs[] += s;
  EndFor
EndFor

// =====================================================
// 4. Extrude Air Layers (No Layers -> Tetrahedral)
// =====================================================
// Bottom Air
outDomBot[] = Extrude {0,0, -zMargin} { Surface{surfs[]}; };
volAirBot[] = {};
For k In {0:#surfs[]-1}
  volAirBot[] += outDomBot[1 + 6*k];
EndFor

// Middle Air (Target for cut)
outDomMid[] = Extrude {0,0, Lz} { Surface{surfs[]}; };
volAirMid[] = {};
topSurfsMid[] = {}; 
For k In {0:#surfs[]-1}
  topSurfsMid[] += outDomMid[0 + 6*k];
  volAirMid[]   += outDomMid[1 + 6*k];
EndFor

// Top Air
outDomTop[] = Extrude {0,0, zMargin} { Surface{topSurfsMid[]}; };
volAirTop[] = {};
For k In {0:#surfs[]-1}
  volAirTop[] += outDomTop[1 + 6*k];
EndFor

// =====================================================
// 5. Coils & Cores Generation
// =====================================================
xC[] = { 0.5*(x1+x2), 0.5*(x2+x3), 0.5*(x3+x4) };
zMid = 0.5*Lz;
rC[] = { 9 * Scale, 9 * Scale, 9 * Scale }; 

rInFactor   = 0.5;
rCoreFactor = 0.5;
ang = Pi/4; ct = Cos(ang); st = Sin(ang);

volCore[] = {};
volCoil[] = {};
coreDnFace[] = {};
coreUpFace[] = {};

For c In {0:2}
  rOut = rC[c]; rIn = rInFactor * rOut;

  // Center & Outer Points
  pc = newp; Point(pc) = {xC[c], yStart, zMid, lcRef};
  pEo = newp; Point(pEo) = {xC[c]+rOut*ct, yStart, zMid+rOut*st, lcRef};
  pNo = newp; Point(pNo) = {xC[c]-rOut*st, yStart, zMid+rOut*ct, lcRef};
  pWo = newp; Point(pWo) = {xC[c]-rOut*ct, yStart, zMid-rOut*st, lcRef};
  pSo = newp; Point(pSo) = {xC[c]+rOut*st, yStart, zMid-rOut*ct, lcRef};

  // Inner Points
  pEi = newp; Point(pEi) = {xC[c]+rIn*ct, yStart, zMid+rIn*st, lcRef};
  pNi = newp; Point(pNi) = {xC[c]-rIn*st, yStart, zMid+rIn*ct, lcRef};
  pWi = newp; Point(pWi) = {xC[c]-rIn*ct, yStart, zMid-rIn*st, lcRef};
  pSi = newp; Point(pSi) = {xC[c]+rIn*st, yStart, zMid-rIn*ct, lcRef};

  // Arcs
  aO0=newl; Circle(aO0)={pEo,pc,pNo}; aO1=newl; Circle(aO1)={pNo,pc,pWo};
  aO2=newl; Circle(aO2)={pWo,pc,pSo}; aO3=newl; Circle(aO3)={pSo,pc,pEo};
  aI0=newl; Circle(aI0)={pEi,pc,pNi}; aI1=newl; Circle(aI1)={pNi,pc,pWi};
  aI2=newl; Circle(aI2)={pWi,pc,pSi}; aI3=newl; Circle(aI3)={pSi,pc,pEi};

  // Radial Lines
  rE=newl; Line(rE)={pEi,pEo}; rN=newl; Line(rN)={pNi,pNo};
  rW=newl; Line(rW)={pWi,pWo}; rS=newl; Line(rS)={pSi,pSo};

  // Core Points
  rD = rCoreFactor * rIn;
  pEd = newp; Point(pEd) = {xC[c]+rD*ct, yStart, zMid+rD*st, lcRef};
  pNd = newp; Point(pNd) = {xC[c]-rD*st, yStart, zMid+rD*ct, lcRef};
  pWd = newp; Point(pWd) = {xC[c]-rD*ct, yStart, zMid-rD*st, lcRef};
  pSd = newp; Point(pSd) = {xC[c]+rD*st, yStart, zMid-rD*ct, lcRef};

  // Core Lines
  d0=newl; Line(d0)={pEd,pNd}; d1=newl; Line(d1)={pNd,pWd};
  d2=newl; Line(d2)={pWd,pSd}; d3=newl; Line(d3)={pSd,pEd};
  cE=newl; Line(cE)={pEd,pEi}; cN=newl; Line(cN)={pNd,pNi};
  cW=newl; Line(cW)={pWd,pWi}; cS=newl; Line(cS)={pSd,pSi};

  // Surfaces
  clCore=newreg; Curve Loop(clCore)={d0,d1,d2,d3};
  sCore=newreg; Plane Surface(sCore)={clCore};

  clIn0=newreg; Curve Loop(clIn0)={cE,aI0,-cN,-d0}; sIn0=newreg; Plane Surface(sIn0)={clIn0};
  clIn1=newreg; Curve Loop(clIn1)={cN,aI1,-cW,-d1}; sIn1=newreg; Plane Surface(sIn1)={clIn1};
  clIn2=newreg; Curve Loop(clIn2)={cW,aI2,-cS,-d2}; sIn2=newreg; Plane Surface(sIn2)={clIn2};
  clIn3=newreg; Curve Loop(clIn3)={cS,aI3,-cE,-d3}; sIn3=newreg; Plane Surface(sIn3)={clIn3};

  cl0=newreg; Curve Loop(cl0)={aI0,rN,-aO0,-rE}; s0=newreg; Plane Surface(s0)={cl0};
  cl1=newreg; Curve Loop(cl1)={aI1,rW,-aO1,-rN}; s1=newreg; Plane Surface(s1)={cl1};
  cl2=newreg; Curve Loop(cl2)={aI2,rS,-aO2,-rW}; s2=newreg; Plane Surface(s2)={cl2};
  cl3=newreg; Curve Loop(cl3)={aI3,rE,-aO3,-rS}; s3=newreg; Plane Surface(s3)={cl3};

  // Extrude Coil (No Layers)
  coilOut[] = Extrude {0, lenC, 0} {
    Surface{sCore,sIn0,sIn1,sIn2,sIn3, s0,s1,s2,s3};
  };
  volCore[] += coilOut[1];
  For k In {1:8}
      volCoil[] += coilOut[1 + 6*k];
  EndFor

  // Extrude Leads (No Layers)
  leadDn[] = Extrude {0, -leadLen, 0} { Surface{sCore}; };
  leadUp[] = Extrude {0, leadLen, 0} { Surface{coilOut[0]}; };

  volCore[] += leadDn[1];
  volCore[] += leadUp[1];
  coreDnFace[] += leadDn[0];
  coreUpFace[] += leadUp[0];
EndFor

// =====================================================
// 6. Yokes Generation
// =====================================================
yokeLenDn = leadLen;
yokeLenUp = leadLen;

ctY = Cos(Pi/4);
aCore = (rCoreFactor * rInFactor * rC[0]) * ctY;
z1c = zMid - aCore; z2c = zMid + aCore;
x0R = xC[0] + aCore; x1L = xC[1] - aCore;
x1R = xC[1] + aCore; x2L = xC[2] - aCore;
epsBB = 1e-4; 

// --- Lower Gaps ---
p0b[] = Point In BoundingBox {x0R-epsBB, y1-epsBB, z1c-epsBB, x0R+epsBB, y1+epsBB, z1c+epsBB};
p1b[] = Point In BoundingBox {x1L-epsBB, y1-epsBB, z1c-epsBB, x1L+epsBB, y1+epsBB, z1c+epsBB};
p1t[] = Point In BoundingBox {x1L-epsBB, y1-epsBB, z2c-epsBB, x1L+epsBB, y1+epsBB, z2c+epsBB};
p0t[] = Point In BoundingBox {x0R-epsBB, y1-epsBB, z2c-epsBB, x0R+epsBB, y1+epsBB, z2c+epsBB};
l0Rdn[] = Line In BoundingBox {x0R-epsBB, y1-epsBB, z1c-epsBB, x0R+epsBB, y1+epsBB, z2c+epsBB};
l1Ldn[] = Line In BoundingBox {x1L-epsBB, y1-epsBB, z1c-epsBB, x1L+epsBB, y1+epsBB, z2c+epsBB};

lBot01=newl; Line(lBot01)={p0b[0], p1b[0]}; lTop01=newl; Line(lTop01)={p1t[0], p0t[0]};
cl01=newreg; Curve Loop(cl01)={lBot01, l1Ldn[0], lTop01, -l0Rdn[0]};
sGapDn01=newreg; Plane Surface(sGapDn01)={cl01};

q1b[] = Point In BoundingBox {x1R-epsBB, y1-epsBB, z1c-epsBB, x1R+epsBB, y1+epsBB, z1c+epsBB};
q2b[] = Point In BoundingBox {x2L-epsBB, y1-epsBB, z1c-epsBB, x2L+epsBB, y1+epsBB, z1c+epsBB};
q2t[] = Point In BoundingBox {x2L-epsBB, y1-epsBB, z2c-epsBB, x2L+epsBB, y1+epsBB, z2c+epsBB};
q1t[] = Point In BoundingBox {x1R-epsBB, y1-epsBB, z2c-epsBB, x1R+epsBB, y1+epsBB, z2c+epsBB};
l1Rdn[] = Line In BoundingBox {x1R-epsBB, y1-epsBB, z1c-epsBB, x1R+epsBB, y1+epsBB, z2c+epsBB};
l2Ldn[] = Line In BoundingBox {x2L-epsBB, y1-epsBB, z1c-epsBB, x2L+epsBB, y1+epsBB, z2c+epsBB};

lBot12=newl; Line(lBot12)={q1b[0], q2b[0]}; lTop12=newl; Line(lTop12)={q2t[0], q1t[0]};
cl12=newreg; Curve Loop(cl12)={lBot12, l2Ldn[0], lTop12, -l1Rdn[0]};
sGapDn12=newreg; Plane Surface(sGapDn12)={cl12};

// --- Upper Gaps ---
u0b[] = Point In BoundingBox {x0R-epsBB, y2-epsBB, z1c-epsBB, x0R+epsBB, y2+epsBB, z1c+epsBB};
u1b[] = Point In BoundingBox {x1L-epsBB, y2-epsBB, z1c-epsBB, x1L+epsBB, y2+epsBB, z1c+epsBB};
u1t[] = Point In BoundingBox {x1L-epsBB, y2-epsBB, z2c-epsBB, x1L+epsBB, y2+epsBB, z2c+epsBB};
u0t[] = Point In BoundingBox {x0R-epsBB, y2-epsBB, z2c-epsBB, x0R+epsBB, y2+epsBB, z2c+epsBB};
l0Rup[] = Line In BoundingBox {x0R-epsBB, y2-epsBB, z1c-epsBB, x0R+epsBB, y2+epsBB, z2c+epsBB};
l1Lup[] = Line In BoundingBox {x1L-epsBB, y2-epsBB, z1c-epsBB, x1L+epsBB, y2+epsBB, z2c+epsBB};

lBotU01=newl; Line(lBotU01)={u0b[0], u1b[0]}; lTopU01=newl; Line(lTopU01)={u1t[0], u0t[0]};
clu01=newreg; Curve Loop(clu01)={lBotU01, l1Lup[0], lTopU01, -l0Rup[0]};
sGapUp01=newreg; Plane Surface(sGapUp01)={clu01};

v1b[] = Point In BoundingBox {x1R-epsBB, y2-epsBB, z1c-epsBB, x1R+epsBB, y2+epsBB, z1c+epsBB};
v2b[] = Point In BoundingBox {x2L-epsBB, y2-epsBB, z1c-epsBB, x2L+epsBB, y2+epsBB, z1c+epsBB};
v2t[] = Point In BoundingBox {x2L-epsBB, y2-epsBB, z2c-epsBB, x2L+epsBB, y2+epsBB, z2c+epsBB};
v1t[] = Point In BoundingBox {x1R-epsBB, y2-epsBB, z2c-epsBB, x1R+epsBB, y2+epsBB, z2c+epsBB};
l1Rup[] = Line In BoundingBox {x1R-epsBB, y2-epsBB, z1c-epsBB, x1R+epsBB, y2+epsBB, z2c+epsBB};
l2Lup[] = Line In BoundingBox {x2L-epsBB, y2-epsBB, z1c-epsBB, x2L+epsBB, y2+epsBB, z2c+epsBB};

lBotU12=newl; Line(lBotU12)={v1b[0], v2b[0]}; lTopU12=newl; Line(lTopU12)={v2t[0], v1t[0]};
clu12=newreg; Curve Loop(clu12)={lBotU12, l2Lup[0], lTopU12, -l1Rup[0]};
sGapUp12=newreg; Plane Surface(sGapUp12)={clu12};

// Extrude Yokes (No Layers)
yokeDnOut[] = Extrude {0, -yokeLenDn, 0} {
  Surface{ coreDnFace[0], sGapDn01, coreDnFace[1], sGapDn12, coreDnFace[2] };
};
yokeUpOut[] = Extrude {0, yokeLenUp, 0} {
  Surface{ coreUpFace[0], sGapUp01, coreUpFace[1], sGapUp12, coreUpFace[2] };
};

yokeVolAll[] = {};
For i In {0:4}
  yokeVolAll[] += yokeDnOut[1 + 6*i];
  yokeVolAll[] += yokeUpOut[1 + 6*i];
EndFor
volCore[] += yokeVolAll[];

// =====================================================
// 7. Boolean Difference & Coherence
// =====================================================
toolVols[] = {volCoil[], volCore[]};

// 空気(Mid)からコイル・鉄心をくり抜く
// 戻り値 airVolsMid[] は「穴の開いた空気」
airVolsMid[] = BooleanDifference{ Volume{volAirMid[]}; Delete; }{ Volume{toolVols[]}; };

finalAirVols[] = {volAirBot[], airVolsMid[], volAirTop[]};

// ★ 重要: 重複頂点を結合し、メッシュの整合性を取る
Coherence;

// =====================================================
// 8. Mesh Field Control (Refinement)
// =====================================================
// コイル・鉄心の全表面を取得
// BooleanDifference後でもtoolVols自体は残っており、
// Coherenceで界面ノードが共有されている状態
coilAndCoreBnd[] = Boundary{ Volume{toolVols[]}; };

Field[1] = Distance;
Field[1].SurfacesList = { coilAndCoreBnd[] };
Field[1].Sampling = 100;

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lcFine;        // 表面の細かさ (1mm)
Field[2].DistMin = 2.0 * Scale; // 2mmまでは細かく維持
Field[2].DistMax = 15.0 * Scale; // 15mmかけて粗くしていく
Field[2].LcMax = lcRef;         // 遠方は基本サイズ

Background Field = 2;

// =====================================================
// 9. Physical Groups
// =====================================================
Physical Volume("DOMAIN") = {finalAirVols[]};

volCoil1[] = {}; volCoil2[] = {}; volCoil3[] = {};
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

// Boundaries
eps = 1e-4; 
sXmin[] = Surface In BoundingBox {xMin-eps, yMin-eps, -zMargin-eps,     xMin+eps, yMax+eps, Lz+zMargin+eps};
sXmax[] = Surface In BoundingBox {xMax-eps, yMin-eps, -zMargin-eps,     xMax+eps, yMax+eps, Lz+zMargin+eps};
sYmin[] = Surface In BoundingBox {xMin-eps, yMin-eps, -zMargin-eps,     xMax+eps, yMin+eps, Lz+zMargin+eps};
sYmax[] = Surface In BoundingBox {xMin-eps, yMax-eps, -zMargin-eps,     xMax+eps, yMax+eps, Lz+zMargin+eps};
sZmin[] = Surface In BoundingBox {xMin-eps, yMin-eps, -zMargin-eps,     xMax+eps, yMax+eps, -zMargin+eps};
sZmax[] = Surface In BoundingBox {xMin-eps, yMin-eps, Lz+zMargin-eps,   xMax+eps, yMax+eps, Lz+zMargin+eps};

allOuterBnd[] = {sXmin[], sXmax[], sYmin[], sYmax[], sZmin[], sZmax[]};
Physical Surface("AIR_OUTER_WALLS") = {allOuterBnd[]};

// 界面の取得 (全空気境界 - 外壁)
allAirBnd[] = Boundary{ Volume{finalAirVols[]}; };
interfaceBnd[] = allAirBnd[];
interfaceBnd[] -= allOuterBnd[];
Physical Surface("INTERFACE_AIR_SOLID") = {interfaceBnd[]};