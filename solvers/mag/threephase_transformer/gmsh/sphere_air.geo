// TEAM11 cubed-sphere HEX: conductive shell (a->b) + outer air (b->R)
// Gmsh 4.12 / Built-in kernel friendly (no deprecated commands)
SetFactory("Built-in");

// ---------------- Parameters ----------------
a  = 0.05;         // inner radius [m]
b  = 0.055;        // outer radius of shell [m]
Ni = 10;           // face divisions (i direction)
Nj = 10;           // face divisions (j direction)
Nr = 16;           // radial divisions (shell: b->a)
gf = 0.98;          // shell radial growth factor (>=1; start at b)
lc = b/20;         // mesh size hint

// --- Outer air layer ---
Rf     = 5.0;      // R = Rf * b (recommend 4~5)
R      = Rf * b;   // far-field artificial boundary radius [m]
Nr_air = 24;       // radial divisions (air: b->R)
gf_air = 1.2;     // air radial growth factor (>=1; start at b)

// cubed-sphere scaling
sA = a / Sqrt(3);
sB = b / Sqrt(3);
sR = R / Sqrt(3);

// center
O = newp; Point(O) = {0,0,0, lc};

// ===========================================================
// ================ Face 0  ( +X face block ) ================
// ===========================================================
pOut0_0 = newp; Point(pOut0_0) = { sB*1, sB*1, sB*1, lc };
pOut0_1 = newp; Point(pOut0_1) = { sB*1, sB*-1, sB*1, lc };
pOut0_2 = newp; Point(pOut0_2) = { sB*1, sB*-1, sB*-1, lc };
pOut0_3 = newp; Point(pOut0_3) = { sB*1, sB*1,  sB*-1, lc };

pIn0_0  = newp; Point(pIn0_0)  = { sA*1, sA*1, sA*1, lc };
pIn0_1  = newp; Point(pIn0_1)  = { sA*1, sA*-1, sA*1, lc };
pIn0_2  = newp; Point(pIn0_2)  = { sA*1, sA*-1, sA*-1, lc };
pIn0_3  = newp; Point(pIn0_3)  = { sA*1, sA*1,  sA*-1, lc };

arcOut0_0 = newl; Circle(arcOut0_0) = { pOut0_0, O, pOut0_1 };
arcOut0_1 = newl; Circle(arcOut0_1) = { pOut0_1, O, pOut0_2 };
arcOut0_2 = newl; Circle(arcOut0_2) = { pOut0_2, O, pOut0_3 };
arcOut0_3 = newl; Circle(arcOut0_3) = { pOut0_3, O, pOut0_0 };

arcIn0_0  = newl; Circle(arcIn0_0)  = { pIn0_0,  O, pIn0_1  };
arcIn0_1  = newl; Circle(arcIn0_1)  = { pIn0_1,  O, pIn0_2  };
arcIn0_2  = newl; Circle(arcIn0_2)  = { pIn0_2,  O, pIn0_3  };
arcIn0_3  = newl; Circle(arcIn0_3)  = { pIn0_3,  O, pIn0_0  };

lRad0_0 = newl; Line(lRad0_0) = { pOut0_0, pIn0_0 };
lRad0_1 = newl; Line(lRad0_1) = { pOut0_1, pIn0_1 };
lRad0_2 = newl; Line(lRad0_2) = { pOut0_2, pIn0_2 };
lRad0_3 = newl; Line(lRad0_3) = { pOut0_3, pIn0_3 };

clOut0 = newll; Curve Loop(clOut0) = { arcOut0_0, arcOut0_1, arcOut0_2, arcOut0_3 };
sOut0 = news; Surface(sOut0) = { clOut0 };
clIn0 = newll; Curve Loop(clIn0) = { arcIn0_0, arcIn0_1, arcIn0_2, arcIn0_3 };
sIn0 = news; Surface(sIn0) = { clIn0 };

clSide0_0 = newll; Curve Loop(clSide0_0) = { arcOut0_0, lRad0_1, -arcIn0_0, -lRad0_0 };
sSide0_0 = news; Surface(sSide0_0) = { clSide0_0 };
clSide0_1 = newll; Curve Loop(clSide0_1) = { arcOut0_1, lRad0_2, -arcIn0_1, -lRad0_1 };
sSide0_1 = news; Surface(sSide0_1) = { clSide0_1 };
clSide0_2 = newll; Curve Loop(clSide0_2) = { arcOut0_2, lRad0_3, -arcIn0_2, -lRad0_2 };
sSide0_2 = news; Surface(sSide0_2) = { clSide0_2 };
clSide0_3 = newll; Curve Loop(clSide0_3) = { arcOut0_3, lRad0_0, -arcIn0_3, -lRad0_3 };
sSide0_3 = news; Surface(sSide0_3) = { clSide0_3 };

sl0 = newsl; Surface Loop(sl0) = { sOut0, sSide0_0, sSide0_1, sSide0_2, sSide0_3, -sIn0 };
v0  = newv; Volume(v0) = { sl0 };

Transfinite Line{ arcOut0_0, arcOut0_2, arcIn0_0, arcIn0_2 } = Ni + 1;
Transfinite Line{ arcOut0_1, arcOut0_3, arcIn0_1, arcIn0_3 } = Nj + 1;
Transfinite Line{ lRad0_0, lRad0_1, lRad0_2, lRad0_3 } = Nr + 1 Using Progression gf;
Transfinite Surface { sOut0, sIn0, sSide0_0, sSide0_1, sSide0_2, sSide0_3 };
Recombine Surface { sOut0, sIn0, sSide0_0, sSide0_1, sSide0_2, sSide0_3 };
Transfinite Volume { v0 } = { pOut0_0, pOut0_1, pOut0_2, pOut0_3, pIn0_0, pIn0_1, pIn0_2, pIn0_3 };
Recombine Volume { v0 };

// ---- Outer air for Face 0 (b->R) ----
pFar0_0 = newp; Point(pFar0_0) = { sR*1, sR*1, sR*1, lc };
pFar0_1 = newp; Point(pFar0_1) = { sR*1, sR*-1, sR*1, lc };
pFar0_2 = newp; Point(pFar0_2) = { sR*1, sR*-1, sR*-1, lc };
pFar0_3 = newp; Point(pFar0_3) = { sR*1, sR*1, sR*-1, lc };

arcFar0_0 = newl; Circle(arcFar0_0) = { pFar0_0, O, pFar0_1 };
arcFar0_1 = newl; Circle(arcFar0_1) = { pFar0_1, O, pFar0_2 };
arcFar0_2 = newl; Circle(arcFar0_2) = { pFar0_2, O, pFar0_3 };
arcFar0_3 = newl; Circle(arcFar0_3) = { pFar0_3, O, pFar0_0 };

lRadBR0_0 = newl; Line(lRadBR0_0) = { pOut0_0, pFar0_0 };
lRadBR0_1 = newl; Line(lRadBR0_1) = { pOut0_1, pFar0_1 };
lRadBR0_2 = newl; Line(lRadBR0_2) = { pOut0_2, pFar0_2 };
lRadBR0_3 = newl; Line(lRadBR0_3) = { pOut0_3, pFar0_3 };

clFar0 = newll; Curve Loop(clFar0) = { arcFar0_0, arcFar0_1, arcFar0_2, arcFar0_3 };
sFar0  = news;  Surface(sFar0)     = { clFar0 };

clSideAir0_0 = newll; Curve Loop(clSideAir0_0) = { arcFar0_0, -lRadBR0_1, -arcOut0_0, lRadBR0_0 };
sSideAir0_0  = news;  Surface(sSideAir0_0)     = { clSideAir0_0 };
clSideAir0_1 = newll; Curve Loop(clSideAir0_1) = { arcFar0_1, -lRadBR0_2, -arcOut0_1, lRadBR0_1 };
sSideAir0_1  = news;  Surface(sSideAir0_1)     = { clSideAir0_1 };
clSideAir0_2 = newll; Curve Loop(clSideAir0_2) = { arcFar0_2, -lRadBR0_3, -arcOut0_2, lRadBR0_2 };
sSideAir0_2  = news;  Surface(sSideAir0_2)     = { clSideAir0_2 };
clSideAir0_3 = newll; Curve Loop(clSideAir0_3) = { arcFar0_3, -lRadBR0_0, -arcOut0_3, lRadBR0_3 };
sSideAir0_3  = news;  Surface(sSideAir0_3)     = { clSideAir0_3 };

slOut0 = newsl; Surface Loop(slOut0) = { sFar0, sSideAir0_0, sSideAir0_1, sSideAir0_2, sSideAir0_3, -sOut0 };
vOut0  = newv;  Volume(vOut0)       = { slOut0 };

Transfinite Line{ arcFar0_0, arcFar0_2 } = Ni + 1;
Transfinite Line{ arcFar0_1, arcFar0_3 } = Nj + 1;
Transfinite Line{ lRadBR0_0, lRadBR0_1, lRadBR0_2, lRadBR0_3 } = Nr_air + 1 Using Progression gf_air;
Transfinite Surface { sFar0, sSideAir0_0, sSideAir0_1, sSideAir0_2, sSideAir0_3 };
Recombine  Surface { sFar0, sSideAir0_0, sSideAir0_1, sSideAir0_2, sSideAir0_3 };
Transfinite Volume  { vOut0 } = { pOut0_0, pOut0_1, pOut0_2, pOut0_3, pFar0_0, pFar0_1, pFar0_2, pFar0_3 };
Recombine  Volume   { vOut0 };

// ===========================================================
// ================ Face 1  ( -X face block ) ================
// ===========================================================
pOut1_0 = newp; Point(pOut1_0) = { sB*-1, sB*1,  sB*1,  lc };
pOut1_1 = newp; Point(pOut1_1) = { sB*-1, sB*1,  sB*-1, lc };
pOut1_2 = newp; Point(pOut1_2) = { sB*-1, sB*-1, sB*-1, lc };
pOut1_3 = newp; Point(pOut1_3) = { sB*-1, sB*-1, sB*1,  lc };

pIn1_0  = newp; Point(pIn1_0)  = { sA*-1, sA*1,  sA*1,  lc };
pIn1_1  = newp; Point(pIn1_1)  = { sA*-1, sA*1,  sA*-1, lc };
pIn1_2  = newp; Point(pIn1_2)  = { sA*-1, sA*-1, sA*-1, lc };
pIn1_3  = newp; Point(pIn1_3)  = { sA*-1, sA*-1, sA*1,  lc };

arcOut1_0 = newl; Circle(arcOut1_0) = { pOut1_0, O, pOut1_1 };
arcOut1_1 = newl; Circle(arcOut1_1) = { pOut1_1, O, pOut1_2 };
arcOut1_2 = newl; Circle(arcOut1_2) = { pOut1_2, O, pOut1_3 };
arcOut1_3 = newl; Circle(arcOut1_3) = { pOut1_3, O, pOut1_0 };

arcIn1_0  = newl; Circle(arcIn1_0)  = { pIn1_0,  O, pIn1_1  };
arcIn1_1  = newl; Circle(arcIn1_1)  = { pIn1_1,  O, pIn1_2  };
arcIn1_2  = newl; Circle(arcIn1_2)  = { pIn1_2,  O, pIn1_3  };
arcIn1_3  = newl; Circle(arcIn1_3)  = { pIn1_3,  O, pIn1_0  };

lRad1_0 = newl; Line(lRad1_0) = { pOut1_0, pIn1_0 };
lRad1_1 = newl; Line(lRad1_1) = { pOut1_1, pIn1_1 };
lRad1_2 = newl; Line(lRad1_2) = { pOut1_2, pIn1_2 };
lRad1_3 = newl; Line(lRad1_3) = { pOut1_3, pIn1_3 };

clOut1 = newll; Curve Loop(clOut1) = { arcOut1_0, arcOut1_1, arcOut1_2, arcOut1_3 };
sOut1 = news; Surface(sOut1) = { clOut1 };
clIn1 = newll; Curve Loop(clIn1) = { arcIn1_0, arcIn1_1, arcIn1_2, arcIn1_3 };
sIn1 = news; Surface(sIn1) = { clIn1 };

clSide1_0 = newll; Curve Loop(clSide1_0) = { arcOut1_0, lRad1_1, -arcIn1_0, -lRad1_0 };
sSide1_0 = news; Surface(sSide1_0) = { clSide1_0 };
clSide1_1 = newll; Curve Loop(clSide1_1) = { arcOut1_1, lRad1_2, -arcIn1_1, -lRad1_1 };
sSide1_1 = news; Surface(sSide1_1) = { clSide1_1 };
clSide1_2 = newll; Curve Loop(clSide1_2) = { arcOut1_2, lRad1_3, -arcIn1_2, -lRad1_2 };
sSide1_2 = news; Surface(sSide1_2) = { clSide1_2 };
clSide1_3 = newll; Curve Loop(clSide1_3) = { arcOut1_3, lRad1_0, -arcIn1_3, -lRad1_3 };
sSide1_3 = news; Surface(sSide1_3) = { clSide1_3 };

sl1 = newsl; Surface Loop(sl1) = { sOut1, sSide1_0, sSide1_1, sSide1_2, sSide1_3, -sIn1 };
v1  = newv; Volume(v1) = { sl1 };

Transfinite Line{ arcOut1_0, arcOut1_2, arcIn1_0, arcIn1_2 } = Ni + 1;
Transfinite Line{ arcOut1_1, arcOut1_3, arcIn1_1, arcIn1_3 } = Nj + 1;
Transfinite Line{ lRad1_0, lRad1_1, lRad1_2, lRad1_3 } = Nr + 1 Using Progression gf;
Transfinite Surface { sOut1, sIn1, sSide1_0, sSide1_1, sSide1_2, sSide1_3 };
Recombine Surface { sOut1, sIn1, sSide1_0, sSide1_1, sSide1_2, sSide1_3 };
Transfinite Volume { v1 } = { pOut1_0, pOut1_1, pOut1_2, pOut1_3, pIn1_0, pIn1_1, pIn1_2, pIn1_3 };
Recombine Volume { v1 };

// ---- Outer air for Face 1 ----
pFar1_0 = newp; Point(pFar1_0) = { sR*-1, sR*1,  sR*1,  lc };
pFar1_1 = newp; Point(pFar1_1) = { sR*-1, sR*1,  sR*-1, lc };
pFar1_2 = newp; Point(pFar1_2) = { sR*-1, sR*-1, sR*-1, lc };
pFar1_3 = newp; Point(pFar1_3) = { sR*-1, sR*-1, sR*1,  lc };

arcFar1_0 = newl; Circle(arcFar1_0) = { pFar1_0, O, pFar1_1 };
arcFar1_1 = newl; Circle(arcFar1_1) = { pFar1_1, O, pFar1_2 };
arcFar1_2 = newl; Circle(arcFar1_2) = { pFar1_2, O, pFar1_3 };
arcFar1_3 = newl; Circle(arcFar1_3) = { pFar1_3, O, pFar1_0 };

lRadBR1_0 = newl; Line(lRadBR1_0) = { pOut1_0, pFar1_0 };
lRadBR1_1 = newl; Line(lRadBR1_1) = { pOut1_1, pFar1_1 };
lRadBR1_2 = newl; Line(lRadBR1_2) = { pOut1_2, pFar1_2 };
lRadBR1_3 = newl; Line(lRadBR1_3) = { pOut1_3, pFar1_3 };

clFar1 = newll; Curve Loop(clFar1) = { arcFar1_0, arcFar1_1, arcFar1_2, arcFar1_3 };
sFar1  = news;  Surface(sFar1)     = { clFar1 };

clSideAir1_0 = newll; Curve Loop(clSideAir1_0) = { arcFar1_0, -lRadBR1_1, -arcOut1_0, lRadBR1_0 };
sSideAir1_0  = news;  Surface(sSideAir1_0)     = { clSideAir1_0 };
clSideAir1_1 = newll; Curve Loop(clSideAir1_1) = { arcFar1_1, -lRadBR1_2, -arcOut1_1, lRadBR1_1 };
sSideAir1_1  = news;  Surface(sSideAir1_1)     = { clSideAir1_1 };
clSideAir1_2 = newll; Curve Loop(clSideAir1_2) = { arcFar1_2, -lRadBR1_3, -arcOut1_2, lRadBR1_2 };
sSideAir1_2  = news;  Surface(sSideAir1_2)     = { clSideAir1_2 };
clSideAir1_3 = newll; Curve Loop(clSideAir1_3) = { arcFar1_3, -lRadBR1_0, -arcOut1_3, lRadBR1_3 };
sSideAir1_3  = news;  Surface(sSideAir1_3)     = { clSideAir1_3 };

slOut1 = newsl; Surface Loop(slOut1) = { sFar1, sSideAir1_0, sSideAir1_1, sSideAir1_2, sSideAir1_3, -sOut1 };
vOut1  = newv;  Volume(vOut1)       = { slOut1 };

Transfinite Line{ arcFar1_0, arcFar1_2 } = Ni + 1;
Transfinite Line{ arcFar1_1, arcFar1_3 } = Nj + 1;
Transfinite Line{ lRadBR1_0, lRadBR1_1, lRadBR1_2, lRadBR1_3 } = Nr_air + 1 Using Progression gf_air;
Transfinite Surface { sFar1, sSideAir1_0, sSideAir1_1, sSideAir1_2, sSideAir1_3 };
Recombine  Surface { sFar1, sSideAir1_0, sSideAir1_1, sSideAir1_2, sSideAir1_3 };
Transfinite Volume  { vOut1 } = { pOut1_0, pOut1_1, pOut1_2, pOut1_3, pFar1_0, pFar1_1, pFar1_2, pFar1_3 };
Recombine  Volume   { vOut1 };

// ===========================================================
// ================ Face 2  ( +Y face block ) ================
// ===========================================================
pOut2_0 = newp; Point(pOut2_0) = { sB*1,  sB*1, sB*1,  lc };
pOut2_1 = newp; Point(pOut2_1) = { sB*-1, sB*1, sB*1,  lc };
pOut2_2 = newp; Point(pOut2_2) = { sB*-1, sB*1, sB*-1, lc };
pOut2_3 = newp; Point(pOut2_3) = { sB*1,  sB*1, sB*-1, lc };

pIn2_0  = newp; Point(pIn2_0)  = { sA*1,  sA*1, sA*1,  lc };
pIn2_1  = newp; Point(pIn2_1)  = { sA*-1, sA*1, sA*1,  lc };
pIn2_2  = newp; Point(pIn2_2)  = { sA*-1, sA*1, sA*-1, lc };
pIn2_3  = newp; Point(pIn2_3)  = { sA*1,  sA*1, sA*-1, lc };

arcOut2_0 = newl; Circle(arcOut2_0) = { pOut2_0, O, pOut2_1 };
arcOut2_1 = newl; Circle(arcOut2_1) = { pOut2_1, O, pOut2_2 };
arcOut2_2 = newl; Circle(arcOut2_2) = { pOut2_2, O, pOut2_3 };
arcOut2_3 = newl; Circle(arcOut2_3) = { pOut2_3, O, pOut2_0 };

arcIn2_0  = newl; Circle(arcIn2_0)  = { pIn2_0,  O, pIn2_1  };
arcIn2_1  = newl; Circle(arcIn2_1)  = { pIn2_1,  O, pIn2_2  };
arcIn2_2  = newl; Circle(arcIn2_2)  = { pIn2_2,  O, pIn2_3  };
arcIn2_3  = newl; Circle(arcIn2_3)  = { pIn2_3,  O, pIn2_0  };

lRad2_0 = newl; Line(lRad2_0) = { pOut2_0, pIn2_0 };
lRad2_1 = newl; Line(lRad2_1) = { pOut2_1, pIn2_1 };
lRad2_2 = newl; Line(lRad2_2) = { pOut2_2, pIn2_2 };
lRad2_3 = newl; Line(lRad2_3) = { pOut2_3, pIn2_3 };

clOut2 = newll; Curve Loop(clOut2) = { arcOut2_0, arcOut2_1, arcOut2_2, arcOut2_3 };
sOut2 = news; Surface(sOut2) = { clOut2 };
clIn2 = newll; Curve Loop(clIn2) = { arcIn2_0, arcIn2_1, arcIn2_2, arcIn2_3 };
sIn2 = news; Surface(sIn2) = { clIn2 };

clSide2_0 = newll; Curve Loop(clSide2_0) = { arcOut2_0, lRad2_1, -arcIn2_0, -lRad2_0 };
sSide2_0 = news; Surface(sSide2_0) = { clSide2_0 };
clSide2_1 = newll; Curve Loop(clSide2_1) = { arcOut2_1, lRad2_2, -arcIn2_1, -lRad2_1 };
sSide2_1 = news; Surface(sSide2_1) = { clSide2_1 };
clSide2_2 = newll; Curve Loop(clSide2_2) = { arcOut2_2, lRad2_3, -arcIn2_2, -lRad2_2 };
sSide2_2 = news; Surface(sSide2_2) = { clSide2_2 };
clSide2_3 = newll; Curve Loop(clSide2_3) = { arcOut2_3, lRad2_0, -arcIn2_3, -lRad2_3 };
sSide2_3 = news; Surface(sSide2_3) = { clSide2_3 };

sl2 = newsl; Surface Loop(sl2) = { sOut2, sSide2_0, sSide2_1, sSide2_2, sSide2_3, -sIn2 };
v2  = newv; Volume(v2) = { sl2 };

Transfinite Line{ arcOut2_0, arcOut2_2, arcIn2_0, arcIn2_2 } = Ni + 1;
Transfinite Line{ arcOut2_1, arcOut2_3, arcIn2_1, arcIn2_3 } = Nj + 1;
Transfinite Line{ lRad2_0, lRad2_1, lRad2_2, lRad2_3 } = Nr + 1 Using Progression gf;
Transfinite Surface { sOut2, sIn2, sSide2_0, sSide2_1, sSide2_2, sSide2_3 };
Recombine Surface { sOut2, sIn2, sSide2_0, sSide2_1, sSide2_2, sSide2_3 };
Transfinite Volume { v2 } = { pOut2_0, pOut2_1, pOut2_2, pOut2_3, pIn2_0, pIn2_1, pIn2_2, pIn2_3 };
Recombine Volume { v2 };

// ---- Outer air for Face 2 ----
pFar2_0 = newp; Point(pFar2_0) = { sR*1,  sR*1, sR*1,  lc };
pFar2_1 = newp; Point(pFar2_1) = { sR*-1, sR*1, sR*1,  lc };
pFar2_2 = newp; Point(pFar2_2) = { sR*-1, sR*1, sR*-1, lc };
pFar2_3 = newp; Point(pFar2_3) = { sR*1,  sR*1, sR*-1, lc };

arcFar2_0 = newl; Circle(arcFar2_0) = { pFar2_0, O, pFar2_1 };
arcFar2_1 = newl; Circle(arcFar2_1) = { pFar2_1, O, pFar2_2 };
arcFar2_2 = newl; Circle(arcFar2_2) = { pFar2_2, O, pFar2_3 };
arcFar2_3 = newl; Circle(arcFar2_3) = { pFar2_3, O, pFar2_0 };

lRadBR2_0 = newl; Line(lRadBR2_0) = { pOut2_0, pFar2_0 };
lRadBR2_1 = newl; Line(lRadBR2_1) = { pOut2_1, pFar2_1 };
lRadBR2_2 = newl; Line(lRadBR2_2) = { pOut2_2, pFar2_2 };
lRadBR2_3 = newl; Line(lRadBR2_3) = { pOut2_3, pFar2_3 };

clFar2 = newll; Curve Loop(clFar2) = { arcFar2_0, arcFar2_1, arcFar2_2, arcFar2_3 };
sFar2  = news;  Surface(sFar2)     = { clFar2 };

clSideAir2_0 = newll; Curve Loop(clSideAir2_0) = { arcFar2_0, -lRadBR2_1, -arcOut2_0, lRadBR2_0 };
sSideAir2_0  = news;  Surface(sSideAir2_0)     = { clSideAir2_0 };
clSideAir2_1 = newll; Curve Loop(clSideAir2_1) = { arcFar2_1, -lRadBR2_2, -arcOut2_1, lRadBR2_1 };
sSideAir2_1  = news;  Surface(sSideAir2_1)     = { clSideAir2_1 };
clSideAir2_2 = newll; Curve Loop(clSideAir2_2) = { arcFar2_2, -lRadBR2_3, -arcOut2_2, lRadBR2_2 };
sSideAir2_2  = news;  Surface(sSideAir2_2)     = { clSideAir2_2 };
clSideAir2_3 = newll; Curve Loop(clSideAir2_3) = { arcFar2_3, -lRadBR2_0, -arcOut2_3, lRadBR2_3 };
sSideAir2_3  = news;  Surface(sSideAir2_3)     = { clSideAir2_3 };

slOut2 = newsl; Surface Loop(slOut2) = { sFar2, sSideAir2_0, sSideAir2_1, sSideAir2_2, sSideAir2_3, -sOut2 };
vOut2  = newv;  Volume(vOut2)       = { slOut2 };

Transfinite Line{ arcFar2_0, arcFar2_2 } = Ni + 1;
Transfinite Line{ arcFar2_1, arcFar2_3 } = Nj + 1;
Transfinite Line{ lRadBR2_0, lRadBR2_1, lRadBR2_2, lRadBR2_3 } = Nr_air + 1 Using Progression gf_air;
Transfinite Surface { sFar2, sSideAir2_0, sSideAir2_1, sSideAir2_2, sSideAir2_3 };
Recombine  Surface { sFar2, sSideAir2_0, sSideAir2_1, sSideAir2_2, sSideAir2_3 };
Transfinite Volume  { vOut2 } = { pOut2_0, pOut2_1, pOut2_2, pOut2_3, pFar2_0, pFar2_1, pFar2_2, pFar2_3 };
Recombine  Volume   { vOut2 };

// ===========================================================
// ================ Face 3  ( -Y face block ) ================
// ===========================================================
pOut3_0 = newp; Point(pOut3_0) = { sB*1,  sB*-1, sB*1,  lc };
pOut3_1 = newp; Point(pOut3_1) = { sB*1,  sB*-1, sB*-1, lc };
pOut3_2 = newp; Point(pOut3_2) = { sB*-1, sB*-1, sB*-1, lc };
pOut3_3 = newp; Point(pOut3_3) = { sB*-1, sB*-1, sB*1,  lc };

pIn3_0  = newp; Point(pIn3_0)  = { sA*1,  sA*-1, sA*1,  lc };
pIn3_1  = newp; Point(pIn3_1)  = { sA*1,  sA*-1, sA*-1, lc };
pIn3_2  = newp; Point(pIn3_2)  = { sA*-1, sA*-1, sA*-1, lc };
pIn3_3  = newp; Point(pIn3_3)  = { sA*-1, sA*-1, sA*1,  lc };

arcOut3_0 = newl; Circle(arcOut3_0) = { pOut3_0, O, pOut3_1 };
arcOut3_1 = newl; Circle(arcOut3_1) = { pOut3_1, O, pOut3_2 };
arcOut3_2 = newl; Circle(arcOut3_2) = { pOut3_2, O, pOut3_3 };
arcOut3_3 = newl; Circle(arcOut3_3) = { pOut3_3, O, pOut3_0 };

arcIn3_0  = newl; Circle(arcIn3_0)  = { pIn3_0,  O, pIn3_1  };
arcIn3_1  = newl; Circle(arcIn3_1)  = { pIn3_1,  O, pIn3_2  };
arcIn3_2  = newl; Circle(arcIn3_2)  = { pIn3_2,  O, pIn3_3  };
arcIn3_3  = newl; Circle(arcIn3_3)  = { pIn3_3,  O, pIn3_0  };

lRad3_0 = newl; Line(lRad3_0) = { pOut3_0, pIn3_0 };
lRad3_1 = newl; Line(lRad3_1) = { pOut3_1, pIn3_1 };
lRad3_2 = newl; Line(lRad3_2) = { pOut3_2, pIn3_2 };
lRad3_3 = newl; Line(lRad3_3) = { pOut3_3, pIn3_3 };

clOut3 = newll; Curve Loop(clOut3) = { arcOut3_0, arcOut3_1, arcOut3_2, arcOut3_3 };
sOut3 = news; Surface(sOut3) = { clOut3 };
clIn3 = newll; Curve Loop(clIn3) = { arcIn3_0, arcIn3_1, arcIn3_2, arcIn3_3 };
sIn3 = news; Surface(sIn3) = { clIn3 };

clSide3_0 = newll; Curve Loop(clSide3_0) = { arcOut3_0, lRad3_1, -arcIn3_0, -lRad3_0 };
sSide3_0 = news; Surface(sSide3_0) = { clSide3_0 };
clSide3_1 = newll; Curve Loop(clSide3_1) = { arcOut3_1, lRad3_2, -arcIn3_1, -lRad3_1 };
sSide3_1 = news; Surface(sSide3_1) = { clSide3_1 };
clSide3_2 = newll; Curve Loop(clSide3_2) = { arcOut3_2, lRad3_3, -arcIn3_2, -lRad3_2 };
sSide3_2 = news; Surface(sSide3_2) = { clSide3_2 };
clSide3_3 = newll; Curve Loop(clSide3_3) = { arcOut3_3, lRad3_0, -arcIn3_3, -lRad3_3 };
sSide3_3 = news; Surface(sSide3_3) = { clSide3_3 };

sl3 = newsl; Surface Loop(sl3) = { sOut3, sSide3_0, sSide3_1, sSide3_2, sSide3_3, -sIn3 };
v3  = newv; Volume(v3) = { sl3 };

Transfinite Line{ arcOut3_0, arcOut3_2, arcIn3_0, arcIn3_2 } = Ni + 1;
Transfinite Line{ arcOut3_1, arcOut3_3, arcIn3_1, arcIn3_3 } = Nj + 1;
Transfinite Line{ lRad3_0, lRad3_1, lRad3_2, lRad3_3 } = Nr + 1 Using Progression gf;
Transfinite Surface { sOut3, sIn3, sSide3_0, sSide3_1, sSide3_2, sSide3_3 };
Recombine Surface { sOut3, sIn3, sSide3_0, sSide3_1, sSide3_2, sSide3_3 };
Transfinite Volume { v3 } = { pOut3_0, pOut3_1, pOut3_2, pOut3_3, pIn3_0, pIn3_1, pIn3_2, pIn3_3 };
Recombine Volume { v3 };

// ---- Outer air for Face 3 ----
pFar3_0 = newp; Point(pFar3_0) = { sR*1,  sR*-1, sR*1,  lc };
pFar3_1 = newp; Point(pFar3_1) = { sR*1,  sR*-1, sR*-1, lc };
pFar3_2 = newp; Point(pFar3_2) = { sR*-1, sR*-1, sR*-1, lc };
pFar3_3 = newp; Point(pFar3_3) = { sR*-1, sR*-1, sR*1,  lc };

arcFar3_0 = newl; Circle(arcFar3_0) = { pFar3_0, O, pFar3_1 };
arcFar3_1 = newl; Circle(arcFar3_1) = { pFar3_1, O, pFar3_2 };
arcFar3_2 = newl; Circle(arcFar3_2) = { pFar3_2, O, pFar3_3 };
arcFar3_3 = newl; Circle(arcFar3_3) = { pFar3_3, O, pFar3_0 };

lRadBR3_0 = newl; Line(lRadBR3_0) = { pOut3_0, pFar3_0 };
lRadBR3_1 = newl; Line(lRadBR3_1) = { pOut3_1, pFar3_1 };
lRadBR3_2 = newl; Line(lRadBR3_2) = { pOut3_2, pFar3_2 };
lRadBR3_3 = newl; Line(lRadBR3_3) = { pOut3_3, pFar3_3 };

clFar3 = newll; Curve Loop(clFar3) = { arcFar3_0, arcFar3_1, arcFar3_2, arcFar3_3 };
sFar3  = news;  Surface(sFar3)     = { clFar3 };

clSideAir3_0 = newll; Curve Loop(clSideAir3_0) = { arcFar3_0, -lRadBR3_1, -arcOut3_0, lRadBR3_0 };
sSideAir3_0  = news;  Surface(sSideAir3_0)     = { clSideAir3_0 };
clSideAir3_1 = newll; Curve Loop(clSideAir3_1) = { arcFar3_1, -lRadBR3_2, -arcOut3_1, lRadBR3_1 };
sSideAir3_1  = news;  Surface(sSideAir3_1)     = { clSideAir3_1 };
clSideAir3_2 = newll; Curve Loop(clSideAir3_2) = { arcFar3_2, -lRadBR3_3, -arcOut3_2, lRadBR3_2 };
sSideAir3_2  = news;  Surface(sSideAir3_2)     = { clSideAir3_2 };
clSideAir3_3 = newll; Curve Loop(clSideAir3_3) = { arcFar3_3, -lRadBR3_0, -arcOut3_3, lRadBR3_3 };
sSideAir3_3  = news;  Surface(sSideAir3_3)     = { clSideAir3_3 };

slOut3 = newsl; Surface Loop(slOut3) = { sFar3, sSideAir3_0, sSideAir3_1, sSideAir3_2, sSideAir3_3, -sOut3 };
vOut3  = newv;  Volume(vOut3)       = { slOut3 };

Transfinite Line{ arcFar3_0, arcFar3_2 } = Ni + 1;
Transfinite Line{ arcFar3_1, arcFar3_3 } = Nj + 1;
Transfinite Line{ lRadBR3_0, lRadBR3_1, lRadBR3_2, lRadBR3_3 } = Nr_air + 1 Using Progression gf_air;
Transfinite Surface { sFar3, sSideAir3_0, sSideAir3_1, sSideAir3_2, sSideAir3_3 };
Recombine  Surface { sFar3, sSideAir3_0, sSideAir3_1, sSideAir3_2, sSideAir3_3 };
Transfinite Volume  { vOut3 } = { pOut3_0, pOut3_1, pOut3_2, pOut3_3, pFar3_0, pFar3_1, pFar3_2, pFar3_3 };
Recombine  Volume   { vOut3 };

// ===========================================================
// ================ Face 4  ( +Z face block ) ================
// ===========================================================
pOut4_0 = newp; Point(pOut4_0) = { sB*1,  sB*1,  sB*1, lc };
pOut4_1 = newp; Point(pOut4_1) = { sB*1,  sB*-1, sB*1, lc };
pOut4_2 = newp; Point(pOut4_2) = { sB*-1, sB*-1, sB*1, lc };
pOut4_3 = newp; Point(pOut4_3) = { sB*-1, sB*1,  sB*1, lc };

pIn4_0  = newp; Point(pIn4_0)  = { sA*1,  sA*1,  sA*1, lc };
pIn4_1  = newp; Point(pIn4_1)  = { sA*1,  sA*-1, sA*1, lc };
pIn4_2  = newp; Point(pIn4_2)  = { sA*-1, sA*-1, sA*1, lc };
pIn4_3  = newp; Point(pIn4_3)  = { sA*-1, sA*1,  sA*1, lc };

arcOut4_0 = newl; Circle(arcOut4_0) = { pOut4_0, O, pOut4_1 };
arcOut4_1 = newl; Circle(arcOut4_1) = { pOut4_1, O, pOut4_2 };
arcOut4_2 = newl; Circle(arcOut4_2) = { pOut4_2, O, pOut4_3 };
arcOut4_3 = newl; Circle(arcOut4_3) = { pOut4_3, O, pOut4_0 };

arcIn4_0  = newl; Circle(arcIn4_0)  = { pIn4_0,  O, pIn4_1  };
arcIn4_1  = newl; Circle(arcIn4_1)  = { pIn4_1,  O, pIn4_2  };
arcIn4_2  = newl; Circle(arcIn4_2)  = { pIn4_2,  O, pIn4_3  };
arcIn4_3  = newl; Circle(arcIn4_3)  = { pIn4_3,  O, pIn4_0  };

lRad4_0 = newl; Line(lRad4_0) = { pOut4_0, pIn4_0 };
lRad4_1 = newl; Line(lRad4_1) = { pOut4_1, pIn4_1 };
lRad4_2 = newl; Line(lRad4_2) = { pOut4_2, pIn4_2 };
lRad4_3 = newl; Line(lRad4_3) = { pOut4_3, pIn4_3 };

clOut4 = newll; Curve Loop(clOut4) = { arcOut4_0, arcOut4_1, arcOut4_2, arcOut4_3 };
sOut4 = news; Surface(sOut4) = { clOut4 };
clIn4 = newll; Curve Loop(clIn4) = { arcIn4_0, arcIn4_1, arcIn4_2, arcIn4_3 };
sIn4 = news; Surface(sIn4) = { clIn4 };

clSide4_0 = newll; Curve Loop(clSide4_0) = { arcOut4_0, lRad4_1, -arcIn4_0, -lRad4_0 };
sSide4_0 = news; Surface(sSide4_0) = { clSide4_0 };
clSide4_1 = newll; Curve Loop(clSide4_1) = { arcOut4_1, lRad4_2, -arcIn4_1, -lRad4_1 };
sSide4_1 = news; Surface(sSide4_1) = { clSide4_1 };
clSide4_2 = newll; Curve Loop(clSide4_2) = { arcOut4_2, lRad4_3, -arcIn4_2, -lRad4_2 };
sSide4_2 = news; Surface(sSide4_2) = { clSide4_2 };
clSide4_3 = newll; Curve Loop(clSide4_3) = { arcOut4_3, lRad4_0, -arcIn4_3, -lRad4_3 };
sSide4_3 = news; Surface(sSide4_3) = { clSide4_3 };

sl4 = newsl; Surface Loop(sl4) = { sOut4, sSide4_0, sSide4_1, sSide4_2, sSide4_3, -sIn4 };
v4  = newv; Volume(v4) = { sl4 };

Transfinite Line{ arcOut4_0, arcOut4_2, arcIn4_0, arcIn4_2 } = Ni + 1;
Transfinite Line{ arcOut4_1, arcOut4_3, arcIn4_1, arcIn4_3 } = Nj + 1;
Transfinite Line{ lRad4_0, lRad4_1, lRad4_2, lRad4_3 } = Nr + 1 Using Progression gf;
Transfinite Surface { sOut4, sIn4, sSide4_0, sSide4_1, sSide4_2, sSide4_3 };
Recombine Surface { sOut4, sIn4, sSide4_0, sSide4_1, sSide4_2, sSide4_3 };
Transfinite Volume { v4 } = { pOut4_0, pOut4_1, pOut4_2, pOut4_3, pIn4_0, pIn4_1, pIn4_2, pIn4_3 };
Recombine Volume { v4 };

// ---- Outer air for Face 4 ----
pFar4_0 = newp; Point(pFar4_0) = { sR*1,  sR*1,  sR*1, lc };
pFar4_1 = newp; Point(pFar4_1) = { sR*1,  sR*-1, sR*1, lc };
pFar4_2 = newp; Point(pFar4_2) = { sR*-1, sR*-1, sR*1, lc };
pFar4_3 = newp; Point(pFar4_3) = { sR*-1, sR*1,  sR*1, lc };

arcFar4_0 = newl; Circle(arcFar4_0) = { pFar4_0, O, pFar4_1 };
arcFar4_1 = newl; Circle(arcFar4_1) = { pFar4_1, O, pFar4_2 };
arcFar4_2 = newl; Circle(arcFar4_2) = { pFar4_2, O, pFar4_3 };
arcFar4_3 = newl; Circle(arcFar4_3) = { pFar4_3, O, pFar4_0 };

lRadBR4_0 = newl; Line(lRadBR4_0) = { pOut4_0, pFar4_0 };
lRadBR4_1 = newl; Line(lRadBR4_1) = { pOut4_1, pFar4_1 };
lRadBR4_2 = newl; Line(lRadBR4_2) = { pOut4_2, pFar4_2 };
lRadBR4_3 = newl; Line(lRadBR4_3) = { pOut4_3, pFar4_3 };

clFar4 = newll; Curve Loop(clFar4) = { arcFar4_0, arcFar4_1, arcFar4_2, arcFar4_3 };
sFar4  = news;  Surface(sFar4)     = { clFar4 };

clSideAir4_0 = newll; Curve Loop(clSideAir4_0) = { arcFar4_0, -lRadBR4_1, -arcOut4_0, lRadBR4_0 };
sSideAir4_0  = news;  Surface(sSideAir4_0)     = { clSideAir4_0 };
clSideAir4_1 = newll; Curve Loop(clSideAir4_1) = { arcFar4_1, -lRadBR4_2, -arcOut4_1, lRadBR4_1 };
sSideAir4_1  = news;  Surface(sSideAir4_1)     = { clSideAir4_1 };
clSideAir4_2 = newll; Curve Loop(clSideAir4_2) = { arcFar4_2, -lRadBR4_3, -arcOut4_2, lRadBR4_2 };
sSideAir4_2  = news;  Surface(sSideAir4_2)     = { clSideAir4_2 };
clSideAir4_3 = newll; Curve Loop(clSideAir4_3) = { arcFar4_3, -lRadBR4_0, -arcOut4_3, lRadBR4_3 };
sSideAir4_3  = news;  Surface(sSideAir4_3)     = { clSideAir4_3 };

slOut4 = newsl; Surface Loop(slOut4) = { sFar4, sSideAir4_0, sSideAir4_1, sSideAir4_2, sSideAir4_3, -sOut4 };
vOut4  = newv;  Volume(vOut4)       = { slOut4 };

Transfinite Line{ arcFar4_0, arcFar4_2 } = Ni + 1;
Transfinite Line{ arcFar4_1, arcFar4_3 } = Nj + 1;
Transfinite Line{ lRadBR4_0, lRadBR4_1, lRadBR4_2, lRadBR4_3 } = Nr_air + 1 Using Progression gf_air;
Transfinite Surface { sFar4, sSideAir4_0, sSideAir4_1, sSideAir4_2, sSideAir4_3 };
Recombine  Surface { sFar4, sSideAir4_0, sSideAir4_1, sSideAir4_2, sSideAir4_3 };
Transfinite Volume  { vOut4 } = { pOut4_0, pOut4_1, pOut4_2, pOut4_3, pFar4_0, pFar4_1, pFar4_2, pFar4_3 };
Recombine  Volume   { vOut4 };

// ===========================================================
// ================ Face 5  ( -Z face block ) ================
// ===========================================================
pOut5_0 = newp; Point(pOut5_0) = { sB*1,  sB*1,  sB*-1, lc };
pOut5_1 = newp; Point(pOut5_1) = { sB*-1, sB*1,  sB*-1, lc };
pOut5_2 = newp; Point(pOut5_2) = { sB*-1, sB*-1, sB*-1, lc };
pOut5_3 = newp; Point(pOut5_3) = { sB*1,  sB*-1, sB*-1, lc };

pIn5_0  = newp; Point(pIn5_0)  = { sA*1,  sA*1,  sA*-1, lc };
pIn5_1  = newp; Point(pIn5_1)  = { sA*-1, sA*1,  sA*-1, lc };
pIn5_2  = newp; Point(pIn5_2)  = { sA*-1, sA*-1, sA*-1, lc };
pIn5_3  = newp; Point(pIn5_3)  = { sA*1,  sA*-1, sA*-1, lc };

arcOut5_0 = newl; Circle(arcOut5_0) = { pOut5_0, O, pOut5_1 };
arcOut5_1 = newl; Circle(arcOut5_1) = { pOut5_1, O, pOut5_2 };
arcOut5_2 = newl; Circle(arcOut5_2) = { pOut5_2, O, pOut5_3 };
arcOut5_3 = newl; Circle(arcOut5_3) = { pOut5_3, O, pOut5_0 };

arcIn5_0  = newl; Circle(arcIn5_0)  = { pIn5_0,  O, pIn5_1  };
arcIn5_1  = newl; Circle(arcIn5_1)  = { pIn5_1,  O, pIn5_2  };
arcIn5_2  = newl; Circle(arcIn5_2)  = { pIn5_2,  O, pIn5_3  };
arcIn5_3  = newl; Circle(arcIn5_3)  = { pIn5_3,  O, pIn5_0  };

lRad5_0 = newl; Line(lRad5_0) = { pOut5_0, pIn5_0 };
lRad5_1 = newl; Line(lRad5_1) = { pOut5_1, pIn5_1 };
lRad5_2 = newl; Line(lRad5_2) = { pOut5_2, pIn5_2 };
lRad5_3 = newl; Line(lRad5_3) = { pOut5_3, pIn5_3 };

clOut5 = newll; Curve Loop(clOut5) = { arcOut5_0, arcOut5_1, arcOut5_2, arcOut5_3 };
sOut5 = news; Surface(sOut5) = { clOut5 };
clIn5 = newll; Curve Loop(clIn5) = { arcIn5_0, arcIn5_1, arcIn5_2, arcIn5_3 };
sIn5 = news; Surface(sIn5) = { clIn5 };

clSide5_0 = newll; Curve Loop(clSide5_0) = { arcOut5_0, lRad5_1, -arcIn5_0, -lRad5_0 };
sSide5_0 = news; Surface(sSide5_0) = { clSide5_0 };
clSide5_1 = newll; Curve Loop(clSide5_1) = { arcOut5_1, lRad5_2, -arcIn5_1, -lRad5_1 };
sSide5_1 = news; Surface(sSide5_1) = { clSide5_1 };
clSide5_2 = newll; Curve Loop(clSide5_2) = { arcOut5_2, lRad5_3, -arcIn5_2, -lRad5_2 };
sSide5_2 = news; Surface(sSide5_2) = { clSide5_2 };
clSide5_3 = newll; Curve Loop(clSide5_3) = { arcOut5_3, lRad5_0, -arcIn5_3, -lRad5_3 };
sSide5_3 = news; Surface(sSide5_3) = { clSide5_3 };

sl5 = newsl; Surface Loop(sl5) = { sOut5, sSide5_0, sSide5_1, sSide5_2, sSide5_3, -sIn5 };
v5  = newv; Volume(v5) = { sl5 };

Transfinite Line{ arcOut5_0, arcOut5_2, arcIn5_0, arcIn5_2 } = Ni + 1;
Transfinite Line{ arcOut5_1, arcOut5_3, arcIn5_1, arcIn5_3 } = Nj + 1;
Transfinite Line{ lRad5_0, lRad5_1, lRad5_2, lRad5_3 } = Nr + 1 Using Progression gf;
Transfinite Surface { sOut5, sIn5, sSide5_0, sSide5_1, sSide5_2, sSide5_3 };
Recombine Surface { sOut5, sIn5, sSide5_0, sSide5_1, sSide5_2, sSide5_3 };
Transfinite Volume { v5 } = { pOut5_0, pOut5_1, pOut5_2, pOut5_3, pIn5_0, pIn5_1, pIn5_2, pIn5_3 };
Recombine Volume { v5 };

// ---- Outer air for Face 5 ----
pFar5_0 = newp; Point(pFar5_0) = { sR*1,  sR*1,  sR*-1, lc };
pFar5_1 = newp; Point(pFar5_1) = { sR*-1, sR*1,  sR*-1, lc };
pFar5_2 = newp; Point(pFar5_2) = { sR*-1, sR*-1, sR*-1, lc };
pFar5_3 = newp; Point(pFar5_3) = { sR*1,  sR*-1, sR*-1, lc };

arcFar5_0 = newl; Circle(arcFar5_0) = { pFar5_0, O, pFar5_1 };
arcFar5_1 = newl; Circle(arcFar5_1) = { pFar5_1, O, pFar5_2 };
arcFar5_2 = newl; Circle(arcFar5_2) = { pFar5_2, O, pFar5_3 };
arcFar5_3 = newl; Circle(arcFar5_3) = { pFar5_3, O, pFar5_0 };

lRadBR5_0 = newl; Line(lRadBR5_0) = { pOut5_0, pFar5_0 };
lRadBR5_1 = newl; Line(lRadBR5_1) = { pOut5_1, pFar5_1 };
lRadBR5_2 = newl; Line(lRadBR5_2) = { pOut5_2, pFar5_2 };
lRadBR5_3 = newl; Line(lRadBR5_3) = { pOut5_3, pFar5_3 };

clFar5 = newll; Curve Loop(clFar5) = { arcFar5_0, arcFar5_1, arcFar5_2, arcFar5_3 };
sFar5  = news;  Surface(sFar5)     = { clFar5 };

clSideAir5_0 = newll; Curve Loop(clSideAir5_0) = { arcFar5_0, -lRadBR5_1, -arcOut5_0, lRadBR5_0 };
sSideAir5_0  = news;  Surface(sSideAir5_0)     = { clSideAir5_0 };
clSideAir5_1 = newll; Curve Loop(clSideAir5_1) = { arcFar5_1, -lRadBR5_2, -arcOut5_1, lRadBR5_1 };
sSideAir5_1  = news;  Surface(sSideAir5_1)     = { clSideAir5_1 };
clSideAir5_2 = newll; Curve Loop(clSideAir5_2) = { arcFar5_2, -lRadBR5_3, -arcOut5_2, lRadBR5_2 };
sSideAir5_2  = news;  Surface(sSideAir5_2)     = { clSideAir5_2 };
clSideAir5_3 = newll; Curve Loop(clSideAir5_3) = { arcFar5_3, -lRadBR5_0, -arcOut5_3, lRadBR5_3 };
sSideAir5_3  = news;  Surface(sSideAir5_3)     = { clSideAir5_3 };

slOut5 = newsl; Surface Loop(slOut5) = { sFar5, sSideAir5_0, sSideAir5_1, sSideAir5_2, sSideAir5_3, -sOut5 };
vOut5  = newv;  Volume(vOut5)       = { slOut5 };

Transfinite Line{ arcFar5_0, arcFar5_2 } = Ni + 1;
Transfinite Line{ arcFar5_1, arcFar5_3 } = Nj + 1;
Transfinite Line{ lRadBR5_0, lRadBR5_1, lRadBR5_2, lRadBR5_3 } = Nr_air + 1 Using Progression gf_air;
Transfinite Surface { sFar5, sSideAir5_0, sSideAir5_1, sSideAir5_2, sSideAir5_3 };
Recombine  Surface { sFar5, sSideAir5_0, sSideAir5_1, sSideAir5_2, sSideAir5_3 };
Transfinite Volume  { vOut5 } = { pOut5_0, pOut5_1, pOut5_2, pOut5_3, pFar5_0, pFar5_1, pFar5_2, pFar5_3 };
Recombine  Volume   { vOut5 };

// ---------------- Post/physical groups & options ----------------
Coherence;

// r = a  (optional, for postprocess)
Physical Surface("InnerSphere")   = { sIn0, sIn1, sIn2, sIn3, sIn4, sIn5 };

// r = b  (interface Shell-OuterAir; keep old name if your pipeline expects it)
Physical Surface("OuterSphere")   = { sOut0, sOut1, sOut2, sOut3, sOut4, sOut5 };

// r = R  (apply background-field Dirichlet for Nedelec edge DOFs here)
Physical Surface("OuterBoundary") = { sFar0, sFar1, sFar2, sFar3, sFar4, sFar5 };

// volumes
Physical Volume("Shell")    = { v0, v1, v2, v3, v4, v5 };
Physical Volume("OuterAir") = { vOut0, vOut1, vOut2, vOut3, vOut4, vOut5 };
Physical Volume("All")    = { v0, v1, v2, v3, v4, v5 , vOut0, vOut1, vOut2, vOut3, vOut4, vOut5};

// Mesh options
Mesh.ElementOrder = 1;
Mesh.RecombineAll = 0;
Mesh.Optimize = 1;
Mesh.Format = 2;           // 2: MSH2; use 4.1 if your solver prefers it
Geometry.Tolerance = 1e-12;

// Usage:
//   gmsh -3 team11_shell_outerair_built_in.geo
