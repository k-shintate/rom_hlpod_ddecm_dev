// ============================================================
// ahmed_ext_stable.geo 
// Target: 粗めのメッシュ + 安定化設定（クラッシュ回避版）
// ============================================================

SetFactory("Built-in");
Mesh.MshFileVersion = 2.2; 

// --- Stability & Algorithm Settings (重要) ---
// HXT(10)は速いですがクラッシュしやすいため、
// 標準Delaunay(1)に変更して計算を完走させます。
Mesh.Algorithm3D = 1; 

// 幾何学的な「ハマり」を防ぐための微小な座標揺らぎ
Mesh.RandomFactor = 1e-6; 
Mesh.RandomFactor3D = 1e-6;

// STLの品質が悪い場合に、より積極的に修正・無視する設定
Mesh.StlRemoveBadTriangles = 2; 

// 3Dメッシュ最適化の回数（品質向上）
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 0; // Netgen最適化は重いので一旦OFF

// --- Performance ---
General.NumThreads = 16; 

gap = 1e-3; 

// ============================================================
// [COARSE SIZES] 計算負荷を下げるための粗い設定
// ============================================================

lcFar     = 0.8;    // 800mm (遠方)
lcWake    = 0.15;   // 150mm (後流)
lcNearBox = 0.08;   // 80mm  (車体周辺)
lcSurf    = 0.02;   // 20mm  (車体表面)

// Threshold settings (表面解像度からの遷移)
dMin   = 0.02;   // 2cmまではlcSurfを維持
dMax   = 0.3;    // 30cmかけてlcNearBoxへ遷移

// ------------------------------------------------------------
// 1) Import STL
Merge "ahmed_25deg_m.stl";
CreateTopology;

stlPts[]   = Point{:};
stlCur[]   = Curve{:};
stlSurfs[] = Surface{:};

// 地面干渉回避のための微小オフセット
Translate {0, 0, +gap} {
  Point{stlPts[]};
  Curve{stlCur[]};
  Surface{stlSurfs[]};
}
Coherence Mesh;

// ------------------------------------------------------------
// 2) Outer box (風洞領域)
z0 = -gap;
H  = 8 + gap;

Point(1) = {-5, -4, z0, lcFar};
Point(2) = {15, -4, z0, lcFar};
Point(3) = {15,  4, z0, lcFar};
Point(4) = {-5,  4, z0, lcFar};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(10) = {1,2,3,4};
Plane Surface(100) = {10}; 

out[] = Extrude {0,0,H} { Surface{100}; };
topS   = out[0];
tmpVol = out[1];
side1  = out[2]; 
side2  = out[3]; 
side3  = out[4]; 
side4  = out[5]; 
Delete { Volume{tmpVol}; }

// ------------------------------------------------------------
// 3) Fluid volume
Surface Loop(200) = {100, topS, side1, side2, side3, side4};
Surface Loop(201) = {stlSurfs[]};
Volume(300) = {200, 201};

// ------------------------------------------------------------
// 4) Size fields (メッシュサイズ制御)

// Field 1: 車体表面からの距離
Field[1] = Distance;
Field[1].SurfacesList = {stlSurfs[]};
Field[1].Sampling = 20; // 軽さを優先して低めに

// Field 2: 表面近傍のグラデーション
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lcSurf;
Field[2].SizeMax = lcNearBox;
Field[2].DistMin = dMin;
Field[2].DistMax = dMax;

// Field 3: 車体周囲ボックス
Field[3] = Box;
Field[3].VIn = lcNearBox;
Field[3].VOut = lcFar;
Field[3].XMin = -1.5; 
Field[3].XMax = 2.5;
Field[3].YMin = -1.0;
Field[3].YMax = 1.0;
Field[3].ZMin = z0;
Field[3].ZMax = 1.5;
Field[3].Thickness = 0.5;

// Field 4: 後流領域ボックス
Field[4] = Box;
Field[4].VIn = lcWake;
Field[4].VOut = lcFar;
Field[4].XMin = 2.5;
Field[4].XMax = 8.0;
Field[4].YMin = -0.8;
Field[4].YMax = 0.8;
Field[4].ZMin = z0;
Field[4].ZMax = 1.2;
Field[4].Thickness = 1.0;

// 統合フィールド
Field[5] = Min;
Field[5].FieldsList = {2, 3, 4};

Background Field = 5;

// ------------------------------------------------------------
// 5) Physical groups
Physical Volume("Fluid") = {300};
Physical Surface("lowerWall") = {100};
Physical Surface("upperWall") = {topS};
Physical Surface("Inlet")  = {side4};           
Physical Surface("outlet") = {side2};           
Physical Surface("frontAndBack") = {side1,side3}; 
Physical Surface("ahmed_body") = {stlSurfs[]};

