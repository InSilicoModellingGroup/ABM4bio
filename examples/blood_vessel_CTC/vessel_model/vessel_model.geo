//+
RAD = 40.0;
HEI = 100.0;
//+
Point(1) = {0.0, 0.0, -HEI};
//+
Point(2) = {RAD, 0.0, -HEI};
Point(3) = {0.0, RAD, -HEI};
Point(4) = {-RAD,0.0, -HEI};
Point(5) = {0.0,-RAD, -HEI};
//+
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
//+
//Extrude {0.0, 0.0, 2*HEI} { Curve{1:4}; Layers{6}; Recombine; }
Extrude {0.0, 0.0, 2*HEI} { Curve{1:4}; }
//+
Physical Surface(2000) = {8:20:4};
//+
Mesh.Algorithm = 2;
Mesh.Algorithm3D = 4;
Mesh.OptimizeNetgen = 1;
Mesh.CharacteristicLengthMin = 20.0;
Mesh.CharacteristicLengthMax = 20.01;
Mesh.Format = 1;
Mesh.MshFileVersion = 2.2;
//+
