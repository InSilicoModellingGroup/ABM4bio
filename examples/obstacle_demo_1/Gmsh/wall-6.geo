// Gmsh project created on Thu Sep  2 14:22:34 2021
//+
Point(1) = {-100, -90, -100};
//+
Extrude {0, 0, +200} {
  Point{1}; 
}
//+
Extrude {+200, 0, 0} {
  Curve{1}; 
}
