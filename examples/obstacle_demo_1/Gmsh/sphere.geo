Point(1) = {45, 0, 0};
//+
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Point{1}; 
}
//+
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Point{2}; 
}
//+
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {
  Curve{1}; Curve{2}; 
}
//+
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {
  Curve{3}; Curve{6}; 
}
//+
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {
  Curve{9}; Curve{12}; 
}
//+
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {
  Curve{15}; Curve{18}; 
}
