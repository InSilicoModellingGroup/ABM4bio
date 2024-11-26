Point(1) = {0,-100,-100};
Extrude {0, 30, 0} {
  Point{1};
}
Extrude {-100, 0, 0} {
  Point{2};
}
Extrude {0, 0, 200} {
  Curve{1}; Curve{2};
}
