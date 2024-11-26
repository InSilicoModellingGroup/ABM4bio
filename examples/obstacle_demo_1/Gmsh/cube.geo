Point(1) = {-45,-45,-45};
Point(2) = { 45,-45,-45};
Point(3) = { 45, 45,-45};
Point(4) = {-45, 45,-45};
Point(5) = {-45,-45, 45};
Point(6) = { 45,-45, 45};
Point(7) = { 45, 45, 45};
Point(8) = {-45, 45, 45};
//+
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 5};
//+
Line(9) = {1, 5};
//+
Line(10) = {2, 6};
//+
Line(11) = {3, 7};
//+
Line(12) = {4, 8};
//+
Curve Loop(1) = {1, 10, -5, -9};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {11, -6, -10, 2};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {-7, 12, 3, -11};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {4, 9, -8, -12};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {-3, -4, -1, -2};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {7, 8, 5, 6};
//+
Plane Surface(6) = {6};
