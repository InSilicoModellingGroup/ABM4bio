R = 480;
//+
Point(1) = {+0,+0,+0};
Point(2) = {+0,+R,+0};
Point(3) = {+R,+0,+0};
Point(4) = {+0,-R,+0};
Point(5) = {-R,+0,+0};
//+
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
//+
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
//+
Mesh.MeshSizeMin = 15;
Mesh.MeshSizeMax = 30;
//+
Mesh.MeshSizeMin = 90;
Mesh.MeshSizeMax =110;
//+
Recursive Delete {
  Point{1}; 
}
//+
Recursive Delete {
  Point{1}; 
}
//+
Recursive Delete {
  Point{1}; 
}
