lc=0.006;
a=1.0;//width
b=1.0;//1.0;//height
Point(1) = {0, 0, 0, lc};
Point(2) = {a, 0, 0, lc};
Point(3) = {a, b, 0, lc};
Point(4) = {0, b, 0, lc};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};
Physical Line(7) = {1};
Physical Line(8) = {2};
Physical Line(9) = {3};
Physical Line(10) = {4};
Physical Surface(11) = {6};
//+
Point(5) = {-0, 0.5, 0, 1.0};
//+
Point(6) = {1, 0.5, -0, 1.0};
//+
Line(5) = {5, 6};
//+
Recursive Delete {
  Line{5}; 
}
