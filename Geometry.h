/*
 * Author: Oscar Martinez
 * Copyright (c) 2012 All rights reserved.
 */

#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include "LPSolver.h"

// structure definitions

typedef struct Point Point;
struct Point {
    float x;
    float y;
};

typedef struct Line Line;
struct Line {
    float m;
    float n;
};

typedef struct Circle Circle;
struct Circle {
  float center_x;
	float center_y;
	float radius;
};

typedef struct Vertice Vertice;
struct Vertice {
    float x;
    float y;
    Vertice *next;
};

// line: y = m*x + n
typedef struct Edge Edge;
struct Edge {
    Line* l; // for inheritance purposes
    Vertice* v1;
    Vertice* v2;
    Edge* next;
};

typedef struct Polygon Polygon;
struct Polygon {
    Edge* first_edge;
    int n;
};

// constant definitions
#define N_CELLS 20
#define M_CELLS 20
#define CONSECUTIVE_MISS	15
#define METHOD_SEQUENTIAL	1
#define METHOD_RANDOMIZED	2
#define METHOD_LP			3


// constructors and destructors

Point* Point_new(float, float);
void Point_destroy(Point*);

Line* Line_new(void);
void Line_destroy(Line*);

Circle* Circle_new(float, float, float);
void Circle_destroy (Circle*);

Vertice* Vertice_new(float, float);
void Vertice_destroy (Vertice*);

Edge* Edge_new(void);
void Edge_destroy (Edge*);

Polygon* Polygon_new(void);
void Polygon_destroy(Polygon*);

// special constructors

Line* Line_new_given_m_and_n(float, float);
Line* Line_new_given_two_points(Point*, Point*);
Line* Line_new_given_coordinates(float, float, float, float);
Line* Line_new_given_point_slope(Point*, float);

Edge* Edge_new_given_coordinates(float, float, float, float);
Edge* Edge_new_given_vertices(Vertice*, Vertice*);

Polygon* Polygon_new_given_first_edge(Edge*);
Polygon* Polygon_new_given_first_vertice(Vertice*, char);
Polygon* Polygon_new_given_matrix(Matrix_2D*, char);


// functions used for debugging

char* Point_to_string(Point*);
char* Line_to_string(Line*);
void Line_print(Line*); //tmp
void Edge_print(Edge*);
char* Circle_to_string(Circle*);
char* Polygon_to_string(Polygon*);


// functions used for visualization

void Line_draw(Line*);
void Edge_draw(Edge*);
void Polygon_draw(Polygon*);


// module-oriented functions

char Geometry_compare_vertices(Vertice*, Vertice*);
char Geometry_is_vertice_in_list(Vertice*, Vertice*);
float Geometry_distance_between_two_points(Point*, Point*);
float Geometry_distance_between_coordinates(float, float, float, float);
float Geometry_distance_between_point_and_vertice(Vertice*, Point*);
float Geometry_distance_between_point_and_line(Line*, Point*);
float Geometry_distance_between_point_and_edge(Edge*, Point*);
float Geometry_distance_between_point_and_polygon(Polygon*, Point*);
Edge* Geometry_closest_polygon_edge_to_point(Polygon*, Point*);
char Geometry_is_point_in_line(Line*, float, float);
char Geometry_is_point_in_edge(Edge*, float, float);
char Geometry_is_point_in_polygon(Polygon*, float, float);
//Point* Geometry_find_point_in_line_given_x(Line*, float);
//Point* Geometry_find_point_in_line_given_y(Line*, float);
Point* Geometry_find_point_in_line_given_point_and_distance(Line*, Point*, float);
Point* Geometry_intersection_between_two_lines(Line*, Line*);
Vertice* Geometry_intersections_between_polygon_and_line(Polygon*, Line*);
Line* Geometry_perpendicular_line_at_point(Line*, float, float);
Polygon* Geometry_generate_random_convex_polygon(float[], int);
Point* Geometry_find_biggest_inscribed_circle_in_polygon(Polygon*, int);


// unit testing

void Geometry_test_line_1(void);
void Geometry_test_polygon_1(void);
void Geometry_test_polygon_2(void);
void Geometry_test_polygon_3(void);
void Geometry_test_find_point_in_line_1(void);
void Geometry_test_point_in_edge(void);
void Geometry_test_two_lines_intersection_1(void);
void Geometry_test_two_lines_intersection_2(void);
void Geometry_test_two_lines_intersection_3(void);
void Geometry_test_line_polygon_intersection_1(void);
void Geometry_test_closest_edge(void);
void Geometry_test_distance_from_point_to_edge_1(void);
void Geometry_test_point_in_polygon(void);

#endif /* GEOMETRY_H_ */