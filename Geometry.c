/*
 * Author: Oscar Martinez
 * Copyright (c) 2012 All rights reserved.
 */

// library includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <sys/time.h>

// global includes
#include "../SVG/Draw.h"

// module includes
#include "Geometry.h"
#include "LPSolver.h"

// point constructor
Point* Point_new (float x, float y) {
    Point *this = malloc(sizeof(Point));
    this->x = x;
    this->y = y;
    return this;
}

void Point_destroy (Point* point) {
  if (point == NULL) return;
	free(point);
}



// line constructor 1
Line* Line_new_given_m_and_n(float m, float n) {
    Line *this = malloc(sizeof (Line));
    this->m = m;
    this->n = n;
    return this;
}

// line constructor 2
Line* Line_new_given_two_points(Point* p1, Point* p2) {
	return Line_new_given_coordinates(p1->x, p1->y, p2->x, p2->y);
}

// line constructor 3
Line* Line_new_given_coordinates(float x0, float y0, float x1, float y1) {

	// special case: vertical line
	if (x0 == x1) return Line_new_given_m_and_n(FLT_MAX, x0);

	// otherwise
	float m = (y1 - y0) / (x1 - x0);
	return Line_new_given_m_and_n(m, -m * x0 + y0);
}

// line constructor 5
Line* Line_new_point_slope(Point* p, float m) {
	float n = p->y - p->x * m;
	return Line_new_given_m_and_n(m, n);
}

void Line_destroy (Line* line) {
	free (line);
}


// circle constructor
Circle* Circle_new(float x, float y, float r) {
    Circle *this = malloc(sizeof (Circle));
    this->center_x = x;
    this->center_y = y;
    this->radius = r;
    return this;
}

void Circle_destroy (Circle* circle) {
	free(circle);
}

// vertice constructor
Vertice* Vertice_new (float x, float y) {
    Vertice *this = malloc(sizeof (Vertice));
    this->x = (float)x;
    this->y = (float)y;
    this->next = NULL;

    return this;
}

void Vertice_destroy (Vertice* v) {
	free(v);
}

// edge constructor 1
Edge* Edge_new_given_coordinates(float x1, float y1, float x2, float y2) {
    Edge *this = malloc(sizeof (Edge));
    this->v1 = Vertice_new(x1, y1);
    this->v2 = Vertice_new(x2, y2);
    this->l = Line_new_given_coordinates(x1, y1, x2, y2);
    this->next = NULL;
    return this;
}

// edge constructor 2
Edge* Edge_new_given_vertices(Vertice* v1, Vertice* v2) {
	Edge* this = malloc(sizeof (Edge));
	this->v1 = Vertice_new(v1->x, v1->y);
	this->v2 = Vertice_new(v2->x, v2->y);
//	this->v1 = v1;
//	this->v2 = v2;
	this->l = Line_new_given_coordinates(v1->x, v1->y, v2->x, v2->y);
	this->next = NULL;
	return this;
}

// TODO: fix edge destroy
void Edge_destroy (Edge* e) {
	if (e->l != NULL) Line_destroy(e->l);
//	if (e->v1 != NULL) Vertice_destroy(e->v1);
//	if (e->v2 != NULL) Vertice_destroy(e->v2);
	free(e);
}

// polygon constructor 1
Polygon* Polygon_new_given_first_edge(Edge* firstedge) {
	Polygon *this = malloc(sizeof (Polygon));
	this->first_edge = firstedge;

	// count edges and link vertices
	int count=0;
	Edge* edge_tmp;
	for (edge_tmp = this->first_edge; edge_tmp != NULL; edge_tmp = edge_tmp->next, count++) {
		edge_tmp->v1->next = edge_tmp->v2;
	}
	this->n = count;
	return this;
}

// polygon constructor 2
Polygon* Polygon_new_given_first_vertice(Vertice* first_vertice, char close_polygon) {
	Polygon *this = malloc(sizeof (Polygon));
	Edge	*first_edge = NULL,
			*edge_prev = NULL;

	int n = 0;
	Vertice* vertice_tmp;
	Edge* edge_tmp = NULL;
	for (vertice_tmp = first_vertice; vertice_tmp != NULL; vertice_tmp = vertice_tmp->next) {

		if (vertice_tmp->next != NULL) {
			edge_tmp = Edge_new_given_vertices(vertice_tmp, vertice_tmp->next);
		} else if (close_polygon) {
			edge_tmp = Edge_new_given_vertices(vertice_tmp, first_vertice);
		}

		// prepare for the next step
		if (edge_prev != NULL) edge_prev->next = edge_tmp;
		else first_edge = edge_tmp;
		edge_prev = edge_tmp;
		edge_tmp->next = NULL;
		n++;
	}

	this->n = n;
	this->first_edge = first_edge;
	return this;
}

// polygon constructor 3
Polygon* Polygon_new_given_matrix(Matrix_2D* matrix, char close_polygon) {
	int i;
	Vertice *vertice_prev = NULL, *vertice_first = NULL;
	for (i = 0; i < matrix->rows; i++) {
		Vertice *vertice_tmp = Vertice_new(
			Matrix_2D_get_cell(matrix, i, 0), Matrix_2D_get_cell(matrix, i, 1));
		if (vertice_prev != NULL) vertice_prev->next = vertice_tmp;
		else vertice_first = vertice_tmp;
		vertice_prev = vertice_tmp;
	}
	Polygon* polygon = Polygon_new_given_first_vertice(vertice_first, close_polygon);
	return polygon;
}

void Polygon_destroy(Polygon* polygon) {

	// destroy vertices
	Vertice* vertice_tmp = NULL, *vertice_prev = NULL;
	for (vertice_tmp = polygon->first_edge->v1; vertice_tmp != NULL; vertice_prev = vertice_tmp, vertice_tmp = vertice_tmp->next) {
		if (vertice_prev != NULL) Vertice_destroy(vertice_prev);
	}

	// destroy last vertice
	if (vertice_prev != NULL) Vertice_destroy(vertice_prev);

	// destroy edges
	Edge* edge_tmp = NULL, *edge_prev = NULL;
	for (edge_tmp = polygon->first_edge;
	edge_tmp != NULL;
	edge_prev = edge_tmp, edge_tmp = edge_tmp->next) {
		if (edge_prev != NULL) {
			edge_prev->v1 = NULL;
			edge_prev->v2 = NULL;
			Edge_destroy(edge_prev);
		}
	}

	// destroy last edge
	if (edge_prev != NULL) {
		edge_prev->v1 = NULL;
		edge_prev->v2 = NULL;
		Edge_destroy(edge_prev);
	}

	free (polygon);
}

char* Polygon_to_string(Polygon* polygon) {
	Vertice* v;
	char buf[128];
	int index = 0;
	for (v = polygon->first_edge->v1; v != NULL; v = v->next) {
		index += sprintf(&buf[index], "(%f,%f) ", v->x, v->y);
	}
	return strdup(buf);
}

void Point_draw(Point* point) {
	Draw_circle(current_canvas, point->x, point->y, 3);
}

void Line_print(Line* l) {
	if (l->m == FLT_MAX) printf("x = %f\n", l->n);
	else if (fabs(l->m) < 0.01f) printf("y = %f\n", l->n);
	else printf("y = %fx + %f\n", l->m, l->n);
}

void Line_draw(Line* l) {
	Draw_line(current_canvas,
			l->m == FLT_MAX ? l->n : -current_canvas->width / 2,
			l->m == FLT_MAX ? 0 : -current_canvas->width / 2 * l->m + l->n,
			l->m == FLT_MAX ? l->n : current_canvas->width / 2,
			l->m == FLT_MAX ? current_canvas->width : current_canvas->width / 2 * l->m + l->n);
}


void Edge_print(Edge* edge) {
	printf("(%f,%f) (%f,%f)\n",
		edge->v1->x, edge->v1->y,
		edge->v2->x, edge->v2->y);
}


void Edge_draw(Edge* edge) {
	Draw_line(current_canvas, edge->v1->x, edge->v1->y, edge->v2->x, edge->v2->y);
	Draw_circle(current_canvas, edge->v1->x, edge->v1->y, 3);
	Draw_circle(current_canvas, edge->v2->x, edge->v2->y, 3);
}




void Polygon_draw(Polygon* polygon) {
	Edge* tmp;
	for (tmp = polygon->first_edge; tmp != NULL; tmp = tmp->next) {
		Draw_line(current_canvas, tmp->v1->x, tmp->v1->y, tmp->v2->x, tmp->v2->y);
		Draw_circle(current_canvas, tmp->v1->x, tmp->v1->y, 3);
		Draw_circle(current_canvas, tmp->v2->x, tmp->v2->y, 3);
	}
}

char Geometry_compare_vertices(Vertice* v1, Vertice* v2) {
	return (fabs(v1->x - v2->x) <= 0.01f && fabs(v1->y - v2->y) <= 0.01f) ? 1 : 0;
}

char Geometry_is_vertice_in_list(Vertice* first_vertice, Vertice* search_vertice) {
	Vertice* v;
	for (v = first_vertice; v != NULL; v = v->next) {
		if (Geometry_compare_vertices(search_vertice, v)) return 1;
	}
	return 0;
}

float Geometry_distance_between_two_points(Point* p1, Point* p2) {
	return Geometry_distance_between_coordinates(p1->x, p1->y, p2->x, p2->y);
}

float Geometry_distance_between_coordinates(float x1, float y1, float x2, float y2) {
	return sqrtf((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

float Geometry_distance_between_point_and_vertice(Vertice* vertice, Point* point) {
	return Geometry_distance_between_coordinates(point->x, point->y, vertice->x, vertice->y);
}


float Geometry_distance_between_point_and_line(Line* line, Point* point) {
	if (line->m == FLT_MAX) return fabs(line->n - point->x); // vertical lines
	float x1 = (point->y - point->x / line->m - line->n) / (line->m - 1 / line->m);
	float y1 = line->m * x1 + line->n;
	return Geometry_distance_between_coordinates(x1, y1, point->x, point->y);
}

// given a line, a point in that line, and a distance find a point in the line
// at distance from the initial point where positive distance is to the right.
// returns NULL if the initial point is not in the line
Point* Geometry_find_point_in_line_given_point_and_distance(Line* line, Point* point,
		float distance) {
	Point* point_tmp;

	// case 1: initial point not in line
	if (!Geometry_is_point_in_line(line, point->x, point->y)) {
		point_tmp = NULL;

	// case 2: vertical line
	} else if (line->m == FLT_MAX) {
		point_tmp = Point_new(point->x, point->y + distance);

	// case 3: any other line
	} else {
		float float_tmp = atanf(line->m);
		point_tmp = Point_new(
				point->x + cosf(float_tmp) * distance, point->y + sinf(float_tmp) * distance);
	}

	return point_tmp;
}

float Geometry_distance_between_point_and_edge(Edge* edge, Point* point) {
	float distance = FLT_MAX;

	// NULL checks
	if (point == NULL || edge == NULL) return distance;

	// find line perpendicular to segment that goes through the point
	Line* t = Geometry_perpendicular_line_at_point(edge->l, point->x, point->y);

	// intersection between that and the edge is the closest point within edge
	Point* closest_pt = Geometry_intersection_between_two_lines(edge->l, t);

	// if point belongs to segment return the euclidean distance
	if (Geometry_is_point_in_edge(edge, closest_pt->x, closest_pt->y)) {
		distance = Geometry_distance_between_two_points(point, closest_pt);

	// otherwise return distance to closest edge
	} else {
		float f1 = Geometry_distance_between_point_and_vertice(edge->v1, point);
		float f2 = Geometry_distance_between_point_and_vertice(edge->v2, point);
		distance = (f1 < f2 ? f1 : f2);
	}

	// cleanup and return
	Point_destroy(closest_pt);
	Line_destroy(t);
	return distance;

}

float Geometry_distance_between_point_and_polygon(Polygon* polygon, Point* point) {
	Edge* edge = Geometry_closest_polygon_edge_to_point(polygon, point);
	return Geometry_distance_between_point_and_edge(edge, point);
}

Edge* Geometry_closest_polygon_edge_to_point(Polygon* polygon, Point* point) {
	Edge* edge_tmp, *closest_edge = NULL;
	float distance_min = FLT_MAX, distance_tmp, distance_comp;
	for (edge_tmp = polygon->first_edge; edge_tmp != NULL; edge_tmp = edge_tmp->next) {
//		if (closest_edge != NULL) {
//			printf("closest edge: ");
//			Edge_print(closest_edge);
//			printf("temp edge: ");
//			Edge_print(edge_tmp);
//		}
		distance_tmp = Geometry_distance_between_point_and_edge(edge_tmp, point);

		// TODO: fix this too!
		if (distance_tmp == distance_min
		&& closest_edge != NULL && Geometry_compare_vertices(closest_edge->v1, edge_tmp->v2)) {
			Point* edge_point_1 = Point_new(closest_edge->v1->x, closest_edge->v1->y);
			Point* edge_point_2 = Point_new(edge_tmp->v2->x, edge_tmp->v2->y);
			Point* p1 = Geometry_find_point_in_line_given_point_and_distance(
					edge_tmp->l, edge_point_1, 0.01f);
			Point* p2 = Geometry_find_point_in_line_given_point_and_distance(
					edge_tmp->l, edge_point_2, 0.01f);
			distance_tmp = Geometry_distance_between_two_points(point, p1);
			distance_comp = Geometry_distance_between_two_points(point, p2);
			Point_destroy(p1);
			Point_destroy(p2);
			Point_destroy(edge_point_1);
			Point_destroy(edge_point_2);

		} else if (distance_tmp == distance_min
		&& closest_edge != NULL && Geometry_compare_vertices(closest_edge->v2, edge_tmp->v1)) {
			Point* edge_point_1 = Point_new(closest_edge->v2->x, closest_edge->v2->y);
			Point* edge_point_2 = Point_new(edge_tmp->v1->x, edge_tmp->v1->y);
			Point* p1 = Geometry_find_point_in_line_given_point_and_distance(
					edge_tmp->l, edge_point_1, 0.01f);
			Point* p2 = Geometry_find_point_in_line_given_point_and_distance(
					edge_tmp->l, edge_point_2, 0.01f);
			distance_tmp = Geometry_distance_between_two_points(point, p1);
			distance_comp = Geometry_distance_between_two_points(point, p2);
			Point_destroy(p1);
			Point_destroy(p2);
			Point_destroy(edge_point_1);
			Point_destroy(edge_point_2);

		} else {
			distance_comp = distance_min;
		}

		if (distance_tmp < distance_comp) {
			distance_min = Geometry_distance_between_point_and_edge(edge_tmp, point);
			closest_edge = edge_tmp;
//			closest_edge = Edge_new_given_vertices(edge_tmp->v1, edge_tmp->v2);
		}
	}
	return closest_edge;
}

char Geometry_is_point_in_line(Line* line, float x, float y) {
	if (line->m == FLT_MAX && line->n == x) return 1;
	if (fabs(line->m * x + line->n - y) < 0.01f) return 1;
	else return 0;
}

char Geometry_is_point_in_edge(Edge* edge, float x, float y) {
	// case 1: point not in line
	if (!Geometry_is_point_in_line(edge->l, x, y)) {
		return 0;

	// case 2: vertical line
	} else if (edge->l->m == FLT_MAX
	&& ((y - 0.01f <= edge->v1->y && y + 0.01f >= edge->v2->y)
	|| (y + 0.01f >= edge->v1->y && y - 0.01f <= edge->v2->y))) {
		return 1;

	// case 3: any other line
	} else if (edge->l->m != FLT_MAX
	&& ((x + 0.01f >= edge->v1->x && x - 0.01f <= edge->v2->x)
	|| (x - 0.01f <= edge->v1->x && x + 0.01f >= edge->v2->x))) {
		return 1;

	// otherwise point not in edge
	} else {
		return 0;
	}
}


/**
 * Function used to test whether a point is inside of a closed polygon or not. Note
 * that both a horizontal and a vertical projection are used since there will be the
 * same number of intersections between the lines and the vertices of the polygon whether
 * the point is inside or outside.
 */
char Geometry_is_point_in_polygon(Polygon* polygon, float x, float y) {
	Line* horizontal_projection = Line_new_given_m_and_n(0.0f, y);
	Line* vertical_projection = Line_new_given_m_and_n(FLT_MAX, x);

	// count the number of intersections between the polygon and the
	// horizontal and vertical projections of the point
	Vertice* i1 = Geometry_intersections_between_polygon_and_line(polygon, horizontal_projection);
	Vertice* i2 = Geometry_intersections_between_polygon_and_line(polygon, vertical_projection);

	int count1 = 0, count2 = 0;
	Vertice* v_tmp;
	for (v_tmp = i1; v_tmp != NULL; v_tmp = v_tmp->next) {
		if (v_tmp->x >= x
		&& Geometry_is_point_in_line(horizontal_projection, v_tmp->x, v_tmp->y)
		&& !Geometry_is_vertice_in_list(polygon->first_edge->v1, v_tmp)) count1++;
	}
	for (v_tmp = i2; v_tmp != NULL; v_tmp = v_tmp->next) {
		if (v_tmp->y >= y
		&& Geometry_is_point_in_line(vertical_projection, v_tmp->x, v_tmp->y)
		&& !Geometry_is_vertice_in_list(polygon->first_edge->v1, v_tmp)) count2++;
	}
	Vertice_destroy(i1);
	Vertice_destroy(i2);

	if (count1 % 2 == 0 || count2 % 2 == 0) return 0;
	else return 1;
}

Line* Geometry_perpendicular_line_at_point(Line* l, float x, float y) {

	// special cases: vertical and horizontal lines
	if (l->m == FLT_MAX) return Line_new_given_m_and_n(0.0, y);
	if (fabs(l->m) < 0.01f) return Line_new_given_m_and_n(FLT_MAX, x);

	// otherwise
	float m = -1/l->m;
	return Line_new_given_m_and_n(m,y-m*x);
}

Point* Geometry_intersection_between_two_lines(Line* l1, Line* l2) {

	// special cases: vertical and horizontal lines
	if (l1->m == FLT_MAX) return Point_new(l1->n, l2->m * l1->n + l2->n);
	if (l2->m == FLT_MAX) return Point_new(l2->n, l1->m * l2->n + l1->n);
	if (l1->m == l2->m && l1->n == l2->n) return Point_new(0, l1->n);

	// paralel lines
	if (fabs(l1->m - l2->m) < 0.01f) return NULL;

	// otherwise
	float x = (l2->n - l1->n) / (l1->m - l2->m);
	return Point_new(x, l1->m * x + l1->n);
}



// find all intersections between a polygon and a line and return a list of points using struct Vertice
Vertice* Geometry_intersections_between_polygon_and_line(Polygon* polygon, Line* line) {
	Vertice* vertice_first = NULL, *vertice_last = NULL;

	Edge* tmp;
	for(tmp = polygon->first_edge; tmp != NULL; tmp = tmp->next) {
		Point* p = Geometry_intersection_between_two_lines(tmp->l, line);
		if (p != NULL && Geometry_is_point_in_edge(tmp, p->x, p->y)) {
			Vertice* vertice_new = Vertice_new(p->x, p->y);
			if (vertice_first == NULL) {
				vertice_last = vertice_first = vertice_new;
			} else if (!Geometry_is_vertice_in_list(vertice_first, vertice_new)) {
				vertice_last = vertice_last->next = Vertice_new(p->x, p->y);
			} else {
//				printf("repeated intersection: (%f,%f) & (%f,%f)\n",
//						p->x, p->y, vertice_last->x, vertice_last->y);
			}
			Point_destroy(p);
		}
	}

	return vertice_first;
}

// given bounds for the x and y coordinates and a number of vertices, generate
// a random monotone polygon
Polygon* Geometry_generate_random_convex_polygon(float bounds[], int vertices) {
	struct timeval time_seed;
	gettimeofday(&time_seed, NULL);
	srand(time_seed.tv_usec);

	// create a triangle as a base for our polygon
	int i;
	Vertice* triangle_vertices[3];
	for (i = 0; i < 3; i++) {
		float x = bounds[0] + (rand() % (int) ((bounds[1] - bounds[0]) * 100000) / 100000.0f);
		float y = bounds[2] + (rand() % (int) ((bounds[3] - bounds[2]) * 100000) / 100000.0f);
		triangle_vertices[i] = Vertice_new(x, y);
		if (i > 0) triangle_vertices[i - 1]->next =  triangle_vertices[i];
	}

	Polygon* polygon = Polygon_new_given_first_vertice(triangle_vertices[0], 1);

	// cleanup the no longer used vertices
	for (i = 0; i < 3; i++)
		Vertice_destroy(triangle_vertices[i]);

	// add vertices to the polygon until the number of vertices required is reached
	float x, y;
	int vertice_count = 3;
	while (vertice_count < vertices) {

		// generate random x and y coordinates
		x = bounds[0] + (rand() % (int) ((bounds[1] - bounds[0]) * 100000) / 100000.0f);
		y = bounds[2] + (rand() % (int) ((bounds[3] - bounds[2]) * 100000) / 100000.0f);

		// if they are outside of the polygon, check intersections with closest edge's vertices
		if (!Geometry_is_point_in_polygon(polygon, x, y)) {

			// point randomly generated
			Point* point_tmp = Point_new(x, y);

			// get the lines between external point vertices of the closest edge
			Edge* closest_edge = Geometry_closest_polygon_edge_to_point(polygon, point_tmp);
			Line* line_1 = Line_new_given_coordinates(
					point_tmp->x, point_tmp->y, closest_edge->v1->x, closest_edge->v1->y);
			Line* line_2 = Line_new_given_coordinates(
					point_tmp->x, point_tmp->y, closest_edge->v2->x, closest_edge->v2->y);

			// make sure it is not part of an edge
			Edge* e_tmp;
			for (e_tmp = polygon->first_edge; e_tmp != NULL; e_tmp = e_tmp->next) {
				if (fabs(e_tmp->l->m - line_1->m) <= 0.01f) {
					goto cleanup;
				}
			}

			// get the intersections between each line and polygon
			Vertice* v1 = Geometry_intersections_between_polygon_and_line(polygon, line_1);
			Vertice* v2 = Geometry_intersections_between_polygon_and_line(polygon, line_2);

			// count the number of intersections and cleanup as we go
			int c1 = 0, c2 = 0;
			Vertice* vertice_tmp, *vertice_last;
			for (vertice_tmp = v1, vertice_last = NULL;
			vertice_tmp != NULL;
			vertice_last = vertice_tmp, vertice_tmp = vertice_tmp->next, c1++)
				if (vertice_last != NULL) Vertice_destroy(vertice_last);
			for (vertice_tmp = v2, vertice_last = NULL;
			vertice_tmp != NULL;
			vertice_last = vertice_tmp, vertice_tmp = vertice_tmp->next, c2++)
				if (vertice_last != NULL) Vertice_destroy(vertice_last);

			// point is only valid if there are no intersections between line and polygon
			if (c1 == 1 && c2 == 1) {
				Vertice* vertice_new = Vertice_new(point_tmp->x, point_tmp->y);

				// create new edges
				Edge* edge_new_1 = Edge_new_given_vertices(closest_edge->v1, vertice_new);
				Edge* edge_new_2 = Edge_new_given_vertices(vertice_new, closest_edge->v2);

				// update closest edge and link edges
				edge_new_2->next = closest_edge->next;
				*closest_edge->v1 = *edge_new_1->v1;
				*closest_edge->v2 = *edge_new_1->v2;
				*closest_edge->l = *edge_new_1->l;
				closest_edge->next = edge_new_2;

				// tmp
				printf("New vertice: (%f,%f)\n", vertice_new->x, vertice_new->y);
				vertice_count++;

				// cleanup edge and vertices
				Edge_destroy(edge_new_1);
				Vertice_destroy(vertice_new);
			}

			// cleanup
			cleanup:
			Point_destroy(point_tmp);
			Line_destroy(line_1);
			Line_destroy(line_2);
		}
//		printf("8\n");
	}

	// link vertices
	Edge* edge_tmp;
	for (edge_tmp = polygon->first_edge; edge_tmp->next != NULL; edge_tmp = edge_tmp->next) {
		edge_tmp->v1->next = edge_tmp->next->v1;
	}

	return polygon;
}

// Inaccessibility of poles problem. Private function using sequential method
Point* Geometry_find_biggest_inscribed_circle_in_polygon_sequential(
		Polygon* polygon, float bounds[]) {
	Point* pia = Point_new((bounds[0] + bounds[1]) / 2, (bounds[2] + bounds[3]) / 2);
	Point* tmp = Point_new(0, 0);

	// calculate the required increment for x and y
	float increment_x = (bounds[1] - bounds[0]) / N_CELLS;
	float increment_y = (bounds[3] - bounds[2]) / M_CELLS;

	// biggest known distance
	float max_distance = 0;

	// loop through x and y coordinates to hit all the nodes
	int i, j;
	float tmp_distance = FLT_MAX;
	for (i = 0; i <= N_CELLS; i++) {

		// find x value
		tmp->x = bounds[0] + i * increment_x;

		for (j = 0; j <= M_CELLS; j++) {


			// find y value
			tmp->y = bounds[2] + j * increment_y;


			// compare with candidate PIA if point is in polygon
			if (Geometry_is_point_in_polygon(polygon, tmp->x, tmp->y)) {
				tmp_distance = Geometry_distance_between_point_and_polygon(polygon, tmp);
//				current_canvas->style = "stroke:blue;stroke-width:2;fill:blue";
//				Draw_circle(current_canvas, tmp->x, tmp->y, 8);
				if (tmp_distance > max_distance) {
					max_distance = tmp_distance;
					pia->x = tmp->x;
					pia->y = tmp->y;
				}
			}
		}
	}

	// cleanup and return
	Point_destroy(tmp);
	return pia;
}

// Inaccessibility of poles problem. Private function using randomized method
Point* Geometry_find_biggest_inscribed_circle_in_polygon_randomized(
		Polygon* polygon, float bounds[]) {
	struct timeval time_seed;
	gettimeofday(&time_seed, NULL);
	srand(time_seed.tv_usec);

	Point* pia = Point_new((bounds[0] + bounds[1]) / 2, (bounds[2] + bounds[3]) / 2);
	Point* tmp = Point_new(0, 0);

	// biggest known distance
	float max_distance = 0;

	// loop through x and y coordinates to hit all the nodes
	int count = 0;
	float tmp_distance = FLT_MAX;
	while(count < CONSECUTIVE_MISS) {

		// find x and y values
		tmp->x = bounds[0] + (rand() % (int) ((bounds[1] - bounds[0]) * 100000) / 100000.0f);
		tmp->y = bounds[2] + (rand() % (int) ((bounds[3] - bounds[2]) * 100000) / 100000.0f);

		// compare with candidate PIA if point is in polygon
		if (Geometry_is_point_in_polygon(polygon, tmp->x, tmp->y)) {
			tmp_distance = Geometry_distance_between_point_and_polygon(polygon, tmp);
//			current_canvas->style = "stroke:blue;stroke-width:2;fill:blue";
//			Draw_circle(current_canvas, tmp->x, tmp->y, 10);
			if (tmp_distance > max_distance) {
				max_distance = tmp_distance;
				pia->x = tmp->x;
				pia->y = tmp->y;
				count = 0;
			} else {
				count++;
			}
		}
	}

	return pia;
}

// Inaccessibility of poles problem. Private function using LP method
Point* Geometry_find_biggest_inscribed_circle_in_polygon_LP(Polygon* polygon) {
	Matrix_2D* m = Matrix_2D_new(polygon->n + 1, 4);

	Edge* edge_tmp;
	int index = 0;
	for (edge_tmp = polygon->first_edge; edge_tmp != NULL; edge_tmp = edge_tmp->next) {

		// get two points slightly away from the edge perpendicularly to test
		// which side of the edge is facing the inside of the polygon
		Point* point_test = Point_new((edge_tmp->v1->x + edge_tmp->v2->x) / 2.0f,
			(edge_tmp->v1->y + edge_tmp->v2->y) / 2);
		Line* line_test = Geometry_perpendicular_line_at_point(
				edge_tmp->l, point_test->x, point_test->y);
		Point* p1 = Geometry_find_point_in_line_given_point_and_distance(line_test, point_test, 1.0f);
		Point* p2 = Geometry_find_point_in_line_given_point_and_distance(line_test, point_test, -1.0f);

		// when the face of the edge facing the inside of the polygon is above the edge's line
		if ((Geometry_is_point_in_polygon(polygon, p1->x, p1->y)
			&& (p1->y - edge_tmp->l->m * p1->x - edge_tmp->l->n) > 0)
		|| (Geometry_is_point_in_polygon(polygon, p2->x, p2->y)
			&& (p2->y - edge_tmp->l->m * p2->x - edge_tmp->l->n) > 0)) {
			Matrix_2D_set_row(m, index,
				edge_tmp->l->m,
				-1.0f,
				0.0f,
				-edge_tmp->l->n);
			Matrix_2D_scale_row(m, index, 1.0f / sqrtf(
					(edge_tmp->l->m * edge_tmp->l->m) + 1));
			Matrix_2D_set_cell(m, index++, 2, 1.0f);

		// when the face of the edge facing the inside of the polygon is below the edge's line
		} else if ((Geometry_is_point_in_polygon(polygon, p1->x, p1->y)
			&& (p1->y - edge_tmp->l->m * p1->x - edge_tmp->l->n) < 0)
		|| (Geometry_is_point_in_polygon(polygon, p2->x, p2->y)
			&& (p2->y - edge_tmp->l->m * p2->x - edge_tmp->l->n) < 0)) {
			Matrix_2D_set_row(m, index,
				edge_tmp->l->m,
				-1.0f,
				0.0f,
				-edge_tmp->l->n);
			Matrix_2D_scale_row(m, index, -1.0f / sqrtf(
					(edge_tmp->l->m * edge_tmp->l->m) + 1));
			Matrix_2D_set_cell(m, index++, 2, 1.0f);
		}

		// cleanup
		Point_destroy(p1);
		Point_destroy(p2);
		Point_destroy(point_test);
		Line_destroy(line_test);
	}


	// the objective function
	Matrix_2D_set_row(m, m->rows - 1, 0.0f, 0.0f, 1.0f, 0.0f);
//	Matrix_2D_print(m);

	// find solution
	Vector* v = LPSolver_simplex_solve(m);
	Point* point_pia = (v == NULL) ? NULL : Point_new(Vector_get(v, 0), Vector_get(v, 1));

	// clean up
	Vector_destroy(v);
	Matrix_2D_destroy(m);

	// return solution
	return point_pia;
}

// Inaccessibility of poles problem
Point* Geometry_find_biggest_inscribed_circle_in_polygon(Polygon* polygon, int method) {

	if (method == METHOD_LP)
		return Geometry_find_biggest_inscribed_circle_in_polygon_LP(polygon);

	float bounds[4] = { polygon->first_edge->v1->x, polygon->first_edge->v1->x,
		polygon->first_edge->v1->y, polygon->first_edge->v1->y };

	Edge* tmp_edge;
	for (tmp_edge = polygon->first_edge; tmp_edge != NULL; tmp_edge = tmp_edge->next) {
		if (tmp_edge->v1->x < bounds[0]) bounds[0] = tmp_edge->v1->x;
		if (tmp_edge->v1->x > bounds[1]) bounds[1] = tmp_edge->v1->x;
		if (tmp_edge->v1->y < bounds[2]) bounds[2] = tmp_edge->v1->y;
		if (tmp_edge->v1->y > bounds[3]) bounds[3] = tmp_edge->v1->y;
	}

	Point* point_pia = Point_new(0, 0);
	Point* point_tmp = Point_new(0, 0);
	float flt_tmp = FLT_MAX;
	int count = 1;
	while (count++) {

		// find new candidate PIA
		if (method == METHOD_SEQUENTIAL) {
			point_tmp =
				Geometry_find_biggest_inscribed_circle_in_polygon_sequential(polygon, bounds);
		} else if (method == METHOD_RANDOMIZED) {
			point_tmp =
				Geometry_find_biggest_inscribed_circle_in_polygon_randomized(polygon, bounds);
		}

		// update current PIA
		point_pia->x = point_tmp->x;
		point_pia->y = point_tmp->y;

		// update the bounds
		flt_tmp = (bounds[1] - bounds[0]) / (sqrtf(2) * 2);
		bounds[0] = point_pia->x - flt_tmp;
		bounds[1] = point_pia->x + flt_tmp;
		flt_tmp = (bounds[3] - bounds[2]) / (sqrtf(2) * 2);
		bounds[2] = point_pia->y - flt_tmp;
		bounds[3] = point_pia->y + flt_tmp;

		// check distance between upper and lower bounds
		if (bounds[1] - bounds[0] < 0.01 || bounds[3] - bounds[2] < 0.01) break;

//		printf("Candidate PIA: (%f,%f)\n", point_pia->x, point_pia->y);
	}

//	printf("Number of iterations: %d, ", count);
	Point_destroy(point_tmp);
	return point_pia;

}


void Geometry_test_line_1() {
	Line* l1 = Line_new_given_coordinates(-250, 0, 250, 250);
//	printf("%s[TEST] Line should be y = 0.5x + 125: ", timestamp());
	Line_print(l1);
	Draw* test_canvas = Draw_new("Geometry_test_line_1", NULL, 500, 500);
	Line_draw(l1);
	Line_destroy(l1);
	Draw_destroy(test_canvas);
}

void Geometry_test_polygon_1() {
	Edge* e0 = Edge_new_given_coordinates(150, 100, 200, 200);
	Edge* e1 = Edge_new_given_coordinates(200, 200, 0, 450);
	Edge* e2 = Edge_new_given_coordinates(0, 450, 150, 100);
	e0->next = e1;
	e1->next = e2;
	Polygon* polygon = Polygon_new_given_first_edge(e0);
	Draw* test_canvas = Draw_new("Geometry_test_polygon_1", NULL, 500, 500);
	Polygon_draw(polygon);
	Draw_destroy(test_canvas);
	Polygon_destroy(polygon);
}

void Geometry_test_polygon_2() {
	Vertice* v1 = Vertice_new(150.0f, 100.0f);
	Vertice* v2 = Vertice_new(200.0f, 200.0f);
	Vertice* v3 = Vertice_new(0.0f, 450.0f);
	v1->next = v2;
	v2->next = v3;
	Polygon* polygon = Polygon_new_given_first_vertice(v1, 1);
	Draw* test_canvas = Draw_new("Geometry_test_polygon_2", NULL, 500, 500);
	Polygon_draw(polygon);
	Draw_destroy(test_canvas);
	Polygon_destroy(polygon);
}

void Geometry_test_polygon_3() {
	Matrix_2D* m = Matrix_2D_new(3, 2);
	Matrix_2D_set_all(m,
			150.0f, 100.0f,
			200.0f, 200.0f,
			0.0f, 450.0f);
	Polygon* polygon = Polygon_new_given_matrix(m, 1);
	Draw* test_canvas = Draw_new("Geometry_test_polygon_3", NULL, 500, 500);
	Polygon_draw(polygon);
	Draw_destroy(test_canvas);
	Polygon_destroy(polygon);
}

void Geometry_test_find_point_in_line_1() {
	struct timeval time_seed;
	gettimeofday(&time_seed, NULL);
	srand(time_seed.tv_usec);

	Draw* test_canvas =
			Draw_new("Geometry_test_find_point_in_line_1", NULL, 1000, 1000);

	Point* p1 = Point_new(
			(rand() % (test_canvas->width - test_canvas->width / 2) *  100000) / 100000.0f,
			(rand() % (test_canvas->height *  100000) / 100000.0f));
	Point* p2 = Point_new(
			(rand() % (test_canvas->width - test_canvas->width / 2) *  100000) / 100000.0f,
			(rand() % (test_canvas->height *  100000) / 100000.0f));

	Line* line = Line_new_given_two_points(p1, p2);
	Point* point = Geometry_find_point_in_line_given_point_and_distance(line, p1, rand() % 50);

	// draw results
	Line_draw(line);
	test_canvas->style = "stroke:blue;fill:blue;stroke-width:1";
	Point_draw(p1);
	test_canvas->style = "stroke:red;fill:red;stroke-width:1";
	Point_draw(point);

	// cleanup
	Line_destroy(line);
	Point_destroy(p1);
	Point_destroy(p2);
	Point_destroy(point);
	Draw_destroy(test_canvas);
}

void Geometry_test_point_in_edge() {
	struct timeval time_seed;
	gettimeofday(&time_seed, NULL);
	srand(time_seed.tv_usec);

	Draw* test_canvas = Draw_new(
			"Geometry_test_point_in_edge", NULL, 1000, 1000);
	float bounds[4] = {	-test_canvas->width / 2.0f,
						 test_canvas->width / 2.0f,
						 0.0f, test_canvas->height };
	Polygon* polygon = Geometry_generate_random_convex_polygon(bounds, 4);
	Polygon_draw(polygon);

	int i;
	test_canvas->style = "stroke:red;fill:red;stroke-width:2";
	for (i = 0; i < 500; i++) {
		Point* p1 = Point_new(
				(rand() % (test_canvas->width - test_canvas->width / 2) *  100000) / 100000.0f,
				(rand() % (test_canvas->height *  100000) / 100000.0f));
		Point* p2 = Point_new(
				(rand() % (test_canvas->width - test_canvas->width / 2) *  100000) / 100000.0f,
				(rand() % (test_canvas->height *  100000) / 100000.0f));
		Line* line = Line_new_given_two_points(p1, p2);

		Vertice* v_tmp;
		Vertice* intersections = Geometry_intersections_between_polygon_and_line(polygon, line);
		for (v_tmp = intersections; v_tmp != NULL; v_tmp = v_tmp->next)
			Draw_circle(test_canvas, v_tmp->x, v_tmp->y, 3);

		// clean up
		Point_destroy(p1);
		Point_destroy(p2);
		Line_destroy(line);
		if (intersections != NULL) Vertice_destroy(intersections);
	}

	Polygon_destroy(polygon);
	Draw_destroy(test_canvas);
}

void Geometry_test_two_lines_intersection_1() {
	Line* l1 = Line_new_given_coordinates(-500, 0, 500, 500);
	Line* l2 = Line_new_given_coordinates(-500, 500, 500, 0);
	Point* i = Geometry_intersection_between_two_lines(l1, l2);
	Draw* test_canvas = Draw_new("Geometry_test_two_lines_intersection_1", NULL, 500, 500);
	Line_draw(l1);
	Line_draw(l2);
	test_canvas->style = "stroke:red;fill:red;stroke-width:2";
	Point_draw(i);

	// clean up
	Draw_destroy(test_canvas);
	Line_destroy(l1);
	Line_destroy(l2);
	Point_destroy(i);
}

void Geometry_test_two_lines_intersection_2() {
	Line* l1 = Line_new_given_m_and_n(1.01, 0);
	Line* l2 = Line_new_given_m_and_n(1.59, -50);
	Point* i = Geometry_intersection_between_two_lines(l1, l2);
//	printf("%s[TEST] Intersection should be approx (86,87): (%f,%f)\n", timestamp(), i->x, i->y);
	Draw* test_canvas = Draw_new("Geometry_test_two_lines_intersection_2", NULL, 500, 500);
	Line_draw(l1);
	Line_draw(l2);
	test_canvas->style = "stroke:red;fill:red;stroke-width:2";
	Point_draw(i);

	// clean up
	Draw_destroy(test_canvas);
	Line_destroy(l1);
	Line_destroy(l2);
	Point_destroy(i);
}

void Geometry_test_two_lines_intersection_3() {
	Line* l1 = Line_new_given_coordinates(100, 100, 100, 300);
	Line* l2 = Line_new_given_m_and_n(1.59, -50);
	Point* i = Geometry_intersection_between_two_lines(l1, l2);
	Draw* test_canvas = Draw_new("Geometry_test_two_lines_intersection_3", NULL, 500, 500);
	Line_draw(l1);
	Line_draw(l2);
	test_canvas->style = "stroke:red;fill:red;stroke-width:2";
	Point_draw(i);

	// clean up
	Draw_destroy(test_canvas);
	Line_destroy(l1);
	Line_destroy(l2);
	Point_destroy(i);
}

void Geometry_test_line_polygon_intersection_1() {
	struct timeval time_seed;
	gettimeofday(&time_seed, NULL);
	srand(time_seed.tv_usec);

	Draw* test_canvas =	Draw_new("Geometry_test_line_polygon_intersection_1", NULL, 1000, 1000);

	float bounds[4] = { -250.0f, 250.0f, 0.0f, 500.0f };
	Polygon* polygon = Geometry_generate_random_convex_polygon(bounds,4);
	Point* polygon_point = Point_new(polygon->first_edge->v1->x, polygon->first_edge->v1->y);

	int i;
	Point* random_points[7];
	for (i = 0; i < 7; i++) {
		random_points[i] = Point_new(
				rand() % test_canvas->width - test_canvas->width / 2,
				rand() % test_canvas->height);
	}

	Line* line1 = Line_new_given_two_points(random_points[0], random_points[1]);
	Line* line2 = Line_new_given_two_points(random_points[2], random_points[3]);
	Line* line3 = Line_new_given_two_points(random_points[4], random_points[5]);
	Line* line4 = Line_new_given_two_points(random_points[6], polygon_point);

	Polygon_draw(polygon);
	test_canvas->style = "stroke:blue;fill:none;stroke-width:2";
	Line_draw(line1);
	Line_draw(line2);
	Line_draw(line3);
	Line_draw(line4);

	// find intersections
	test_canvas->style = "stroke:red;fill:red;stroke-width:2";
	Vertice* i1 = Geometry_intersections_between_polygon_and_line(polygon, line1);
	Vertice* i2 = Geometry_intersections_between_polygon_and_line(polygon, line2);
	Vertice* i3 = Geometry_intersections_between_polygon_and_line(polygon, line3);
	Vertice* i4 = Geometry_intersections_between_polygon_and_line(polygon, line4);
	Vertice* i_list[4] = { i1, i2, i3, i4 };
	Vertice* vertice_tmp;

	int count[4] = {0,0,0,0};
	for (i = 0; i < 4; i++) {
		for (vertice_tmp = i_list[i]; vertice_tmp != NULL; vertice_tmp = vertice_tmp->next, count[i]++) {
			Point* point_tmp = Point_new(vertice_tmp->x, vertice_tmp->y);
			Point_draw(point_tmp);
			Point_destroy(point_tmp);
		}
	}
	printf("The counts are: %d, %d, %d and %d\n", count[0], count[1], count[2], count[3]);

	// clean up
	Draw_destroy(test_canvas);
	Polygon_destroy(polygon);
	Line_destroy(line1);
	Line_destroy(line2);
	Point_destroy(random_points[0]);
	Point_destroy(random_points[1]);
	Point_destroy(random_points[2]);
	Point_destroy(random_points[3]);
	Point_destroy(random_points[4]);
	Point_destroy(random_points[5]);
	Point_destroy(random_points[6]);
}

void Geometry_test_closest_edge() {
	struct timeval time_seed;
	gettimeofday(&time_seed, NULL);
	srand(time_seed.tv_usec);

	Draw* test_canvas =
		Draw_new("Geometry_test_closest_edge", NULL, 1000, 1000);

	float bounds[4] = {
		-test_canvas->width / 2.0f, test_canvas->width / 2.0f, 0, test_canvas->height };

	int i;
	Point* point_tmp = NULL;
	Vertice* triangle_vertices[3];
	for (i = 0; i < 4; i++) {
		float x = bounds[0] + (rand() % (int) ((bounds[1] - bounds[0]) * 100000) / 100000.0f);
		float y = bounds[2] + (rand() % (int) ((bounds[3] - bounds[2]) * 100000) / 100000.0f);
		if (i < 3) triangle_vertices[i] = Vertice_new(x, y);
		if (i < 3 && i > 0) triangle_vertices[i - 1]->next =  triangle_vertices[i];
		if (i == 3) {
			point_tmp = Point_new(x, y);
		}
	}
	Polygon* polygon = Polygon_new_given_first_vertice(triangle_vertices[0], 1);
	Polygon_draw(polygon);

	test_canvas->style = "stroke:blue;fill:blue;stroke-width:2";
	Draw_circle(test_canvas, point_tmp->x, point_tmp->y, 3);
	Edge* edge_tmp = Geometry_closest_polygon_edge_to_point(polygon, point_tmp);
	test_canvas->style = "stroke:red;fill:red;stroke-width:2";
	Draw_line(test_canvas, edge_tmp->v1->x, edge_tmp->v1->y, edge_tmp->v2->x, edge_tmp->v2->y);

	// clean up
	Draw_destroy(test_canvas);
	Polygon_destroy(polygon);
}

void Geometry_test_distance_from_point_to_edge_1() {
	struct timeval time_seed;
	gettimeofday(&time_seed, NULL);
	srand(time_seed.tv_usec);

	Draw* test_canvas =
		Draw_new("Geometry_test_distance_from_point_to_edge_1", NULL, 500, 500);
	Edge* edge_tmp = Edge_new_given_coordinates(
		rand() % test_canvas->width - test_canvas->width / 2,
		rand() % test_canvas->height,
		rand() % test_canvas->width - test_canvas->width / 2,
		rand() % test_canvas->height);
	Point* point_tmp =
		Point_new(rand() % test_canvas->width - 50, rand() % test_canvas->height);
	float distance = Geometry_distance_between_point_and_edge(edge_tmp, point_tmp);
	char label[22];
	sprintf(label, "Distance: %f", distance);

	// draw everything
	Edge_draw(edge_tmp);
	Point_draw(point_tmp);
	Draw_text(test_canvas, -test_canvas->width / 2, test_canvas->height - 20, 14, label);

	// clean up
	Draw_destroy(test_canvas);
	Edge_destroy(edge_tmp);
	Point_destroy(point_tmp);
}

Polygon* get_polygon() {
	Matrix_2D* m = Matrix_2D_new(4, 2);
	Matrix_2D_set_all(m,
			-150.915985,484.934998,
			4.669006,392.812012,
			76.716980,713.557983,
			-105.128998,534.853027);
	Polygon* polygon = Polygon_new_given_matrix(m, 1);
	return polygon;
}

void Geometry_test_point_in_polygon() {
	struct timeval time_seed;
	gettimeofday(&time_seed, NULL);
	srand(time_seed.tv_usec);

	Draw* test_canvas = Draw_new("Geometry_test_point_in_polygon_1", NULL, 1000, 1000);
	float bounds[4] = { -500.0f, 500.0f, 0.0f, 1000.0f };
	Polygon* polygon = Geometry_generate_random_convex_polygon(bounds, 4);
	Polygon_draw(polygon);

	Point* point_tmp = Point_new(0, 0);
	int i;
	for (i = 0; i < 10000; i++) {
		point_tmp->x = (rand() % test_canvas->width) - test_canvas->width / 2;
		point_tmp->y = rand() % test_canvas->height;
		test_canvas->style = Geometry_is_point_in_polygon(polygon, point_tmp->x, point_tmp->y) ?
				"stroke:blue;stroke-width:2;fill:blue" : "stroke:red;stroke-width:2;fill:red";
		Draw_circle(test_canvas, point_tmp->x, point_tmp->y, 2);
	}

	Point* point_pia = Geometry_find_biggest_inscribed_circle_in_polygon(polygon, METHOD_SEQUENTIAL);
	test_canvas->style = "stroke:yellow;stroke-width:10;fill:yellow";
	Draw_circle(test_canvas, point_pia->x, point_pia->y, 3);
	test_canvas->style = "stroke:yellow;stroke-width:10;fill:none";
	Draw_circle(test_canvas, point_pia->x, point_pia->y,
			Geometry_distance_between_point_and_polygon(polygon, point_pia));

	// clean up
	Draw_destroy(test_canvas);
	Polygon_destroy(polygon);
	Point_destroy(point_tmp);
}

