/*
 * Author: Oscar Martinez
 * Copyright (c) 2012 All rights reserved.
 */

// library includes
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <stdint.h>

// module includes
#include "../Math/LPSolver.h"
#include "../Math/Geometry.h"
#include "../SVG/Draw.h"




Polygon* get_polygon_1() {
	Matrix_2D* m = Matrix_2D_new(3, 2);
	Matrix_2D_set_all(m,
			-50.0,	0.0,
			0.0,	100.0,
			50.0,	0.0);
	Polygon* polygon = Polygon_new_given_matrix(m, 1);
	return polygon;
}

Polygon* get_polygon_2() {
	Matrix_2D* m = Matrix_2D_new(6, 2);
	Matrix_2D_set_all(m,
			-2500.0,	0.0,
			2500.0,	1000.0,
			5000.0,	5000.0,
			2500.0,	10000.0,
			-2500.0,	9000.0,
			-5000.0,	5000.0);
	Polygon* polygon = Polygon_new_given_matrix(m, 1);
	return polygon;
}

Polygon* get_polygon_3() {
	Matrix_2D* m = Matrix_2D_new(3, 2);
	Matrix_2D_set_all(m,
			-25.0,	30.0,
			25.0,	10.0,
			20.0,	75.0);
	Polygon* polygon = Polygon_new_given_matrix(m, 1);
	return polygon;
}

Polygon* get_polygon_4() {
	Matrix_2D* m = Matrix_2D_new(6, 2);
	Matrix_2D_set_all(m,
			-25.0,	0.0,
			25.0,	10.0,
			50.0,	50.0,
			25.0,	100.0,
			-25.0,	80.0,
			-45.0,	30.0);
	Polygon* polygon = Polygon_new_given_matrix(m, 1);
	return polygon;
}

Polygon* get_polygon_5() {
	Matrix_2D* m = Matrix_2D_new(4, 2);
	Matrix_2D_set_all(m,
			-150.915985,484.934998,
			4.669006,392.812012,
			76.716980,713.557983,
			-105.128998,534.853027);
	Polygon* polygon = Polygon_new_given_matrix(m, 1);
	return polygon;
}

int main(int argc, char *argv[]) {

//	Geometry_test_point_in_polygon();
//	Geometry_test_two_lines_intersection_3();
//	Geometry_test_point_in_edge();
//	Geometry_test_line_polygon_intersection_1();
//	Geometry_test_closest_edge();
//	return 0;

	int method = 0;
	sscanf(argv[1], "%d", &method);

	// generate numbers
	if (method == -1) {
		char buf1[128], buf2[128];
		FILE* fd1 = fopen(argv[2],"r");
		FILE* fd2 = fopen(argv[3],"r");
		float x1, x2, y1, y2;
		int good = 0, error_01 = 0, error_10 = 0, error_x = 0;
		while (fgets(buf1, 128, fd1) != NULL && fgets(buf2, 128, fd2) != NULL) {
			sscanf(buf1, "(%f,%f)", &x1, &y1);
			sscanf(buf2, "(%f,%f)", &x2, &y2);
			float distance = Geometry_distance_between_coordinates(x1, y1, x2, y2);
			if (distance == 0.0f) good++;
			else if (distance <= 0.1f) error_01++;
			else if (distance <= 1.0f) error_10++;
			else if (distance > 1.0f) error_x++;
		}
//		printf("Good: %d\nError <= 0.01%%: %d\nError <= 0.1%%: %d\nError > 0.1%%: %d\n",
//				good, error_01, error_10, error_x);
		printf("%d,%d,%d,%d\n",
				good, error_01, error_10, error_x);

		return 0;
	}

	// generate polygons and write them into a file
	if (method == 0) {
		int i;
		char buf[128];
		FILE* outfile = fopen("polyfile","w+");
		float bounds[4] = { 0.0f, 5000.0f, 0.0f, 5000.0f };
		for (i = 0; i < 500; i++) {
			Polygon* polygon = Geometry_generate_random_convex_polygon(bounds,3);
			sprintf(buf, "%s\n", Polygon_to_string(polygon));
			fputs(buf, outfile);
		}
		return 0;
	}

	// read polygon from file
//	int count = 1;
	Vertice* v1 = Vertice_new(0,0);
	Vertice* v2 = Vertice_new(0,0);
	Vertice* v3 = Vertice_new(0,0);
	v1->next = v2;
	v2->next = v3;
	char buf[128];
//	char canvasname[64];
	FILE* fd = fopen(argv[2],"r");
	while (fgets(buf, 128, fd) != NULL) {
//		sprintf(canvasname, "biggest_circle_%03d_%d.svg", count++, method);
//		Draw* canvas = Draw_new(canvasname, NULL, 10000, 10000);
		sscanf(buf, "(%f,%f) (%f,%f) (%f,%f)",
			&v1->x, &v1->y, &v2->x, &v2->y, &v3->x, &v3->y);
		Polygon* polygon = Polygon_new_given_first_vertice(v1, 1);
		Point* point_pia = Geometry_find_biggest_inscribed_circle_in_polygon(polygon, method);

		if (point_pia == NULL) {
			printf("Error finding solution\n");

		} else {

			// output
			printf("(%.2f,%.2f)\n", point_pia->x, point_pia->y);
//			canvas->style = "stroke:black;stroke-width:10;fill:none";
//			Polygon_draw(polygon);
//			canvas->style = "stroke:blue;stroke-width:10;fill:blue";
//			Draw_circle(canvas, point_pia->x, point_pia->y, 3);
//			canvas->style = "stroke:blue;stroke-width:10;fill:none";
//			Draw_circle(canvas, point_pia->x, point_pia->y,
//					Geometry_distance_between_point_and_polygon(polygon, point_pia));
		}

		// cleanup
//		Polygon_destroy(polygon);
		Point_destroy(point_pia);
//		Draw_destroy(canvas);
	}

	return 0;

//	Draw* canvas = Draw_new("biggest_circle.svg", "stroke:black;stroke-width:10;fill:none", 5000, 5000);
//	float bounds[4] = { -2500.0f, 2500.0f, 0.0f, 5000.0f };
//	Polygon* polygon = Geometry_generate_random_convex_polygon(bounds,3);
//	Polygon_draw(polygon);
//
//	Point* point_pia = Geometry_find_biggest_inscribed_circle_in_polygon(polygon, method);
//
//	canvas->style = "stroke:yellow;stroke-width:10;fill:yellow";
//	Draw_circle(canvas, point_pia->x, point_pia->y, 3);
//	canvas->style = "stroke:yellow;stroke-width:10;fill:none";
//	Draw_circle(canvas, point_pia->x, point_pia->y,
//			Geometry_distance_between_point_and_polygon(polygon, point_pia));
//
//	printf("The center is: (%f, %f) and it's %f from the closest edge\n", point_pia->x, point_pia->y,
//			Geometry_distance_between_point_and_polygon(polygon, point_pia));
//
//	Draw_destroy(canvas);
//	printf("(%f,%f)\n", point_pia->x, point_pia->y);
	return 0;
}
