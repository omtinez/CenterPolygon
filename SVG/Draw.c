/*
 * Author: Oscar Martinez
 * Copyright (c) 2012 All rights reserved.
 */

// library includes
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// module includes
#include "Draw.h"

char drawbuffer[200];

Draw* Draw_new(char* filepath, char* style, int width, int height) {
    Draw* this = malloc(sizeof(Draw));
    this->svgfile = fopen(filepath,"w+");
    this->style = (style == NULL) ? "stroke:black;stroke-width:2" : style;
    this->height = height;
    this->width = width;

    sprintf(drawbuffer,"<svg width=\"%d\" height=\"%d\">\n", width, height);
    fputs(drawbuffer, this->svgfile);
    current_canvas = this;
    return this;
}

void Draw_destroy(Draw* draw) {
	fputs("</svg>", draw->svgfile);
    fclose(draw->svgfile);
    free (draw);
}

void Draw_text(Draw* draw, float x, float y, int size, char* text) {
    sprintf(drawbuffer,"<text x=\"%f\" y=\"%f\" font-size=\"%dpx\" style=\"%s\">%s</text>\n",
		x + draw->width / 2, draw->width - y, size, draw->style, text);
    fputs(drawbuffer, draw->svgfile);
}

void Draw_circle(Draw* draw, float cx, float cy, float r) {
    sprintf(drawbuffer,"<circle cx=\"%f\" cy=\"%f\" r=\"%f\" style=\"%s\"/>",
    		cx + draw->width/2, draw->width-cy, r, draw->style);
	fputs(drawbuffer, draw->svgfile);
}

void Draw_arc(Draw* draw, float x0, float y0, float x1, float y1, float r, int flag) {
    sprintf(drawbuffer,"<path d=\"M%f,%f A%f,%f 0 0,%d %f,%f\" style=\"%s\"/>",
    		x0 + draw->width/2, draw->width - y0,
    		r, r, flag,
    		x1 + draw->width/2, draw->width - y1,
    		draw->style);
	fputs(drawbuffer, draw->svgfile);
}

void Draw_line(Draw* draw, float x1, float y1, float x2, float y2) {
    sprintf(drawbuffer,"<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" style=\"%s\"/>",
    		x1 + draw->width/2, draw->width - y1, x2 + draw->width/2, draw->width - y2, draw->style);
    fputs(drawbuffer, draw->svgfile);
}

void Draw_bezier_curve(Draw* draw, int x0, int y0, int x1, int y1, int x2, int y2, int x, int y) {
	sprintf(drawbuffer,"<path d=\"M%d,%d C%d,%d %d,%d %d,%d\" style=\"%s\"/>",
			x0,500-y0,x1,500-y1,x2,500-y2,x,500-y,draw->style);
	fputs(drawbuffer, draw->svgfile);
}

void Draw_test_arc_1() {
	Draw* draw = Draw_new("Draw_test_arc_1", NULL, 500, 500);
	Draw_text(draw, -50, 400, 24, "sweep-flag: 1");
	draw->style = "stroke:black;stroke-width:5;fill:none";
	Draw_arc(draw, 0, 50, 100, 200, 300, 1);
	Draw_destroy(draw);
}

void Draw_test_arc_2() {
	Draw* draw = Draw_new("Draw_test_arc_2", NULL, 500, 500);
	Draw_text(draw, -50, 400, 24, "sweep-flag: 0");
	draw->style = "stroke:black;stroke-width:5;fill:none";
	Draw_arc(draw, 0, 50, -100, 200, 300, 0);
	Draw_destroy(draw);
}

