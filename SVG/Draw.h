/*
 * Author: Oscar Martinez
 * Copyright (c) 2012 All rights reserved.
 */

#ifndef DRAW_H_
#define DRAW_H_

// structure definitions

typedef struct Draw Draw;
struct Draw {
    FILE* svgfile;
    char* style;
    int width;
    int height;
};


// constructors and destructors

Draw* Draw_new(char* filepath, char* style, int width, int height);
void Draw_destroy(Draw* draw);


// module-oriented functions

void Draw_text(Draw* draw, float x, float y, int size, char* text);
void Draw_circle(Draw* draw, float cx, float cy, float r);
void Draw_arc(Draw* draw, float x0, float y0, float x1, float y1, float r, int flag);
void Draw_line(Draw* draw, float x1, float y1, float x2, float y2);
void Draw_bezier_curve(Draw* draw, int x0, int y0, int x1, int y1, int x2, int y2, int x, int y);

// unit testing
void Draw_test_arc_1(void);
void Draw_test_arc_2(void);



// global variables

Draw* current_canvas;

#endif /* DRAW_H_ */
