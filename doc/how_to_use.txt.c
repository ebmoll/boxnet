/*
	Copyright 2012 Samuel Moll
	License: AGPLv3
*/

// Everything you need to know to use boxnet in your C application
// ===============================================================

#include "boxnet.h"
#include "stdio.h"
int main(void) {

// Let's say, you want to collide Circles defined by:

    typedef struct Circle {
    	char* name;   //name of the circle
        double x,y,r; //position and radius!
        Box* box;     //boxnet bounding box!
    } Circle;

    Circle circle1 = { .name="c1", .x=-4, .y=0, .r=1.5 };
    Circle circle2 = { .name="c2", .x=4, .y=0, .r=3   };
    Circle circle3 = { .name="c3", .x=0, .y=6, .r=5   };

// Now that you have 3 circles at different positions, you want
// to efficiently find out which ones are overlapping using boxnet!
// First, create a new 2D collision space by calling

    Boxnet* my_space = Boxnet_new();

// Then, start adding colliding boxes for every circle:

    circle1.box = Boxnet_addbox(my_space,
        circle1.x - circle1.r, circle1.y - circle1.r,
        circle1.x + circle1.r, circle1.y + circle1.r,
        NULL, (void*) &circle1);
    circle2.box = Boxnet_addbox(my_space,
        circle2.x - circle2.r, circle2.y - circle2.r,
        circle2.x + circle2.r, circle2.y + circle2.r,
        NULL, (void*) &circle2);
    circle3.box = Boxnet_addbox(my_space,
        circle3.x - circle3.r, circle3.y - circle3.r,
        circle3.x + circle3.r, circle3.y + circle3.r,
        NULL, (void*) &circle3);

// As you see, we pass the extents of the bounding box of the circle
// to Boxnet_addbox() (x-r to x+r and y-r to y+r), then a NULL argument
// (see the full docs for what it does, but you can always use NULL if
// you want) and finally a user pointer, so boxnet can make the connection
// to your circles.

// To start colliding, you first need to have a collision handling
// function. This function gets called by boxnet when it detects
// overlapping bounding boxes. This function will have to check
// if the circles really do overlap (since the bounding box is
// larger than the circles), and then do some cool collision
// response.
// If you need to delete a colliding object in response to a
// collision, just flag it for deletion and delete it later.
// Never call Boxnet_delbox() or Boxnet_delbox_byusrdata() inside
// the callback!
// The callback looks like this:

    void collision(void* obj1, void* obj2, void* data) {
        Circle c1 = *(Circle*)obj1;
        Circle c2 = *(Circle*)obj2;
        // check if the two circles do really overlap
        double RR = (c1.x-c2.x)*(c1.x-c2.x) + (c1.y-c2.y)*(c1.y-c2.y);
        if ( (c1.r+c2.r)*(c1.r+c2.r) < RR )
            return; // no collision
        // Now do something in response to the collision
        printf( "%s collided with %s!\n", c1.name, c2.name);
        printf( "additional data was \"%s\"\n", (char*) data);
    }

// the *data pointer can be used to have access to additional
// data structures that you might need for collision response.

// Now you can call the broadphase collision detection algorithm!
// This will call the collision callback with all circles
// that have overlapping bounding boxes, so in this case it will
// do a bunch of "printf"'s telling you about the collisions.

    Boxnet_collide(my_space, collision, "nothing");

// Now move your objects again, update the bounding boxes
// like above, and call Boxnet_collide() again, and so on.
// If you move your circles around, make sure to update the bounding box
// information like so:

// move the circle 4 units to the right
	printf("\nmoving circle1...\n\n");
    circle1.x += 4;
// now update bounding box information
    circle1.box->posx  = circle1.x - circle1.r;
    circle1.box->posy  = circle1.y - circle1.r;
    circle1.box->right = circle1.x + circle1.r;
    circle1.box->top   = circle1.y + circle1.r;
// and collide again, this will be very fast if all objects
// only moved a little since the last call to Boxnet_collide()
    Boxnet_collide(my_space, collision, "second time step!");

// To delete objects from the boxnet, do
    Boxnet_delbox(my_space,circle1.box);
// or, alternatively delete by the usrdata pointer:
    Boxnet_delbox_byusrdata(my_space,&circle2);

// To free all memory that was allocated by the boxnet algorithm, call:
    Boxnet_free( my_space );
	
	return 0;
}



