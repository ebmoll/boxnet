/*
	Copyright 2012 Samuel Moll
	License: LGPLv3

	Broadphase 2D collision detection using the boxnet algorithm
*/


// TODO: tune these numbers correctly, or replace the
// simplistic memory allocation with something better.
//
// These are only numbers for incremental memory allocation
// ("vector_append" macro in boxnet.c)
// e.g. BOXES_SIZE_INIT 100 means initially there will be
// 100 boxes and when they are used up, another 100 will be
// allocated, and so on.
#define BOXES_SIZE_INIT 100
#define COLLISIONS_SIZE_INIT 200
#define REPAIR_QUEUE_INIT 100
#define BC_QUEUE_SIZE_INIT 40



struct Box;
struct Junction;
struct RepairQueue;

typedef struct Junction {
	struct Junction*	nb[4];		// neighbors; can be Null
	struct Box*			pos[2];		// posx and posy; points to struct Box
	unsigned char		dir;		// 0=up, 1=left, 2=down, 3=right,
									// 4=cross, 5=infinity/disconnected
	unsigned char		beamdir;	// direction of non-terminating beam
	unsigned char		enqueued;	// marker for the repair queue
} Junction;

typedef struct Box {
	struct Junction		jnc;
	struct Junction		rayend[4];
	double				posx;		// == left
	double				posy;		// == bottom
	double				right;		// bounding box goes from posx to right
	double				top;		// and from posy to top
	void*				usrdata;	// user pointer; normally points
									// to user-defined object
	struct Box*			marked;
} Box;

typedef struct Boxnet {
	struct Box**		boxes;
	int					boxes_size;
	int					boxes_size_max;
} Boxnet;

typedef void (*collisionCallback)(void* obj1, void* obj2, void* data);


Boxnet* Boxnet_new();
void Boxnet_free(Boxnet* net);
Box* Boxnet_addbox(Boxnet* net, double x, double y,
							double right, double top,
							Box* near, void* usrdata);
void Boxnet_delbox(Boxnet* net, Box* box);
void Boxnet_delbox_byusrdata(Boxnet* net, void* usrdata);
void Boxnet_collide(Boxnet* net, collisionCallback func, void* data);

