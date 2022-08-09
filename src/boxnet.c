/*
	Copyright 2012 Samuel Moll
	License: AGPLv3
	
	Broadphase 2D collision detection using the boxnet algorithm.
*/

/*****************
	TODO: strictly seperate the net structure from the spatial
	information. May simplify the code quite a bit.
	 - in Junction_flip(), there are still spatial comparisions.
	   this is necessary, because Junction_flip() is used after
	   Boxnet_repair() and before collisions are found to ensure
	   the lower ray of each box ends after the right edge...
	
	TODO: add support for repairing after only a few bodies have
	moved; this would probably imply a box_update() function or similar
	that also adds the box to the repair queue or so.
	
	TODO: implement a segment query, a proper raycast and a BB query.
	this might require reactivating (and probably rewriting) the
	navigate() function :|
	
	TODO: devise a good solution for static geometry...
	
	TODO:
	reduce memory footprint by only saving one pos[] per Junction
	(for dir=0 and dir=2 the x-Position)
	maybe remove beamdir / merge into dir
	
	TODO: profiling and performance optimizations
	
	TODO: analyze repair algorithm and apply heuristics to reduce
	total number of slide and flip operations.
	
	TODO: maybe really interpret Junction.dir==5 as "infinity"
	and also connect the ends of rays to those; might simplify
	code by getting rid of all those "if(next==NULL)"s.
	Maybe split dir=5 into 4 different values; maybe introduce
	four special junctions for the infinity corners?
	
*****************/

/* comment out for debugging */
//#define NDEBUG
#define NOTEST


#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "boxnet.h"

// debugging functions
static void validate(Boxnet* net);
static int repair_check(Boxnet* net);


struct Connection {
	struct Junction*	jnc;
	unsigned char		tdir;
};

typedef struct RepairQueue {
	struct Connection*	queue;
	int					size;
	int					size_max;
} RepairQueue;



static Junction* Junction_flip(Junction* jnc, struct RepairQueue* queue);
static void detach(Junction* jnc);


/*
	vector handling function...
*/
#define vector_append(v, append, v_size,\
					 v_size_max, v_inc) {\
	if((v_size) == (v_size_max)) {\
		/* TODO: add error checking (NULL pointer)*/\
		v_size_max+=v_inc;\
		v = realloc(v, (v_size_max)*sizeof *(v));\
		assert((v)!=NULL);\
	}\
	(v)[v_size]=append;\
	(v_size)++;\
}



// TODO: maybe allocate boxes in a continuous array for better
//       cache coherence
//       write a box relocation function (swap memory location
//       of two boxes; all links to Junction members of the box
//       must be changed too then
static Box* Box_new() {
	Box* new = malloc(sizeof *new);
	assert(new!=NULL);
	new->jnc.dir = 4; // those never change...
	new->jnc.pos[0] = new;
	new->jnc.pos[1] = new;
	new->jnc.enqueued = 0;
	for(int d=0;d<4;d++) {
		new->rayend[d].pos[d%2] = new;
		new->rayend[d].enqueued = 0;
		new->rayend[d].dir = 5;
	}
	return new;
}


/*
	removes all associated junctions from the net
	and frees memory for box.
*/
static void Box_free(Box* box) {
	assert(box!=NULL);
	/* disconnect the associated junction from the net */
	for(int d=0;d<4;d++) {
		Junction* jnc = box->jnc.nb[d];
		if(jnc!=NULL && jnc->dir!=(d^2))
			Junction_flip(jnc,NULL);
	}
	for(int d=0;d<4;d++) {
		Junction* jnc = &box->rayend[d];
		if(jnc->dir!=5)
			detach(jnc);
	}
	free(box);
}


static RepairQueue RepairQueue_new() {
	RepairQueue	q;
	q.queue = malloc(REPAIR_QUEUE_INIT * sizeof *q.queue);
	assert(q.queue!=NULL);
	q.size_max=REPAIR_QUEUE_INIT;
	return q;
}

static void RepairQueue_append(Junction* jnc, unsigned char tdir, RepairQueue* q) {
	if((jnc->enqueued & (1<<tdir)) == 0) {
		struct Connection conn;
		conn.jnc = jnc;
		conn.tdir = tdir;
		vector_append(q->queue, conn, q->size, q->size_max, REPAIR_QUEUE_INIT);
		jnc->enqueued ^= (1<<tdir);
	}
}



Boxnet* Boxnet_new() {
	Boxnet* new = malloc(sizeof *new);
	assert(new!=NULL);
	new->boxes = malloc(BOXES_SIZE_INIT * sizeof *new->boxes);
	assert(new->boxes!=NULL);
	new->boxes_size_max = BOXES_SIZE_INIT;
	new->boxes_size = 0;
	return new;
}

/*
	recursively deletes the boxnet and all associated
	Box and Junction structures.
*/
void Boxnet_free(Boxnet* net) {
	for(int i=0;i<net->boxes_size;i++) {
		Box_free(net->boxes[i]);
	}
	free(net->boxes);
	free(net);
}

#ifdef COMMENT_THIS_OUT
/*
	returns the direction of (px,py) in relation to jnc.
	0=upper left, 1=lower left, 2=lower right, 3=upper right
	(with positive x right, positive y up)
*/
static unsigned char quadrant(Junction* jnc, double px, double py) {
	return ((px > jnc->pos[0]->posx)<<1) ^ ((py < jnc->pos[1]->posy) ^
	(px > jnc->pos[0]->posx));
}


/*
	find a Junction that is guaranteed to be a neighbor of (x,y),
	that is a junction that lies on the edges of the surronding
	rectangle.
	Is more or less only needed when inserting new Junctions.
	
	this function was deactivated and is not currently used;
	new boxes are just inserted wherever and the boxnet repair
	function takes care of the rest.
	
	TODO:
	maybe rewrite this function and use it again as an
	optimization for inserting points.
	
	this function has a bug (occuring rarely) that might have
	something to do with junctions having the same coordinates.
	
	since "navigate" is currently not used, the bug is not showing.
*/
static Junction* navigate(Junction* start, double x, double y) {
	Junction* cur = start;
	// choose initial direction
	unsigned char quadr = quadrant(cur,x,y);
	unsigned char curdir = quadr;
	if((curdir^2) == cur->dir)
		curdir = (curdir+1)%4;
	unsigned char found = 0; // 1=pass, 2=pass+change
	while(1) {
		Junction* next = cur->nb[curdir];
		// are we going to infinity?
		if(next==NULL) {
			// counts as pass
			if(found==2)
				break;
			found=1;
			curdir=curdir^2; // backtrack...
		} else {
			cur = next;
		}
		// try to change direction
		unsigned char quadrnew = quadrant(cur,x,y);
		if(quadrnew!=quadr) {
			if(found==2)
				break;
			found=1;
			quadr = quadrnew;
		}
		// turn to point if possible
		unsigned char newdir = (4+quadr-curdir)%4 < 2 ? (curdir+1)%4 : (curdir+3)%4;
		if((newdir^2) != cur->dir) {
			found = found==1 ? 2 : 0; // reset if change+change
			curdir = newdir;
		}
	}
	return cur;
}
#endif


/*
	Inserts a Junction jnc, that must be of type dir=4,
	adjacent to start into the Junction network. start has to be a member of
	the network.
*/
static void Junction_insert(Junction* jnc, Junction* start) {
	assert(start!=NULL);
	assert(jnc->dir==4);
	assert(jnc->pos[0]==jnc->pos[1]);
	assert(start->dir<=4);
	int inserted = 0;
	// initialize nb's...
	jnc->nb[0] = jnc->nb[1] = jnc->nb[2] = jnc->nb[3] = NULL;
	unsigned char initdir = start->dir==4 ? 0 : start->dir;
	// try clockwise and counterclockwise
	for(unsigned char cwccw=1; cwccw<4; cwccw+=2) {
		Junction* cur = start;
		// choose initial direction
		unsigned char curdir;
		if(cwccw==1)
			curdir = initdir;
		else
			curdir = (initdir+1)%4;
		while(inserted!=4) {
			Junction* next = cur->nb[curdir];
			void insert() {
				// determine new beamdir
				Junction* newjnc = &jnc->pos[0]->rayend[(curdir+cwccw)%4];
				newjnc->dir = (curdir+cwccw)%4;
				if(cur->dir==4) {
					newjnc->beamdir = curdir;
				} else {
					if(curdir==cur->dir)
						newjnc->beamdir = curdir^2;
					else
						newjnc->beamdir = cur->beamdir;
				}
				if(curdir%2==0)
					newjnc->pos[0] = cur->pos[0];
				else
					newjnc->pos[1] = cur->pos[1];
				newjnc->nb[curdir^2] = cur;
				cur->nb[curdir] = newjnc;
				newjnc->nb[curdir] = next;
				if(next!=NULL)
					next->nb[curdir^2] = newjnc;
				newjnc->nb[(curdir+cwccw)%4] = jnc;
				jnc->nb[(curdir+cwccw+2)%4] = newjnc;
				inserted++;
			}
			if(next==NULL) {
				// insert and continue with other cwccw direction
				insert();
				break;
			}
			// change direction if possible
			if((curdir+cwccw+2)%4 != next->dir) {
				insert();
				curdir = (curdir+cwccw)%4;
			}
			cur=next;
		}
	}
}


/*
	removes a T-Junction from the net; acts as if the
	ending ray was removed
	the caller is responsible for clearing up the net...
*/
static void detach(Junction* jnc) {
	assert(jnc->dir<4);
	Junction* next = jnc->nb[jnc->beamdir];
	Junction* prev = jnc->nb[jnc->beamdir^2];
	prev->nb[jnc->beamdir] = next;
	if(next!=NULL)
		next->nb[jnc->beamdir^2] = prev;
	jnc->dir=5;
}

/*
	flips this Junction. jnc->dir must be != 4.
	if queue==NULL, enqueuing is not used
	
	flipping means that at the T-junction, the ray that
	previously stopped now continues, and the other ray
	stops there.
	
	^
	|
	|----      ==>    <-------
	|                     |
*/
static Junction* Junction_flipone(Junction* jnc, RepairQueue* queue);
static Junction* Junction_flip(Junction* jnc, RepairQueue* queue) {
	Junction* cur = jnc;
	Junction* next = cur->nb[jnc->beamdir];
	while(next != NULL && next->dir != (jnc->beamdir^2)) {
		cur = next;
		next = cur->nb[jnc->beamdir];
	}
	while(1) {
		if(cur==jnc)
			return Junction_flipone(cur,queue);
		cur = Junction_flipone(cur,queue);
		cur = cur->nb[(jnc->beamdir)^2];
	}
}
/*
static Junction* Junction_flipone_rewrite(Junction* jnc, RepairQueue* queue) {
	Junction* detach = jnc->nb[jnc->beamdir];
	// find position of new Junction
	Junction* insert = jnc->nb[jnc->beamdir^2];
	while(insert->dir==jnc->dir)
		insert = insert->nb[jnc->beamdir^2];
	insert = insert->nb[jnc->dir^2];
	while(insert!=NULL && insert->dir==(jnc->beamdir^2))
		insert = insert->nb[jnc->dir^2];
	//while(insert!=NULL && insert
		
}*/

static Junction* Junction_flipone(Junction* jnc, RepairQueue* queue) {
	Junction* next = jnc->nb[jnc->beamdir];
	if(next != NULL) {
		if(queue!=NULL) {
			RepairQueue_append(next->nb[next->beamdir^2], next->beamdir, queue);
		}
		detach(next);
	}
	// flip
	Junction* flipped;
	if(jnc->dir%2==0) {
		flipped = &jnc->pos[1]->rayend[jnc->beamdir^2];
		assert(flipped->pos[1]==jnc->pos[1]);
		flipped->pos[0] = jnc->pos[0];
	} else {
		flipped = &jnc->pos[0]->rayend[jnc->beamdir^2];
		assert(flipped->pos[0]==jnc->pos[0]);
		flipped->pos[1] = jnc->pos[1];
	}
	flipped->dir = jnc->beamdir^2;
	flipped->beamdir = jnc->dir^2;
	flipped->nb[flipped->dir] = jnc->nb[flipped->dir];
	jnc->nb[flipped->dir]->nb[flipped->dir^2] = flipped;
	flipped->nb[jnc->dir] = jnc->nb[jnc->dir];
	jnc->nb[jnc->dir]->nb[flipped->beamdir] = flipped;
	jnc->dir = 5;
	// reconnect loose end
	Junction* cur = flipped->nb[flipped->dir];
	while(cur->dir==(flipped->beamdir^2))
		cur = cur->nb[flipped->dir];
	// search in direction of new beamdir
	do {
		cur = cur->nb[flipped->beamdir];
		if(cur==NULL){
			flipped->nb[flipped->beamdir] = NULL;
			return flipped;
		}
	} while(cur->dir == flipped->dir);
	// find new Junction
	next = cur->nb[flipped->dir^2];
	while(next != NULL && next->dir == flipped->beamdir) {
		if(flipped->dir==0) {
			if(flipped->pos[1]->posy > next->pos[1]->posy) break;
		} else if (flipped->dir==1) {
			if(flipped->pos[0]->posx < next->pos[0]->posx) break;
		} else if (flipped->dir==2) {
			if(flipped->pos[1]->posy < next->pos[1]->posy) break;
		} else if (flipped->dir==3) {
			if(flipped->pos[0]->posx > next->pos[0]->posx) break;
		}
		cur = next;
		next = cur->nb[flipped->dir^2];
	}
	// insert new intersection and connect it
	Junction* newjnc;
	if(flipped->beamdir%2==0) { // up or down
		newjnc = &flipped->pos[0]->rayend[flipped->beamdir^2];
		assert(newjnc->pos[0]==flipped->pos[0]);
		newjnc->pos[1] = cur->pos[1];
	} else {
		newjnc = &flipped->pos[1]->rayend[flipped->beamdir^2];
		assert(newjnc->pos[1]==flipped->pos[1]);
		newjnc->pos[0] = cur->pos[0];
	}
	newjnc->dir = flipped->beamdir^2;
	if(cur->dir==4)
		newjnc->beamdir = flipped->dir^2;
	else if(cur->dir == (flipped->dir^2))
		newjnc->beamdir = flipped->dir;
	else
		newjnc->beamdir = cur->beamdir;
	
	newjnc->nb[flipped->beamdir^2] = flipped;
	flipped->nb[flipped->beamdir] = newjnc;
	newjnc->nb[flipped->dir^2] = next;
	if(next!=NULL)
		next->nb[flipped->dir] = newjnc;
	newjnc->nb[flipped->dir] = cur;
	cur->nb[flipped->dir^2] = newjnc;
	if(queue!=NULL) {
		RepairQueue_append(newjnc, newjnc->beamdir, queue);
		RepairQueue_append(flipped, flipped->beamdir, queue);
		RepairQueue_append(newjnc->nb[newjnc->beamdir^2], newjnc->beamdir, queue);
	}
	return flipped;
}

/*
	returns True if self needs a flip with the junction
	in direction d
*/
static int needsflip(Junction* jnc, unsigned char d) {
	double nbpos;
	double jncpos;
	if( (d+1)%2 ) {
		nbpos = jnc->nb[d]->pos[1]->posy;
		jncpos = jnc->pos[1]->posy;
	} else {
		nbpos = jnc->nb[d]->pos[0]->posx;
		jncpos = jnc->pos[0]->posx;
	}
	// needed for some reason to avoid wrong conclusions
	if(nbpos==jncpos)
		return 0;
	return (nbpos<jncpos) != ((d+1)%4)/2;
}

#ifdef COMMENT_THIS_OUT
/*
	returns True if self needs a flip with the junction
	in direction d
	alternative version; returns the same.
*/
static int needsflip_alt(Junction* jnc, unsigned char d) {
	Junction* nb = jnc->nb[d];
	assert(nb!=NULL);
	if(d%2==0) {
		if(nb->pos[1]->posy==jnc->pos[1]->posy) return 0;
		return (nb->pos[1]->posy > jnc->pos[1]->posy) != (d/2==0);
	} else {
		if(nb->pos[0]->posx==jnc->pos[0]->posx) return 0;
		return (nb->pos[0]->posx < jnc->pos[0]->posx) != (d/2==0);
	}
}
#endif



static void reconnect_linear(Junction* start, Junction* next, unsigned char d) {
	start->nb[d] = next->nb[d];
	if(next->nb[d]!=NULL)
		next->nb[d]->nb[d^2] = start;
	next->nb[d^2] = start->nb[d^2];
	if(start->nb[d^2]!=NULL)
		start->nb[d^2]->nb[d] = next;
	start->nb[d^2] = next;
	next->nb[d] = start;
}

/*
	slides an intersection of dir 4 into the direction tdir.
	this is done by removing one of the rays perpendicular to
	tdir and reinserting it after the next intersection in tdir.
*/
static void Junction_slide(Junction* jnc, unsigned char tdir, RepairQueue* queue) {
	assert(jnc->dir==4);
	assert(needsflip(jnc,tdir));
	Junction* bar = jnc->nb[tdir];
	if(bar->dir==(tdir^2))
		bar = Junction_flip(bar,queue);
	unsigned char ndir = bar->dir;
	Junction* next = jnc->nb[ndir];
	// Note: next can't be NULL
	if(next->dir!=(ndir^2))
		next = Junction_flip(next,queue);
	// remove ray
	RepairQueue_append(next->nb[next->beamdir^2], next->beamdir, queue);
	detach(next);
	// swap isc and bar
	reconnect_linear(jnc,bar,tdir);
	bar->beamdir = tdir^2;
	// append to queue...
	RepairQueue_append(jnc, tdir, queue);
	RepairQueue_append(bar, tdir^2, queue);
	// reinsert ray
	next = bar->nb[ndir];
	while(next->dir==(tdir^2))
		next = next->nb[ndir];
	Junction* newjnc = &jnc->pos[0]->rayend[ndir^2];
	if(next->dir==4) newjnc->beamdir=tdir;
	else if(next->dir==tdir) newjnc->beamdir=tdir^2;
	else newjnc->beamdir=next->beamdir;
	if(tdir%2==0) { // tdir is up or down
		assert(newjnc->pos[1]==jnc->pos[1]);
		newjnc->pos[0] = next->pos[0];
	} else {
		assert(newjnc->pos[0]==jnc->pos[0]);
		newjnc->pos[1] = next->pos[1];
	}
	newjnc->dir = ndir^2;
	// set neighbors
	newjnc->nb[ndir^2] = jnc;
	jnc->nb[ndir] = newjnc;
	newjnc->nb[tdir] = next->nb[tdir];
	if(newjnc->nb[tdir]!=NULL)
		newjnc->nb[tdir]->nb[tdir^2] = newjnc;
	newjnc->nb[tdir^2] = next;
	next->nb[tdir] = newjnc;
	// append new connections to queue
	RepairQueue_append(newjnc, newjnc->beamdir, queue);
	RepairQueue_append(newjnc->nb[newjnc->beamdir^2], newjnc->beamdir, queue);
	RepairQueue_append(jnc, ndir,queue);
}

/*
	slides two T-Junctions past each other, if possible.
	takes jnc and its neighbor in beamdir.
*/
static void Junction_slide_T(Junction* jnc, RepairQueue* queue) {
	assert(jnc->dir<4);
	assert(needsflip(jnc,jnc->beamdir));
	Junction* next = jnc->nb[jnc->beamdir];
	if(next->dir==jnc->dir || next->beamdir==(jnc->dir^2)) return;
	if(jnc->beamdir!=next->beamdir)
		next = Junction_flip(next,queue);
	assert(jnc->beamdir==next->beamdir);
	assert(jnc->dir==(next->dir^2));
	reconnect_linear(jnc,next,jnc->beamdir);
	// append to queue
	RepairQueue_append(jnc, jnc->beamdir, queue);
	next = next->nb[jnc->beamdir^2];
	if(next!=NULL)
		RepairQueue_append(next, jnc->beamdir, queue);
}

Box* Boxnet_addbox(Boxnet* net, double x, double y,
							double right, double top,
							Box* near, void* usrdata) {
	Box* new = Box_new();
	new->usrdata = usrdata;
	assert(right>=x && top>=y);
	new->posx = x;		new->posy = y;
	new->right = right;	new->top = top;
	if(near==NULL && net->boxes_size!=0)
		near = net->boxes[0];
	if(near!=NULL) {
		Junction_insert(&new->jnc, &near->jnc);
	} else {
		for(int d=0;d<4;d++)
			new->jnc.nb[d] = NULL;
	}
	vector_append(net->boxes, new, net->boxes_size,
					net->boxes_size_max, BOXES_SIZE_INIT);
	return new;
}

void Boxnet_delbox(Boxnet* net, Box* box) {
	// TODO: this is stupid, going through the whole array
	//       just to remove one element...
	for(int n=0;n<net->boxes_size;n++) {
		if(net->boxes[n]==box) {
			Box_free(box);
			net->boxes_size--;
			net->boxes[n] = net->boxes[net->boxes_size];
			return;
		}
	}
	assert(0);
}

// CAUTION: only removes the FIRST element that matches usrdata...
void Boxnet_delbox_byusrdata(Boxnet* net, void* usrdata) {
	// TODO: this is stupid, going through the whole array
	//       just to remove one element...
	for(int n=0;n<net->boxes_size;n++) {
		if(net->boxes[n]->usrdata==usrdata) {
			Box_free(net->boxes[n]);
			net->boxes_size--;
			net->boxes[n] = net->boxes[net->boxes_size];
			return;
		}
	}
	assert(0); // should never be reached
}

static double bnabs(double a) {
	return a > 0 ? a : -a;
}

/*
	TODO: This function could also try to optimize the
	memory layout...
*/
static void Boxnet_optimize(Boxnet* net) {
	static int n=0;
	if(n>=net->boxes_size) n=0;
	
	Box* box = net->boxes[n];
	for(int i=0;i<4;i++) {
		Junction* jnc = &box->rayend[i];
		assert(jnc->dir!=4);
		if(jnc->dir!=5) {
			if(jnc->dir%2==0) {
				if(bnabs( jnc->pos[1]->posx - box->posx ) >
					bnabs( jnc->pos[1]->posy - box->posy ))
					Junction_flip(jnc,NULL);
			} else {
				if(bnabs( jnc->pos[0]->posy - box->posy ) >
					bnabs( jnc->pos[0]->posx - box->posx ))
					Junction_flip(jnc,NULL);
			}
		}
	}
	n++;
}

/*
	Ensures that the coordinate relations implied by the
	boxnet structure are consistent with the explicit
	coordinates of the boxes.
	Gets called by collide(), so the user should never
	need to call this for normal usage...
*/
void Boxnet_repair(Boxnet* net) {
	static RepairQueue	queue1 = {NULL,0,0};
	static RepairQueue	queue2 = {NULL,0,0};
	// queue can not get longer than 4*(number of boxes)
	if(queue1.queue==NULL) {
		queue1 = RepairQueue_new();
		queue2 = RepairQueue_new();
	}
	void solve_conn(Junction* jnc, unsigned char tdir, RepairQueue* q) {
		assert(jnc->enqueued!=0);
		jnc->enqueued ^= (1<<tdir);
		assert((jnc->enqueued & (1<<tdir))==0);
		if(jnc->dir==5) return;
		Junction* next = jnc->nb[tdir];
		if(next==NULL) return;
		if(!needsflip(jnc,tdir)) return;
		if(jnc->dir==4) {
			Junction_slide(jnc,tdir,q);
		} else {
			if(jnc->beamdir!=tdir) return;
			Junction_slide_T(jnc,q);
		}
	}
	queue1.size=0;
	queue2.size=0;
	//int max_queue_size = 0;
	for(int i=0;i<net->boxes_size;i++) {
		for(unsigned char tdir=0;tdir<4;tdir++) {
			RepairQueue_append(&net->boxes[i]->jnc,tdir, &queue1);
			Junction* jnc = &net->boxes[i]->rayend[tdir];
			if(jnc->dir!=5)
				RepairQueue_append(jnc,jnc->beamdir, &queue1);
		}
		while(queue1.size>0) {
			while(queue1.size>0) {
				//if(queue1.size>max_queue_size)
				//	max_queue_size=queue1.size;
				queue1.size--;
				solve_conn(queue1.queue[queue1.size].jnc, queue1.queue[queue1.size].tdir, &queue2);
			}
			while(queue2.size>0) {
				//if(queue2.size>max_queue_size)
				//	max_queue_size=queue2.size;
				queue2.size--;
				solve_conn(queue2.queue[queue2.size].jnc, queue2.queue[queue2.size].tdir, &queue1);
			}
		}
	}
	//printf("max queue size: %i\n",max_queue_size);
	//Boxnet_optimize(net);
}

/*
	finds collisions for this bounding box; this does
	not find all collisions of box, in the sense that each collision
	pair is reported only by one of the two involved physpoints.
	The collision net needs to be valid.
	CAUTION: this function assumes that the lower edge of
	the BB "stands" completely on rays, i.e. if you go
	from self.intersection to the right, you encounter
	the lower right vertex of the BB before you encounter
	an intersection of type -| (dir=1)
	also, rasterization results are meaningless if the
	net changes between two rasterizations of BBs that
	should be tested for overlap.
	
	CAUTION: do NOT call Boxnet_delbox from the collision callback
	function "func", or else you will have buggy behaviour!
*/
static void boxcollisions(Box* box, Boxnet* net, collisionCallback func, void* data) {
	//Box_overlap_right_append(Box* box, Box* append)
	static Box**		queue = NULL;
	static int			queue_size_max;
	int					queue_size;
	if(queue==NULL) {
		// TODO: do error checking... (NULL pointer)
		queue = malloc(BC_QUEUE_SIZE_INIT * sizeof *queue);
		assert(queue!=NULL);
		queue_size_max=BC_QUEUE_SIZE_INIT;
	}
	void queue_append(Box* append) {
		if(append->marked==box)
			return;
		append->marked=box;
		// add to overlap regions
		assert(append->posy <= box->top);
		assert(append->top >= box->posy);
		assert(box!=append);
		if(append->posx <= box->right &&
					append->right >= box->posx) {
				func(box->usrdata,append->usrdata,data);
		}
		vector_append(queue, append, queue_size, queue_size_max, BC_QUEUE_SIZE_INIT);
	}
	queue[0] = box;
	for(queue_size = 1;queue_size>0;) {
		// go left
		// BEWARE: nearly duplicated code below...
		queue_size--;
		Junction* jnc = &queue[queue_size]->jnc;
		Junction* root = jnc;
		while(root!=NULL && root->dir!=3 && root->pos[0]->posx > box->posx) {
			if(root->dir != 2) {
				Junction* next = root->nb[0];
				// go upwards until we can go forward
				while(next!=NULL && next->pos[1]->posy <= box->top) {
					if(next->dir!=3) {
						queue_append(next->pos[1]);
						break;
					}
					next = next->nb[0];
				}
			}
			root = root->nb[1];
		}
		// go right
		// BEWARE: nearly duplicated code above...
		root = jnc;
		while(root!=NULL && root->dir!=1 &&
						root->pos[0]->posx <= box->right) {
			if(root->dir != 2) {
				Junction* next = root->nb[0];
				// go upwards until we can go forward
				while(next!=NULL && next->pos[1]->posy <= box->top) {
					if(next->dir!=1) {
						queue_append(next->pos[1]);
						break;
					}
					next = next->nb[0];
				}
			}
			root = root->nb[3];
		}
	}
}

/*
	find all collisions between BBs.
	repairs the net before finding collisions.
*/
void Boxnet_collide(Boxnet* net, collisionCallback func, void* data) {
	Boxnet_repair(net);
	// prepare net for collisions
	for(int i=0;i<net->boxes_size;i++) {
		Box* box = net->boxes[i];
		box->marked=NULL;
		for(Junction* next = box->jnc.nb[3];
				next != NULL && next->pos[0]->posx <= box->right;
				next = next->nb[3]) {
			if(next->dir==1)
				next = Junction_flip(next,NULL);
		}
	}
	// find collisions
	for(int i=0;i<net->boxes_size;i++)
		boxcollisions(net->boxes[i], net, func, data);
}












/*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *    <------------ DEBUGGING FUNCTIONS ---------------->
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
*/




/*
	checks the boxnet for flipped connections
	returns 1 if all spatial relationships are
	as in the absolute box positions, 0 if
	boxnet is not correctly repaired
*/
static int repair_check(Boxnet* net) {
	// control results
	for(int i=0;i<net->boxes_size;i++) {
		for(unsigned char tdir=0;tdir<4;tdir++) {
			Junction* next = net->boxes[i]->jnc.nb[tdir];
			if(next!=NULL)
				if(needsflip(next,tdir^2))
					return 0;
			Junction* jnc = &net->boxes[i]->rayend[tdir];
			if(jnc->dir!=5)
				if(jnc->nb[jnc->beamdir]!=NULL)
					if(needsflip(jnc,jnc->beamdir))
						return 0;
		}
	}
	return 1;
}


/*
	finds inconsistencies by trying to deduce possible
	positions for all points from the boxnet structure.
	if this is not possible, the boxnet is invalid...
*/
static void find_inconsistencies(Boxnet* net) {
	int move(Junction* jnc, unsigned char d) {
		if((d^2)==jnc->dir) return 1;
		Junction* nb = jnc->nb[d];
		if(nb==NULL) return 1;
		if(d%2) {
			if( jnc->pos[0]->posx - nb->pos[0]->posx >= 1.)
				return 1;
			jnc->pos[0]->posx = nb->pos[0]->posx + 1;
		} else {
			if( jnc->pos[1]->posy - nb->pos[1]->posy >= 1.)
				return 1;
			jnc->pos[1]->posy = nb->pos[1]->posy + 1;
		}
		return 0;
	}
	// for restoring the positions later
	double posx[net->boxes_size];
	double posy[net->boxes_size];
	for(int i=0;i<net->boxes_size;i++) {
		posx[i] = net->boxes[i]->posx;
		posy[i] = net->boxes[i]->posy;
		net->boxes[i]->posx = 0.;
		net->boxes[i]->posy = 0.;
	}
	void restore() {
		for(int i=0;i<net->boxes_size;i++) {
			net->boxes[i]->posx = posx[i];
			net->boxes[i]->posy = posy[i];
		}
	}
	
	double maxsize = net->boxes_size;
	int done = 0;
	while(!done) {
		done=1;
		for(int i=0;i<net->boxes_size;i++) {
			Box* box = net->boxes[i];
			if(box->posx > maxsize || box->posy > maxsize) {
				restore();
				assert(0); // net is invalid
			}
			done = done & move(&box->jnc,1) & move(&box->jnc,2);
			for( Junction* cur = box->jnc.nb[3]; cur != NULL; cur=cur->nb[3]) {
				done = done & move(cur,1) & move(cur,2);
				if(cur->dir==1) break;
			}
			for( Junction* cur = box->jnc.nb[1]; cur != NULL; cur=cur->nb[1]) {
				done = done & move(cur,1) & move(cur,2);
				if(cur->dir==3) break;
			}
		}
	}
	restore();
}

/*
	tries to find every possible error that could be present in
	the structure of the boxnet;
	returns 1 if no error could be found, 0 if an error was found
*/
static void validate(Boxnet* net) {
	// find simple errors like wrong beamdirs
	// and wrong links
	for(int i=0;i<net->boxes_size;i++) {
		Box* box = net->boxes[i];
		for(unsigned char tdir=0;tdir<4;tdir++) {
			Junction* prev = &box->jnc;
			Junction* next = box->jnc.nb[tdir];
			while(next!=NULL) {
				// check connection to previous
				assert(next->nb[tdir^2]==prev);
				// check enqueued flag
				assert(next->enqueued==0);
				// dir < 4
				assert(next->dir<4);
				// posx and posy need to be from different boxes
				assert(next->pos[0]!=next->pos[1]);
				// main position needs to be from originating box
				assert(next->pos[tdir%2]==box);
				// end junction needs to be one of rayend[]
				if(next->dir==(tdir^2)) {
					assert(next == &box->rayend[tdir^2]);
					break;
				}
				// check beamdir
				assert(next->beamdir==tdir);
				prev = next;
				next = next->nb[tdir];
			}
			// check if unused rayends have dir==5 set
			if(next==NULL)
				assert(box->rayend[tdir^2].dir==5);
		}
	}
	find_inconsistencies(net);
}


/*
	prints the net; only works if connected boxes don't have the
	same x or y value.
*/
void print_net(Boxnet* net) {
	// check net for degenerate positions (in that case we can't print
	// the net because of the simple save format)
	for(int i=0;i<net->boxes_size;i++) {
		for(int j=0;j<net->boxes_size;j++) {
			if(net->boxes[i]->posx == net->boxes[j]->posx ||
						net->boxes[i]->posy == net->boxes[j]->posy)
				printf("can't save net (duplicate positions).\n");
		}
	}
	printf("boxnet dump:\n");
	for(int i=0;i<net->boxes_size;i++) {
		Box* b = net->boxes[i];
		printf("P:%f,%f,%f,%f:",b->posx,b->posy,b->right,b->top);
		int c[4]={0,0,0,0};
		for(unsigned char d=0;d<4;d++) {
			Junction* next = b->jnc.nb[d];
			while(next!=NULL && next->dir!=(d^2)) {
				c[d]++;
				next=next->nb[d];
			}
		}
		printf("%i,%i,%i,%i\n",c[0],c[1],c[2],c[3]);
	}
}

#ifndef NOTEST

#include <time.h>
#include <math.h>
//#include <time, random, ... .h>


// helper stuff for storing collisions...
typedef struct Collision {
	Box*		box1;
	Box*		box2;
} Collision;

typedef struct Collisions {
	struct Collision*	cols;
	int					size;
	int					size_max;
} Collisions;

Collisions* Collisions_new() {
	Collisions* cols = malloc(sizeof *cols);
	cols->cols = NULL;
	cols->size_max = 0;
	cols->size = 0;
	return cols;
}

void col_callback(void* obj1, void* obj2, void* data) {
	Collisions* cols = (Collisions*) data;
	Collision col;
	col.box1 = (Box*) obj1;
	col.box2 = (Box*) obj2;
	vector_append(cols->cols, col, cols->size, cols->size_max, 256);
}

void Boxnet_collide_store(Boxnet* net, Collisions* cols) {
	cols->size = 0;
	Boxnet_collide(net,col_callback,cols);
}


/*
	brute-force test collision results for correctness;
	call this after Boxnet_collide() to check the results
	stored in boxnet.collisions[].
	
	returns 1 if all collisions are already in the buffer,
	0 if new collisions were found or false positives
	were in the buffer.
*/
int collide_control(Boxnet* net, Collisions* cols) {
	//test for false negatives
	for(int i=0;i<net->boxes_size;i++) {
		Box* box1 = net->boxes[i];
		for(int j=i+1;j<net->boxes_size;j++) {
			Box* box2 = net->boxes[j];
			if(			box1->posx <= box2->right &&
						box1->right >= box2->posx &&
						box1->posy <= box2->top &&
						box1->top >= box2->posy ) {
				int found=0;
				for(int n=0;n<cols->size;n++) {
					if( (cols->cols[n].box1==box1 &&
							cols->cols[n].box2==box2) ||
							(cols->cols[n].box1==box2 &&
							cols->cols[n].box2==box1) ) {
						found=1;
						break;
					}
				}
				if(!found) {
					printf("false negative: %i  %i\n",(int)box1, (int)box2);
					/*printf("pair not found: %i  %i\n",(int)box1, (int)box2);
					printf("(%.3f %.3f %.3f %.3f) (%.3f %.3f %.3f %.3f)\n",
							box1->posx,box1->posy,box1->right,box1->top,
							box2->posx,box2->posy,box2->right,box2->top);*/
					return 0;
				}
			}
		}
	}
	// test for false positives
	for(int n=0;n<cols->size;n++) {
		Box* box1 = cols->cols[n].box1;
		Box* box2 = cols->cols[n].box2;
		if(	!(box1->posx <= box2->right &&
				box1->right >= box2->posx &&
				box1->posy <= box2->top &&
				box1->top >= box2->posy) ) {
			printf("false positive: %i  %i\n",(int)box1, (int)box2);
			return 0;
		}
	}
	return 1;
}


double random_d() {
	return (double)rand()/RAND_MAX;
}


/*
	tries to uncover bugs by shuffling boxes randomly around.
	if assert() is disabled, acts as a benchmark
	nbox		number of boxes
	ncycl		number of move/repair iterations
	ndelete		number of boxes deleted and recreated each iteration
	discrete	1 if only discrete positions and sizes should be used,
				0 otherwise
*/
void stresstest(int nbox, int ncycl, int ndelete, int discrete, double stepcoeff) {
	Collisions* cols = Collisions_new();
	assert(ndelete<=nbox);
	int Ndis = (int)(0.1*sqrt(nbox)+1);
	void quantize(Box* box) {
		box->posx = ((int)(box->posx*Ndis))/(double)Ndis;
		box->posy = ((int)(box->posy*Ndis))/(double)Ndis;
	}
	void resize(Box* box) {
		if(discrete) {
			double unit = 1/(double)Ndis;
			if(random_d()<0.8)
				box->right = box->posx + unit;
			else
				box->right = box->posx;
			for(int i=0;i<Ndis && random_d()<0.2;i++)
				box->right += unit;
			if(random_d()<0.8)
				box->top = box->posy + unit;
			else
				box->top = box->posy;
			for(int i=0;i<Ndis && random_d()<0.2;i++)
				box->top += unit;
		} else {
			box->right = box->posx + random_d()*sqrt(1/(double)nbox);
			box->top = box->posy + random_d()*sqrt(1/(double)nbox);
		}
	}
	void create(Boxnet* net) {
		double x,y;
		x = random_d();
		y = random_d();
		Box* box = Boxnet_addbox(net, x,y,x,y,NULL,NULL,1);
		box->usrdata = box;
		if(discrete)
			quantize(box);
		resize(box);
	}
	void move(Box* box, double step) {
		double triangle(double x) {
			return 2*fabs(0.5*x-floor(0.5*x+0.5));
		}
		box->posx = triangle((box->posx+step*(0.5-random_d())));
		box->posy = triangle((box->posy+step*(0.5-random_d())));
		if(discrete)
			quantize(box);
		resize(box);
		assert(box->right>=box->posx && box->top>=box->posy);
	}
	printf("creating %i boxes...\n",nbox);
	Boxnet* net = Boxnet_new();
	for(int n=0;n<nbox;n++) {
		create(net);
	}
	printf("shuffling boxes wildly...\n");
	for(int n=0;n<ncycl;n++) {
		// deleting and creating boxes
		for(int i=0;i<ndelete;i++) {
			assert(net->boxes_size>0);
			Boxnet_delbox(net, net->boxes[rand()%net->boxes_size]);
		}
		for(int i=0;i<ndelete;i++)
			create(net);
		double step = random_d();
		step = stepcoeff*2*(step*step*step*step);
		for(int i=0;i<net->boxes_size;i++) {
			Box* box = net->boxes[i];
			move(box,step);
		}
		assert(nbox==net->boxes_size);
		
#ifndef NDEBUG
		Boxnet_repair(net);
		assert(repair_check(net));
#endif
		
		Boxnet_collide_store(net,cols);
		printf("n =%7i/%i, step=%.3f, collisions: %i\n",n+1,ncycl,step,cols->size);
		assert(collide_control(net,cols));
		
#ifndef NDEBUG
		validate(net);
#endif
	}
	Boxnet_free(net);	
}


int main() {
	// reproducability
	srand(10389);
	//stresstest(10000,10000,100,1);
#ifndef NDEBUG
	stresstest(1000,200,100,0,1.0);
	stresstest(200,70,10,1,1.0);
#endif
#ifdef NDEBUG
	stresstest(10000,1000,10,0,0.003);
#endif
	
}

#endif





