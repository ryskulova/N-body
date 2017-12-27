#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <mpi.h>
#define PI 3.14159265358979323846

typedef struct _body {
	float x, y; // position
	float ax, ay; // acceleration
	float vx, vy; // velocity
	float mass; // mass
} body;

body  *initBodies (int nBodies);
void integrate(body *body, float deltaTime);
void calculateNewtonGravityAcceleration(body *a, body *b, float *ax, float *ay);
void simulateWithBruteforce(int rank, int totalBodies, int *bodies, float dt);
float randerValue();

int main(int argc, char **argv) {
	
	int nBodies = 100;
	float simulationTime = 1.0;
	float dt = 0.1;
	if(argc > 3){
		nBodies = atof(argv[1]);
		simulationTime = atoi(argv[2]);
		dt = atoi(argv[3]);
	}
	double parallel_average_time = 0.0;
    MPI_Init(&argc, &argv);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


	MPI_Datatype dt_body;
  	MPI_Type_contiguous(7, MPI_FLOAT, &dt_body);
  	MPI_Type_commit(&dt_body);
	body *bodies = initBodies(nBodies);
    if (rank == 0) {
        parallel_average_time -= MPI_Wtime();
	 
    }
    size_t items_per_process = nBodies / world_size;
    
	body *local_bodies =
        (body *) malloc(sizeof(*local_bodies) * items_per_process);
    MPI_Scatter(
        bodies,
        items_per_process,
        dt_body,
        local_bodies,
        items_per_process,
        dt_body,
        0,
        MPI_COMM_WORLD
    );
	simulateWithBruteforce(rank, nBodies, items_per_process, bodies, local_bodies, dt);
	body *gathered_bodies = NULL;
	if (rank == 0) {
		gathered_bodies =
        (body *) malloc(sizeof(*gathered_bodies) * nBodies);
	}
	MPI_Gather(
        local_bodies,
        items_per_process,
        dt_body,
        gathered_bodies,
        items_per_process,
        dt_body,
        0,
        MPI_COMM_WORLD
        );
	if (rank == 0) {
        parallel_average_time += MPI_Wtime();
          printf("%d\n", nBodies);
         printf("%.2f\n", simulationTime);
        printf("%.2f\n", dt);
    
        for (size_t i = 0; i < nBodies; ++i){
		printf("%.2f %.2f\n", bodies[i].x, bodies[i].y);
         	printf("%.2f %.2f\n", bodies[i].ax, bodies[i].ay);
         	printf("%.2f %.2f\n", bodies[i].vx, bodies[i].vy);
         	printf("%.2f\n", bodies[i].mass);
        for (float j = 0.0; j < simulationTime; j += dt)
    	{
 	    	for (size_t i = 0; i < nBodies; ++i){
 
 	    		printf("%.2f %.2f\n", gathered_bodies[i].ax, gathered_bodies[i].ay);
 	    		integrate(&gathered_bodies[i], dt);
 	    	}
 	    }
        	
        }
    }
	if (bodies != NULL) {
        free(bodies);
    }
    free(local_bodies);
    MPI_Finalize();
	return 0;
}
void integrate(body *body, float deltaTime) {
        body->x += body->vx * deltaTime;
	body->y += body->vy * deltaTime;
	body->vx += body->ax * deltaTime;
	body->vy += body->ay * deltaTime;
}
void calculateNewtonGravityAcceleration(body *a, body *b, float *ax, float *ay) {
	float distanceX = b->x - a->x
	float distanceY = b->y - a->y
	float vectorDistance = a->x * a->x + a->y * a->y
	float vectorDistanceCubed = vectorDistance * vectorDistance * vectorDistance;
    float inverse = 1.0 / sqrt( vectorDistanceCubed);
	float soft = 2000;
    float scale = b->mass * inverse;
	*ax = (distanceX * scale);
	*ay = (distanceY * scale);
}
void simulateWithBruteforce(int rank, int totalBodies, int *bodies, float dt) {
	for(size_t i = 0; i < nBodies; i++) {
		float total_ax = 0, total_ay = 0;
		for (size_t j = 0; j < totalBodies; j++) {
			if (j == nBodies * rank + i) {
				continue;
			}
			float ax, ay;
			calculateNewtonGravityAcceleration(&local_bodies[i], &bodies[j], &ax, &ay);
			total_ax += ax;
			total_ay += ay;
		}
		local_bodies[i].ax = total_ax;
		local_bodies[i].ay = total_ay;
		integrate(&local_bodies[i], dt);
	}
}

float randerValue(){
	return ((float) rand() / RAND_MAX);
}
body  *initBodies (int nBodies) {
	srand(time(NULL));
	const float accelerationScale = 100.0;
	body *bodies = (body *) malloc(sizeof(*bodies) * nBodies);
	for (int i = 0; i < nBodies; i++) {
		float angle = 	((float) i / nBodies) * 2.0 * PI + 
						((randerValue() - 0.5) * 0.5);
		float initialMass = 200;
		body object = {	
			.x = randerValue(), .y = randerValue(),
			.vx = cos(angle) * accelerationScale * randerValue(), 
			.vy = sin(angle) * accelerationScale * randerValue(), 
			.mass = initialMass * randValue() + initialMass * 0.5
		};
        bodies[i] = object;
    }
    return bodies;
}
