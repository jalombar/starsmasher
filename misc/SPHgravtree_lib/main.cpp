/*
 * 
 * 
 * Test program for the Bonsai library
 * 
 * 
 * Author: Jeroen Bédorf
 * Date: 19 - Jan - 2012
 * 
 */

#include <iostream>
#include <stdlib.h>
#include <vector>
#include <fstream>

#include "../lib/include/sequoiaInterface.h"
#include "../lib/include/my_cuda.h"


#include "tipsydefs.h"

void read_tipsy_file_parallel(vector<real4> &bodyPositions, vector<real4> &bodyVelocities,  vector<int> &bodiesIDs,  float eps2,
    string fileName, int rank, int procs, int &NTotal2, int &NFirst, int &NSecond, int &NThird)  
{
  //Process 0 does the file reading and sends the data
  //to the other processes

  //Now we have different types of files, try to determine which one is used
  /*****
    If individual softening is on there is only one option:
    Header is formatted as follows: 
    N     #       #       #
    so read the first number and compute how particles should be distributed

    If individual softening is NOT enabled, i can be anything, but for ease I assume standard dumbp files:
    no Header
    ID mass x y z vx vy vz
    now the next step is risky, we assume mass adds up to 1, so number of particles will be : 1 / mass
    use this as initial particle distribution

*/

  char fullFileName[256];
  sprintf(fullFileName, "%s", fileName.c_str());

  cout << "Trying to read file: " << fullFileName << endl;



  ifstream inputFile(fullFileName, ios::in | ios::binary);
  if(!inputFile.is_open())
  {
    cout << "Can't open input file \n";
    exit(0);
  }

  dump  h;
  inputFile.read((char*)&h, sizeof(h));  

  int NTotal;
  int idummy;
  real4 positions;
  real4 velocity;


  //Read tipsy header  
  NTotal        = h.nbodies;
  NFirst        = h.ndark;
  NSecond       = h.nstar;
  NThird        = h.nsph;

  //Rough divide
  uint perProc = NTotal / procs;
  bodyPositions.reserve(perProc+10);
  bodyVelocities.reserve(perProc+10);
  bodiesIDs.reserve(perProc+10);
  perProc -= 1;

  //Start reading
  int particleCount = 0;
  //   int procCntr = 1;

  dark_particle d;
  star_particle s;

  for(int i=0; i < NTotal; i++)
  {
    if(i < NFirst)
    {
      inputFile.read((char*)&d, sizeof(d));
      velocity.w        = d.eps;
      positions.w       = d.mass;
      positions.x       = d.pos[0];
      positions.y       = d.pos[1];
      positions.z       = d.pos[2];
      velocity.x        = d.vel[0];
      velocity.y        = d.vel[1];
      velocity.z        = d.vel[2];
      idummy            = d.phi;
    }
    else
    {
      inputFile.read((char*)&s, sizeof(s));
      velocity.w        = s.eps;
      positions.w       = s.mass;
      positions.x       = s.pos[0];
      positions.y       = s.pos[1];
      positions.z       = s.pos[2];
      velocity.x        = s.vel[0];
      velocity.y        = s.vel[1];
      velocity.z        = s.vel[2];
      idummy            = s.phi;
    }

    bodyPositions.push_back(positions);
    bodyVelocities.push_back(velocity);
    bodiesIDs.push_back(idummy);  

    particleCount++;


    //     if(bodyPositions.size() > perProc && procCntr != procs)
    //     { 
    //       tree->ICSend(procCntr,  &bodyPositions[0], &bodyVelocities[0],  &bodiesIDs[0], bodyPositions.size());
    //       procCntr++;
    //       
    //       bodyPositions.clear();
    //       bodyVelocities.clear();
    //       bodiesIDs.clear();
    //     }
  }//end while

  inputFile.close();

  //Clear the last one since its double
  //   bodyPositions.resize(bodyPositions.size()-1);  
  //   NTotal2 = particleCount-1;
  NTotal2 = particleCount;
  cerr << "NTotal: " << NTotal << "\tper proc: " << perProc << "\tFor ourself:" << bodiesIDs.size() << endl;
}




int CPUNeighbourCount(my_dev::dev_mem<real4> &j_bodies_pos,                           
    my_dev::dev_mem<int>   &j_bodies_ids,
    int                     n_jbodies,
    real4                   pos_i, 
    int                     i_bodies_ids,
    float                   h_i)

{

  int ngbCount = 0;
  for(int i=0; i < n_jbodies; i++)
  {
    real4 pos_j = j_bodies_pos[i];
    float3 dr  = {pos_i.x - pos_j.x,
      pos_i.y - pos_j.y,
      pos_i.z - pos_j.z};

    float ds2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;



    ngbCount += (ds2 <= h_i*h_i) * (i_bodies_ids !=  j_bodies_ids[i]);
    if(ds2 <= h_i*h_i)
    {
      //fprintf(stderr,"ds2: %f  h_i*h_i: %f  id: %d %d \n", ds2,  h_i*h_i, i_bodies_ids, j_bodies_ids[i]);
    }
  }  

  return ngbCount;



}



using namespace std;


int main(int argc, char** argv)
{

  //Initialize the code

  float eps = 0.05;
  my_dev::context devContext =  sequoia_init(argv, 0, 0.75, eps);


  //Get the particles
  vector<real4> bodyPositions;
  vector<real4> bodyVelocities;
  vector<int> bodyIDs;

  string fileName = "model3_child_compact.tipsy";

  if(argc > 1)
    fileName.assign(argv[1]);

  int NTotal, NFirst, NSecond, NThird;
  //Read particles
  read_tipsy_file_parallel(bodyPositions, bodyVelocities, bodyIDs, eps, fileName,
      0, 1, NTotal, NFirst, NSecond, NThird) ;

  int n_bodies = bodyPositions.size();                      
  my_dev::dev_mem<real4> bodies_pos(devContext, n_bodies);    //Bodies positions
  my_dev::dev_mem<real4> bodies_vel(devContext, n_bodies);    //Bodies velocities
  my_dev::dev_mem<int>   bodies_ids(devContext, n_bodies);      //Bodies ids

  //Output
  my_dev::dev_mem<real4> bodies_acc(devContext, n_bodies);    //Bodies Accelerations
  my_dev::dev_mem<real>  bodies_ds2(devContext, n_bodies);    //Bodies distance to nearest neighbour squared
  my_dev::dev_mem<int>   bodies_ngb(devContext, n_bodies);    //Bodies nearest neighbour



  //    my_dev::dev_mem<uint> sortPermutation(devContext, n_bodies);    //Bodies

  //Copy data in the GPU Host and devices buffers
  memcpy(&bodies_pos[0], &bodyPositions[0],  n_bodies * sizeof(real4));
  bodies_pos.h2d();
  memcpy(&bodies_vel[0], &bodyVelocities[0], n_bodies * sizeof(real4));
  bodies_vel.h2d();
  memcpy(&bodies_ids[0], &bodyIDs[0], n_bodies * sizeof(int));
  bodies_ids.h2d();

  //Sort and reorder the bodies
  //    sequoia_sortBodies(bodies_pos, sortPermutation, n_bodies);   
  //    sequoia_reorderReal4(bodies_pos, sortPermutation, n_bodies);
  //    sequoia_reorderReal4(bodies_vel, sortPermutation, n_bodies);
  //    sequoia_reorderInt1(bodies_ids, sortPermutation, n_bodies);   
  //    
  //    sequoia_buildTreeStructure(bodies_pos, n_bodies);
  //    sequoia_createGroups(bodies_pos, n_bodies);
  //       
  //    sequoia_computeGravity(bodies_pos, bodies_acc, n_bodies);

  for(int i=0; i < 1; i++)
  {
    sequoia_setParticlesAndGetGravity(bodies_pos,      //Positions J-particles
        bodies_ids,      //Particle IDs J-particles
        n_bodies,         //Number of J-particles
        bodies_pos,      //Positions I-particles (Can be the same or different than J-particles)
        bodies_ids,      //Particle IDs J-particles (Can be the same or different than J-particles)
        n_bodies,         //Number of I-particles (Can be the same or different than J-particles)
        true,        //Do we need to sort J-particles?
        false,        //Do we need to sort I-particles? Can be false if i-bodies are the same as j-bodies           
        bodies_acc,      //OUT  Accelerations for I-particles
        bodies_ds2,      //OUT  min distance squared for I-particles
        bodies_ngb);      //OUT  J-ID of the nearest neighbour for I-particles

  }

  /*for(int i=0; i < 10; i++)
    {

    fprintf(stderr, "%d \t acc: %f %f %f %f \t\tds2: %f \tngb-idx: %d \t ngb-id: %d\n", i*10,
    bodies_acc[i*10].x,bodies_acc[i*10].y,bodies_acc[i*10].z,bodies_acc[i*10].w,
    bodies_ds2[i*10],bodies_ngb[i*10],
    -1);      
  //  (bodies_ngb[i*10] >=0) ? bodies_ids[bodies_ngb[i*10]] : -1);
  }*/





  sequoia_get_ngb(bodies_pos,                           
      bodies_ids, 
      bodies_pos, 
      bodies_ids, 
      bodies_ngb);




  //Copy results to the host
  bodies_acc.d2h();
  bodies_ds2.d2h();
  bodies_ngb.d2h();
  bodies_ids.d2h();

  bodies_pos.d2h();//Needed for CPU test

  fprintf(stderr, "\n");

  int jump = 1;
  long long sum = 0;
  for(int i=0; i < n_bodies; i++)
  {
    /*  int cpuN = CPUNeighbourCount(bodies_pos, bodies_ids,
        n_bodies, bodies_pos[i*jump],
        bodies_ids[i*jump], 0.1);*/
    int cpuN = 0;
    if(cpuN >= 0)
    {
      // if(bodies_ngb[i*jump] > 0)
      {
        /*  cpuN = CPUNeighbourCount(bodies_pos, bodies_ids,
            n_bodies, bodies_pos[i*jump],
            bodies_ids[i*jump], 0.1);

*/
        sum += bodies_ngb[i*jump];

        /*    fprintf(stderr, "TESTSTAT %d \t acc: %f %f %f %f \t\tds2: %f \tgpu-nngb: %d cpu-nngb: %d  \n", i*jump,
              bodies_acc[i*jump].x,bodies_acc[i*jump].y,bodies_acc[i*jump].z,bodies_acc[i*jump].w,
              bodies_ds2[i*jump],bodies_ngb[i*jump], cpuN);*/

        if(bodies_ngb[i*jump] != cpuN)
        {
          //   fprintf(stderr, "ERROR\n");
        }
      }
    }
  }

  fprintf(stderr, "%ld %d   %f \n", sum, n_bodies, sum / (float)n_bodies);


  //Make sure to free memory before cleanup
  //    sortPermutation.free_mem();
  bodies_acc.free_mem();
  bodies_pos.free_mem();
  bodies_vel.free_mem();
  bodies_ids.free_mem();
  bodies_ds2.free_mem();
  bodies_ngb.free_mem();

  //Clean up
  sequoia_cleanup();

  return 0;  
}


#if 0

int main(int argc, char** argv)
{

  //Initialize the code

  float eps = 0.05;
  my_dev::context devContext =  sequoia_init(argv, 0, 0.75, eps);


  //Get the particles
  vector<real4> bodyPositions;
  vector<real4> bodyVelocities;
  vector<int> bodyIDs;

  string fileName = "model3_child_compact.tipsy";

  if(argc > 1)
    fileName.assign(argv[1]);

  int NTotal, NFirst, NSecond, NThird;
  //Read particles
  read_tipsy_file_parallel(bodyPositions, bodyVelocities, bodyIDs, eps, fileName,
      0, 1, NTotal, NFirst, NSecond, NThird) ;

  int n_bodies = bodyPositions.size();                      
  my_dev::dev_mem<real4> bodies_pos(devContext, n_bodies);    //Bodies positions
  my_dev::dev_mem<real4> bodies_acc(devContext, n_bodies);    //Bodies positions
  my_dev::dev_mem<real4> bodies_vel(devContext, n_bodies);    //Bodies velocities
  my_dev::dev_mem<int> bodies_ids(devContext, n_bodies);      //Bodies ids

  my_dev::dev_mem<uint> sortPermutation(devContext, n_bodies);    //Bodies

  memcpy(&bodies_pos[0], &bodyPositions[0],  n_bodies * sizeof(real4));
  bodies_pos.h2d();
  memcpy(&bodies_vel[0], &bodyVelocities[0], n_bodies * sizeof(real4));
  bodies_vel.h2d();
  memcpy(&bodies_ids[0], &bodyIDs[0], n_bodies * sizeof(int));
  bodies_ids.h2d();

  //Sort and reorder the bodies
  sequoia_sortBodies(bodies_pos, sortPermutation, n_bodies);   
  sequoia_reorderReal4(bodies_pos, sortPermutation, n_bodies);
  sequoia_reorderReal4(bodies_vel, sortPermutation, n_bodies);
  sequoia_reorderInt1(bodies_ids, sortPermutation, n_bodies);


  sequoia_buildTreeStructure(bodies_pos, n_bodies);
  sequoia_createGroups(bodies_pos, n_bodies);


  sequoia_computeGravity(bodies_pos, bodies_acc, n_bodies);


  //Make sure to free memory before cleanup
  sortPermutation.free_mem();
  bodies_acc.free_mem();
  bodies_pos.free_mem();
  bodies_vel.free_mem();
  bodies_ids.free_mem();
  //Clean up
  sequoia_cleanup();

  return 0;  
}
#endif
