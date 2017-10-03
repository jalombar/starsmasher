#include "support_kernels.cu"
#include "dev_shared_traverse_functions.cu"

#include <stdio.h>

//Some settings

#define TEXTURES
#define OLDPREFIX
#define DOGRAV

#if 1      /* use this if you want the potential to be in double precision */
#define _GPOTDP_
#endif


struct ds64
{
  float2 val;
  __host__ __device__ ds64() {}
  __host__ __device__ ds64(float x) : val(make_float2(x, x)) {}
  __host__ __device__ ds64 operator+=(const float x) 
  {
    const float vx = val.x + x;
    const float vy = val.y - ((vx - val.x) - x);
    val = make_float2(vx, vy);
    return *this;
  }
  __host__ __device__ double to_double() const { return (double)val.x + (double)val.y; }
  __host__ __device__ float to_float() const { return (float)((double)val.x + (double)val.y);}
};

struct devForce
{
  ds64 ax, ay, az;   // 6
#ifdef _GPOTDP_
  ds64 pot;          // 8
#else
  float pot;         // 7
  int  iPad;         // 8
#endif
  __host__ __device__ devForce() {}
  __device__ devForce(const float v) : ax(v), ay(v), az(v), pot(v) {}
  __device__ float4 to_float4() 
  {
#ifdef _GPOTDP_
    return (float4){ax.to_float(), ay.to_float(), az.to_float(), pot.to_float()};
#else
    return (float4){ax.to_float(), ay.to_float(), az.to_float(), pot};
#endif
  }
};


texture<float4, 1, cudaReadModeElementType> texNodeSize;
texture<float4, 1, cudaReadModeElementType> texNodeCenter;
texture<float4, 1, cudaReadModeElementType> texMultipole;
texture<float4, 1, cudaReadModeElementType> texBody;


__device__ devForce body_body(
    devForce iforce,
    const float4 posi,
    const float4 _hi,
    const float  massj,
    const float3 posj,
    const float3 _hj,
    const bool   not_selfGravity)
{
  const float3 dr = {posj.x - posi.x, posj.y - posi.y, posj.z - posi.z};
  const float  r2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;

  const float  rinv  = rsqrtf(r2);
  const float  rinv2 = rinv   * rinv;
  const float mrinv1 = rinv   * massj;
  const float mrinv3 = rinv2  * mrinv1;
  
  const float hi    = _hi.x;
  const float invhi = _hi.y;
  const float hj    = _hj.x;
  const float invhj = _hj.y;
  const float hflag = _hj.z;

  const float h2i = hi*hi;
  const float h2j = hj*hj;

  if (r2 > fmaxf(h2i, h2j))
  {
    iforce.ax  +=   mrinv3 * dr.x;
    iforce.ay  +=   mrinv3 * dr.y;
    iforce.az  +=   mrinv3 * dr.z;
    iforce.pot += (-mrinv1);
  } 
  else
  {
    const float2 invq  = not_selfGravity ? (float2){rinv * hi,   rinv * hj  } : (float2){0.0f, 0.0f};
    const float2 q     = not_selfGravity ? (float2){1.0f/invq.x, 1.0f/invq.y} : (float2){0.0f, 0.0f};
    const float2 q2    = {q.x*q .x, q.y*q. y};
    const float2 q3    = {q.x*q2.x, q.y*q2.y};
    const float2 invq2 = {invq.x*invq. x, invq.y*invq .y};
    const float2 invq3 = {invq.x*invq2.x, invq.y*invq2.y};
    
    const float2 f   = {q.x < 0.5f, q.y < 0.5f};
    const float2 acc = 
    {
        (       f.x) * (10.666666666667f + q2.x * (32.0f * q.x + (-38.4f))) + 
        (1.0f - f.x) * (21.333333333333f + (-48.0f)*q.x + 38.4f*q2.x + (-10.666666666667f)*q3.x + (-0.066666666667f) * invq3.x),
        (       f.y) * (10.666666666667f + q2.y * (32.0f * q.y + (-38.4f))) + 
        (1.0f - f.y) * (21.333333333333f + (-48.0f)*q.y + 38.4f*q2.y + (-10.666666666667f)*q3.y + (-0.066666666667f) * invq3.y)
    };
    
    const float2 pot =
    {
        (       f.x) * ((-2.8f) + q2.x * (5.333333333333f  + q2.x * (6.4f * q.x + (-9.6f)))) + 
        (1.0f - f.x) * ((-3.2f) + 0.066666666667f * invq.x + q2.x * (10.666666666667f + q.x * ((-16.0f) + q.x * (9.6f + (-2.133333333333f) * q.x)))),
        (       f.y) * ((-2.8f) + q2.y * (5.333333333333f  + q2.y * (6.4f * q.y + (-9.6f)))) + 
        (1.0f - f.y) * ((-3.2f) + 0.066666666667f * invq.y + q2.y * (10.666666666667f + q.y * ((-16.0f) + q.y * (9.6f + (-2.133333333333f) * q.y)))),
    };
     
    const float2 mj1 = {massj * invhi,       massj * invhj };
    const float2 mj2 = {mj1.x * invhi*invhi, mj1.y * invhj*invhj};

    const float2 g = {r2 <= h2i, r2 <= h2j};
    const float gacc = not_selfGravity ?
      0.5f*(g.x*mj2.x*acc.x + (1.0f - g.x)*mrinv3 + 
            g.y*mj2.y*acc.y + (1.0f - g.y)*mrinv3) : 0.0f;
    const float gpot = not_selfGravity ?
      0.5f*(g.x*mj1.x*pot.x + (g.x - 1.0f)*mrinv1 + 
            g.y*mj1.y*pot.y + (g.y - 1.0f)*mrinv1) : (-1.4f)*hflag*(mj1.x + mj1.y);


    iforce.ax  += gacc * dr.x;
    iforce.ay  += gacc * dr.y;
    iforce.az  += gacc * dr.z;
    iforce.pot += gpot;
  }

  return iforce;
}

__device__ float4 get_D04(const float ds2)
{
  const float ids = rsqrt(ds2);
  const float ids2 = ids *ids;
  const float ids3 = ids2*ids;  
  const float ids5 = ids3*ids2;
  const float ids7 = ids5*ids2;
  return (float4){+ids, -ids3, +3.0f*ids5, -15.0f*ids7};
} 

__device__ devForce body_node(
    devForce acc, 
    const float4 pos,
    const float mass, 
    const float3 com,
    const float3 Q0,  
    const float3 Q1)
{
  const float3  dr = {pos.x - com.x, pos.y - com.y, pos.z - com.z};
  const float   r2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
  const float4 D04 = get_D04(r2);

  const float  D0  = D04.x*mass;
  const float  D1  = D04.y*mass;
  const float  D2  = D04.z*mass;
  const float  D3  = D04.w*mass;

  const float oct_q11 = Q0.x;
  const float oct_q22 = Q0.y;
  const float oct_q33 = Q0.z;
  const float oct_q12 = Q1.x;
  const float oct_q13 = Q1.y;
  const float oct_q23 = Q1.z;

  const float Qii = oct_q11 + oct_q22 + oct_q33;
  const float QijRiRj =
    (oct_q11*dr.x*dr.x + oct_q22*dr.y*dr.y + oct_q33*dr.z*dr.z) +
    2.0f*(oct_q12*dr.y*dr.x + oct_q13*dr.z*dr.x + oct_q23*dr.y*dr.z);

  const float C01a = D1 + 0.5f*D2*Qii + 0.5f*D3*QijRiRj;
  acc.pot   +=     -(D0 + 0.5f*D1*Qii + 0.5f*D2*QijRiRj);
  acc.ax    += C01a*dr.x + D2*(oct_q11*dr.x + oct_q12*dr.y + oct_q13*dr.z);
  acc.ay    += C01a*dr.y + D2*(oct_q12*dr.x + oct_q22*dr.y + oct_q23*dr.z);
  acc.az    += C01a*dr.z + D2*(oct_q13*dr.x + oct_q23*dr.y + oct_q33*dr.z);

  return acc;
}

template<int DIM2, int SHIFT>
__device__ float4 approximate_gravity(
    int DIM2x, int DIM2y,
    int tid, 
    int tx, int ty,
    int    body_i, 
    float4 pos_i,
    float4 h_i,
    real4 group_pos,
    float eps2,
    uint2 node_begend,
    real4 *multipole_data,
    real4 *body_pos,
    int *shmem,
    int *lmem,

    int &apprCount, int &direCount,
    volatile float4 *boxSizeInfo,
    float4 groupSize,
    volatile float4 *boxCenterInfo,
    float group_eps,

    real4 *body_h,
    int *bodies_IDs) 
{
  devForce acc_i(0.0f);

  /*********** set necessary thread constants **********/

  const int DIMx = 1  << DIM2x;
  const int DIMy = 1  << DIM2y;
  const int DIM  = 1  << DIM2;
  const int offs = ty << DIM2x;

  /*********** shared memory distribution **********/

  //  begin,    end,   size
  // -----------------------
  int *approx = (int*)&shmem [     0];            //  0*DIM,  2*DIM,  2*DIM
  int *direct = (int*)&approx[ 2*DIM];            //  2*DIM,  3*DIM,  1*DIM
  int *nodes  = (int*)&direct[   DIM];            //  3*DIM, 13*DIM, 10*DIM
  int *prefix = (int*)&nodes [10*DIM];            // 13*DIM, 15*DIM,  2*DIM

  float  *node_mon0 = (float* )&nodes    [DIM];   //  4*DIM,  5*DIM,  1*DIM
  float3 *node_mon1 = (float3*)&node_mon0[DIM];   //  5*DIM,  8*DIM,  3*DIM
  float3 *node_oct0 = (float3*)&node_mon1[DIM];   //  8*DIM, 11*DIM,  3*DIM
  float3 *node_oct1 = (float3*)&node_oct0[DIM];   // 11*DIM, 14*DIM,  3*DIM

  int    *body_list = (int*   )&nodes    [  DIM]; //  4*DIM,  8*DIM,  4*DIM
  float  *sh_mass   = (float* )&body_list[4*DIM]; //  8*DIM,  9*DIM,  1*DIM
  float3 *sh_pos    = (float3*)&sh_mass  [  DIM]; //  9*DIM, 12*DIM   3*DIM
  int    *sh_jid    = (int*   )&sh_pos   [  DIM]; // 12*DIM, 13*DIM,  1*DIM
  float3 *sh_h      = (float3*)&sh_jid   [3*DIM]; // 15*DIM, 18*DIM,  3*DIM

  devForce *shForce = (devForce*)shmem;

  /*********** stack **********/

  int *nstack = lmem;

  /*********** begin tree-walk **********/

  int n_approx = 0;
  int n_direct = 0;

  for (int root_node = node_begend.x; root_node < node_begend.y; root_node += DIM) 
  {
    int n_nodes0 = min(node_begend.y - root_node, DIM);
    int n_stack0 = 0;
    int n_stack_pre = 0;

    { 
      nstack[ACCS<SHIFT>(n_stack0)] = root_node + tid;   
      n_stack0++; 
    }

    /*********** walk each level **********/
    while (n_nodes0 > 0) 
    {
      int n_nodes1 = 0;
      int n_offset = 0;

      int n_stack1 = n_stack0;
      int c_stack0 = n_stack_pre;

      /*********** walk a level **********/
      while(c_stack0 < n_stack0) 
      {
        /***
         **** --> fetch the list of nodes rom LMEM
         ***/
        bool use_node = tid <  n_nodes0;
        { 
          prefix[tid] = nstack[ACCS<SHIFT>(c_stack0)];   
          c_stack0++; 
        }
        __syncthreads();
        int node  = prefix[min(tid, n_nodes0 - 1)];

        if(n_nodes0 > 0)       //Work around pre 4.1 compiler bug
          n_nodes0 -= DIM;

        /***
         **** --> process each of the nodes in the list in parallel
         ***/

#ifndef TEXTURES
        float4 nodeSize = get_float4(boxSizeInfo[node]);                   //Fetch the size of the box. Size.w = child info
        float4 node_pos = get_float4(boxCenterInfo[node]);                 //Fetch the center of the box. center.w = opening info
#else
        float4 nodeSize =  tex1Dfetch(texNodeSize, node);
        float4 node_pos =  tex1Dfetch(texNodeCenter, node);
#endif

        int node_data = __float_as_int(nodeSize.w);

        //Check if a cell has to be opened
#if 0 // def IMPBH
        //Improved barnes hut method

#ifndef TEXTURES
        float4 nodeCOM = multipole_data[node*3];
#else
        float4 nodeCOM = tex1Dfetch(texMultipole,node*3);
#endif  

        nodeCOM.w      = node_pos.w;
        bool   split   = split_node_grav_impbh(nodeCOM, group_pos, groupSize);
#else   /* IMPBH */
        //Minimum distance method
        bool   split   = split_node_grav_md(node_pos, nodeSize, group_pos, groupSize);  //Check if node should be split
#endif /* IMPBH */

        bool leaf       = node_pos.w <= 0;  //Small AND equal incase of a 1 particle cell       //Check if it is a leaf
        //         split = true;


        uint mask    = BTEST((split && !leaf) && use_node);               // mask = #FFFFFFFF if use_node+split+not_a_leaf==true, otherwise zero
        int child    =    node_data & 0x0FFFFFFF;                         //Index to the first child of the node
        int nchild   = (((node_data & 0xF0000000) >> 28)) & mask;         //The number of children this node has


        /***
         **** --> calculate prefix
         ***/

        int *prefix0 = &prefix[  0];
        int *prefix1 = &prefix[DIM];

#ifdef OLDPREFIX
        int n_total = calc_prefix<DIM2>(prefix, tid,  nchild);
        prefix[tid] += n_offset - nchild;
        __syncthreads();
#else
        inclusive_scan_block<ADDOP<int>, int>(prefix, nchild, tid);        // inclusive scan to compute memory offset of each child
        int n_total = prefix[blockDim.x - 1];                              // fetch total number of children, i.e. offset of the last child -1
        __syncthreads();                                                   // thread barrier to make sure that warps completed their jobs
        prefix[tid] += n_offset - nchild;                                  // convert inclusive into exclusive scan for referencing purpose
        __syncthreads();                                                   // thread barrier
#endif

        for (int i = n_offset; i < n_offset + n_total; i += DIM)         //nullify part of the array that will be filled with children
          nodes[tid + i] = 0;                                          //but do not touch those parts which has already been filled
        __syncthreads();                                                 //Thread barrier to make sure all warps finished writing data

        bool flag = (split && !leaf) && use_node;                        //Flag = use_node + split + not_a_leaf;Use only non_leaf nodes that are to be split
        if (flag) nodes[prefix[tid]] = child;                            //Thread with the node that is about to be split
        __syncthreads();                                                 //writes the first child in the array of nodes

        /*** in the following 8 lines, we calculate indexes of all the children that have to be walked from the index of the first child***/
        if (flag && nodes[prefix[tid] + 1] == 0) nodes[prefix[tid] + 1] = child + 1; __syncthreads();
        if (flag && nodes[prefix[tid] + 2] == 0) nodes[prefix[tid] + 2] = child + 2; __syncthreads();
        if (flag && nodes[prefix[tid] + 3] == 0) nodes[prefix[tid] + 3] = child + 3; __syncthreads();
        if (flag && nodes[prefix[tid] + 4] == 0) nodes[prefix[tid] + 4] = child + 4; __syncthreads();
        if (flag && nodes[prefix[tid] + 5] == 0) nodes[prefix[tid] + 5] = child + 5; __syncthreads();
        if (flag && nodes[prefix[tid] + 6] == 0) nodes[prefix[tid] + 6] = child + 6; __syncthreads();
        if (flag && nodes[prefix[tid] + 7] == 0) nodes[prefix[tid] + 7] = child + 7; __syncthreads();

        n_offset += n_total;    //Increase the offset in the array by the number of newly added nodes


        /***
         **** --> save list of nodes to LMEM
         ***/

        /*** if half of shared memory or more is filled with the the nodes, dump these into slowmem stack ***/
        while(n_offset >= DIM) 
        {
          n_offset -= DIM;
          const int offs1 = ACCS<SHIFT>(n_stack1);
          nstack[offs1] = nodes[n_offset + tid];   n_stack1++;
          n_nodes1 += DIM;

          if((n_stack1 - c_stack0) >= (LMEM_STACK_SIZE << SHIFT))
          {
            apprCount = -1; 
            return acc_i.to_float4();
          }
        }

        __syncthreads();

        /******************************/
        /******************************/
        /*****     EVALUATION     *****/
        /******************************/
        /******************************/
#if 1
        /***********************************/
        /******       APPROX          ******/
        /***********************************/

#ifdef OLDPREFIX
        n_total = calc_prefix<DIM2>(prefix, tid,  1 - (split || !use_node));
#else
        inclusive_scan_block<ADDOP<int>, int>(prefix, 1 - (split || !use_node), tid);
        n_total = prefix[blockDim.x - 1];
#endif


        // 	n_total = calc_prefix<DIM2>(prefix, tid,  !split && use_node);         // for some unkown reason this does not work right on the GPU
        if (!split && use_node) approx[n_approx + prefix[tid] - 1] = node;
        __syncthreads();
        n_approx += n_total;

        while (n_approx >= DIM) 
        {
          n_approx -= DIM;
          int address      = (approx[n_approx + tid] << 1) + approx[n_approx + tid];
#ifndef TEXTURES
          float4 monopole  = multipole_data[address    ];
          float4 octopole0 = multipole_data[address + 1];
          float4 octopole1 = multipole_data[address + 2];
#else
          float4 monopole  = tex1Dfetch(texMultipole, address);
          float4 octopole0 = tex1Dfetch(texMultipole, address + 1);
          float4 octopole1 = tex1Dfetch(texMultipole, address + 2);
#endif

          node_mon0[tid] = monopole.w;
          node_mon1[tid] = (float3){monopole.x,  monopole.y,  monopole.z};
          node_oct0[tid] = (float3){octopole0.x, octopole0.y, octopole0.z};
          node_oct1[tid] = (float3){octopole1.x, octopole1.y, octopole1.z};

          __syncthreads();

#pragma unroll
          for (int i = 0; i < DIMx; i++)
          {
            apprCount++;
            acc_i = body_node(
                acc_i, pos_i,
                node_mon0[offs + i], node_mon1[offs + i],
                node_oct0[offs + i], node_oct1[offs + i]);
          }
          __syncthreads();
        }
        __syncthreads();
#endif

#if 1  /* DIRECT */
        /***********************************/
        /******       DIRECT          ******/
        /***********************************/

        int *sh_body = &approx[DIM];

        flag         = split && leaf && use_node;                                //flag = split + leaf + use_node
        int  jbody   = node_data & BODYMASK;                                     //the first body in the leaf
        int  nbody   = (((node_data & INVBMASK) >> LEAFBIT)+1) & BTEST(flag);    //number of bodies in the leaf masked with the flag

        body_list[tid] = direct[tid];                                            //copy list of bodies from previous pass to body_list
        sh_body  [tid] = jbody;                                                  //store the leafs first body id into shared memory

        // step 1
#ifdef OLDPREFIX
        calc_prefix<DIM2>(prefix0, tid, flag);
#else
        inclusive_scan_block<ADDOP<int>, int>(prefix0, (int)flag, tid);       // inclusive scan on flags to construct array
#endif

        if (flag) prefix1[prefix0[tid] - 1] = tid;                             //with tidś whose leaves have to be opened
        __syncthreads();                                                      //thread barrier, make sure all warps completed the job

        // step 2
#ifdef OLDPREFIX
        int n_bodies  = calc_prefix<DIM2>(prefix0, tid, nbody);
#else
        inclusive_scan_block<ADDOP<int>, int>(prefix0, nbody, tid);        // inclusive scan to compute memory offset for each body
        int n_bodies = prefix0[blockDim.x - 1];                            //Total number of bides extract from the leaves
        __syncthreads();                                                   // thread barrier to make sure that warps completed their jobs
#endif

        direct [tid]  = prefix0[tid];                                       //Store a copy of inclusive scan in direct
        prefix0[tid] -= nbody;                                              //convert inclusive int oexclusive scan
        prefix0[tid] += 1;                                                  //add unity, since later prefix0[tid] == 0 used to check barrier

        int nl_pre = 0;                                                     //Number of leaves that have already been processed

#define NJMAX (DIM*4)
        while (n_bodies > 0) 
        {
          int nb    = min(n_bodies, NJMAX - n_direct);                    //Make sure number of bides to be extracted does not exceed
          //the amount of allocated shared memory

          // step 0                                                      //nullify part of the body_list that will be filled with bodies
          for (int i = n_direct; i < n_direct + nb; i += DIM)            //from the leaves that are being processed
            body_list[i + tid] = 0;
          __syncthreads();

          //step 1:
          if (flag && (direct[tid] <= nb) && (prefix0[tid] > 0))        //make sure that the thread indeed carries a leaf
            body_list[n_direct + prefix0[tid] - 1] = 1;                 //whose bodies will be extracted
          __syncthreads();

          //step 2:
#ifdef OLDPREFIX
          int nl = calc_prefix<DIM2>(nb, &body_list[n_direct], tid);
#else
          int nl = inclusive_scan_array<ADDOP<int>, int>              // inclusive scan to compute number of leaves to process
            (&body_list[n_direct], nb, tid);                          // to make sure that there is enough shared memory for bodies
#endif
          nb = direct[prefix1[nl_pre + nl - 1]];                        // number of bodies stored in these leaves

          // step 3:
          for (int i = n_direct; i < n_direct + nb; i += DIM) 
          {                                                              //segmented fill of the body_list
            int j = prefix1[nl_pre + body_list[i + tid] - 1];            // compute the first body in shared j-body array
            body_list[i + tid] = (i + tid - n_direct) -                 //add to the index of the first j-body in a child
              (prefix0[j] - 1) + sh_body[j];         //the index of the first child in body_list array
          }
          __syncthreads();


          /**************************************************
           *  example of what is accomplished in steps 0-4   *
           *       ---------------------------               *
           * step 0: body_list = 000000000000000000000       *
           * step 1: body_list = 100010001000000100100       *
           * step 2: body_list = 111122223333333444555       *
           * step 3: body_list = 012301230123456012012       *
           *         assuming that sh_body[j] = 0            *
           ***************************************************/

          n_bodies     -= nb;                                   //subtract from n_bodies number of bodies that have been extracted
          nl_pre       += nl;                                   //increase the number of leaves that where processed
          direct [tid] -= nb;                                   //subtract the number of extracted bodies in this pass
          prefix0[tid] = max(prefix0[tid] - nb, 0);             //same here, but do not let the number be negative (GT200 bug!?)
          n_direct     += nb;                                  //increase the number of bodies to be procssed

          while(n_direct >= DIM) 
          {
            n_direct -= DIM;

            const float4 posj  = body_pos[body_list[n_direct + tid]];
            const float4 hj    = body_h  [body_list[n_direct + tid]];
            sh_mass[tid] = posj.w;
            sh_pos [tid] = (float3){posj.x, posj.y, posj.z};
            sh_h   [tid] = (float3){hj  .x, hj  .y, hj  .z};  /* h, hinv, hflag */
            sh_jid [tid] = bodies_IDs[body_list[n_direct + tid]];
            __syncthreads();

#pragma unroll
            for (int j = 0; j < DIMx; j++)
            {
              direCount++;
              acc_i = body_body(
                  acc_i, 
                  pos_i, 
                  h_i,
                  sh_mass[offs + j], 
                  sh_pos [offs + j], 
                  sh_h   [offs + j],
                  body_i != sh_jid[offs + j]);
            }
            __syncthreads();
          }
        }
        direct[tid] = body_list[tid];
        __syncthreads();
#endif /* DIRECT */
      } //end of the tree level, proceed to the next one


      n_nodes1 += n_offset;
      if (n_offset > 0)
      { 
        nstack[ACCS<SHIFT>(n_stack1)] = nodes[tid];   n_stack1++; 
        if((n_stack1 - c_stack0) >= (LMEM_STACK_SIZE << SHIFT))
        {
          //We overwrote our current stack
          apprCount = -1; 
          return acc_i.to_float4();
        }
      }
      __syncthreads();


      /***
       **** --> copy nodes1 to nodes0: done by reassigning the pointers
       ***/
      n_nodes0    = n_nodes1;

      n_stack_pre = n_stack0;
      n_stack0    = n_stack1;

    }//end while   levels
  }//end for


  if(n_approx > 0)
  {

    if (tid < n_approx) {
      int address      = (approx[tid] << 1) + approx[tid];
#ifndef TEXTURES
      float4 monopole  = multipole_data[address    ];
      float4 octopole0 = multipole_data[address + 1];
      float4 octopole1 = multipole_data[address + 2];
#else
      float4 monopole  = tex1Dfetch(texMultipole, address);
      float4 octopole0 = tex1Dfetch(texMultipole, address + 1);
      float4 octopole1 = tex1Dfetch(texMultipole, address + 2);
#endif

      node_mon0[tid] = monopole.w;
      node_mon1[tid] = (float3){monopole.x,  monopole.y,  monopole.z};
      node_oct0[tid] = (float3){octopole0.x, octopole0.y, octopole0.z};
      node_oct1[tid] = (float3){octopole1.x, octopole1.y, octopole1.z};

    } else {
    
      //Set non-active memory locations to zero
      node_mon0[tid] = 0.0f;
      node_mon1[tid] = (float3){1.0e10f, 1.0e10f, 1.0e10f};
      node_oct0[tid] = (float3){0.0f, 0.0f, 0.0f};
      node_oct1[tid] = (float3){0.0f, 0.0f, 0.0f};
    }
    __syncthreads();

#pragma unroll
    for (int i = 0; i < DIMx; i++)
    {
      apprCount++;
      acc_i = body_node(acc_i, pos_i,
          node_mon0[offs + i], node_mon1[offs + i],
          node_oct0[offs + i], node_oct1[offs + i]);
    }

    __syncthreads();
  } //if n_approx > 0

  if(n_direct > 0)
  {
    if (tid < n_direct) 
    {
      const float4 posj = body_pos[direct[tid]];
      const float4 hj   = body_h  [direct[tid]];
      sh_mass[tid] = posj.w;
      sh_pos [tid] = (float3){posj.x, posj.y, posj.z};
      sh_h   [tid] = (float3){hj  .x, hj  .y, hj  .z};  /* h, hinv, hflag */
      sh_jid [tid] = bodies_IDs[direct[tid]];

    } 
    else 
    {
      sh_mass[tid] = 0.0f;
      sh_pos [tid] = (float3){1.0e10f, 1.0e10f, 1.0e10f};
      sh_jid [tid] = -1;
      sh_h   [tid] = (float3){0.0f, 0.0f, 0.0f};
    }
    __syncthreads();

#pragma unroll
    for (int j = 0; j < DIMx; j++) 
      if ((sh_jid[offs + j] >= 0)) 
      {
        direCount++;
        acc_i = body_body(
            acc_i, 
            pos_i, 
            h_i,
            sh_mass[offs + j], 
            sh_pos [offs + j], 
            sh_h   [offs + j],
            body_i != sh_jid[offs + j]);
      }
    __syncthreads();
  }

  shForce[tid] = acc_i;
  __syncthreads();
  if (ty == 0)
    for (int i = 1; i < DIMy; i++) 
    {
      const int idx = (i << DIM2x) + tx;
      acc_i.ax  += shForce[idx].ax .to_float();
      acc_i.ay  += shForce[idx].ay .to_float();
      acc_i.az  += shForce[idx].az .to_float();
#ifdef _GPOTDP_
      acc_i.pot += shForce[idx].pot.to_float();
#else
      acc_i.pot += shForce[idx].pot;
#endif
    }
  __syncthreads();

  return acc_i.to_float4();
}


  extern "C" __global__ void
  __launch_bounds__(NTHREAD)
dev_approximate_gravity(
    const int n_active_groups,
    float eps2,
    uint2 node_begend,   
    int    *atomicValues,
    real4  *body_pos,                        
    real4  *body_h,
    float4 *acc_out,
    real4  *group_body_pos, 
    real4  *group_body_h,
    real   *ds2_out,
    int    *ngb_out,
    int    *active_inout,
    int2   *interactions,
    uint2  *group_list,
    real4  *multipole_data,
    float4  *boxSizeInfo,                        
    float4  *boxCenterInfo,                        
    int     *MEM_BUF,
    int     *bodies_IDs,
    int     *group_bodies_IDs) 
{


  const int blockDim2 = NTHREAD2;
  __shared__ int shmem[18*(1 << blockDim2)];
  //    __shared__ int shmem[24*(1 << blockDim2)]; is possible on FERMI
  //    int             lmem[LMEM_STACK_SIZE];



  /*********** check if this block is linked to a leaf **********/

  int bid = gridDim.x * blockIdx.y + blockIdx.x;

  while(true)
  {

    if(threadIdx.x == 0)
    {
      bid         = atomicAdd(&atomicValues[0], 1);
      shmem[0]    = bid;
    }
    __syncthreads();

    bid   = shmem[0];

    if (bid >= n_active_groups) return;


    int tid = threadIdx.y * blockDim.x + threadIdx.x;

    //   volatile int *lmem = &MEM_BUF[blockIdx.x*LMEM_STACK_SIZE*blockDim.x + threadIdx.x*LMEM_STACK_SIZE];
    //   int *lmem = &MEM_BUF[blockIdx.x*LMEM_STACK_SIZE*blockDim.x + threadIdx.x*LMEM_STACK_SIZE];
    int *lmem = &MEM_BUF[blockIdx.x*LMEM_STACK_SIZE*blockDim.x];


    /*********** set necessary thread constants **********/

    //   real4 curGroupSize    = groupSizeInfo[active_groups[bid + grpOffset]];
    //   int   groupData       = __float_as_int(curGroupSize.w);
    //   uint body_i           =   groupData & CRITMASK;
    //   uint nb_i             = ((groupData & INVCMASK) >> CRITBIT) + 1;
    // 
    //   real4 group_pos       = groupCenterInfo[active_groups[bid + grpOffset]];

    //   if(tid == 0)
    //   printf("[%f %f %f %f ] \n [%f %f %f %f ] %d %d \n",
    //           curGroupSize.x, curGroupSize.y, curGroupSize.z, curGroupSize.w,
    //           group_pos.x, group_pos.y, group_pos.z, group_pos.w, body_i, nb_i);
    uint2 grpInfo = group_list[bid];
    uint body_i = grpInfo.x;
    uint nb_i   = (grpInfo.y - grpInfo.x) + 1;

    int DIM2x = 0;
    while (((nb_i - 1) >> DIM2x) > 0) DIM2x++;

    DIM2x     = max(DIM2x,4);
    int DIM2y = blockDim2 - DIM2x;

    int tx = tid & ((1 << DIM2x) - 1);
    int ty = tid >> DIM2x;

    body_i += tx%nb_i;
   
    float4 acc_i = {0.0f, 0.0f, 0.0f, 0.0f};

    const float4 pos_i = group_body_pos[body_i];
    const float4 h_i   = group_body_h  [body_i];

    real4 group_pos;
    real4 curGroupSize;

    computeGroupProps(group_pos, curGroupSize, pos_i, shmem);
    float group_eps = 0;        //This is disabled for the moment

    int ngb_i = -1;
    float ds2 = -1;

    int apprCount = 0;
    int direCount = 0;

    const int body_id = group_bodies_IDs[body_i];

#if 1
    acc_i = approximate_gravity<blockDim2, 0>(
        DIM2x, DIM2y, tid, tx, ty,
        body_id, pos_i, h_i, group_pos,
        eps2, node_begend,
        multipole_data, body_pos,
        shmem, lmem, 
        apprCount, direCount, boxSizeInfo, curGroupSize, boxCenterInfo,
        group_eps, 
        body_h,
        bodies_IDs);
#endif
#if 0
    if (apprCount < 0)
    {
      printf(" --------- Error:: body_i= %d \n", body_i);
    }
#else
    if(apprCount < 0)
    {

      //Try to get access to the big stack, only one block per time is allowed
      if(threadIdx.x == 0)
      {
        int res = atomicExch(&atomicValues[1], 1); //If the old value (res) is 0 we can go otherwise sleep
        int waitCounter  = 0;
        while(res != 0)
        {
          //Sleep
          for(int i=0; i < (1024); i++)
          {
            waitCounter += 1;
          }
          //Test again
          shmem[0] = waitCounter;
          res = atomicExch(&atomicValues[1], 1); 
        }
      }

      __syncthreads();

      lmem = &MEM_BUF[gridDim.x*LMEM_STACK_SIZE*blockDim.x];    //Use the extra large buffer
      apprCount = direCount = 0;
      acc_i = approximate_gravity<blockDim2, 8>( DIM2x, DIM2y, tid, tx, ty,
          body_id, pos_i, h_i, group_pos,
          eps2, node_begend,
          multipole_data, body_pos,
          shmem, lmem, 
          apprCount, direCount, boxSizeInfo, curGroupSize, boxCenterInfo,
          group_eps, 
          body_h,
          bodies_IDs);

      lmem = &MEM_BUF[blockIdx.x*LMEM_STACK_SIZE*blockDim.x]; //Back to normal location

      if(threadIdx.x == 0)
      {
        atomicExch(&atomicValues[1], 0); //Release the lock
      }
    }//end if apprCount < 0
#endif


    if (tid < nb_i) 
    {
      acc_out     [body_i] = acc_i;
      ngb_out     [body_i] = ngb_i;
      ds2_out     [body_i] = ds2;
      active_inout[body_i] = 1;
      interactions[body_i].x = apprCount;
      interactions[body_i].y = direCount ;
    }
  }     //end while
}


