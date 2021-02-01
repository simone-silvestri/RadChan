/*******************************************************************//**
*
*     NARROWBAND MONTE CARLO OERM GPU by Simone Silvestri
*          
*     runs on p_col*num_gpu GPUs in p_col nodes
*     
*     imax*jmax*kmax --> domain size (in x-, y- and z-direction, respectively)
*  
*     num_gpu    --> number of GPUs available in one node 
*
*     p_col      --> column discretization of the DNS code,
*                    must match with the number of nodes employed!
*
*     Lx, Ly, Lz --> length of the domain 
*
*     fact_mesh  --> parameter for grid streching in x-direction
*
*     nVar       --> number of independent MC runs for variance calculation
*                    (set to 1 if variance on a particular grid has already been assessed)
*
*     photmax    --> number of statistical sampled per cell in one MC run
*
*     grid_num   --> number of coarsened grids employed for MC simulations
*                    The mesh dimension are specified in maxi[],maxj[],maxk[]
*                    each mesh should be an integer multiple of its coarser version
*  
*     maxs[]     --> number of steps employed in each grid 
*                    the last one should be large enought to ensure total absorption
*
*     num_steams --> number of streams used for the domain MC routine
*                    should divide exactly imax
*
*     num_steams --> number of streams used for the surface MC routine
*                    should divide exactly jmax
*
*     blockdim   --> thread number in a block.
*                    keep it an integer multiple of the warpsize (32)
*
*     nblocks    --> minimum number of blocks employed. 
*     	             The actual number of blocks used is calculated 
*     	             as MIN(nblocks, Domain/blockdim) 
*
*     srt        --> if equal to 1 it triggers the narrow-band sorting 
*                    procedure
*
*     grey       --> if equal to 1 ignores the narrow-band formulation and
*                    performs a grey gas MC simulation
*
*     calcwalls  --> if equal to 1 (with at least one of the adiab option enabled)
*                    calculates the radiative heat flux at the enables wall with
*                    a grey/narrow-band OERMC formulation
*
*     floatpoint --> if equal to 1 single precision, if equal to 0 double precision
*                    texture memory is not yet implemented in double precision,
*                    thus always keep this equal to 1 (if want to change to 0 the 
*                    memory_copy.cu file needs an update). Advisable anyway equal to 1
*                    for speed related issues. If a higher accuracy is required, it is 
*                    advised to increase photmax.
*
*     epsx       --> emissivity of wall x (where x=w,e,n,s,t,b)
*                    w,e: west , east   (i=0,i=imax+1)
*                    t,b: top  , bottom (k=0,k=kmax+1)
*                    n,s: north, south   (j=0,j=jmax+1)
*
*     Twx        --> temperature of wall x
*                    Used only if adiabx == 0               
*
*     adiabx     --> if equal to 1, then Twx is taken from the temperature field and
*                    is not a parameter. At least one among w,e should be 1 to enable
*                    the calculation of radiative heat flux at the walls
*
*     toll,tollS --> tollerances for ray marching for domain and surface, respectively
*	             (higher than 1e-3 is not advisable)
*
*     bdx        --> flag for boundary condition for wall x
*                    1 : periodic wall
*                    2 : black wall
*                    3 : specularly reflecting wall
*                    4 : adiabatic black wall conditions (implemented only for outlet)
*                    5 : lost ray (temperature equal to the emitter) 
*
*     The constraints on the parameters are:
*
*     - imax must be an integer multiple of num_streams
*     - jmax must be a multiple integer of num_gpu
*     - kmax must be a multiple integer of p_col
*     - if wall radiative flux is to be calculated (calcwalls==1) then
*       jmax must be also an integer multiple of num_strS
*     - num_strS must be smaller than num_streams
*
**********************************************************************/
