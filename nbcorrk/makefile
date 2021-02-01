FLAGS  = -c -mcmodel=large #-CB -fpe0 -traceback#O5   
PROF =  
LIB =  -lm 
COMP  = g++  
ckdist:   ckdist_hitran_part.o correlatedK.o read_files.o calc_prob.o alloc_bands.o write_tables.o
	$(COMP) $(DBG) ckdist_hitran_part.o correlatedK.o read_files.o calc_prob.o alloc_bands.o write_tables.o $(LIB) -o ckdist 
ckdist_hitran_part.o: ckdist_hitran_part.cpp 
	$(COMP) $(DBG) $(FLAGS) ckdist_hitran_part.cpp
write_tables.o: write_tables.cpp 
	$(COMP) $(DBG) $(FLAGS) write_tables.cpp
calc_prob.o: calc_prob.cpp 
	$(COMP) $(DBG) $(FLAGS) calc_prob.cpp
correlatedK.o: correlatedK.cpp 
	$(COMP) $(DBG) $(FLAGS) correlatedK.cpp
read_files.o: read_files.cpp 
	$(COMP) $(DBG) $(FLAGS) read_files.cpp
alloc_bands.o: alloc_bands.cpp 
	$(COMP) $(DBG) $(FLAGS) alloc_bands.cpp
clean:
	rm *.o ckdist
