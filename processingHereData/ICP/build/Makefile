
NVCC        = nvcc

NVCC_FLAGS  = -I/usr/local/cuda/include -gencode=arch=compute_50,code=\"sm_50,compute_50\"
ifdef dbg
	NVCC_FLAGS  += -g 
else
	NVCC_FLAGS  += -O2
endif

LD_FLAGS    = -lcudart -L/usr/local/cuda/lib64
EXE	        = icp
OBJ	        = icp_cu.o icp_cpp.o

default: $(EXE)

icp_cu.o: icp_main.cu icp_kernel.cu globals.h
	$(NVCC) -c -o $@ matrixmul.cu $(NVCC_FLAGS)

icp_cpp.o: icp_gold.cpp globals.h
	$(NVCC) -c -o $@ icp_gold.cpp $(NVCC_FLAGS) 

$(EXE): $(OBJ)
	$(NVCC) $(OBJ) -o $(EXE) $(LD_FLAGS) $(NVCC_FLAGS)

clean:
	rm -rf *.o $(EXE)
