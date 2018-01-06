################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../util/recorder.cpp 

CU_SRCS += \
../util/cudaVector.cu \
../util/pinnedVector.cu 

CU_DEPS += \
./util/cudaVector.d \
./util/pinnedVector.d 

OBJS += \
./util/cudaVector.o \
./util/pinnedVector.o \
./util/recorder.o 

CPP_DEPS += \
./util/recorder.d 


# Each subdirectory must supply rules for building sources it contributes
util/%.o: ../util/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-8.0/bin/nvcc -I/home/omega/Dev/external_sources/cub-1.5.2 -I/usr/local/cuda/lib64/ -G -g -lineinfo -O0 -optf /home/omega/Dev/dSimL/nvcc_options -std=c++11 -gencode arch=compute_60,code=sm_60  -odir "util" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-8.0/bin/nvcc -I/home/omega/Dev/external_sources/cub-1.5.2 -I/usr/local/cuda/lib64/ -G -g -lineinfo -O0 -optf /home/omega/Dev/dSimL/nvcc_options -std=c++11 --compile --relocatable-device-code=true -gencode arch=compute_60,code=compute_60 -gencode arch=compute_60,code=sm_60  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

util/%.o: ../util/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-8.0/bin/nvcc -I/home/omega/Dev/external_sources/cub-1.5.2 -I/usr/local/cuda/lib64/ -G -g -lineinfo -O0 -optf /home/omega/Dev/dSimL/nvcc_options -std=c++11 -gencode arch=compute_60,code=sm_60  -odir "util" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-8.0/bin/nvcc -I/home/omega/Dev/external_sources/cub-1.5.2 -I/usr/local/cuda/lib64/ -G -g -lineinfo -O0 -optf /home/omega/Dev/dSimL/nvcc_options -std=c++11 --compile  -x c++ -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


