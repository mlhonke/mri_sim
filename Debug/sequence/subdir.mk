################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CU_SRCS += \
../sequence/GRE.cu \
../sequence/PGSE_D.cu \
../sequence/pulses.cu \
../sequence/sequence.cu 

CU_DEPS += \
./sequence/GRE.d \
./sequence/PGSE_D.d \
./sequence/pulses.d \
./sequence/sequence.d 

OBJS += \
./sequence/GRE.o \
./sequence/PGSE_D.o \
./sequence/pulses.o \
./sequence/sequence.o 


# Each subdirectory must supply rules for building sources it contributes
sequence/%.o: ../sequence/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-8.0/bin/nvcc -I/home/omega/Dev/external_sources/cub-1.5.2 -I/usr/local/cuda/lib64/ -G -g -lineinfo -O0 -optf /home/omega/Dev/dSimL/nvcc_options -std=c++11 -gencode arch=compute_60,code=sm_60  -odir "sequence" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-8.0/bin/nvcc -I/home/omega/Dev/external_sources/cub-1.5.2 -I/usr/local/cuda/lib64/ -G -g -lineinfo -O0 -optf /home/omega/Dev/dSimL/nvcc_options -std=c++11 --compile --relocatable-device-code=true -gencode arch=compute_60,code=compute_60 -gencode arch=compute_60,code=sm_60  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


