################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CU_SRCS += \
../scanner/k_space.cu \
../scanner/scanner.cu 

CU_DEPS += \
./scanner/k_space.d \
./scanner/scanner.d 

OBJS += \
./scanner/k_space.o \
./scanner/scanner.o 


# Each subdirectory must supply rules for building sources it contributes
scanner/%.o: ../scanner/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-8.0/bin/nvcc -I/home/omega/Dev/external_sources/cub-1.5.2 -I/usr/local/cuda/lib64/ -G -g -lineinfo -O0 -optf /home/omega/Dev/dSimL/nvcc_options -std=c++11 -gencode arch=compute_60,code=sm_60  -odir "scanner" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-8.0/bin/nvcc -I/home/omega/Dev/external_sources/cub-1.5.2 -I/usr/local/cuda/lib64/ -G -g -lineinfo -O0 -optf /home/omega/Dev/dSimL/nvcc_options -std=c++11 --compile --relocatable-device-code=true -gencode arch=compute_60,code=compute_60 -gencode arch=compute_60,code=sm_60  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


