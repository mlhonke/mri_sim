################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CU_SRCS += \
../acquisition/magAcquisitionStream.cu \
../acquisition/magAcquisitionStreamCPU.cu 

CU_DEPS += \
./acquisition/magAcquisitionStream.d \
./acquisition/magAcquisitionStreamCPU.d 

OBJS += \
./acquisition/magAcquisitionStream.o \
./acquisition/magAcquisitionStreamCPU.o 


# Each subdirectory must supply rules for building sources it contributes
acquisition/%.o: ../acquisition/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-8.0/bin/nvcc -I/home/omega/Dev/external_sources/cub-1.5.2 -I/usr/local/cuda/lib64/ -G -g -lineinfo -O0 -optf /home/omega/Dev/dSimL/nvcc_options -std=c++11 -gencode arch=compute_60,code=sm_60  -odir "acquisition" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-8.0/bin/nvcc -I/home/omega/Dev/external_sources/cub-1.5.2 -I/usr/local/cuda/lib64/ -G -g -lineinfo -O0 -optf /home/omega/Dev/dSimL/nvcc_options -std=c++11 --compile --relocatable-device-code=true -gencode arch=compute_60,code=compute_60 -gencode arch=compute_60,code=sm_60  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


