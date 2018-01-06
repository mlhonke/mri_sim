################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CU_SRCS += \
../demo_scan_full.cu 

CU_DEPS += \
./demo_scan_full.d 

OBJS += \
./demo_scan_full.o 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-8.0/bin/nvcc -I/home/omega/Dev/external_sources/cub-1.5.2 -I"/home/omega/Dev/dSimL/acquisition" -I"/home/omega/Dev/dSimL/primitives" -I"/home/omega/Dev/dSimL/sequence" -I/usr/local/cuda/lib64/ -O3 -optf /home/omega/Dev/dSimL/nvcc_options -std=c++11 -gencode arch=compute_60,code=sm_60  -odir "." -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-8.0/bin/nvcc -I/home/omega/Dev/external_sources/cub-1.5.2 -I"/home/omega/Dev/dSimL/acquisition" -I"/home/omega/Dev/dSimL/primitives" -I"/home/omega/Dev/dSimL/sequence" -I/usr/local/cuda/lib64/ -O3 -optf /home/omega/Dev/dSimL/nvcc_options -std=c++11 --compile --relocatable-device-code=true -gencode arch=compute_60,code=compute_60 -gencode arch=compute_60,code=sm_60  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


