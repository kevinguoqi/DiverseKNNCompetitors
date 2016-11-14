################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Dataset.cpp \
../src/DistanceFunction.cpp \
../src/DiverseKNNComparison.cpp \
../src/Evaluation.cpp \
../src/IDDistance.cpp \
../src/KNDN.cpp \
../src/Point.cpp \
../src/WeightedFunction.cpp 

OBJS += \
./src/Dataset.o \
./src/DistanceFunction.o \
./src/DiverseKNNComparison.o \
./src/Evaluation.o \
./src/IDDistance.o \
./src/KNDN.o \
./src/Point.o \
./src/WeightedFunction.o 

CPP_DEPS += \
./src/Dataset.d \
./src/DistanceFunction.d \
./src/DiverseKNNComparison.d \
./src/Evaluation.d \
./src/IDDistance.d \
./src/KNDN.d \
./src/Point.d \
./src/WeightedFunction.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -D__GXX_EXPERIMENTAL_CXX0X__ -O3 -Wall -c -fmessage-length=0 -std=c++0x -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


