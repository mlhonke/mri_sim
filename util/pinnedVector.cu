#ifndef _PINNED_VECTOR_CU
#define _PINNED_VECTOR_CU

#include "pinnedVector.cuh"
	
template <class T>
pinnedVector<T>::pinnedVector(std::size_t size){
  _size = size; cudaHostAlloc( (void**)&pinned_ptr, _size*sizeof(T),cudaHostAllocDefault );
}
	
template <class T>
pinnedVector<T>::pinnedVector(std::size_t size, T & val){
  _size = size; 
  cudaHostAlloc( (void**)&pinned_ptr, _size*sizeof(T),cudaHostAllocDefault );
  for (int i = 0; i < _size; i++){pinned_ptr[i] = val;}
}
	
template <class T>
pinnedVector<T>::~pinnedVector(){
  cudaFreeHost(pinned_ptr);
}
	
template <class T>	
std::size_t pinnedVector<T>::size(){
  return _size;
}

template <class T>
void pinnedVector<T>::alloc(std::size_t size){
_size = size; cudaHostAlloc( (void**)&pinned_ptr, _size*sizeof(T),cudaHostAllocDefault );
}

template <class T>
void pinnedVector<T>::alloc(std::size_t size, T & val){
  _size = size; 
  cudaHostAlloc( (void**)&pinned_ptr, _size*sizeof(T),cudaHostAllocDefault );
  for (int i = 0; i < _size; i++){pinned_ptr[i] = val;}
}
	
template <class T>
T& pinnedVector<T>::operator [] (int i){
  return pinned_ptr[i];
}
		
template <class T>
void pinnedVector<T>::copyToDevice(cudaVector<T> & dev, cudaStream_t & stream){
	cudaMemcpyAsync( dev.getPointer(), pinned_ptr, _size*sizeof(T), cudaMemcpyHostToDevice, stream );
}
	
template <class T>
void pinnedVector<T>::copyFromDevice(cudaVector<T> & dev, cudaStream_t & stream){
  cudaMemcpyAsync( pinned_ptr, dev.getPointer(), _size*sizeof(T), cudaMemcpyDeviceToHost, stream );
}
		
template <class T>	
void pinnedVector<T>::operator= (std::vector<T> & v){
  int iter;
  if (v.size() >= _size){ iter = _size; }
  else {iter = v.size();}
		
  for (int i = 0; i < iter; i++){
    pinned_ptr[i] = v[i];
  }

}
	
template <class T>	
void pinnedVector<T>::operator= (cudaVector<T> & v){
	
  int iter;
  if (v.size() >= _size){ iter = _size; }
  else {iter = v.size();}
	
  for (int i = 0; i < iter; i++){
    pinned_ptr[i] = v[i];
  }
}

template <class T>
void pinnedVector<T>::copyTo(std::vector<T> & v){

  int iter;
  if (v.size() >= _size){ iter = _size; }
  else {iter = v.size();}
	
  for (int i = 0; i < iter; i++){
    v[i] = pinned_ptr[i];
  }
}
	
template <class T>	
T* pinnedVector<T>::getPointer(){
  return pinned_ptr;
}
		
#endif
