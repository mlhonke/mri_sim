#ifndef _PINNEDVECTOR_H_
#define _PINNEDVECTOR_H_

#include <vector>

template <class X> class cudaVector;
template <class x> class cudaScalar;

  template <class T>
  class pinnedVector {

  private:
    std::size_t _size;
    T* pinned_ptr;
		
  public:

    pinnedVector(){}
	pinnedVector(std::size_t);
	pinnedVector(std::size_t, T&);
    ~pinnedVector();
	void alloc(std::size_t);
	void alloc(std::size_t, T&);
    std::size_t size();
    T& operator [] (int);
    void copyToDevice(cudaVector<T> &, cudaStream_t &);
    void copyFromDevice(cudaVector<T> &, cudaStream_t &);
    void operator= (std::vector<T> &);
    void operator= (cudaVector<T> &);
	void copyTo(std::vector<T> &);
    T* getPointer();
		
  };
  
  
 template <class T>	
class pinnedScalar {

    T* pinned_ptr;

	public:
	
	pinnedScalar(T & val){ 
		 cudaHostAlloc( (void**)&pinned_ptr, sizeof(T),cudaHostAllocDefault );
		 *pinned_ptr = val;
	}	
	
	pinnedScalar(){ 
		 cudaHostAlloc( (void**)&pinned_ptr, sizeof(T),cudaHostAllocDefault );
	}
		

	void copyToDevice(cudaScalar<T> & dev, cudaStream_t & stream){
	  cudaMemcpyAsync( dev.getPointer(), pinned_ptr, sizeof(T), cudaMemcpyHostToDevice, stream );
	}
		

	void copyFromDevice(cudaScalar<T> & dev, cudaStream_t & stream){
	  cudaMemcpyAsync( pinned_ptr, dev.getPointer(), sizeof(T), cudaMemcpyDeviceToHost, stream );
	}
		

void operator= (T v){

	*pinned_ptr = v;

}



void copyTo(T & v){

	v = *pinned_ptr;
	
}
	

T* getPointer(){
  return pinned_ptr;
}

}; 

//MH: Include the implementation of the above methods, since the compiler doesn't instantiate until it knows the type.
#include "pinnedVector.cu"

#endif
