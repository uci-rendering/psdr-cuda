#pragma once

#include <optix.h>
#include <stdexcept>
#include <sstream>

//------------------------------------------------------------------------------
//
// OptiX error-checking
//
//------------------------------------------------------------------------------

#define OPTIX_CHECK( call )                                                    \
    do                                                                         \
    {                                                                          \
        OptixResult res = call;                                                \
        if( res != OPTIX_SUCCESS )                                             \
        {                                                                      \
            std::stringstream ss;                                              \
            ss << "Optix call '" << #call << "' failed: " __FILE__ ":"         \
               << __LINE__ << ")\n";                                           \
            throw psdr::Exception( res, ss.str().c_str() );                    \
        }                                                                      \
    } while( 0 )


#define OPTIX_CHECK_LOG( call )                                                \
    do                                                                         \
    {                                                                          \
        OptixResult res = call;                                                \
        if( res != OPTIX_SUCCESS )                                             \
        {                                                                      \
            std::stringstream ss;                                              \
            ss << "Optix call '" << #call << "' failed: " __FILE__ ":"         \
               << __LINE__ << ")\nLog:\n" << log                               \
               << ( sizeof_log > sizeof( log ) ? "<TRUNCATED>" : "" )          \
               << "\n";                                                        \
            throw psdr::Exception( res, ss.str().c_str() );                    \
        }                                                                      \
    } while( 0 )


//------------------------------------------------------------------------------
//
// CUDA error-checking
//
//------------------------------------------------------------------------------

#define CUDA_CHECK( call )                                                     \
    do                                                                         \
    {                                                                          \
        cudaError_t error = call;                                              \
        if( error != cudaSuccess )                                             \
        {                                                                      \
            std::stringstream ss;                                              \
            ss << "CUDA call (" << #call << " ) failed with error: '"          \
               << cudaGetErrorString( error )                                  \
               << "' (" __FILE__ << ":" << __LINE__ << ")\n";                  \
            throw psdr::Exception( ss.str().c_str() );                         \
        }                                                                      \
    } while( 0 )


#define CUDA_SYNC_CHECK()                                                      \
    do                                                                         \
    {                                                                          \
        cudaDeviceSynchronize();                                               \
        cudaError_t error = cudaGetLastError();                                \
        if( error != cudaSuccess )                                             \
        {                                                                      \
            std::stringstream ss;                                              \
            ss << "CUDA error on synchronize with error '"                     \
               << cudaGetErrorString( error )                                  \
               << "' (" __FILE__ << ":" << __LINE__ << ")\n";                  \
            throw psdr::Exception( ss.str().c_str() );                         \
        }                                                                      \
    } while( 0 )


//------------------------------------------------------------------------------
//
// Assertions
//
//------------------------------------------------------------------------------

#define PSDR_ASSERT( cond )                                                    \
    do                                                                         \
    {                                                                          \
        if( !(cond) )                                                          \
        {                                                                      \
            std::stringstream ss;                                              \
            ss << __FILE__ << " (" << __LINE__ << "): " << #cond;              \
            throw psdr::Exception( ss.str().c_str() );                         \
        }                                                                      \
    } while( 0 )


#define PSDR_ASSERT_MSG( cond, msg )                                           \
    do                                                                         \
    {                                                                          \
        if( !(cond) )                                                          \
        {                                                                      \
            std::stringstream ss;                                              \
            ss << "\n  File \"" << __FILE__ << "\", line " << __LINE__;        \
            throw psdr::Exception( ( std::string(msg) + ss.str() ).c_str() );  \
        }                                                                      \
    } while( 0 )



namespace psdr
{

class Exception : public std::runtime_error
{
 public:
     Exception( const char* msg )
         : std::runtime_error( msg )
     { }

     Exception( OptixResult res, const char* msg )
         : std::runtime_error( createMessage( res, msg ).c_str() )
     { }

 private:
     std::string createMessage( OptixResult res, const char* msg )
     {
         std::ostringstream out;
         out << optixGetErrorName( res ) << ": " << msg;
         return out.str();
     }
};

} // namespace psdr
