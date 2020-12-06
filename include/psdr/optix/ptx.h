#pragma once

#if defined(_WIN32) || defined(_WIN64)
#	define PSDRAPI __declspec(dllimport)
#   define PSDRCLASSAPI
#else
#	define PSDRAPI __attribute__ ((visibility ("default")))
#	define PSDRCLASSAPI PSDRAPI
#endif

#define OPTIX_SAMPLE_NAME_STRINGIFY2(name) #name
#define OPTIX_SAMPLE_NAME_STRINGIFY(name) OPTIX_SAMPLE_NAME_STRINGIFY2(name)
#define OPTIX_SAMPLE_NAME OPTIX_SAMPLE_NAME_STRINGIFY(OPTIX_SAMPLE_NAME_DEFINE)

namespace psdr
{

PSDRAPI const char* getPtxString(const char* sample, const char* filename, const char** log = NULL );

} // end namespace PSDR
