#include <fstream>
#include <iterator>
#include <map>
#ifdef _WIN32
#   include <windows.h>
#else
#   include <dirent.h>
#endif
#include <misc/Exception.h>
#include <psdr/optix/ptx.h>

namespace psdr
{

static bool fileExists( const char* path )
{
    std::ifstream str( path );
    return static_cast<bool>( str );
}

static bool fileExists( const std::string& path )
{
    return fileExists( path.c_str() );
}

static bool readSourceFile( std::string& str, const std::string& filename )
{
    // Try to open file
    std::ifstream file( filename.c_str() );
    if( file.good() )
    {
        // Found usable source file
        std::stringstream source_buffer;
        source_buffer << file.rdbuf();
        str = source_buffer.str();
        return true;
    }
    return false;
}

static std::string samplePTXFilePath( const char* sampleName, const char* fileName )
{
    // Allow for overrides.
    static const char* directories[] =
    {
        PTX_OUTPUT_DIR
    };
    for( const char* directory : directories )
    {
        if( directory )
        {
            std::string path = directory;
            path += '/';
            path += "ptx";
            path += "_generated_";
            path += fileName;
            path += ".ptx";
            if( fileExists( path ) )
                return path;
        }
    }

    std::string error = "samplePTXFilePath couldn't locate ";
    error += fileName;
    error += " for sample ";
    error += sampleName;
    throw Exception( error.c_str() );
}

static void getPtxStringFromFile( std::string& ptx, const char* sample_name, const char* filename )
{
    const std::string sourceFilePath = samplePTXFilePath( sample_name, filename );

    // Try to open source PTX file
    if( !readSourceFile( ptx, sourceFilePath ) )
    {
        std::string err = "Couldn't open source file " + sourceFilePath;
        throw std::runtime_error( err.c_str() );
    }
}

struct PtxSourceCache
{
    std::map<std::string, std::string*> map;
    ~PtxSourceCache()
    {
        for( std::map<std::string, std::string*>::const_iterator it = map.begin(); it != map.end(); ++it )
            delete it->second;
    }
};
static PtxSourceCache g_ptxSourceCache;

const char* getPtxString( const char* sample, const char* filename, const char** log )
{
    if( log )
        *log = NULL;

    std::string *                                 ptx, cu;
    std::string                                   key  = std::string( filename ) + ";" + ( sample ? sample : "" );
    std::map<std::string, std::string*>::iterator elem = g_ptxSourceCache.map.find( key );

    if( elem == g_ptxSourceCache.map.end() )
    {
        ptx = new std::string();
        getPtxStringFromFile( *ptx, sample, filename );
        g_ptxSourceCache.map[key] = ptx;
    }
    else
    {
        ptx = elem->second;
    }

    return ptx->c_str();
}

} // namespace sutil
