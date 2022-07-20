#pragma once

#include <cuda_runtime.h>
#include <optix.h>
#include <optix_function_table_definition.h>
#include <optix_stubs.h>

#include <misc/Exception.h>
#include <psdr/optix/ptx.h>
#include <cuda/psdr_cuda.h>

#include <stdint.h>
#include <iomanip>
#include <psdr/macros.h>

template <typename IntegerType>
static inline IntegerType roundUp(IntegerType x, IntegerType y) {
    return ( ( x + y - 1 ) / y ) * y;
}

struct PathTracerState
{
    OptixDeviceContext             context                  = 0;
    OptixTraversableHandle         gas_handle               = 0;  // Traversable handle for triangle AS
    CUdeviceptr                    d_gas_output_buffer      = 0;  // Triangle AS memory
    CUdeviceptr                    d_vertices               = 0;
    OptixModule                    ptx_module               = 0;
    OptixPipelineCompileOptions    pipeline_compile_options = {};
    OptixPipeline                  pipeline                 = 0;
    OptixProgramGroup              raygen_prog_group        = 0;
    OptixProgramGroup              radiance_miss_group      = 0;
    OptixProgramGroup              radiance_hit_group       = 0;
    CUstream                       stream                   = 0;
    Params                         params;
    Params*                        d_params                 = nullptr;
    OptixShaderBindingTable        sbt                      = {};
};

template <typename T>
struct Record
{
    __align__( OPTIX_SBT_RECORD_ALIGNMENT ) char header[OPTIX_SBT_RECORD_HEADER_SIZE];
    T data;
};

typedef Record<RayGenData>   RayGenRecord;
typedef Record<MissData>     MissRecord;
typedef Record<HitGroupData> HitGroupRecord;

static void context_log_cb( unsigned int level, const char* tag, const char* message, void* /*cbdata */ )
{
    std::cerr << "[" << std::setw( 2 ) << level << "][" << std::setw( 12 ) << tag << "]: " << message << "\n";
}

void optix_release( PathTracerState& state )
{
    OPTIX_CHECK( optixPipelineDestroy( state.pipeline ) );
    OPTIX_CHECK( optixProgramGroupDestroy( state.raygen_prog_group ) );
    OPTIX_CHECK( optixProgramGroupDestroy( state.radiance_miss_group ) );
    OPTIX_CHECK( optixProgramGroupDestroy( state.radiance_hit_group ) );
    OPTIX_CHECK( optixModuleDestroy( state.ptx_module ) );
    OPTIX_CHECK( optixDeviceContextDestroy( state.context ) );

    CUDA_CHECK( cudaFree( reinterpret_cast<void*>( state.sbt.raygenRecord ) ) );
    CUDA_CHECK( cudaFree( reinterpret_cast<void*>( state.sbt.missRecordBase ) ) );
    CUDA_CHECK( cudaFree( reinterpret_cast<void*>( state.sbt.hitgroupRecordBase ) ) );
    // CUDA_CHECK( cudaFree( reinterpret_cast<void*>( state.d_gas_output_buffer ) ) );
    cuda_free(reinterpret_cast<void*>(state.d_gas_output_buffer));
    CUDA_CHECK( cudaFree( reinterpret_cast<void*>( state.d_params ) ) );
}

void optix_config(PathTracerState& state, const std::vector<int> &face_offset) {

    // ------------------------
    //  OptiX context creation
    // ------------------------
    CUDA_CHECK( cudaFree( 0 ) );
    CUcontext cu_ctx                  = 0;
    OPTIX_CHECK( optixInit() );
    OptixDeviceContextOptions options = {};
    options.logCallbackFunction       = &context_log_cb;
#ifdef PSDR_OPTIX_DEBUG
    options.logCallbackLevel          = 4;
#else
    options.logCallbackLevel          = 0;
#endif
    OPTIX_CHECK( optixDeviceContextCreate( cu_ctx, &options, &state.context ) );
    char   log[2048];
    size_t sizeof_log = sizeof( log );


    // ----------------------------------------------
    //  Pipeline generation - Create Module from PTX
    // ----------------------------------------------
    OptixModuleCompileOptions module_compile_options = {};
    module_compile_options.maxRegisterCount  = OPTIX_COMPILE_DEFAULT_MAX_REGISTER_COUNT;

#if defined(PSDR_OPTIX_DEBUG)
    module_compile_options.optLevel         = OPTIX_COMPILE_OPTIMIZATION_LEVEL_0;
    module_compile_options.debugLevel       = OPTIX_COMPILE_DEBUG_LEVEL_FULL;
#else
    module_compile_options.optLevel         = OPTIX_COMPILE_OPTIMIZATION_DEFAULT;
    module_compile_options.debugLevel       = OPTIX_COMPILE_DEBUG_LEVEL_NONE;
#endif

    state.pipeline_compile_options.usesMotionBlur        = false;
    state.pipeline_compile_options.traversableGraphFlags = OPTIX_TRAVERSABLE_GRAPH_FLAG_ALLOW_SINGLE_GAS;
    state.pipeline_compile_options.numPayloadValues      = 2;
    state.pipeline_compile_options.numAttributeValues    = 2;
    state.pipeline_compile_options.exceptionFlags        = OPTIX_EXCEPTION_FLAG_NONE;
    state.pipeline_compile_options.pipelineLaunchParamsVariableName = "params";

    const std::string ptx = psdr::getPtxString( OPTIX_SAMPLE_NAME, (std::string(PSDR_CUDA_FILE) + ".cu").c_str() );

    OPTIX_CHECK_LOG( optixModuleCreateFromPTX(
                state.context,
                &module_compile_options,
                &state.pipeline_compile_options,
                ptx.c_str(),
                ptx.size(),
                log,
                &sizeof_log,
                &state.ptx_module
                ) );


    // ---------------------------------------------
    //  Pipeline generation - Create program groups
    // ---------------------------------------------
    OptixProgramGroupOptions  program_group_options = {};
    {
        OptixProgramGroupDesc raygen_prog_group_desc    = {};
        raygen_prog_group_desc.kind                     = OPTIX_PROGRAM_GROUP_KIND_RAYGEN;
        raygen_prog_group_desc.raygen.module            = state.ptx_module;
        raygen_prog_group_desc.raygen.entryFunctionName = "__raygen__psdr_rg";

        OPTIX_CHECK_LOG( optixProgramGroupCreate(
                    state.context, &raygen_prog_group_desc,
                    1,  // num program groups
                    &program_group_options,
                    log,
                    &sizeof_log,
                    &state.raygen_prog_group
                    ) );
    }

    {
        OptixProgramGroupDesc moe_miss_prog_group_desc  = {};
        moe_miss_prog_group_desc.kind                   = OPTIX_PROGRAM_GROUP_KIND_MISS;
        moe_miss_prog_group_desc.miss.module            = state.ptx_module;
        moe_miss_prog_group_desc.miss.entryFunctionName = "__miss__psdr_ms";
        OPTIX_CHECK_LOG( optixProgramGroupCreate(
                    state.context, &moe_miss_prog_group_desc,
                    1,  // num program groups
                    &program_group_options,
                    log, &sizeof_log,
                    &state.radiance_miss_group
                    ) );
    }

    {
        OptixProgramGroupDesc hit_prog_group_desc        = {};
        hit_prog_group_desc.kind                         = OPTIX_PROGRAM_GROUP_KIND_HITGROUP;
        hit_prog_group_desc.hitgroup.moduleCH            = state.ptx_module;
        hit_prog_group_desc.hitgroup.entryFunctionNameCH = "__closesthit__psdr_ch";
        sizeof_log                                       = sizeof( log );
        OPTIX_CHECK_LOG( optixProgramGroupCreate(
                    state.context,
                    &hit_prog_group_desc,
                    1,  // num program groups
                    &program_group_options,
                    log,
                    &sizeof_log,
                    &state.radiance_hit_group
                    ) );
    }


    // ---------------------------------------
    //  Pipeline generation - Create pipeline
    // ---------------------------------------
    OptixProgramGroup program_groups[] =
    {
        state.raygen_prog_group,
        state.radiance_miss_group,
        state.radiance_hit_group,
    };

    OptixPipelineLinkOptions pipeline_link_options  = {};
    pipeline_link_options.maxTraceDepth             = 1;
#if defined(PSDR_OPTIX_DEBUG)
    pipeline_link_options.debugLevel                = OPTIX_COMPILE_DEBUG_LEVEL_FULL;
#else
    pipeline_link_options.debugLevel                = OPTIX_COMPILE_DEBUG_LEVEL_NONE;
#endif

    OPTIX_CHECK_LOG( optixPipelineCreate(
                state.context,
                &state.pipeline_compile_options,
                &pipeline_link_options,
                program_groups,
                sizeof( program_groups ) / sizeof( program_groups[0] ),
                log,
                &sizeof_log,
                &state.pipeline
                ) );


    // ---------------------------------
    //  Shader Binding Table generation
    // ---------------------------------
    CUdeviceptr  d_raygen_record;
    const size_t raygen_record_size = sizeof( RayGenRecord );
    CUDA_CHECK( cudaMalloc( reinterpret_cast<void**>( &d_raygen_record ), raygen_record_size ) );

    RayGenRecord rg_sbt = {};
    OPTIX_CHECK( optixSbtRecordPackHeader( state.raygen_prog_group, &rg_sbt ) );

    CUDA_CHECK( cudaMemcpy(
                reinterpret_cast<void*>( d_raygen_record ),
                &rg_sbt,
                raygen_record_size,
                cudaMemcpyHostToDevice
                ) );


    CUdeviceptr  d_miss_records;
    const size_t miss_record_size = sizeof( MissRecord );
    CUDA_CHECK( cudaMalloc( reinterpret_cast<void**>( &d_miss_records ), miss_record_size ) );

    MissRecord ms_sbt;
    OPTIX_CHECK( optixSbtRecordPackHeader( state.radiance_miss_group,  &ms_sbt ) );

    CUDA_CHECK( cudaMemcpy(
                reinterpret_cast<void*>( d_miss_records ),
                &ms_sbt,
                miss_record_size,
                cudaMemcpyHostToDevice
                ) );

    std::vector<HitGroupRecord> hg_sbts;
    for( int i = 0; i < static_cast<int>(face_offset.size()); ++i )
    {
        HitGroupRecord rec = {};
        rec.data.shape_offset = face_offset[i];
        rec.data.shape_id = i;
        OPTIX_CHECK( optixSbtRecordPackHeader( state.radiance_hit_group, &rec ) );
        hg_sbts.push_back( rec );
    }

    CUdeviceptr  d_hitgroup_records;
    const size_t hitgroup_record_size = sizeof( HitGroupRecord );
    CUDA_CHECK( cudaMalloc(
                reinterpret_cast<void**>( &d_hitgroup_records ),
                hitgroup_record_size*hg_sbts.size()
                ) );

    CUDA_CHECK( cudaMemcpy(
                reinterpret_cast<void*>( d_hitgroup_records ),
                hg_sbts.data(),
                hitgroup_record_size*hg_sbts.size(),
                cudaMemcpyHostToDevice
                ) );

    state.sbt.raygenRecord                = d_raygen_record;
    state.sbt.missRecordBase              = d_miss_records;
    state.sbt.missRecordStrideInBytes     = static_cast<uint32_t>( miss_record_size );
    state.sbt.missRecordCount             = 1;
    state.sbt.hitgroupRecordBase          = d_hitgroup_records;
    state.sbt.hitgroupRecordStrideInBytes = sizeof(HitGroupRecord);
    state.sbt.hitgroupRecordCount         = hg_sbts.size();

    CUDA_CHECK( cudaStreamCreate( &state.stream ) );
    CUDA_CHECK( cudaMalloc( reinterpret_cast<void**>( &state.d_params ), sizeof( Params ) ) );
}

void build_accel(PathTracerState& state, const std::vector<OptixBuildInput>& triangle_input) {
    int num_build_inputs = static_cast<int>(triangle_input.size());
    OptixAccelBuildOptions accel_options = {};
    accel_options.buildFlags             = OPTIX_BUILD_FLAG_ALLOW_COMPACTION;
    accel_options.operation              = OPTIX_BUILD_OPERATION_BUILD;
    accel_options.motionOptions.numKeys = 0;

    OptixAccelBufferSizes buffer_sizes;

    OPTIX_CHECK( optixAccelComputeMemoryUsage(
                state.context,
                &accel_options,
                triangle_input.data(),
                num_build_inputs,
                &buffer_sizes
                ) );

    void* d_temp_buffer = cuda_malloc(buffer_sizes.tempSizeInBytes);
    void* output_buffer = cuda_malloc(buffer_sizes.outputSizeInBytes + 8);

    OptixAccelEmitDesc emit_property = {};
    emit_property.type   = OPTIX_PROPERTY_TYPE_COMPACTED_SIZE;
    emit_property.result = (CUdeviceptr)((char*)output_buffer + buffer_sizes.outputSizeInBytes);

    OPTIX_CHECK( optixAccelBuild(
                state.context,
                0,
                &accel_options,
                triangle_input.data(),
                num_build_inputs,
                (CUdeviceptr)d_temp_buffer,
                buffer_sizes.tempSizeInBytes,
                (CUdeviceptr)output_buffer,
                buffer_sizes.outputSizeInBytes,
                &state.gas_handle,
                &emit_property,
                1
                ) );

    cuda_free(d_temp_buffer);

    size_t compact_size;
    CUDA_CHECK( cudaMemcpy( &compact_size, reinterpret_cast<void*>(emit_property.result), sizeof(size_t), cudaMemcpyDeviceToHost ) );

        if (compact_size < buffer_sizes.outputSizeInBytes) {
            void* compact_buffer = cuda_malloc(compact_size);
            OPTIX_CHECK(optixAccelCompact(
                state.context,
                0, // CUDA stream
                state.gas_handle,
                (CUdeviceptr)compact_buffer,
                compact_size,
                &state.gas_handle
            ));
            cuda_free(output_buffer);
            output_buffer = compact_buffer;
        }

    // build CUDA stream
    state.params.handle         = state.gas_handle;

    if (state.d_gas_output_buffer)
        cuda_free(reinterpret_cast<void*>(state.d_gas_output_buffer));
    state.d_gas_output_buffer = (CUdeviceptr)output_buffer;
}
