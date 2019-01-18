#include "__cf_sldemo_mdlref_counter.h"
#include "rtw_capi.h"
#ifdef HOST_CAPI_BUILD
#include "sldemo_mdlref_counter_capi_host.h"
#define sizeof(s) ((size_t)(0xFFFF))
#undef rt_offsetof
#define rt_offsetof(s,el) ((uint16_T)(0xFFFF))
#define TARGET_CONST
#define TARGET_STRING(s) (s)    
#else
#include "builtin_typeid_types.h"
#include "sldemo_mdlref_counter.h"
#include "sldemo_mdlref_counter_capi.h"
#include "sldemo_mdlref_counter_private.h"
#ifdef LIGHT_WEIGHT_CAPI
#define TARGET_CONST                  
#define TARGET_STRING(s)               (NULL)                    
#else
#define TARGET_CONST                   const
#define TARGET_STRING(s)               (s)
#endif
#endif
static rtwCAPI_Signals rtBlockSignals [ ] = { { 0 , 0 , ( NULL ) , ( NULL ) ,
0 , 0 , 0 , 0 , 0 } } ; static rtwCAPI_States rtBlockStates [ ] = { { 0 , - 1
, TARGET_STRING ( "sldemo_mdlref_counter/Previous Output" ) , TARGET_STRING (
"DSTATE" ) , "" , 0 , 0 , 0 , 0 , 0 , 0 , - 1 , 0 } , { 0 , - 1 , ( NULL ) ,
( NULL ) , ( NULL ) , 0 , 0 , 0 , 0 , 0 , 0 , - 1 , 0 } } ;
#ifndef HOST_CAPI_BUILD
static void sldemo_mdlref_counter_InitializeDataAddr ( void * dataAddr [ ] ,
bc1s13imfa * localDW ) { dataAddr [ 0 ] = ( void * ) ( & localDW ->
eantbbong1 ) ; }
#endif
#ifndef HOST_CAPI_BUILD
static void sldemo_mdlref_counter_InitializeVarDimsAddr ( int32_T *
vardimsAddr [ ] ) { vardimsAddr [ 0 ] = ( NULL ) ; }
#endif
static TARGET_CONST rtwCAPI_DataTypeMap rtDataTypeMap [ ] = { { "double" ,
"real_T" , 0 , 0 , sizeof ( real_T ) , SS_DOUBLE , 0 , 0 } } ;
#ifdef HOST_CAPI_BUILD
#undef sizeof
#endif
static TARGET_CONST rtwCAPI_ElementMap rtElementMap [ ] = { { ( NULL ) , 0 ,
0 , 0 , 0 } , } ; static rtwCAPI_DimensionMap rtDimensionMap [ ] = { {
rtwCAPI_SCALAR , 0 , 2 , 0 } } ; static uint_T rtDimensionArray [ ] = { 1 , 1
} ; static const real_T rtcapiStoredFloats [ ] = { 0.2 , 0.0 } ; static
rtwCAPI_FixPtMap rtFixPtMap [ ] = { { ( NULL ) , ( NULL ) ,
rtwCAPI_FIX_RESERVED , 0 , 0 , 0 } , } ; static rtwCAPI_SampleTimeMap
rtSampleTimeMap [ ] = { { ( const void * ) & rtcapiStoredFloats [ 0 ] , (
const void * ) & rtcapiStoredFloats [ 1 ] , 0 , 0 } } ; static int_T
rtContextSystems [ 3 ] ; static rtwCAPI_LoggingMetaInfo loggingMetaInfo [ ] =
{ { 0 , 0 , "" , 0 } } ; static rtwCAPI_ModelMapLoggingStaticInfo
mmiStaticInfoLogging = { 3 , rtContextSystems , loggingMetaInfo , 0 , NULL ,
{ 0 , NULL , NULL } , 0 , ( NULL ) } ; static rtwCAPI_ModelMappingStaticInfo
mmiStatic = { { rtBlockSignals , 0 , ( NULL ) , 0 , ( NULL ) , 0 } , { ( NULL
) , 0 , ( NULL ) , 0 } , { rtBlockStates , 1 } , { rtDataTypeMap ,
rtDimensionMap , rtFixPtMap , rtElementMap , rtSampleTimeMap ,
rtDimensionArray } , "float" , { 1647128993U , 362101413U , 1586741126U ,
3924366976U } , & mmiStaticInfoLogging , 0 , 0 } ; const
rtwCAPI_ModelMappingStaticInfo * sldemo_mdlref_counter_GetCAPIStaticMap (
void ) { return & mmiStatic ; }
#ifndef HOST_CAPI_BUILD
static void sldemo_mdlref_counter_InitializeSystemRan ( oi2ouxwxvu * const
ilxj1cme0p , sysRanDType * systemRan [ ] , bc1s13imfa * localDW , int_T
systemTid [ ] , void * rootSysRanPtr , int rootTid ) { UNUSED_PARAMETER (
ilxj1cme0p ) ; UNUSED_PARAMETER ( localDW ) ; systemRan [ 0 ] = ( sysRanDType
* ) rootSysRanPtr ; systemRan [ 1 ] = ( NULL ) ; systemRan [ 2 ] = ( NULL ) ;
systemTid [ 1 ] = rootTid ; systemTid [ 2 ] = rootTid ; systemTid [ 0 ] =
rootTid ; rtContextSystems [ 0 ] = 0 ; rtContextSystems [ 1 ] = 0 ;
rtContextSystems [ 2 ] = 0 ; }
#endif
#ifndef HOST_CAPI_BUILD
void sldemo_mdlref_counter_InitializeDataMapInfo ( oi2ouxwxvu * const
ilxj1cme0p , bc1s13imfa * localDW , void * sysRanPtr , int contextTid ) {
rtwCAPI_SetVersion ( ilxj1cme0p -> DataMapInfo . mmi , 1 ) ;
rtwCAPI_SetStaticMap ( ilxj1cme0p -> DataMapInfo . mmi , & mmiStatic ) ;
rtwCAPI_SetLoggingStaticMap ( ilxj1cme0p -> DataMapInfo . mmi , &
mmiStaticInfoLogging ) ; sldemo_mdlref_counter_InitializeDataAddr (
ilxj1cme0p -> DataMapInfo . dataAddress , localDW ) ;
rtwCAPI_SetDataAddressMap ( ilxj1cme0p -> DataMapInfo . mmi , ilxj1cme0p ->
DataMapInfo . dataAddress ) ; sldemo_mdlref_counter_InitializeVarDimsAddr (
ilxj1cme0p -> DataMapInfo . vardimsAddress ) ; rtwCAPI_SetVarDimsAddressMap (
ilxj1cme0p -> DataMapInfo . mmi , ilxj1cme0p -> DataMapInfo . vardimsAddress
) ; rtwCAPI_SetPath ( ilxj1cme0p -> DataMapInfo . mmi , ( NULL ) ) ;
rtwCAPI_SetFullPath ( ilxj1cme0p -> DataMapInfo . mmi , ( NULL ) ) ;
rtwCAPI_SetInstanceLoggingInfo ( ilxj1cme0p -> DataMapInfo . mmi , &
ilxj1cme0p -> DataMapInfo . mmiLogInstanceInfo ) ; rtwCAPI_SetChildMMIArray (
ilxj1cme0p -> DataMapInfo . mmi , ( NULL ) ) ; rtwCAPI_SetChildMMIArrayLen (
ilxj1cme0p -> DataMapInfo . mmi , 0 ) ;
sldemo_mdlref_counter_InitializeSystemRan ( ilxj1cme0p , ilxj1cme0p ->
DataMapInfo . systemRan , localDW , ilxj1cme0p -> DataMapInfo . systemTid ,
sysRanPtr , contextTid ) ; rtwCAPI_SetSystemRan ( ilxj1cme0p -> DataMapInfo .
mmi , ilxj1cme0p -> DataMapInfo . systemRan ) ; rtwCAPI_SetSystemTid (
ilxj1cme0p -> DataMapInfo . mmi , ilxj1cme0p -> DataMapInfo . systemTid ) ;
rtwCAPI_SetGlobalTIDMap ( ilxj1cme0p -> DataMapInfo . mmi , & ilxj1cme0p ->
Timing . mdlref_GlobalTID [ 0 ] ) ; }
#else
#ifdef __cplusplus
extern "C" {
#endif
void sldemo_mdlref_counter_host_InitializeDataMapInfo (
sldemo_mdlref_counter_host_DataMapInfo_T * dataMap , const char * path ) {
rtwCAPI_SetVersion ( dataMap -> mmi , 1 ) ; rtwCAPI_SetStaticMap ( dataMap ->
mmi , & mmiStatic ) ; rtwCAPI_SetDataAddressMap ( dataMap -> mmi , NULL ) ;
rtwCAPI_SetVarDimsAddressMap ( dataMap -> mmi , NULL ) ; rtwCAPI_SetPath (
dataMap -> mmi , path ) ; rtwCAPI_SetFullPath ( dataMap -> mmi , NULL ) ;
rtwCAPI_SetChildMMIArray ( dataMap -> mmi , ( NULL ) ) ;
rtwCAPI_SetChildMMIArrayLen ( dataMap -> mmi , 0 ) ; }
#ifdef __cplusplus
}
#endif
#endif
