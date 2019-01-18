#include "__cf_sldemo_mdlref_counter.h"
#ifndef RTW_HEADER_sldemo_mdlref_counter_h_
#define RTW_HEADER_sldemo_mdlref_counter_h_
#include <string.h>
#include <stddef.h>
#include "rtw_modelmap.h"
#ifndef sldemo_mdlref_counter_COMMON_INCLUDES_
#define sldemo_mdlref_counter_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "simstruc.h"
#include "fixedpoint.h"
#endif
#include "sldemo_mdlref_counter_types.h"
#include "multiword_types.h"
#include "rt_nonfinite.h"
typedef struct { real_T eantbbong1 ; } bc1s13imfa ; struct j43akkmfxbg_ {
real_T P_0 ; } ; struct epvdpbhwrw { struct SimStruct_tag * _mdlRefSfcnS ;
struct { rtwCAPI_ModelMappingInfo mmi ; rtwCAPI_ModelMapLoggingInstanceInfo
mmiLogInstanceInfo ; void * dataAddress [ 1 ] ; int32_T * vardimsAddress [ 1
] ; sysRanDType * systemRan [ 3 ] ; int_T systemTid [ 3 ] ; } DataMapInfo ;
struct { int_T mdlref_GlobalTID [ 1 ] ; } Timing ; } ; typedef struct {
bc1s13imfa rtdw ; oi2ouxwxvu rtm ; } j0v54a0zxtl ; extern void b51pkukdad (
SimStruct * _mdlRefSfcnS , int_T mdlref_TID0 , oi2ouxwxvu * const ilxj1cme0p
, bc1s13imfa * localDW , void * sysRanPtr , int contextTid ,
rtwCAPI_ModelMappingInfo * rt_ParentMMI , const char_T * rt_ChildPath , int_T
rt_ChildMMIIdx , int_T rt_CSTATEIdx ) ; extern void
mr_sldemo_mdlref_counter_MdlInfoRegFcn ( SimStruct * mdlRefSfcnS , char_T *
modelName , int_T * retVal ) ; extern mxArray *
mr_sldemo_mdlref_counter_GetDWork ( const j0v54a0zxtl * mdlrefDW ) ; extern
void mr_sldemo_mdlref_counter_SetDWork ( j0v54a0zxtl * mdlrefDW , const
mxArray * ssDW ) ; extern void
mr_sldemo_mdlref_counter_RegisterSimStateChecksum ( SimStruct * S ) ; extern
mxArray * mr_sldemo_mdlref_counter_GetSimStateDisallowedBlocks ( ) ; extern
const rtwCAPI_ModelMappingStaticInfo * sldemo_mdlref_counter_GetCAPIStaticMap
( void ) ; extern void i0gdxss3tt ( bc1s13imfa * localDW ) ; extern void
f3gnqdmgzm ( bc1s13imfa * localDW ) ; extern void hdnimqzru5 ( real_T *
p5szrpovik , bc1s13imfa * localDW ) ; extern void sldemo_mdlref_counter (
const real_T * podswudgkn , const real_T * cz0atirrtv , const real_T *
jhx1i2w1lb , real_T * p5szrpovik , bc1s13imfa * localDW ) ; extern void
hagwau5kpb ( oi2ouxwxvu * const ilxj1cme0p ) ;
#endif
