#include "__cf_sldemo_mdlref_counter.h"
#ifndef RTW_HEADER_sldemo_mdlref_counter_cap_host_h_
#define RTW_HEADER_sldemo_mdlref_counter_cap_host_h_
#ifdef HOST_CAPI_BUILD
#include "rtw_capi.h"
#include "rtw_modelmap.h"
typedef struct { rtwCAPI_ModelMappingInfo mmi ; }
sldemo_mdlref_counter_host_DataMapInfo_T ;
#ifdef __cplusplus
extern "C" {
#endif
void sldemo_mdlref_counter_host_InitializeDataMapInfo (
sldemo_mdlref_counter_host_DataMapInfo_T * dataMap , const char * path ) ;
#ifdef __cplusplus
}
#endif
#endif
#endif
