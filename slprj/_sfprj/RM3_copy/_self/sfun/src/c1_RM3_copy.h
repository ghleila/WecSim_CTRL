#ifndef __c1_RM3_copy_h__
#define __c1_RM3_copy_h__

/* Type Definitions */
#ifndef struct_tag_sI3an6DNyfQAOCkXY61B1nE
#define struct_tag_sI3an6DNyfQAOCkXY61B1nE

struct tag_sI3an6DNyfQAOCkXY61B1nE
{
  real_T StepRatio;
};

#endif                                 /*struct_tag_sI3an6DNyfQAOCkXY61B1nE*/

#ifndef typedef_c1_sI3an6DNyfQAOCkXY61B1nE
#define typedef_c1_sI3an6DNyfQAOCkXY61B1nE

typedef struct tag_sI3an6DNyfQAOCkXY61B1nE c1_sI3an6DNyfQAOCkXY61B1nE;

#endif                                 /*typedef_c1_sI3an6DNyfQAOCkXY61B1nE*/

#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T

struct emxArray_real_T
{
  real_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray_real_T*/

#ifndef typedef_c1_emxArray_real_T
#define typedef_c1_emxArray_real_T

typedef struct emxArray_real_T c1_emxArray_real_T;

#endif                                 /*typedef_c1_emxArray_real_T*/

#ifndef typedef_SFc1_RM3_copyInstanceStruct
#define typedef_SFc1_RM3_copyInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c1_sfEvent;
  boolean_T c1_doneDoubleBufferReInit;
  uint8_T c1_is_active_c1_RM3_copy;
  real_T c1_y[150000];
  real_T c1_dv13[150000];
  real_T c1_u[150000];
  real_T c1_m_data[5000];
  real_T c1_range_data[5000];
  real_T c1_Out_data[5000];
  real_T c1_tmp_data[5000];
  real_T c1_b_y[15000];
  real_T *c1_B2_address;
  int32_T c1__index;
  boolean_T c1_dsmdiag_B2;
  real_T (*c1__indexFinal_Rew_matt_address)[1000];
  int32_T c1__indexRew_matt;
  boolean_T c1_dsmdiag_Final_Rew_matt;
  real_T *c1__indexK2_address;
  int32_T c1_b__index;
  boolean_T c1_dsmdiag_K2;
  real_T (*c1__indexMeanR_address)[10000];
  int32_T c1_c__index;
  boolean_T c1_dsmdiag_MeanR;
  real_T *c1__indexMemory_D_size_address;
  int32_T c1__index_D_size;
  boolean_T c1_dsmdiag_Memory_D_size;
  real_T (*c1__indexMemory_Stor_address)[150000];
  int32_T c1__index_Stor;
  boolean_T c1_dsmdiag_Memory_Stor;
  real_T *c1__indexNN_dd_address;
  int32_T c1_d__index;
  boolean_T c1_dsmdiag_NN_dd;
  real_T (*c1__indexQ_address)[864];
  int32_T c1_e__index;
  boolean_T c1_dsmdiag_Q;
  real_T (*c1__indexQ_prime_address)[864];
  int32_T c1__indexe;
  boolean_T c1_dsmdiag_Q_prime;
  real_T (*c1__indexQ_targ_address)[864];
  int32_T c1_f__index;
  boolean_T c1_dsmdiag_Q_targ;
  real_T (*c1__indexR_address)[3600];
  int32_T c1_g__index;
  boolean_T c1_dsmdiag_R;
  real_T (*c1__indexStor_eta_address)[15000];
  int32_T c1__indexta;
  boolean_T c1_dsmdiag_Stor_eta;
  real_T (*c1__indexaction_address)[24];
  int32_T c1_h__index;
  boolean_T c1_dsmdiag_action;
  int16_T (*c1__indexnum_R_address)[10000];
  int32_T c1_i__index;
  boolean_T c1_dsmdiag_num_R;
  real_T *c1__indexove_count_address;
  int32_T c1__indexunt;
  boolean_T c1_dsmdiag_ove_count;
  int16_T *c1__indexstate_idx_address;
  int32_T c1__indexidx;
  boolean_T c1_dsmdiag_state_idx;
  real_T (*c1__indexstates_space_address)[288];
  int32_T c1__index_space;
  boolean_T c1_dsmdiag_states_space;
  real_T (*c1__indextarget_address)[864];
  int32_T c1_j__index;
  boolean_T c1_dsmdiag_target;
  void *c1_fEmlrtCtx;
  real_T *c1_c_y;
} SFc1_RM3_copyInstanceStruct;

#endif                                 /*typedef_SFc1_RM3_copyInstanceStruct*/

/* Named Constants */

/* Variable Declarations */
extern struct SfDebugInstanceStruct *sfGlobalDebugInstanceStruct;

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c1_RM3_copy_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c1_RM3_copy_get_check_sum(mxArray *plhs[]);
extern void c1_RM3_copy_method_dispatcher(SimStruct *S, int_T method, void *data);

#endif
