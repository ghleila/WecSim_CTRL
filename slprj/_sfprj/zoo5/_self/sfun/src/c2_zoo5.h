#ifndef __c2_zoo5_h__
#define __c2_zoo5_h__

/* Type Definitions */
#ifndef typedef_SFc2_zoo5InstanceStruct
#define typedef_SFc2_zoo5InstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c2_sfEvent;
  boolean_T c2_doneDoubleBufferReInit;
  uint8_T c2_is_active_c2_zoo5;
  void *c2_fEmlrtCtx;
  real_T *c2_u;
  real_T (*c2_y)[2];
} SFc2_zoo5InstanceStruct;

#endif                                 /*typedef_SFc2_zoo5InstanceStruct*/

/* Named Constants */

/* Variable Declarations */
extern struct SfDebugInstanceStruct *sfGlobalDebugInstanceStruct;

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c2_zoo5_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c2_zoo5_get_check_sum(mxArray *plhs[]);
extern void c2_zoo5_method_dispatcher(SimStruct *S, int_T method, void *data);

#endif
