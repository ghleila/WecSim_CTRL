/* Include files */

#include "RM3_copy_sfun.h"
#include "c1_RM3_copy.h"
#include <string.h>
#include "mwmathutil.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "RM3_copy_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c_with_debugger(S, sfGlobalDebugInstanceStruct);

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization);
static void chart_debug_initialize_data_addresses(SimStruct *S);
static const mxArray* sf_opaque_get_hover_data_for_msg(void *chartInstance,
  int32_T msgSSID);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static real_T _sfTime_;
static const char * c1_debug_family_names[37] = { "current_action1", "PptoAvg",
  "gamma", "opts", "alpha", "current_action", "action_idx", "action_idx_pr",
  "Height", "Ts", "next_stPrimee_idx", "st_chosen", "Height_prime", "Ts_prime",
  "next_state_idx", "next_reward", "Period", "R_idx", "MaxR", "Final_Rew", "In",
  "Output", "nodes", "Out", "m", "range", "Input", "m11", "net", "tr", "y_prime",
  "c", "Q_prev", "k", "nargin", "nargout", "y" };

static emlrtRTEInfo c1_emlrtRTEI = { 44,/* lineNo */
  10,                                  /* colNo */
  "find",                              /* fName */
  "/Applications/MATLAB_R2017b.app/toolbox/eml/lib/matlab/elmat/find.m"/* pName */
};

static emlrtRTEInfo c1_b_emlrtRTEI = { 44,/* lineNo */
  5,                                   /* colNo */
  "find",                              /* fName */
  "/Applications/MATLAB_R2017b.app/toolbox/eml/lib/matlab/elmat/find.m"/* pName */
};

static emlrtRTEInfo c1_c_emlrtRTEI = { 39,/* lineNo */
  9,                                   /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_d_emlrtRTEI = { 56,/* lineNo */
  9,                                   /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_e_emlrtRTEI = { 67,/* lineNo */
  9,                                   /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_f_emlrtRTEI = { 70,/* lineNo */
  15,                                  /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_g_emlrtRTEI = { 70,/* lineNo */
  34,                                  /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_h_emlrtRTEI = { 74,/* lineNo */
  12,                                  /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_i_emlrtRTEI = { 74,/* lineNo */
  30,                                  /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_j_emlrtRTEI = { 75,/* lineNo */
  9,                                   /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_k_emlrtRTEI = { 80,/* lineNo */
  9,                                   /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_l_emlrtRTEI = { 83,/* lineNo */
  9,                                   /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_m_emlrtRTEI = { 94,/* lineNo */
  9,                                   /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_n_emlrtRTEI = { 96,/* lineNo */
  20,                                  /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_o_emlrtRTEI = { 96,/* lineNo */
  9,                                   /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_p_emlrtRTEI = { 124,/* lineNo */
  17,                                  /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_q_emlrtRTEI = { 116,/* lineNo */
  17,                                  /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_r_emlrtRTEI = { 126,/* lineNo */
  17,                                  /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_s_emlrtRTEI = { 118,/* lineNo */
  17,                                  /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_t_emlrtRTEI = { 140,/* lineNo */
  17,                                  /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_u_emlrtRTEI = { 141,/* lineNo */
  25,                                  /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_v_emlrtRTEI = { 19,/* lineNo */
  16,                                  /* colNo */
  "minOrMax",                          /* fName */
  "/Applications/MATLAB_R2017b.app/toolbox/eml/eml/+coder/+internal/minOrMax.m"/* pName */
};

static emlrtRTEInfo c1_w_emlrtRTEI = { 142,/* lineNo */
  29,                                  /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_x_emlrtRTEI = { 142,/* lineNo */
  17,                                  /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_y_emlrtRTEI = { 143,/* lineNo */
  23,                                  /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_ab_emlrtRTEI = { 143,/* lineNo */
  35,                                  /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_bb_emlrtRTEI = { 143,/* lineNo */
  17,                                  /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_cb_emlrtRTEI = { 146,/* lineNo */
  17,                                  /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_db_emlrtRTEI = { 147,/* lineNo */
  27,                                  /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_eb_emlrtRTEI = { 147,/* lineNo */
  17,                                  /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_fb_emlrtRTEI = { 148,/* lineNo */
  29,                                  /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_gb_emlrtRTEI = { 148,/* lineNo */
  25,                                  /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_hb_emlrtRTEI = { 148,/* lineNo */
  17,                                  /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_ib_emlrtRTEI = { 149,/* lineNo */
  26,                                  /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_jb_emlrtRTEI = { 149,/* lineNo */
  42,                                  /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_kb_emlrtRTEI = { 149,/* lineNo */
  17,                                  /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_lb_emlrtRTEI = { 153,/* lineNo */
  34,                                  /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_mb_emlrtRTEI = { 153,/* lineNo */
  41,                                  /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_nb_emlrtRTEI = { 156,/* lineNo */
  26,                                  /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_ob_emlrtRTEI = { 157,/* lineNo */
  21,                                  /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_pb_emlrtRTEI = { 139,/* lineNo */
  17,                                  /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_qb_emlrtRTEI = { 145,/* lineNo */
  17,                                  /* colNo */
  "MATLAB Function1",                  /* fName */
  "#RM3_copy:893"                      /* pName */
};

static emlrtRTEInfo c1_rb_emlrtRTEI = { 44,/* lineNo */
  32,                                  /* colNo */
  "mpower",                            /* fName */
  "/Applications/MATLAB_R2017b.app/toolbox/eml/lib/matlab/ops/mpower.m"/* pName */
};

static emlrtRTEInfo c1_sb_emlrtRTEI = { 205,/* lineNo */
  42,                                  /* colNo */
  "matrix_to_integer_power",           /* fName */
  "/Applications/MATLAB_R2017b.app/toolbox/eml/lib/matlab/ops/private/matrix_to_integer_power.m"/* pName */
};

static emlrtRTEInfo c1_tb_emlrtRTEI = { 210,/* lineNo */
  5,                                   /* colNo */
  "matrix_to_integer_power",           /* fName */
  "/Applications/MATLAB_R2017b.app/toolbox/eml/lib/matlab/ops/private/matrix_to_integer_power.m"/* pName */
};

static emlrtRTEInfo c1_ub_emlrtRTEI = { 73,/* lineNo */
  5,                                   /* colNo */
  "mpower",                            /* fName */
  "/Applications/MATLAB_R2017b.app/toolbox/eml/lib/matlab/ops/mpower.m"/* pName */
};

static emlrtRTEInfo c1_vb_emlrtRTEI = { 26,/* lineNo */
  1,                                   /* colNo */
  "rdivide",                           /* fName */
  "/Applications/MATLAB_R2017b.app/toolbox/eml/lib/matlab/ops/rdivide.m"/* pName */
};

static const char_T c1_cv0[15] = { 'M', 'A', 'T', 'L', 'A', 'B', ':', 'd', 'i',
  'm', 'a', 'g', 'r', 'e', 'e' };

/* Function Declarations */
static void initialize_c1_RM3_copy(SFc1_RM3_copyInstanceStruct *chartInstance);
static void initialize_params_c1_RM3_copy(SFc1_RM3_copyInstanceStruct
  *chartInstance);
static void enable_c1_RM3_copy(SFc1_RM3_copyInstanceStruct *chartInstance);
static void disable_c1_RM3_copy(SFc1_RM3_copyInstanceStruct *chartInstance);
static void c1_update_debugger_state_c1_RM3_copy(SFc1_RM3_copyInstanceStruct
  *chartInstance);
static const mxArray *get_sim_state_c1_RM3_copy(SFc1_RM3_copyInstanceStruct
  *chartInstance);
static void set_sim_state_c1_RM3_copy(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_st);
static void finalize_c1_RM3_copy(SFc1_RM3_copyInstanceStruct *chartInstance);
static void sf_gateway_c1_RM3_copy(SFc1_RM3_copyInstanceStruct *chartInstance);
static void mdl_start_c1_RM3_copy(SFc1_RM3_copyInstanceStruct *chartInstance);
static void c1_chartstep_c1_RM3_copy(SFc1_RM3_copyInstanceStruct *chartInstance);
static void initSimStructsc1_RM3_copy(SFc1_RM3_copyInstanceStruct *chartInstance);
static void init_script_number_translation(uint32_T c1_machineNumber, uint32_T
  c1_chartNumber, uint32_T c1_instanceNumber);
static const mxArray *c1_sf_marshallOut(void *chartInstanceVoid, void *c1_inData);
static void c1_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_b_sf_marshallOut(void *chartInstanceVoid, real_T
  c1_inData_data[], int32_T c1_inData_size[2]);
static void c1_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T
  c1_y_data[], int32_T c1_y_size[2]);
static void c1_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, real_T c1_outData_data[], int32_T
  c1_outData_size[2]);
static const mxArray *c1_c_sf_marshallOut(void *chartInstanceVoid,
  c1_emxArray_real_T *c1_inData);
static void c1_b_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId,
  c1_emxArray_real_T *c1_d_y);
static void c1_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, c1_emxArray_real_T *c1_outData);
static const mxArray *c1_d_sf_marshallOut(void *chartInstanceVoid, real_T
  c1_inData_data[], int32_T c1_inData_size[1]);
static void c1_c_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T
  c1_y_data[], int32_T c1_y_size[1]);
static void c1_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, real_T c1_outData_data[], int32_T
  c1_outData_size[1]);
static const mxArray *c1_e_sf_marshallOut(void *chartInstanceVoid,
  c1_emxArray_real_T *c1_inData);
static void c1_d_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId,
  c1_emxArray_real_T *c1_d_y);
static void c1_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, c1_emxArray_real_T *c1_outData);
static const mxArray *c1_f_sf_marshallOut(void *chartInstanceVoid, real_T
  c1_inData_data[], int32_T c1_inData_size[1]);
static void c1_e_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T
  c1_y_data[], int32_T c1_y_size[1]);
static void c1_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, real_T c1_outData_data[], int32_T
  c1_outData_size[1]);
static const mxArray *c1_g_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_f_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_d_y
  [864]);
static void c1_g_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_h_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static const mxArray *c1_i_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_g_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_h_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_j_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_i_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_k_sf_marshallOut(void *chartInstanceVoid, real_T
  c1_inData_data[], int32_T c1_inData_size[1]);
static void c1_h_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T
  c1_y_data[], int32_T c1_y_size[1]);
static void c1_j_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, real_T c1_outData_data[], int32_T
  c1_outData_size[1]);
static const mxArray *c1_l_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static c1_sI3an6DNyfQAOCkXY61B1nE c1_i_emlrt_marshallIn
  (SFc1_RM3_copyInstanceStruct *chartInstance, const mxArray *c1_b_u, const
   emlrtMsgIdentifier *c1_parentId);
static void c1_k_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static real_T c1_sum(SFc1_RM3_copyInstanceStruct *chartInstance, real_T c1_x[50]);
static int16_T c1_rdivide(SFc1_RM3_copyInstanceStruct *chartInstance, real_T
  c1_x, int16_T c1_d_y);
static void c1_mpower(SFc1_RM3_copyInstanceStruct *chartInstance, real_T
                      c1_a_data[], int32_T c1_a_size[1], real_T c1_c_data[],
                      int32_T c1_c_size[1]);
static void c1_matrix_to_scalar_power(SFc1_RM3_copyInstanceStruct *chartInstance,
  real_T c1_a_data[], int32_T c1_a_size[1], real_T c1_c_data[], int32_T
  c1_c_size[1]);
static void c1_b_rdivide(SFc1_RM3_copyInstanceStruct *chartInstance, real_T
  c1_x_data[], int32_T c1_x_size[1], real_T c1_y_data[], int32_T c1_y_size[1],
  real_T c1_z_data[], int32_T c1_z_size[1]);
static void c1_c_rdivide(SFc1_RM3_copyInstanceStruct *chartInstance, real_T
  c1_x_data[], int32_T c1_x_size[2], real_T c1_y_data[], int32_T c1_y_size[2],
  real_T c1_z_data[], int32_T c1_z_size[2]);
static const mxArray *c1_emlrt_marshallOut(SFc1_RM3_copyInstanceStruct
  *chartInstance, const real_T c1_b_u[24]);
static const mxArray *c1_b_emlrt_marshallOut(SFc1_RM3_copyInstanceStruct
  *chartInstance, const real_T c1_b_u);
static void c1_j_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_datasample, const char_T *c1_identifier, real_T c1_d_y[2]);
static void c1_k_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_d_y[2]);
static const mxArray *c1_c_emlrt_marshallOut(SFc1_RM3_copyInstanceStruct
  *chartInstance, const real_T c1_b_u[72]);
static real_T c1_l_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_datasample, const char_T *c1_identifier);
static real_T c1_m_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId);
static const mxArray *c1_d_emlrt_marshallOut(SFc1_RM3_copyInstanceStruct
  *chartInstance, const char_T c1_b_u[16]);
static const mxArray *c1_e_emlrt_marshallOut(SFc1_RM3_copyInstanceStruct
  *chartInstance, const int16_T c1_b_u);
static const mxArray *c1_f_emlrt_marshallOut(SFc1_RM3_copyInstanceStruct
  *chartInstance, const char_T c1_b_u[14]);
static const mxArray *c1_g_emlrt_marshallOut(SFc1_RM3_copyInstanceStruct
  *chartInstance, const real_T c1_u_data[], const int32_T c1_u_size[1]);
static const mxArray *c1_h_emlrt_marshallOut(SFc1_RM3_copyInstanceStruct
  *chartInstance, const real_T c1_b_u[2]);
static const mxArray *c1_i_emlrt_marshallOut(SFc1_RM3_copyInstanceStruct
  *chartInstance, const char_T c1_b_u[15]);
static const mxArray *c1_j_emlrt_marshallOut(SFc1_RM3_copyInstanceStruct
  *chartInstance, const real_T c1_u_data[], const int32_T c1_u_size[1]);
static const mxArray *c1_k_emlrt_marshallOut(SFc1_RM3_copyInstanceStruct
  *chartInstance, const real_T c1_u_data[], const int32_T c1_u_size[2]);
static const mxArray *c1_l_emlrt_marshallOut(SFc1_RM3_copyInstanceStruct
  *chartInstance, const real_T c1_u_data[], const int32_T c1_u_size[2]);
static void c1_n_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_transpose, const char_T *c1_identifier, real_T c1_y_data[],
  int32_T c1_y_size[2]);
static void c1_o_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T
  c1_y_data[], int32_T c1_y_size[2]);
static const mxArray *c1_m_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static int32_T c1_p_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_l_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_n_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_q_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_d_y[24]);
static void c1_m_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_o_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_r_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_d_y
  [15000]);
static void c1_n_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_p_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static int16_T c1_s_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_o_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_q_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_t_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_d_y
  [288]);
static void c1_p_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_r_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_u_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, int16_T c1_d_y
  [10000]);
static void c1_q_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_s_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_v_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_d_y
  [3600]);
static void c1_r_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_t_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_w_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_d_y
  [10000]);
static void c1_s_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_u_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_x_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_d_y
  [1000]);
static void c1_t_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_v_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_y_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_d_y
  [150000]);
static void c1_u_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static uint8_T c1_ab_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_is_active_c1_RM3_copy, const char_T *c1_identifier);
static uint8_T c1_bb_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_emxEnsureCapacity_real_T(SFc1_RM3_copyInstanceStruct
  *chartInstance, c1_emxArray_real_T *c1_emxArray, int32_T c1_oldNumel, const
  emlrtRTEInfo *c1_srcLocation);
static void c1_emxInit_real_T(SFc1_RM3_copyInstanceStruct *chartInstance,
  c1_emxArray_real_T **c1_pEmxArray, int32_T c1_numDimensions, const
  emlrtRTEInfo *c1_srcLocation);
static void c1_emxFree_real_T(SFc1_RM3_copyInstanceStruct *chartInstance,
  c1_emxArray_real_T **c1_pEmxArray);
static int32_T c1__s32_d_(SFc1_RM3_copyInstanceStruct *chartInstance, real_T
  c1_b, uint32_T c1_ssid_src_loc, int32_T c1_offset_src_loc, int32_T
  c1_length_src_loc);
static real_T c1_get_B2(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex);
static void c1_set_B2(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
                      c1_elementIndex, real_T c1_elementValue);
static real_T *c1_access_B2(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_rdOnly);
static real_T c1_get_Final_Rew_matt(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_elementIndex);
static void c1_set_Final_Rew_matt(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_elementIndex, real_T c1_elementValue);
static real_T *c1_access_Final_Rew_matt(SFc1_RM3_copyInstanceStruct
  *chartInstance, uint32_T c1_rdOnly);
static real_T c1_get_K2(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex);
static void c1_set_K2(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
                      c1_elementIndex, real_T c1_elementValue);
static real_T *c1_access_K2(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_rdOnly);
static real_T c1_get_MeanR(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex);
static void c1_set_MeanR(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex, real_T c1_elementValue);
static real_T *c1_access_MeanR(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_rdOnly);
static real_T c1_get_Memory_D_size(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_elementIndex);
static void c1_set_Memory_D_size(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_elementIndex, real_T c1_elementValue);
static real_T *c1_access_Memory_D_size(SFc1_RM3_copyInstanceStruct
  *chartInstance, uint32_T c1_rdOnly);
static real_T c1_get_Memory_Stor(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_elementIndex);
static void c1_set_Memory_Stor(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_elementIndex, real_T c1_elementValue);
static real_T *c1_access_Memory_Stor(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_rdOnly);
static real_T c1_get_NN_dd(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex);
static void c1_set_NN_dd(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex, real_T c1_elementValue);
static real_T *c1_access_NN_dd(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_rdOnly);
static real_T c1_get_Q(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex);
static void c1_set_Q(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
                     c1_elementIndex, real_T c1_elementValue);
static real_T *c1_access_Q(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_rdOnly);
static real_T c1_get_Q_prime(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_elementIndex);
static void c1_set_Q_prime(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex, real_T c1_elementValue);
static real_T *c1_access_Q_prime(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_rdOnly);
static real_T c1_get_Q_targ(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex);
static void c1_set_Q_targ(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex, real_T c1_elementValue);
static real_T *c1_access_Q_targ(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_rdOnly);
static real_T c1_get_R(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex);
static void c1_set_R(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
                     c1_elementIndex, real_T c1_elementValue);
static real_T *c1_access_R(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_rdOnly);
static real_T c1_get_Stor_eta(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_elementIndex);
static void c1_set_Stor_eta(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex, real_T c1_elementValue);
static real_T *c1_access_Stor_eta(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_rdOnly);
static real_T c1_get_action(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex);
static void c1_set_action(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex, real_T c1_elementValue);
static real_T *c1_access_action(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_rdOnly);
static int16_T c1_get_num_R(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex);
static void c1_set_num_R(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex, int16_T c1_elementValue);
static int16_T *c1_access_num_R(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_rdOnly);
static real_T c1_get_ove_count(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_elementIndex);
static void c1_set_ove_count(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_elementIndex, real_T c1_elementValue);
static real_T *c1_access_ove_count(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_rdOnly);
static int16_T c1_get_state_idx(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_elementIndex);
static void c1_set_state_idx(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_elementIndex, int16_T c1_elementValue);
static int16_T *c1_access_state_idx(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_rdOnly);
static real_T c1_get_states_space(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_elementIndex);
static void c1_set_states_space(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_elementIndex, real_T c1_elementValue);
static real_T *c1_access_states_space(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_rdOnly);
static real_T c1_get_target(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex);
static void c1_set_target(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex, real_T c1_elementValue);
static real_T *c1_access_target(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_rdOnly);
static void init_dsm_address_info(SFc1_RM3_copyInstanceStruct *chartInstance);
static void init_simulink_io_address(SFc1_RM3_copyInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c1_RM3_copy(SFc1_RM3_copyInstanceStruct *chartInstance)
{
  if (sf_is_first_init_cond(chartInstance->S)) {
    initSimStructsc1_RM3_copy(chartInstance);
    chart_debug_initialize_data_addresses(chartInstance->S);
  }

  chartInstance->c1_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c1_is_active_c1_RM3_copy = 0U;
}

static void initialize_params_c1_RM3_copy(SFc1_RM3_copyInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void enable_c1_RM3_copy(SFc1_RM3_copyInstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c1_RM3_copy(SFc1_RM3_copyInstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c1_update_debugger_state_c1_RM3_copy(SFc1_RM3_copyInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c1_RM3_copy(SFc1_RM3_copyInstanceStruct
  *chartInstance)
{
  const mxArray *c1_st;
  const mxArray *c1_d_y = NULL;
  real_T c1_hoistedGlobal;
  const mxArray *c1_e_y = NULL;
  uint8_T c1_b_hoistedGlobal;
  const mxArray *c1_f_y = NULL;
  c1_st = NULL;
  c1_st = NULL;
  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_createcellmatrix(2, 1), false);
  c1_hoistedGlobal = *chartInstance->c1_c_y;
  c1_e_y = NULL;
  sf_mex_assign(&c1_e_y, sf_mex_create("y", &c1_hoistedGlobal, 0, 0U, 0U, 0U, 0),
                false);
  sf_mex_setcell(c1_d_y, 0, c1_e_y);
  c1_b_hoistedGlobal = chartInstance->c1_is_active_c1_RM3_copy;
  c1_f_y = NULL;
  sf_mex_assign(&c1_f_y, sf_mex_create("y", &c1_b_hoistedGlobal, 3, 0U, 0U, 0U,
    0), false);
  sf_mex_setcell(c1_d_y, 1, c1_f_y);
  sf_mex_assign(&c1_st, c1_d_y, false);
  return c1_st;
}

static void set_sim_state_c1_RM3_copy(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_st)
{
  const mxArray *c1_b_u;
  chartInstance->c1_doneDoubleBufferReInit = true;
  c1_b_u = sf_mex_dup(c1_st);
  *chartInstance->c1_c_y = c1_l_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getcell(c1_b_u, 0)), "y");
  chartInstance->c1_is_active_c1_RM3_copy = c1_ab_emlrt_marshallIn(chartInstance,
    sf_mex_dup(sf_mex_getcell(c1_b_u, 1)), "is_active_c1_RM3_copy");
  sf_mex_destroy(&c1_b_u);
  c1_update_debugger_state_c1_RM3_copy(chartInstance);
  sf_mex_destroy(&c1_st);
}

static void finalize_c1_RM3_copy(SFc1_RM3_copyInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c1_RM3_copy(SFc1_RM3_copyInstanceStruct *chartInstance)
{
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 0U, chartInstance->c1_sfEvent);
  chartInstance->c1_sfEvent = CALL_EVENT;
  c1_chartstep_c1_RM3_copy(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_RM3_copyMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  _SFD_DATA_RANGE_CHECK(*chartInstance->c1_c_y, 0U);
}

static void mdl_start_c1_RM3_copy(SFc1_RM3_copyInstanceStruct *chartInstance)
{
  sim_mode_is_external(chartInstance->S);
}

static void c1_chartstep_c1_RM3_copy(SFc1_RM3_copyInstanceStruct *chartInstance)
{
  c1_emxArray_real_T *c1_In;
  c1_emxArray_real_T *c1_Output;
  c1_emxArray_real_T *c1_Out;
  c1_emxArray_real_T *c1_Input;
  boolean_T c1_b0;
  uint32_T c1_debug_family_var_map[37];
  real_T c1_current_action1[2];
  real_T c1_PptoAvg;
  real_T c1_gamma;
  c1_sI3an6DNyfQAOCkXY61B1nE c1_opts;
  real_T c1_alpha;
  real_T c1_current_action[2];
  real_T c1_action_idx;
  real_T c1_action_idx_pr_data[12];
  int32_T c1_action_idx_pr_size[1];
  real_T c1_Height;
  real_T c1_Ts;
  real_T c1_next_stPrimee_idx_data[72];
  int32_T c1_next_stPrimee_idx_size[1];
  real_T c1_st_chosen;
  real_T c1_Height_prime;
  real_T c1_Ts_prime;
  real_T c1_next_state_idx_data[72];
  int32_T c1_next_state_idx_size[1];
  real_T c1_next_reward;
  real_T c1_Period_data[72];
  int32_T c1_Period_size[1];
  real_T c1_R_idx_data[72];
  int32_T c1_R_idx_size[1];
  real_T c1_MaxR;
  real_T c1_Final_Rew_data[72];
  int32_T c1_Final_Rew_size[1];
  real_T c1_nodes[2];
  int32_T c1_m_size[1];
  int32_T c1_range_size[1];
  real_T c1_m11_data[4];
  int32_T c1_m11_size[2];
  const mxArray *c1_net = NULL;
  const mxArray *c1_tr = NULL;
  const mxArray *c1_y_prime = NULL;
  real_T c1_c;
  real_T c1_Q_prev[864];
  real_T c1_k;
  real_T c1_Height_data[72];
  int32_T c1_Height_size[1];
  real_T c1_Period;
  real_T c1_R_idx;
  int32_T c1_Out_size[1];
  real_T c1_b_range_data[4];
  int32_T c1_b_range_size[2];
  real_T c1_Input_data[4];
  int32_T c1_Input_size[2];
  real_T c1_next_state_idx;
  real_T c1_Final_Rew;
  real_T c1_nargin = 0.0;
  real_T c1_nargout = 1.0;
  real_T c1_d_y;
  int32_T c1_i0;
  int32_T c1_i1;
  int32_T c1_i2;
  real_T c1_dv0[24];
  real_T c1_dv1[2];
  int32_T c1_i3;
  real_T c1_b_current_action;
  int32_T c1_i4;
  real_T c1_c_current_action;
  boolean_T c1_bv0[12];
  int32_T c1_i5;
  int32_T c1_idx;
  boolean_T c1_bv1[12];
  int32_T c1_ii_size[1];
  int32_T c1_ii;
  int32_T c1_i6;
  int32_T c1_ii_data[12];
  int32_T c1_loop_ub;
  int32_T c1_i7;
  int32_T c1_action_idx_pr[2];
  int32_T c1_i8;
  int32_T c1_i9;
  boolean_T c1_x[72];
  int32_T c1_i10;
  boolean_T c1_bv2[72];
  int32_T c1_b_idx;
  int32_T c1_b_ii_size[1];
  int32_T c1_b_ii;
  int32_T c1_i11;
  int32_T c1_b_ii_data[72];
  int32_T c1_b_loop_ub;
  int32_T c1_i12;
  int32_T c1_i13;
  real_T c1_dv2[72];
  int32_T c1_i14;
  int32_T c1_i15;
  int32_T c1_i16;
  int32_T c1_c_idx;
  int32_T c1_c_ii;
  int32_T c1_i17;
  int32_T c1_c_loop_ub;
  int32_T c1_i18;
  int32_T c1_d_loop_ub;
  int32_T c1_i19;
  int32_T c1_tmp_size[1];
  int32_T c1_e_loop_ub;
  int32_T c1_i20;
  int16_T c1_b_tmp_data[72];
  int32_T c1_f_loop_ub;
  int32_T c1_i21;
  int32_T c1_i22;
  real_T c1_a;
  int32_T c1_g_loop_ub;
  int32_T c1_i23;
  int32_T c1_b_tmp_size[1];
  int32_T c1_h_loop_ub;
  int32_T c1_i24;
  int32_T c1_d_ii[2];
  int32_T c1_c_tmp_data[72];
  int32_T c1_i_loop_ub;
  int32_T c1_i25;
  int32_T c1_j_loop_ub;
  int32_T c1_i26;
  int32_T c1_k_loop_ub;
  int32_T c1_i27;
  int32_T c1_b_Height[2];
  int32_T c1_l_loop_ub;
  int32_T c1_i28;
  int32_T c1_b_Period[2];
  int32_T c1_i29;
  int32_T c1_i30;
  int32_T c1_d_idx;
  int32_T c1_e_ii;
  int32_T c1_i31;
  int32_T c1_m_loop_ub;
  int32_T c1_i32;
  int32_T c1_b_R_idx[2];
  int32_T c1_c_R_idx;
  int32_T c1_i33;
  real_T c1_dv3[50];
  real_T c1_varargin_1;
  real_T c1_A;
  real_T c1_B;
  int32_T c1_n_loop_ub;
  int32_T c1_i34;
  int32_T c1_b_Final_Rew_size[1];
  int32_T c1_o_loop_ub;
  int32_T c1_i35;
  real_T c1_b_Final_Rew_data[72];
  real_T c1_d_tmp_data[72];
  int32_T c1_c_tmp_size[1];
  int32_T c1_p_loop_ub;
  int32_T c1_i36;
  int32_T c1_b_action_idx;
  int32_T c1_i37;
  static char_T c1_cv1[16] = { 'c', 'u', 'r', 'r', 'e', 'n', 't', ' ', 's', 't',
    'a', 't', 'e', ' ', ':', ' ' };

  static char_T c1_cv2[14] = { ' ', 'n', 'e', 'x', 't', ' ', 's', 't', 'a', 't',
    'e', ' ', ':', ' ' };

  static char_T c1_cv3[16] = { ' ', 't', 'a', 'k', 'e', 'n', ' ', 'a', 'c', 't',
    'i', 'o', 'n', ' ', ':', ' ' };

  real_T c1_dv4[2];
  static char_T c1_cv4[15] = { ' ', 'n', 'e', 'x', 't', ' ', 'r', 'e', 'w', 'a',
    'r', 'd', ' ', ':', ' ' };

  int32_T c1_i38;
  int32_T c1_i39;
  int32_T c1_i40;
  int32_T c1_i41;
  int32_T c1_i42;
  int32_T c1_i43;
  int32_T c1_i44;
  int32_T c1_i45;
  int32_T c1_i46;
  int32_T c1_q_loop_ub;
  int32_T c1_i47;
  int32_T c1_i48;
  int32_T c1_i49;
  int32_T c1_i50;
  int32_T c1_i51;
  int32_T c1_i52;
  int32_T c1_i53;
  const mxArray *c1_b_net = NULL;
  const mxArray *c1_b_tr = NULL;
  int32_T c1_i54;
  int32_T c1_i55;
  int32_T c1_r_loop_ub;
  int32_T c1_i56;
  int32_T c1_i57;
  int32_T c1_b_Out;
  int32_T c1_c_Out;
  int32_T c1_s_loop_ub;
  int32_T c1_i58;
  c1_emxArray_real_T *c1_b_varargin_1;
  int32_T c1_i59;
  int32_T c1_c_varargin_1;
  int32_T c1_d_varargin_1;
  int32_T c1_t_loop_ub;
  int32_T c1_i60;
  int32_T c1_vstride;
  int32_T c1_j;
  int32_T c1_ixstart;
  int32_T c1_i61;
  int32_T c1_ixstop;
  real_T c1_mtmp;
  int32_T c1_e_varargin_1;
  const mxArray *c1_e_y = NULL;
  int32_T c1_f_varargin_1;
  static char_T c1_cv5[27] = { 'C', 'o', 'd', 'e', 'r', ':', 'b', 'u', 'i', 'l',
    't', 'i', 'n', 's', ':', 'V', 'e', 'c', 't', 'o', 'r', 'S', 't', 'r', 'i',
    'd', 'e' };

  int32_T c1_ix;
  int32_T c1_u_loop_ub;
  const mxArray *c1_f_y = NULL;
  int32_T c1_i62;
  static char_T c1_cv6[27] = { 'C', 'o', 'd', 'e', 'r', ':', 'b', 'u', 'i', 'l',
    't', 'i', 'n', 's', ':', 'V', 'e', 'c', 't', 'o', 'r', 'S', 't', 'r', 'i',
    'd', 'e' };

  int32_T c1_b_ix;
  int32_T c1_d_tmp_size[1];
  int32_T c1_b_vstride;
  int32_T c1_b_j;
  int32_T c1_b_ixstart;
  int32_T c1_b_ixstop;
  real_T c1_b_mtmp;
  int32_T c1_v_loop_ub;
  int32_T c1_i63;
  const mxArray *c1_g_y = NULL;
  static char_T c1_cv7[27] = { 'C', 'o', 'd', 'e', 'r', ':', 'b', 'u', 'i', 'l',
    't', 'i', 'n', 's', ':', 'V', 'e', 'c', 't', 'o', 'r', 'S', 't', 'r', 'i',
    'd', 'e' };

  int32_T c1_c_ix;
  const mxArray *c1_h_y = NULL;
  static char_T c1_cv8[27] = { 'C', 'o', 'd', 'e', 'r', ':', 'b', 'u', 'i', 'l',
    't', 'i', 'n', 's', ':', 'V', 'e', 'c', 't', 'o', 'r', 'S', 't', 'r', 'i',
    'd', 'e' };

  int32_T c1_d_ix;
  int32_T c1_w_loop_ub;
  int32_T c1_i64;
  int32_T c1_e_tmp_size[1];
  int32_T c1_x_loop_ub;
  int32_T c1_i65;
  int32_T c1_c_range_size[1];
  real_T c1_e_tmp_data[5000];
  int32_T c1_y_loop_ub;
  int32_T c1_i66;
  real_T c1_c_range_data[5000];
  real_T c1_f_tmp_data[5000];
  int32_T c1_f_tmp_size[1];
  int32_T c1_ab_loop_ub;
  int32_T c1_i67;
  int32_T c1_i68;
  int32_T c1_b_Input;
  int32_T c1_c_Input;
  int32_T c1_bb_loop_ub;
  int32_T c1_i69;
  c1_emxArray_real_T *c1_g_varargin_1;
  int32_T c1_i70;
  int32_T c1_h_varargin_1;
  int32_T c1_i_varargin_1;
  int32_T c1_cb_loop_ub;
  int32_T c1_i71;
  boolean_T c1_b1;
  const mxArray *c1_i_y = NULL;
  static char_T c1_cv9[36] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'a', 'u', 't', 'o', 'D', 'i', 'm', 'I', 'n', 'c', 'o',
    'm', 'p', 'a', 't', 'i', 'b', 'i', 'l', 'i', 't', 'y' };

  const mxArray *c1_j_y = NULL;
  static char_T c1_cv10[39] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'e', 'm', 'l', '_', 'm', 'i', 'n', '_', 'o', 'r', '_',
    'm', 'a', 'x', '_', 'v', 'a', 'r', 'D', 'i', 'm', 'Z', 'e', 'r', 'o' };

  int32_T c1_n;
  int32_T c1_i;
  int32_T c1_e_ix;
  int32_T c1_i72;
  int32_T c1_c_ixstart;
  int32_T c1_c_ixstop;
  real_T c1_c_mtmp;
  int32_T c1_j_varargin_1;
  int32_T c1_k_varargin_1;
  int32_T c1_f_ix;
  int32_T c1_db_loop_ub;
  int32_T c1_i73;
  int32_T c1_g_ix;
  boolean_T c1_b2;
  const mxArray *c1_k_y = NULL;
  static char_T c1_cv11[36] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'a', 'u', 't', 'o', 'D', 'i', 'm', 'I', 'n', 'c', 'o',
    'm', 'p', 'a', 't', 'i', 'b', 'i', 'l', 'i', 't', 'y' };

  const mxArray *c1_l_y = NULL;
  static char_T c1_cv12[39] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'e', 'm', 'l', '_', 'm', 'i', 'n', '_', 'o', 'r', '_',
    'm', 'a', 'x', '_', 'v', 'a', 'r', 'D', 'i', 'm', 'Z', 'e', 'r', 'o' };

  int32_T c1_g_tmp_size[2];
  int32_T c1_b_n;
  int32_T c1_b_i;
  int32_T c1_h_ix;
  int32_T c1_d_ixstart;
  int32_T c1_d_ixstop;
  int32_T c1_range;
  real_T c1_d_mtmp;
  int32_T c1_b_range;
  int32_T c1_i74;
  real_T c1_g_tmp_data[4];
  int32_T c1_i_ix;
  int32_T c1_j_ix;
  int32_T c1_i75;
  int32_T c1_i76;
  int32_T c1_d_Input[2];
  int32_T c1_m11[2];
  int32_T c1_b_Input_size[2];
  int32_T c1_e_Input;
  int32_T c1_f_Input;
  int32_T c1_eb_loop_ub;
  int32_T c1_i77;
  int32_T c1_d_range_size[2];
  real_T c1_b_Input_data[4];
  int32_T c1_c_range;
  int32_T c1_d_range;
  int32_T c1_i78;
  real_T c1_d_range_data[4];
  real_T c1_h_tmp_data[4];
  int32_T c1_h_tmp_size[2];
  int32_T c1_g_Input;
  int32_T c1_h_Input;
  int32_T c1_i79;
  int32_T c1_c_Input_size[1];
  int32_T c1_i80;
  int32_T c1_b_Out_size[2];
  real_T c1_c_Input_data[4];
  int32_T c1_fb_loop_ub;
  int32_T c1_i81;
  real_T c1_b_Out_data[5000];
  int32_T c1_i_tmp_size[2];
  int32_T c1_gb_loop_ub;
  int32_T c1_i82;
  real_T c1_i_tmp_data[288];
  int32_T c1_i83;
  int32_T c1_hb_loop_ub;
  int32_T c1_i84;
  real_T c1_j_tmp_data[864];
  int32_T c1_j_tmp_size[2];
  int32_T c1_i85;
  int32_T c1_ib_loop_ub;
  int32_T c1_i86;
  int32_T c1_b_next_state_idx[2];
  int32_T c1_b_Final_Rew[2];
  int32_T c1_c_next_state_idx;
  int32_T c1_i87;
  int32_T c1_e_ixstart;
  real_T c1_l_varargin_1[12];
  real_T c1_e_mtmp;
  int32_T c1_k_ix;
  int32_T c1_l_ix;
  int32_T c1_i88;
  int32_T c1_i89;
  boolean_T exitg1;
  c1_emxInit_real_T(chartInstance, &c1_In, 2, &c1_q_emlrtRTEI);
  c1_emxInit_real_T(chartInstance, &c1_Output, 2, &c1_s_emlrtRTEI);
  c1_emxInit_real_T(chartInstance, &c1_Out, 2, &c1_pb_emlrtRTEI);
  c1_emxInit_real_T(chartInstance, &c1_Input, 2, &c1_qb_emlrtRTEI);
  c1_b0 = false;
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 0U, chartInstance->c1_sfEvent);
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 37U, 47U, c1_debug_family_names,
    c1_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_current_action1, 0U,
    c1_j_sf_marshallOut, c1_i_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c1_PptoAvg, 1U, c1_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_gamma, 2U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_opts, 3U, c1_l_sf_marshallOut,
    c1_k_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_alpha, 4U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_current_action, 5U,
    c1_j_sf_marshallOut, c1_i_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_action_idx, 6U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_DYN_EMX_IMPORTABLE((void *)&c1_action_idx_pr_data, (
    const int32_T *)&c1_action_idx_pr_size, NULL, 0, 7, (void *)
    c1_k_sf_marshallOut, (void *)c1_j_sf_marshallIn, (void *)
    &c1_action_idx_pr_data, false);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_Height, MAX_uint32_T,
    c1_sf_marshallOut, c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_Ts, 9U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_DYN_EMX_IMPORTABLE((void *)
    &c1_next_stPrimee_idx_data, (const int32_T *)&c1_next_stPrimee_idx_size,
    NULL, 0, 10, (void *)c1_f_sf_marshallOut, (void *)c1_f_sf_marshallIn, (void *)
    &c1_next_stPrimee_idx_data, false);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_st_chosen, 11U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_Height_prime, 12U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_Ts_prime, 13U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_DYN_EMX_IMPORTABLE((void *)&c1_next_state_idx_data,
    (const int32_T *)&c1_next_state_idx_size, NULL, 0, -1, (void *)
    c1_f_sf_marshallOut, (void *)c1_f_sf_marshallIn, (void *)
    &c1_next_state_idx_data, false);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_next_reward, 15U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_DYN_EMX_IMPORTABLE((void *)&c1_Period_data, (const
    int32_T *)&c1_Period_size, NULL, 0, -1, (void *)c1_f_sf_marshallOut, (void *)
    c1_f_sf_marshallIn, (void *)&c1_Period_data, false);
  _SFD_SYMBOL_SCOPE_ADD_EML_DYN_EMX_IMPORTABLE((void *)&c1_R_idx_data, (const
    int32_T *)&c1_R_idx_size, NULL, 0, -1, (void *)c1_f_sf_marshallOut, (void *)
    c1_f_sf_marshallIn, (void *)&c1_R_idx_data, false);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_MaxR, 18U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_DYN_EMX_IMPORTABLE((void *)&c1_Final_Rew_data, (
    const int32_T *)&c1_Final_Rew_size, NULL, 0, -1, (void *)c1_f_sf_marshallOut,
    (void *)c1_f_sf_marshallIn, (void *)&c1_Final_Rew_data, false);
  _SFD_SYMBOL_SCOPE_ADD_EML_DYN_EMX_IMPORTABLE(c1_In->data, (const int32_T *)
    c1_In->size, NULL, 0, 20, (void *)c1_c_sf_marshallOut, (void *)
    c1_c_sf_marshallIn, (void *)c1_In, true);
  _SFD_SYMBOL_SCOPE_ADD_EML_DYN_EMX_IMPORTABLE(c1_Output->data, (const int32_T *)
    c1_Output->size, NULL, 0, 21, (void *)c1_e_sf_marshallOut, (void *)
    c1_e_sf_marshallIn, (void *)c1_Output, true);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_nodes, 22U, c1_j_sf_marshallOut,
    c1_i_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(NULL, MAX_uint32_T, c1_i_sf_marshallOut,
    c1_h_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_DYN_EMX_IMPORTABLE((void *)&chartInstance->c1_m_data,
    (const int32_T *)&c1_m_size, NULL, 0, 24, (void *)c1_d_sf_marshallOut, (void
    *)c1_d_sf_marshallIn, (void *)&chartInstance->c1_m_data, false);
  _SFD_SYMBOL_SCOPE_ADD_EML_DYN_EMX_IMPORTABLE((void *)
    &chartInstance->c1_range_data, (const int32_T *)&c1_range_size, NULL, 0, -1,
    (void *)c1_d_sf_marshallOut, (void *)c1_d_sf_marshallIn, (void *)
    &chartInstance->c1_range_data, false);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(NULL, MAX_uint32_T, c1_i_sf_marshallOut,
    c1_h_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_DYN_EMX_IMPORTABLE((void *)&c1_m11_data, (const
    int32_T *)&c1_m11_size, NULL, 0, 27, (void *)c1_b_sf_marshallOut, (void *)
    c1_b_sf_marshallIn, (void *)&c1_m11_data, false);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c1_net, 28U, c1_h_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c1_tr, 29U, c1_h_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c1_y_prime, 30U, c1_h_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_c, 31U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_Q_prev, 32U, c1_g_sf_marshallOut,
    c1_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_k, 33U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_DYN_EMX_IMPORTABLE((void *)&c1_Height_data, (const
    int32_T *)&c1_Height_size, NULL, 0, -1, (void *)c1_f_sf_marshallOut, (void *)
    c1_f_sf_marshallIn, (void *)&c1_Height_data, false);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_Period, MAX_uint32_T,
    c1_sf_marshallOut, c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_R_idx, MAX_uint32_T,
    c1_sf_marshallOut, c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_DYN_EMX_IMPORTABLE(c1_Out->data, (const int32_T *)
    c1_Out->size, NULL, 0, -1, (void *)c1_e_sf_marshallOut, (void *)
    c1_e_sf_marshallIn, (void *)c1_Out, true);
  _SFD_SYMBOL_SCOPE_ADD_EML_DYN_EMX_IMPORTABLE((void *)
    &chartInstance->c1_Out_data, (const int32_T *)&c1_Out_size, NULL, 0, -1,
    (void *)c1_d_sf_marshallOut, (void *)c1_d_sf_marshallIn, (void *)
    &chartInstance->c1_Out_data, false);
  _SFD_SYMBOL_SCOPE_ADD_EML_DYN_EMX_IMPORTABLE(c1_Input->data, (const int32_T *)
    c1_Input->size, NULL, 0, -1, (void *)c1_c_sf_marshallOut, (void *)
    c1_c_sf_marshallIn, (void *)c1_Input, true);
  _SFD_SYMBOL_SCOPE_ADD_EML_DYN_EMX_IMPORTABLE((void *)&c1_b_range_data, (const
    int32_T *)&c1_b_range_size, NULL, 0, -1, (void *)c1_b_sf_marshallOut, (void *)
    c1_b_sf_marshallIn, (void *)&c1_b_range_data, false);
  _SFD_SYMBOL_SCOPE_ADD_EML_DYN_EMX_IMPORTABLE((void *)&c1_Input_data, (const
    int32_T *)&c1_Input_size, NULL, 0, -1, (void *)c1_b_sf_marshallOut, (void *)
    c1_b_sf_marshallIn, (void *)&c1_Input_data, false);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_next_state_idx, MAX_uint32_T,
    c1_sf_marshallOut, c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_Final_Rew, MAX_uint32_T,
    c1_sf_marshallOut, c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargin, 34U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargout, 35U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_d_y, 36U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 3);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 5);
  for (c1_i0 = 0; c1_i0 < 2; c1_i0++) {
    c1_current_action1[c1_i0] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 6);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 8);
  c1_PptoAvg = 50000.0;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 10);
  if (CV_EML_IF(0, 1, 0, CV_RELATIONAL_EVAL(4U, 0U, 0, c1_get_ove_count
        (chartInstance, 0), 400.0, -1, 5U, c1_get_ove_count(chartInstance, 0) >=
        400.0))) {
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 12);
    c1_gamma = 0.6;
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 13);
    c1_opts.StepRatio = 0.01;
  } else {
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 16);
    c1_gamma = 0.4;
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 17);
    c1_opts.StepRatio = 0.1;
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 23);
  if (CV_EML_IF(0, 1, 1, CV_RELATIONAL_EVAL(4U, 0U, 1, c1_get_ove_count
        (chartInstance, 0), 200.0, -1, 5U, c1_get_ove_count(chartInstance, 0) >=
        200.0))) {
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 25);
    c1_alpha = 0.6;
  } else {
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 29);
    c1_alpha = 0.4;
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 35);
  for (c1_i1 = 0; c1_i1 < 2; c1_i1++) {
    c1_current_action[c1_i1] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 36);
  c1_action_idx = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 38);
  for (c1_i2 = 0; c1_i2 < 24; c1_i2++) {
    c1_dv0[c1_i2] = c1_get_action(chartInstance, c1_i2);
  }

  c1_j_emlrt_marshallIn(chartInstance, sf_mex_call_debug
                        (sfGlobalDebugInstanceStruct, "datasample", 1U, 2U, 14,
    c1_emlrt_marshallOut(chartInstance, c1_dv0), 14, c1_b_emlrt_marshallOut
    (chartInstance, 1.0)), "datasample", c1_dv1);
  for (c1_i3 = 0; c1_i3 < 2; c1_i3++) {
    c1_current_action[c1_i3] = c1_dv1[c1_i3];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 39);
  c1_b_current_action = c1_current_action[0];
  for (c1_i4 = 0; c1_i4 < 12; c1_i4++) {
    c1_bv0[c1_i4] = (c1_get_action(chartInstance, c1_i4) == c1_b_current_action);
  }

  c1_c_current_action = c1_current_action[1];
  for (c1_i5 = 0; c1_i5 < 12; c1_i5++) {
    c1_bv1[c1_i5] = (c1_get_action(chartInstance, c1_i5 + 12) ==
                     c1_c_current_action);
  }

  c1_idx = 0;
  c1_ii_size[0] = 12;
  c1_ii = 1;
  exitg1 = false;
  while ((!exitg1) && (c1_ii < 13)) {
    if (c1_bv0[c1_ii - 1] && c1_bv1[c1_ii - 1]) {
      c1_idx++;
      c1_ii_data[c1_idx - 1] = c1_ii;
      if (c1_idx >= 12) {
        exitg1 = true;
      } else {
        c1_ii++;
      }
    } else {
      c1_ii++;
    }
  }

  if (1 > c1_idx) {
    c1_i6 = 0;
  } else {
    c1_i6 = c1_idx;
  }

  c1_ii_size[0] = c1_i6;
  c1_action_idx_pr_size[0] = c1_ii_size[0];
  c1_loop_ub = c1_ii_size[0] - 1;
  for (c1_i7 = 0; c1_i7 <= c1_loop_ub; c1_i7++) {
    c1_action_idx_pr_data[c1_i7] = (real_T)c1_ii_data[c1_i7];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 40);
  (real_T)sf_eml_array_bounds_check(sfGlobalDebugInstanceStruct,
    chartInstance->S, 1U, 985, 1, MAX_uint32_T, 1, 1, c1_action_idx_pr_size[0]);
  c1_action_idx_pr[0] = c1_action_idx_pr_size[0];
  c1_action_idx_pr[1] = 1;
  c1_action_idx = c1_action_idx_pr_data[0];
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 42);
  c1_set_B2(chartInstance, 0, c1_get_action(chartInstance,
             sf_eml_array_bounds_check(sfGlobalDebugInstanceStruct,
              chartInstance->S, 1U, 1010, 20, 13U, (int32_T)c1_action_idx, 1, 12)
             - 1));
  ssUpdateDataStoreLog_wrapper(chartInstance->S, 0);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 45);
  c1_set_K2(chartInstance, 0, c1_get_action(chartInstance, c1__s32_d_
             (chartInstance, c1_action_idx, 1U, 1080U, 20U) + 11));
  ssUpdateDataStoreLog_wrapper(chartInstance->S, 2);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 48);
  c1_set_Stor_eta(chartInstance, sf_eml_array_bounds_check
                  (sfGlobalDebugInstanceStruct, chartInstance->S, 1U, 1137, 9,
                   12U, (int32_T)sf_integer_check(chartInstance->S, 1U, 1137U,
    9U, c1_get_ove_count(chartInstance, 0)), 1, 3000) - 1, (real_T)
                  c1_get_state_idx(chartInstance, 0));
  ssUpdateDataStoreLog_wrapper(chartInstance->S, 11);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 49);
  c1_set_Stor_eta(chartInstance, sf_eml_array_bounds_check
                  (sfGlobalDebugInstanceStruct, chartInstance->S, 1U, 1178, 9,
                   12U, (int32_T)sf_integer_check(chartInstance->S, 1U, 1178U,
    9U, c1_get_ove_count(chartInstance, 0)), 1, 3000) + 2999, c1_action_idx);
  ssUpdateDataStoreLog_wrapper(chartInstance->S, 11);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 51);
  c1_Height = c1_get_states_space(chartInstance, (int16_T)
    sf_eml_array_bounds_check(sfGlobalDebugInstanceStruct, chartInstance->S, 1U,
    1228, 26, 17U, (int32_T)c1_get_state_idx(chartInstance, 0), 1, 72) - 1);
  _SFD_SYMBOL_SWITCH(8U, 8U);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 53);
  c1_Ts = c1_get_states_space(chartInstance, (int16_T)sf_eml_array_bounds_check
    (sfGlobalDebugInstanceStruct, chartInstance->S, 1U, 1322, 26, 17U, (int32_T)
     c1_get_state_idx(chartInstance, 0), 1, 72) + 71);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 56);
  for (c1_i8 = 0; c1_i8 < 72; c1_i8++) {
    c1_x[c1_i8] = (c1_get_states_space(chartInstance, c1_i8) == c1_Height);
  }

  for (c1_i9 = 0; c1_i9 < 72; c1_i9++) {
    c1_bv2[c1_i9] = (c1_get_states_space(chartInstance, c1_i9 + 72) == c1_Ts);
  }

  for (c1_i10 = 0; c1_i10 < 72; c1_i10++) {
    c1_x[c1_i10] = (c1_x[c1_i10] && c1_bv2[c1_i10] && (c1_get_states_space
      (chartInstance, c1_i10 + 144) == c1_get_B2(chartInstance, 0)) &&
                    (c1_get_states_space(chartInstance, c1_i10 + 216) ==
                     c1_get_K2(chartInstance, 0)));
  }

  c1_b_idx = 0;
  c1_b_ii_size[0] = 72;
  c1_b_ii = 1;
  exitg1 = false;
  while ((!exitg1) && (c1_b_ii < 73)) {
    if (c1_x[c1_b_ii - 1]) {
      c1_b_idx++;
      c1_b_ii_data[c1_b_idx - 1] = c1_b_ii;
      if (c1_b_idx >= 72) {
        exitg1 = true;
      } else {
        c1_b_ii++;
      }
    } else {
      c1_b_ii++;
    }
  }

  if (1 > c1_b_idx) {
    c1_i11 = 0;
  } else {
    c1_i11 = c1_b_idx;
  }

  c1_b_ii_size[0] = c1_i11;
  c1_next_stPrimee_idx_size[0] = c1_b_ii_size[0];
  c1_b_loop_ub = c1_b_ii_size[0] - 1;
  for (c1_i12 = 0; c1_i12 <= c1_b_loop_ub; c1_i12++) {
    c1_next_stPrimee_idx_data[c1_i12] = (real_T)c1_b_ii_data[c1_i12];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 60);
  c1_st_chosen = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 61);
  if (!c1_b0) {
    for (c1_i13 = 0; c1_i13 < 72; c1_i13++) {
      c1_dv2[c1_i13] = 1.0 + (real_T)c1_i13;
    }

    c1_b0 = true;
  }

  c1_st_chosen = c1_l_emlrt_marshallIn(chartInstance, sf_mex_call_debug
    (sfGlobalDebugInstanceStruct, "datasample", 1U, 2U, 14,
     c1_c_emlrt_marshallOut(chartInstance, c1_dv2), 14, c1_b_emlrt_marshallOut
     (chartInstance, 1.0)), "datasample");
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 64);
  c1_Height_prime = c1_get_states_space(chartInstance, sf_eml_array_bounds_check
    (sfGlobalDebugInstanceStruct, chartInstance->S, 1U, 1718, 26, 17U, (int32_T)
     sf_integer_check(chartInstance->S, 1U, 1718U, 26U, c1_st_chosen), 1, 72) -
    1);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 65);
  c1_Ts_prime = c1_get_states_space(chartInstance, c1__s32_d_(chartInstance,
    c1_st_chosen, 1U, 1809U, 26U) + 71);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 67);
  for (c1_i14 = 0; c1_i14 < 72; c1_i14++) {
    c1_x[c1_i14] = (c1_get_states_space(chartInstance, c1_i14) ==
                    c1_Height_prime);
  }

  for (c1_i15 = 0; c1_i15 < 72; c1_i15++) {
    c1_bv2[c1_i15] = (c1_get_states_space(chartInstance, c1_i15 + 72) ==
                      c1_Ts_prime);
  }

  for (c1_i16 = 0; c1_i16 < 72; c1_i16++) {
    c1_x[c1_i16] = (c1_x[c1_i16] && c1_bv2[c1_i16] && (c1_get_states_space
      (chartInstance, c1_i16 + 144) == c1_get_B2(chartInstance, 0)) &&
                    (c1_get_states_space(chartInstance, c1_i16 + 216) ==
                     c1_get_K2(chartInstance, 0)));
  }

  c1_c_idx = 0;
  c1_b_ii_size[0] = 72;
  c1_c_ii = 1;
  exitg1 = false;
  while ((!exitg1) && (c1_c_ii < 73)) {
    if (c1_x[c1_c_ii - 1]) {
      c1_c_idx++;
      c1_b_ii_data[c1_c_idx - 1] = c1_c_ii;
      if (c1_c_idx >= 72) {
        exitg1 = true;
      } else {
        c1_c_ii++;
      }
    } else {
      c1_c_ii++;
    }
  }

  if (1 > c1_c_idx) {
    c1_i17 = 0;
  } else {
    c1_i17 = c1_c_idx;
  }

  c1_b_ii_size[0] = c1_i17;
  c1_next_state_idx_size[0] = c1_b_ii_size[0];
  c1_c_loop_ub = c1_b_ii_size[0] - 1;
  for (c1_i18 = 0; c1_i18 <= c1_c_loop_ub; c1_i18++) {
    c1_next_state_idx_data[c1_i18] = (real_T)c1_b_ii_data[c1_i18];
  }

  _SFD_SYMBOL_SWITCH(14U, 14U);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 70);
  c1_b_ii_size[0] = c1_next_stPrimee_idx_size[0];
  c1_d_loop_ub = c1_next_stPrimee_idx_size[0] - 1;
  for (c1_i19 = 0; c1_i19 <= c1_d_loop_ub; c1_i19++) {
    c1_b_ii_data[c1_i19] = sf_eml_array_bounds_check(sfGlobalDebugInstanceStruct,
      chartInstance->S, 1U, 2092, 17, MAX_uint32_T, (int32_T)
      c1_next_stPrimee_idx_data[c1_i19], 1, 10000);
  }

  c1_tmp_size[0] = c1_next_stPrimee_idx_size[0];
  c1_e_loop_ub = c1_next_stPrimee_idx_size[0] - 1;
  for (c1_i20 = 0; c1_i20 <= c1_e_loop_ub; c1_i20++) {
    c1_b_tmp_data[c1_i20] = c1_get_num_R(chartInstance, c1__s32_d_(chartInstance,
      c1_next_stPrimee_idx_data[c1_i20], 1U, 2117U, 17U) - 1);
  }

  _SFD_SUB_ASSIGN_SIZE_CHECK_1D(c1_b_ii_size[0], c1_tmp_size[0]);
  c1_f_loop_ub = c1_tmp_size[0] - 1;
  for (c1_i21 = 0; c1_i21 <= c1_f_loop_ub; c1_i21++) {
    c1_i22 = c1_b_tmp_data[c1_i21] + 1;
    if (c1_i22 > 32767) {
      CV_SATURATION_EVAL(4, 0, 0, 0, 1);
      c1_i22 = 32767;
      _SFD_OVERFLOW_DETECTION(SFDB_SATURATE, 1U, 2111U, 26U);
    } else {
      if (CV_SATURATION_EVAL(4, 0, 0, 0, c1_i22 < -32768)) {
        c1_i22 = -32768;
        _SFD_OVERFLOW_DETECTION(SFDB_SATURATE, 1U, 2111U, 26U);
      }
    }

    c1_set_num_R(chartInstance, c1_b_ii_data[c1_i21] - 1, (int16_T)c1_i22);
  }

  ssUpdateDataStoreLog_wrapper(chartInstance->S, 13);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 72);
  c1_a = c1_Height;
  c1_next_reward = 50000.0 / (c1_a * c1_a);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 74);
  c1_b_ii_size[0] = c1_next_stPrimee_idx_size[0];
  c1_g_loop_ub = c1_next_stPrimee_idx_size[0] - 1;
  for (c1_i23 = 0; c1_i23 <= c1_g_loop_ub; c1_i23++) {
    c1_b_ii_data[c1_i23] = c1__s32_d_(chartInstance,
      c1_next_stPrimee_idx_data[c1_i23], 1U, 2230U, 17U);
  }

  c1_b_tmp_size[0] = c1_next_stPrimee_idx_size[0];
  c1_h_loop_ub = c1_next_stPrimee_idx_size[0] - 1;
  for (c1_i24 = 0; c1_i24 <= c1_h_loop_ub; c1_i24++) {
    c1_c_tmp_data[c1_i24] = (int16_T)sf_eml_array_bounds_check
      (sfGlobalDebugInstanceStruct, chartInstance->S, 1U, 2248, 24, MAX_uint32_T,
       (int32_T)c1_get_num_R(chartInstance, c1__s32_d_(chartInstance,
         c1_next_stPrimee_idx_data[c1_i24], 1U, 2254U, 17U) - 1), 1, 50);
  }

  c1_d_ii[0] = c1_b_ii_size[0];
  c1_d_ii[1] = c1_b_tmp_size[0];
  c1_i_loop_ub = c1_d_ii[1] - 1;
  for (c1_i25 = 0; c1_i25 <= c1_i_loop_ub; c1_i25++) {
    c1_j_loop_ub = c1_d_ii[0] - 1;
    for (c1_i26 = 0; c1_i26 <= c1_j_loop_ub; c1_i26++) {
      c1_set_R(chartInstance, (c1_b_ii_data[c1_i26] + 72 * (c1_c_tmp_data[c1_i25]
                 - 1)) - 1, c1_next_reward);
    }
  }

  ssUpdateDataStoreLog_wrapper(chartInstance->S, 10);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 75);
  c1_Height_size[0] = c1_next_stPrimee_idx_size[0];
  c1_k_loop_ub = c1_next_stPrimee_idx_size[0] - 1;
  for (c1_i27 = 0; c1_i27 <= c1_k_loop_ub; c1_i27++) {
    c1_Height_data[c1_i27] = c1_get_states_space(chartInstance, c1__s32_d_
      (chartInstance, c1_next_stPrimee_idx_data[c1_i27], 1U, 2315U, 17U) - 1);
  }

  _SFD_SYMBOL_SWITCH(8U, 34U);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 76);
  (real_T)sf_eml_array_bounds_check(sfGlobalDebugInstanceStruct,
    chartInstance->S, 1U, 2359, 1, MAX_uint32_T, 1, 1, c1_Height_size[0]);
  c1_b_Height[0] = c1_Height_size[0];
  c1_b_Height[1] = 1;
  c1_Height = c1_Height_data[0];
  _SFD_SYMBOL_SWITCH(8U, 8U);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 78);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 80);
  c1_Period_size[0] = c1_next_stPrimee_idx_size[0];
  c1_l_loop_ub = c1_next_stPrimee_idx_size[0] - 1;
  for (c1_i28 = 0; c1_i28 <= c1_l_loop_ub; c1_i28++) {
    c1_Period_data[c1_i28] = c1_get_states_space(chartInstance, c1__s32_d_
      (chartInstance, c1_next_stPrimee_idx_data[c1_i28], 1U, 2435U, 17U) + 71);
  }

  _SFD_SYMBOL_SWITCH(16U, 16U);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 81);
  (real_T)sf_eml_array_bounds_check(sfGlobalDebugInstanceStruct,
    chartInstance->S, 1U, 2479, 1, MAX_uint32_T, 1, 1, c1_Period_size[0]);
  c1_b_Period[0] = c1_Period_size[0];
  c1_b_Period[1] = 1;
  c1_Period = c1_Period_data[0];
  _SFD_SYMBOL_SWITCH(16U, 35U);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 83);
  for (c1_i29 = 0; c1_i29 < 72; c1_i29++) {
    c1_bv2[c1_i29] = (c1_get_states_space(chartInstance, c1_i29) == c1_Height);
  }

  for (c1_i30 = 0; c1_i30 < 72; c1_i30++) {
    c1_x[c1_i30] = (c1_get_states_space(chartInstance, c1_i30 + 72) == c1_Period);
  }

  c1_d_idx = 0;
  c1_b_ii_size[0] = 72;
  c1_e_ii = 1;
  exitg1 = false;
  while ((!exitg1) && (c1_e_ii < 73)) {
    if (c1_bv2[c1_e_ii - 1] && c1_x[c1_e_ii - 1]) {
      c1_d_idx++;
      c1_b_ii_data[c1_d_idx - 1] = c1_e_ii;
      if (c1_d_idx >= 72) {
        exitg1 = true;
      } else {
        c1_e_ii++;
      }
    } else {
      c1_e_ii++;
    }
  }

  if (1 > c1_d_idx) {
    c1_i31 = 0;
  } else {
    c1_i31 = c1_d_idx;
  }

  c1_b_ii_size[0] = c1_i31;
  c1_R_idx_size[0] = c1_b_ii_size[0];
  c1_m_loop_ub = c1_b_ii_size[0] - 1;
  for (c1_i32 = 0; c1_i32 <= c1_m_loop_ub; c1_i32++) {
    c1_R_idx_data[c1_i32] = (real_T)c1_b_ii_data[c1_i32];
  }

  _SFD_SYMBOL_SWITCH(17U, 17U);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 84);
  (real_T)sf_eml_array_bounds_check(sfGlobalDebugInstanceStruct,
    chartInstance->S, 1U, 2590, 1, MAX_uint32_T, 1, 1, c1_R_idx_size[0]);
  c1_b_R_idx[0] = c1_R_idx_size[0];
  c1_b_R_idx[1] = 1;
  c1_R_idx = c1_R_idx_data[0];
  _SFD_SYMBOL_SWITCH(17U, 36U);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 86);
  c1_c_R_idx = sf_eml_array_bounds_check(sfGlobalDebugInstanceStruct,
    chartInstance->S, 1U, 2635, 5, MAX_uint32_T, (int32_T)c1_R_idx, 1, 72) - 1;
  for (c1_i33 = 0; c1_i33 < 50; c1_i33++) {
    c1_dv3[c1_i33] = c1_get_R(chartInstance, c1_c_R_idx + 72 * c1_i33);
  }

  c1_set_MeanR(chartInstance, sf_eml_array_bounds_check
               (sfGlobalDebugInstanceStruct, chartInstance->S, 1U, 2614, 14, 4U,
                (int32_T)c1_R_idx, 1, 100) - 1, (real_T)c1_rdivide(chartInstance,
    c1_sum(chartInstance, c1_dv3), c1_get_num_R(chartInstance,
    sf_eml_array_bounds_check(sfGlobalDebugInstanceStruct, chartInstance->S, 1U,
    2648, 12, 14U, (int32_T)c1_R_idx, 1, 10000) - 1)));
  ssUpdateDataStoreLog_wrapper(chartInstance->S, 3);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 90);
  c1_varargin_1 = c1_get_MeanR(chartInstance, c1__s32_d_(chartInstance, c1_R_idx,
    1U, 2739U, 12U) - 1);
  c1_MaxR = c1_varargin_1;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 92);
  c1_A = c1_get_MeanR(chartInstance, c1__s32_d_(chartInstance, c1_R_idx, 1U,
    2797U, 14U) - 1);
  c1_B = c1_MaxR;
  c1_set_Final_Rew_matt(chartInstance, c1__s32_d_(chartInstance, c1_R_idx, 1U,
    2790U, 5U) - 1, c1_A / c1_B);
  ssUpdateDataStoreLog_wrapper(chartInstance->S, 1);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 94);
  c1_Final_Rew_size[0] = c1_next_stPrimee_idx_size[0];
  c1_n_loop_ub = c1_next_stPrimee_idx_size[0] - 1;
  for (c1_i34 = 0; c1_i34 <= c1_n_loop_ub; c1_i34++) {
    c1_Final_Rew_data[c1_i34] = c1_get_Final_Rew_matt(chartInstance, c1__s32_d_
      (chartInstance, c1_next_stPrimee_idx_data[c1_i34], 1U, 2860U, 17U) - 1);
  }

  _SFD_SYMBOL_SWITCH(19U, 19U);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 96);
  c1_b_Final_Rew_size[0] = c1_Final_Rew_size[0];
  c1_o_loop_ub = c1_Final_Rew_size[0] - 1;
  for (c1_i35 = 0; c1_i35 <= c1_o_loop_ub; c1_i35++) {
    c1_b_Final_Rew_data[c1_i35] = c1_Final_Rew_data[c1_i35];
  }

  c1_mpower(chartInstance, c1_b_Final_Rew_data, c1_b_Final_Rew_size,
            c1_d_tmp_data, c1_c_tmp_size);
  c1_Final_Rew_size[0] = c1_c_tmp_size[0];
  c1_p_loop_ub = c1_c_tmp_size[0] - 1;
  for (c1_i36 = 0; c1_i36 <= c1_p_loop_ub; c1_i36++) {
    c1_Final_Rew_data[c1_i36] = c1_d_tmp_data[c1_i36];
  }

  _SFD_SYMBOL_SWITCH(19U, 19U);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 98);
  c1_set_Stor_eta(chartInstance, sf_eml_array_bounds_check
                  (sfGlobalDebugInstanceStruct, chartInstance->S, 1U, 2948, 9,
                   12U, (int32_T)sf_integer_check(chartInstance->S, 1U, 2948U,
    9U, c1_get_ove_count(chartInstance, 0)), 1, 3000) + 5999, c1_get_B2
                  (chartInstance, 0));
  ssUpdateDataStoreLog_wrapper(chartInstance->S, 11);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 99);
  c1_set_Stor_eta(chartInstance, sf_eml_array_bounds_check
                  (sfGlobalDebugInstanceStruct, chartInstance->S, 1U, 2982, 9,
                   12U, (int32_T)sf_integer_check(chartInstance->S, 1U, 2982U,
    9U, c1_get_ove_count(chartInstance, 0)), 1, 3000) + 8999, c1_get_K2
                  (chartInstance, 0));
  ssUpdateDataStoreLog_wrapper(chartInstance->S, 11);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 99);
  c1_set_Stor_eta(chartInstance, sf_eml_array_bounds_check
                  (sfGlobalDebugInstanceStruct, chartInstance->S, 1U, 3009, 9,
                   12U, (int32_T)sf_integer_check(chartInstance->S, 1U, 3009U,
    9U, c1_get_ove_count(chartInstance, 0)), 1, 3000) + 11999, 50000.0);
  ssUpdateDataStoreLog_wrapper(chartInstance->S, 11);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 105);
  c1_b_action_idx = c1__s32_d_(chartInstance, c1_action_idx, 1U, 3313U, 10U) - 1;
  for (c1_i37 = 0; c1_i37 < 2; c1_i37++) {
    c1_dv4[c1_i37] = c1_get_action(chartInstance, c1_b_action_idx + 12 * c1_i37);
  }

  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "disp", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "horzcat", 1U,
    6U, 14, c1_d_emlrt_marshallOut(chartInstance, c1_cv1), 14, sf_mex_call_debug
    (sfGlobalDebugInstanceStruct, "num2str", 1U, 1U, 14, c1_e_emlrt_marshallOut
     (chartInstance, c1_get_state_idx(chartInstance, 0))), 14,
    c1_f_emlrt_marshallOut(chartInstance, c1_cv2), 14, sf_mex_call_debug
    (sfGlobalDebugInstanceStruct, "num2str", 1U, 1U, 14, c1_g_emlrt_marshallOut
     (chartInstance, c1_next_state_idx_data, c1_next_state_idx_size)), 14,
    c1_d_emlrt_marshallOut(chartInstance, c1_cv3), 14, sf_mex_call_debug
    (sfGlobalDebugInstanceStruct, "num2str", 1U, 1U, 14, c1_h_emlrt_marshallOut
     (chartInstance, c1_dv4))));
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 107);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "disp", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "horzcat", 1U,
    2U, 14, c1_i_emlrt_marshallOut(chartInstance, c1_cv4), 14, sf_mex_call_debug
    (sfGlobalDebugInstanceStruct, "num2str", 1U, 1U, 14, c1_g_emlrt_marshallOut
     (chartInstance, c1_Final_Rew_data, c1_Final_Rew_size))));
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 110);
  if (CV_EML_IF(0, 1, 2, CV_RELATIONAL_EVAL(4U, 0U, 2, c1_get_ove_count
        (chartInstance, 0), c1_get_NN_dd(chartInstance, 0), -1, 4U,
        c1_get_ove_count(chartInstance, 0) > c1_get_NN_dd(chartInstance, 0)))) {
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 114);
    if (CV_EML_IF(0, 1, 3, CV_RELATIONAL_EVAL(4U, 0U, 3, c1_get_ove_count
          (chartInstance, 0), c1_get_Memory_D_size(chartInstance, 0), -1, 4U,
          c1_get_ove_count(chartInstance, 0) > c1_get_Memory_D_size
          (chartInstance, 0)))) {
      _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 116);
      if ((5000.0 - c1_get_Memory_D_size(chartInstance, 0)) + 1.0 > 5000.0) {
        c1_i40 = 1;
        c1_i41 = 0;
      } else {
        c1_i40 = sf_eml_array_bounds_check(sfGlobalDebugInstanceStruct,
          chartInstance->S, 1U, 3642, 19, MAX_uint32_T, (int32_T)
          sf_integer_check(chartInstance->S, 1U, 3642U, 19U, (5000.0 -
          c1_get_Memory_D_size(chartInstance, 0)) + 1.0), 1, 5000);
        c1_i41 = 5000;
      }

      c1_i43 = c1_In->size[0] * c1_In->size[1];
      c1_In->size[0] = (c1_i41 - c1_i40) + 1;
      c1_In->size[1] = 4;
      c1_emxEnsureCapacity_real_T(chartInstance, c1_In, c1_i43, &c1_q_emlrtRTEI);
      for (c1_i46 = 0; c1_i46 < 4; c1_i46++) {
        c1_q_loop_ub = c1_i41 - c1_i40;
        for (c1_i48 = 0; c1_i48 <= c1_q_loop_ub; c1_i48++) {
          c1_In->data[c1_i48 + c1_In->size[0] * c1_i46] = c1_get_Memory_Stor
            (chartInstance, ((c1_i40 + c1_i48) + 5000 * (4 + c1_i46)) - 1);
        }
      }

      _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 118);
      if ((5000.0 - c1_get_Memory_D_size(chartInstance, 0)) + 1.0 > 5000.0) {
        c1_i49 = 1;
        c1_i50 = 0;
      } else {
        c1_i49 = sf_eml_array_bounds_check(sfGlobalDebugInstanceStruct,
          chartInstance->S, 1U, 3724, 19, MAX_uint32_T, (int32_T)
          sf_integer_check(chartInstance->S, 1U, 3724U, 19U, (5000.0 -
          c1_get_Memory_D_size(chartInstance, 0)) + 1.0), 1, 5000);
        c1_i50 = 5000;
      }

      c1_i53 = c1_Output->size[0] * c1_Output->size[1];
      c1_Output->size[0] = (c1_i50 - c1_i49) + 1;
      c1_Output->size[1] = 22;
      c1_emxEnsureCapacity_real_T(chartInstance, c1_Output, c1_i53,
        &c1_s_emlrtRTEI);
      for (c1_i54 = 0; c1_i54 < 22; c1_i54++) {
        c1_r_loop_ub = c1_i50 - c1_i49;
        for (c1_i57 = 0; c1_i57 <= c1_r_loop_ub; c1_i57++) {
          c1_Output->data[c1_i57 + c1_Output->size[0] * c1_i54] =
            c1_get_Memory_Stor(chartInstance, ((c1_i49 + c1_i57) + 5000 * (8 +
            c1_i54)) - 1);
        }
      }

      _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 120);
      for (c1_i56 = 0; c1_i56 < 2; c1_i56++) {
        c1_nodes[c1_i56] = 4.0 + 8.0 * (real_T)c1_i56;
      }
    } else {
      _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 124);
      c1_i39 = c1_In->size[0] * c1_In->size[1];
      c1_In->size[0] = 5000;
      c1_In->size[1] = 4;
      c1_emxEnsureCapacity_real_T(chartInstance, c1_In, c1_i39, &c1_p_emlrtRTEI);
      for (c1_i42 = 0; c1_i42 < 4; c1_i42++) {
        for (c1_i44 = 0; c1_i44 < 5000; c1_i44++) {
          c1_In->data[c1_i44 + c1_In->size[0] * c1_i42] = c1_get_Memory_Stor
            (chartInstance, c1_i44 + 5000 * (4 + c1_i42));
        }
      }

      _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 126);
      c1_i45 = c1_Output->size[0] * c1_Output->size[1];
      c1_Output->size[0] = 5000;
      c1_Output->size[1] = 22;
      c1_emxEnsureCapacity_real_T(chartInstance, c1_Output, c1_i45,
        &c1_r_emlrtRTEI);
      for (c1_i47 = 0; c1_i47 < 22; c1_i47++) {
        for (c1_i51 = 0; c1_i51 < 5000; c1_i51++) {
          c1_Output->data[c1_i51 + c1_Output->size[0] * c1_i47] =
            c1_get_Memory_Stor(chartInstance, c1_i51 + 5000 * (8 + c1_i47));
        }
      }

      _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 128U);
      for (c1_i52 = 0; c1_i52 < 2; c1_i52++) {
        c1_nodes[c1_i52] = 4.0 + 8.0 * (real_T)c1_i52;
      }
    }

    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 139U);
    _SFD_SYMBOL_SWITCH(23U, 23U);
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 140U);
    c1_i55 = c1_Out->size[0] * c1_Out->size[1];
    c1_Out->size[0] = c1_Output->size[0];
    c1_Out->size[1] = 22;
    c1_emxEnsureCapacity_real_T(chartInstance, c1_Out, c1_i55, &c1_t_emlrtRTEI);
    c1_b_Out = c1_Out->size[0];
    c1_c_Out = c1_Out->size[1];
    c1_s_loop_ub = c1_Output->size[0] * c1_Output->size[1] - 1;
    for (c1_i58 = 0; c1_i58 <= c1_s_loop_ub; c1_i58++) {
      c1_Out->data[c1_i58] = c1_Output->data[c1_i58];
    }

    c1_emxInit_real_T(chartInstance, &c1_b_varargin_1, 2, &c1_u_emlrtRTEI);
    _SFD_SYMBOL_SWITCH(23U, 37U);
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 141U);
    c1_i59 = c1_b_varargin_1->size[0] * c1_b_varargin_1->size[1];
    c1_b_varargin_1->size[0] = c1_Out->size[0];
    c1_b_varargin_1->size[1] = 22;
    c1_emxEnsureCapacity_real_T(chartInstance, c1_b_varargin_1, c1_i59,
      &c1_u_emlrtRTEI);
    c1_c_varargin_1 = c1_b_varargin_1->size[0];
    c1_d_varargin_1 = c1_b_varargin_1->size[1];
    c1_t_loop_ub = c1_Out->size[0] * c1_Out->size[1] - 1;
    for (c1_i60 = 0; c1_i60 <= c1_t_loop_ub; c1_i60++) {
      c1_b_varargin_1->data[c1_i60] = c1_Out->data[c1_i60];
    }

    c1_m_size[0] = c1_b_varargin_1->size[0];
    c1_vstride = c1_b_varargin_1->size[0];
    for (c1_j = 1; c1_j <= c1_vstride; c1_j++) {
      c1_ixstart = c1_j;
      c1_ixstop = c1_j + 21 * c1_vstride;
      c1_mtmp = c1_b_varargin_1->data[c1_j - 1];
      if (muDoubleScalarIsNaN(c1_mtmp)) {
        if (c1_vstride != 0) {
        } else {
          c1_e_y = NULL;
          sf_mex_assign(&c1_e_y, sf_mex_create("y", c1_cv5, 10, 0U, 1U, 0U, 2, 1,
            27), false);
          sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                            sf_mex_call_debug(sfGlobalDebugInstanceStruct,
            "message", 1U, 1U, 14, c1_e_y));
        }

        c1_ix = c1_j + c1_vstride;
        exitg1 = false;
        while ((!exitg1) && ((c1_vstride > 0) && (c1_ix <= c1_ixstop))) {
          c1_ixstart = c1_ix;
          if (!muDoubleScalarIsNaN(c1_b_varargin_1->data[c1_ix - 1])) {
            c1_mtmp = c1_b_varargin_1->data[c1_ix - 1];
            exitg1 = true;
          } else {
            c1_ix += c1_vstride;
          }
        }
      }

      if (c1_ixstart < c1_ixstop) {
        if (c1_vstride != 0) {
        } else {
          c1_f_y = NULL;
          sf_mex_assign(&c1_f_y, sf_mex_create("y", c1_cv6, 10, 0U, 1U, 0U, 2, 1,
            27), false);
          sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                            sf_mex_call_debug(sfGlobalDebugInstanceStruct,
            "message", 1U, 1U, 14, c1_f_y));
        }

        c1_b_ix = (c1_ixstart + c1_vstride) - 1;
        while ((c1_vstride > 0) && (c1_b_ix + 1 <= c1_ixstop)) {
          if (c1_b_varargin_1->data[c1_b_ix] < c1_mtmp) {
            c1_mtmp = c1_b_varargin_1->data[c1_b_ix];
          }

          c1_b_ix += c1_vstride;
        }
      }

      chartInstance->c1_m_data[c1_j - 1] = c1_mtmp;
    }

    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 142U);
    c1_i61 = c1_b_varargin_1->size[0] * c1_b_varargin_1->size[1];
    c1_b_varargin_1->size[0] = c1_Out->size[0];
    c1_b_varargin_1->size[1] = 22;
    c1_emxEnsureCapacity_real_T(chartInstance, c1_b_varargin_1, c1_i61,
      &c1_w_emlrtRTEI);
    c1_e_varargin_1 = c1_b_varargin_1->size[0];
    c1_f_varargin_1 = c1_b_varargin_1->size[1];
    c1_u_loop_ub = c1_Out->size[0] * c1_Out->size[1] - 1;
    for (c1_i62 = 0; c1_i62 <= c1_u_loop_ub; c1_i62++) {
      c1_b_varargin_1->data[c1_i62] = c1_Out->data[c1_i62];
    }

    c1_d_tmp_size[0] = c1_b_varargin_1->size[0];
    c1_b_vstride = c1_b_varargin_1->size[0];
    for (c1_b_j = 1; c1_b_j <= c1_b_vstride; c1_b_j++) {
      c1_b_ixstart = c1_b_j;
      c1_b_ixstop = c1_b_j + 21 * c1_b_vstride;
      c1_b_mtmp = c1_b_varargin_1->data[c1_b_j - 1];
      if (muDoubleScalarIsNaN(c1_b_mtmp)) {
        if (c1_b_vstride != 0) {
        } else {
          c1_g_y = NULL;
          sf_mex_assign(&c1_g_y, sf_mex_create("y", c1_cv7, 10, 0U, 1U, 0U, 2, 1,
            27), false);
          sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                            sf_mex_call_debug(sfGlobalDebugInstanceStruct,
            "message", 1U, 1U, 14, c1_g_y));
        }

        c1_c_ix = c1_b_j + c1_b_vstride;
        exitg1 = false;
        while ((!exitg1) && ((c1_b_vstride > 0) && (c1_c_ix <= c1_b_ixstop))) {
          c1_b_ixstart = c1_c_ix;
          if (!muDoubleScalarIsNaN(c1_b_varargin_1->data[c1_c_ix - 1])) {
            c1_b_mtmp = c1_b_varargin_1->data[c1_c_ix - 1];
            exitg1 = true;
          } else {
            c1_c_ix += c1_b_vstride;
          }
        }
      }

      if (c1_b_ixstart < c1_b_ixstop) {
        if (c1_b_vstride != 0) {
        } else {
          c1_h_y = NULL;
          sf_mex_assign(&c1_h_y, sf_mex_create("y", c1_cv8, 10, 0U, 1U, 0U, 2, 1,
            27), false);
          sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                            sf_mex_call_debug(sfGlobalDebugInstanceStruct,
            "message", 1U, 1U, 14, c1_h_y));
        }

        c1_d_ix = (c1_b_ixstart + c1_b_vstride) - 1;
        while ((c1_b_vstride > 0) && (c1_d_ix + 1 <= c1_b_ixstop)) {
          if (c1_b_varargin_1->data[c1_d_ix] > c1_b_mtmp) {
            c1_b_mtmp = c1_b_varargin_1->data[c1_d_ix];
          }

          c1_d_ix += c1_b_vstride;
        }
      }

      chartInstance->c1_tmp_data[c1_b_j - 1] = c1_b_mtmp;
    }

    c1_emxFree_real_T(chartInstance, &c1_b_varargin_1);
    _SFD_SIZE_EQ_CHECK_1D(c1_d_tmp_size[0], c1_m_size[0]);
    c1_range_size[0] = c1_d_tmp_size[0];
    c1_v_loop_ub = c1_d_tmp_size[0] - 1;
    for (c1_i63 = 0; c1_i63 <= c1_v_loop_ub; c1_i63++) {
      chartInstance->c1_range_data[c1_i63] = chartInstance->c1_tmp_data[c1_i63]
        - chartInstance->c1_m_data[c1_i63];
    }

    _SFD_SYMBOL_SWITCH(25U, 25U);
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 143U);
    _SFD_DIM_SIZE_EQ_CHECK(1, 22, 2);
    c1_d_tmp_size[0] = c1_Out->size[0];
    c1_w_loop_ub = c1_Out->size[0] - 1;
    for (c1_i64 = 0; c1_i64 <= c1_w_loop_ub; c1_i64++) {
      chartInstance->c1_tmp_data[c1_i64] = c1_Out->data[c1_i64];
    }

    _SFD_SIZE_EQ_CHECK_1D(c1_d_tmp_size[0], c1_m_size[0]);
    c1_e_tmp_size[0] = c1_d_tmp_size[0];
    c1_x_loop_ub = c1_d_tmp_size[0] - 1;
    for (c1_i65 = 0; c1_i65 <= c1_x_loop_ub; c1_i65++) {
      c1_e_tmp_data[c1_i65] = chartInstance->c1_tmp_data[c1_i65] -
        chartInstance->c1_m_data[c1_i65];
    }

    c1_c_range_size[0] = c1_range_size[0];
    c1_y_loop_ub = c1_range_size[0] - 1;
    for (c1_i66 = 0; c1_i66 <= c1_y_loop_ub; c1_i66++) {
      c1_c_range_data[c1_i66] = chartInstance->c1_range_data[c1_i66];
    }

    c1_b_rdivide(chartInstance, c1_e_tmp_data, c1_e_tmp_size, c1_c_range_data,
                 c1_c_range_size, c1_f_tmp_data, c1_f_tmp_size);
    c1_Out_size[0] = c1_f_tmp_size[0];
    c1_ab_loop_ub = c1_f_tmp_size[0] - 1;
    for (c1_i67 = 0; c1_i67 <= c1_ab_loop_ub; c1_i67++) {
      chartInstance->c1_Out_data[c1_i67] = c1_f_tmp_data[c1_i67];
    }

    _SFD_SYMBOL_SWITCH(23U, 38U);
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 145U);
    _SFD_SYMBOL_SWITCH(26U, 26U);
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 146U);
    c1_i68 = c1_Input->size[0] * c1_Input->size[1];
    c1_Input->size[0] = c1_In->size[0];
    c1_Input->size[1] = 4;
    c1_emxEnsureCapacity_real_T(chartInstance, c1_Input, c1_i68,
      &c1_cb_emlrtRTEI);
    c1_b_Input = c1_Input->size[0];
    c1_c_Input = c1_Input->size[1];
    c1_bb_loop_ub = c1_In->size[0] * c1_In->size[1] - 1;
    for (c1_i69 = 0; c1_i69 <= c1_bb_loop_ub; c1_i69++) {
      c1_Input->data[c1_i69] = c1_In->data[c1_i69];
    }

    c1_emxInit_real_T(chartInstance, &c1_g_varargin_1, 2, &c1_db_emlrtRTEI);
    _SFD_SYMBOL_SWITCH(26U, 39U);
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 147U);
    c1_i70 = c1_g_varargin_1->size[0] * c1_g_varargin_1->size[1];
    c1_g_varargin_1->size[0] = c1_Input->size[0];
    c1_g_varargin_1->size[1] = 4;
    c1_emxEnsureCapacity_real_T(chartInstance, c1_g_varargin_1, c1_i70,
      &c1_db_emlrtRTEI);
    c1_h_varargin_1 = c1_g_varargin_1->size[0];
    c1_i_varargin_1 = c1_g_varargin_1->size[1];
    c1_cb_loop_ub = c1_Input->size[0] * c1_Input->size[1] - 1;
    for (c1_i71 = 0; c1_i71 <= c1_cb_loop_ub; c1_i71++) {
      c1_g_varargin_1->data[c1_i71] = c1_Input->data[c1_i71];
    }

    c1_b1 = ((real_T)c1_g_varargin_1->size[0] != 1.0);
    if (c1_b1) {
    } else {
      c1_i_y = NULL;
      sf_mex_assign(&c1_i_y, sf_mex_create("y", c1_cv9, 10, 0U, 1U, 0U, 2, 1, 36),
                    false);
      sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                        sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message",
        1U, 1U, 14, c1_i_y));
    }

    if ((real_T)c1_g_varargin_1->size[0] > 0.0) {
    } else {
      c1_j_y = NULL;
      sf_mex_assign(&c1_j_y, sf_mex_create("y", c1_cv10, 10, 0U, 1U, 0U, 2, 1,
        39), false);
      sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                        sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message",
        1U, 1U, 14, c1_j_y));
    }

    c1_m11_size[0] = 1;
    c1_m11_size[1] = 4;
    c1_n = c1_g_varargin_1->size[0];
    for (c1_i = 0; c1_i < 4; c1_i++) {
      c1_e_ix = c1_i * c1_n;
      c1_c_ixstart = c1_i * c1_n + 2;
      c1_c_ixstop = c1_e_ix + c1_n;
      c1_c_mtmp = c1_g_varargin_1->data[c1_e_ix];
      if (c1_n > 1) {
        if (muDoubleScalarIsNaN(c1_c_mtmp)) {
          c1_f_ix = c1_c_ixstart;
          exitg1 = false;
          while ((!exitg1) && (c1_f_ix <= c1_c_ixstop)) {
            c1_c_ixstart = c1_f_ix + 1;
            if (!muDoubleScalarIsNaN(c1_g_varargin_1->data[c1_f_ix - 1])) {
              c1_c_mtmp = c1_g_varargin_1->data[c1_f_ix - 1];
              exitg1 = true;
            } else {
              c1_f_ix++;
            }
          }
        }

        if (c1_c_ixstart - 1 < c1_c_ixstop) {
          for (c1_g_ix = c1_c_ixstart - 1; c1_g_ix + 1 <= c1_c_ixstop; c1_g_ix++)
          {
            if (c1_g_varargin_1->data[c1_g_ix] < c1_c_mtmp) {
              c1_c_mtmp = c1_g_varargin_1->data[c1_g_ix];
            }
          }
        }
      }

      c1_m11_data[c1_i] = c1_c_mtmp;
    }

    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 148U);
    c1_i72 = c1_g_varargin_1->size[0] * c1_g_varargin_1->size[1];
    c1_g_varargin_1->size[0] = c1_Input->size[0];
    c1_g_varargin_1->size[1] = 4;
    c1_emxEnsureCapacity_real_T(chartInstance, c1_g_varargin_1, c1_i72,
      &c1_fb_emlrtRTEI);
    c1_j_varargin_1 = c1_g_varargin_1->size[0];
    c1_k_varargin_1 = c1_g_varargin_1->size[1];
    c1_db_loop_ub = c1_Input->size[0] * c1_Input->size[1] - 1;
    for (c1_i73 = 0; c1_i73 <= c1_db_loop_ub; c1_i73++) {
      c1_g_varargin_1->data[c1_i73] = c1_Input->data[c1_i73];
    }

    c1_b2 = ((real_T)c1_g_varargin_1->size[0] != 1.0);
    if (c1_b2) {
    } else {
      c1_k_y = NULL;
      sf_mex_assign(&c1_k_y, sf_mex_create("y", c1_cv11, 10, 0U, 1U, 0U, 2, 1,
        36), false);
      sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                        sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message",
        1U, 1U, 14, c1_k_y));
    }

    if ((real_T)c1_g_varargin_1->size[0] > 0.0) {
    } else {
      c1_l_y = NULL;
      sf_mex_assign(&c1_l_y, sf_mex_create("y", c1_cv12, 10, 0U, 1U, 0U, 2, 1,
        39), false);
      sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                        sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message",
        1U, 1U, 14, c1_l_y));
    }

    c1_g_tmp_size[0] = 1;
    c1_g_tmp_size[1] = 4;
    c1_b_n = c1_g_varargin_1->size[0];
    for (c1_b_i = 0; c1_b_i < 4; c1_b_i++) {
      c1_h_ix = c1_b_i * c1_b_n;
      c1_d_ixstart = c1_b_i * c1_b_n + 2;
      c1_d_ixstop = c1_h_ix + c1_b_n;
      c1_d_mtmp = c1_g_varargin_1->data[c1_h_ix];
      if (c1_b_n > 1) {
        if (muDoubleScalarIsNaN(c1_d_mtmp)) {
          c1_i_ix = c1_d_ixstart;
          exitg1 = false;
          while ((!exitg1) && (c1_i_ix <= c1_d_ixstop)) {
            c1_d_ixstart = c1_i_ix + 1;
            if (!muDoubleScalarIsNaN(c1_g_varargin_1->data[c1_i_ix - 1])) {
              c1_d_mtmp = c1_g_varargin_1->data[c1_i_ix - 1];
              exitg1 = true;
            } else {
              c1_i_ix++;
            }
          }
        }

        if (c1_d_ixstart - 1 < c1_d_ixstop) {
          for (c1_j_ix = c1_d_ixstart - 1; c1_j_ix + 1 <= c1_d_ixstop; c1_j_ix++)
          {
            if (c1_g_varargin_1->data[c1_j_ix] > c1_d_mtmp) {
              c1_d_mtmp = c1_g_varargin_1->data[c1_j_ix];
            }
          }
        }
      }

      c1_g_tmp_data[c1_b_i] = c1_d_mtmp;
    }

    c1_emxFree_real_T(chartInstance, &c1_g_varargin_1);
    c1_b_range_size[0] = 1;
    c1_b_range_size[1] = 4;
    c1_range = c1_b_range_size[0];
    c1_b_range = c1_b_range_size[1];
    for (c1_i74 = 0; c1_i74 < 4; c1_i74++) {
      c1_b_range_data[c1_i74] = c1_g_tmp_data[c1_i74] - c1_m11_data[c1_i74];
    }

    _SFD_SYMBOL_SWITCH(25U, 40U);
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 149U);
    for (c1_i75 = 0; c1_i75 < 2; c1_i75++) {
      c1_d_Input[c1_i75] = c1_Input->size[c1_i75];
    }

    for (c1_i76 = 0; c1_i76 < 2; c1_i76++) {
      c1_m11[c1_i76] = c1_m11_size[c1_i76];
    }

    _SFD_SIZE_EQ_CHECK_ND(c1_d_Input, c1_m11, 2);
    c1_b_Input_size[0] = c1_Input->size[0];
    c1_b_Input_size[1] = 4;
    c1_e_Input = c1_b_Input_size[0];
    c1_f_Input = c1_b_Input_size[1];
    c1_eb_loop_ub = c1_Input->size[0] * c1_Input->size[1] - 1;
    for (c1_i77 = 0; c1_i77 <= c1_eb_loop_ub; c1_i77++) {
      c1_b_Input_data[c1_i77] = c1_Input->data[c1_i77] - c1_m11_data[c1_i77];
    }

    c1_d_range_size[0] = 1;
    c1_d_range_size[1] = 4;
    c1_c_range = c1_d_range_size[0];
    c1_d_range = c1_d_range_size[1];
    for (c1_i78 = 0; c1_i78 < 4; c1_i78++) {
      c1_d_range_data[c1_i78] = c1_b_range_data[c1_i78];
    }

    c1_c_rdivide(chartInstance, c1_b_Input_data, c1_b_Input_size,
                 c1_d_range_data, c1_d_range_size, c1_h_tmp_data, c1_h_tmp_size);
    c1_Input_size[0] = 1;
    c1_Input_size[1] = 4;
    c1_g_Input = c1_Input_size[0];
    c1_h_Input = c1_Input_size[1];
    for (c1_i79 = 0; c1_i79 < 4; c1_i79++) {
      c1_Input_data[c1_i79] = c1_h_tmp_data[c1_i79];
    }

    _SFD_SYMBOL_SWITCH(26U, 41U);
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 152U);
    sf_mex_assign(&c1_net, sf_mex_call_debug(sfGlobalDebugInstanceStruct,
      "feedforwardnet", 1U, 1U, 14, c1_b_emlrt_marshallOut(chartInstance, 5.0)),
                  false);
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 153U);
    c1_c_Input_size[0] = 4;
    for (c1_i80 = 0; c1_i80 < 4; c1_i80++) {
      c1_c_Input_data[c1_i80] = c1_Input_data[c1_Input_size[0] * c1_i80];
    }

    c1_b_Out_size[0] = 1;
    c1_b_Out_size[1] = c1_Out_size[0];
    c1_fb_loop_ub = c1_Out_size[0] - 1;
    for (c1_i81 = 0; c1_i81 <= c1_fb_loop_ub; c1_i81++) {
      c1_b_Out_data[c1_b_Out_size[0] * c1_i81] = chartInstance->
        c1_Out_data[c1_i81];
    }

    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "train", 2U, 3U, 14,
                      sf_mex_dup(c1_net), 14, c1_j_emlrt_marshallOut
                      (chartInstance, c1_c_Input_data, c1_c_Input_size), 14,
                      c1_k_emlrt_marshallOut(chartInstance, c1_b_Out_data,
      c1_b_Out_size), &c1_b_net, &c1_b_tr);
    sf_mex_assign(&c1_net, sf_mex_dup(c1_b_net), false);
    sf_mex_assign(&c1_tr, sf_mex_dup(c1_b_tr), false);
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 156U);
    c1_i_tmp_size[0] = 4;
    c1_i_tmp_size[1] = c1_next_state_idx_size[0];
    c1_gb_loop_ub = c1_next_state_idx_size[0] - 1;
    for (c1_i82 = 0; c1_i82 <= c1_gb_loop_ub; c1_i82++) {
      for (c1_i83 = 0; c1_i83 < 4; c1_i83++) {
        c1_i_tmp_data[c1_i83 + c1_i_tmp_size[0] * c1_i82] = c1_get_states_space
          (chartInstance, (sf_eml_array_bounds_check(sfGlobalDebugInstanceStruct,
             chartInstance->S, 1U, 4815, 14, 17U, (int32_T)
             c1_next_state_idx_data[c1_i82], 1, 72) + 72 * c1_i83) - 1);
      }
    }

    sf_mex_assign(&c1_y_prime, sf_mex_call_debug(sfGlobalDebugInstanceStruct,
      "coder.internal.mxSubscript", 1U, 2U, 14, sf_mex_dup(c1_net), 14,
      c1_l_emlrt_marshallOut(chartInstance, c1_i_tmp_data, c1_i_tmp_size)),
                  false);
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 157U);
    c1_b_ii_size[0] = c1_next_state_idx_size[0];
    c1_hb_loop_ub = c1_next_state_idx_size[0] - 1;
    for (c1_i84 = 0; c1_i84 <= c1_hb_loop_ub; c1_i84++) {
      c1_b_ii_data[c1_i84] = c1__s32_d_(chartInstance,
        c1_next_state_idx_data[c1_i84], 1U, 4856U, 14U) - 1;
    }

    c1_n_emlrt_marshallIn(chartInstance, sf_mex_emlrtCoerceToClass
                          (sf_mex_call_debug(sfGlobalDebugInstanceStruct,
      "transpose", 1U, 1U, 14, sf_mex_dup(c1_y_prime)), "double"), "transpose",
                          c1_j_tmp_data, c1_j_tmp_size);
    for (c1_i85 = 0; c1_i85 < 12; c1_i85++) {
      c1_ib_loop_ub = c1_j_tmp_size[0] - 1;
      for (c1_i86 = 0; c1_i86 <= c1_ib_loop_ub; c1_i86++) {
        c1_set_Q_prime(chartInstance, c1_b_ii_data[c1_i86] + 72 * c1_i85,
                       c1_j_tmp_data[c1_i86 + c1_j_tmp_size[0] * c1_i85]);
      }
    }

    ssUpdateDataStoreLog_wrapper(chartInstance->S, 8);
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 171U);
    (real_T)sf_eml_array_bounds_check(sfGlobalDebugInstanceStruct,
      chartInstance->S, 1U, 5284, 1, MAX_uint32_T, 1, 1, c1_next_state_idx_size
      [0]);
    c1_b_next_state_idx[0] = c1_next_state_idx_size[0];
    c1_b_next_state_idx[1] = 1;
    c1_next_state_idx = c1_next_state_idx_data[0];
    _SFD_SYMBOL_SWITCH(14U, 42U);
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 172U);
    (real_T)sf_eml_array_bounds_check(sfGlobalDebugInstanceStruct,
      chartInstance->S, 1U, 5322, 1, MAX_uint32_T, 1, 1, c1_Final_Rew_size[0]);
    c1_b_Final_Rew[0] = c1_Final_Rew_size[0];
    c1_b_Final_Rew[1] = 1;
    c1_Final_Rew = c1_Final_Rew_data[0];
    _SFD_SYMBOL_SWITCH(19U, 43U);
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 173U);
    c1_c_next_state_idx = c1__s32_d_(chartInstance, c1_next_state_idx, 1U, 5371U,
      14U) - 1;
    for (c1_i87 = 0; c1_i87 < 12; c1_i87++) {
      c1_l_varargin_1[c1_i87] = c1_get_Q_prime(chartInstance,
        c1_c_next_state_idx + 72 * c1_i87);
    }

    c1_e_ixstart = 1;
    c1_e_mtmp = c1_l_varargin_1[0];
    if (muDoubleScalarIsNaN(c1_l_varargin_1[0])) {
      c1_k_ix = 1;
      exitg1 = false;
      while ((!exitg1) && (c1_k_ix + 1 < 13)) {
        c1_e_ixstart = c1_k_ix + 1;
        if (!muDoubleScalarIsNaN(c1_l_varargin_1[c1_k_ix])) {
          c1_e_mtmp = c1_l_varargin_1[c1_k_ix];
          exitg1 = true;
        } else {
          c1_k_ix++;
        }
      }
    }

    if (c1_e_ixstart < 12) {
      for (c1_l_ix = c1_e_ixstart; c1_l_ix + 1 < 13; c1_l_ix++) {
        if (c1_l_varargin_1[c1_l_ix] > c1_e_mtmp) {
          c1_e_mtmp = c1_l_varargin_1[c1_l_ix];
        }
      }
    }

    c1_c = c1_Final_Rew + c1_gamma * c1_e_mtmp;
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 174U);
    sf_mex_printf("%s =\\n", "c");
    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "disp", 0U, 1U, 14,
                      c1_b_emlrt_marshallOut(chartInstance, c1_c));
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 175U);
    c1_set_Q_targ(chartInstance, ((int16_T)sf_eml_array_bounds_check
      (sfGlobalDebugInstanceStruct, chartInstance->S, 1U, 5451, 9, 10U, (int32_T)
       c1_get_state_idx(chartInstance, 0), 1, 72) + 72 * (c1__s32_d_
      (chartInstance, c1_action_idx, 1U, 5461U, 10U) - 1)) - 1, c1_c);
    ssUpdateDataStoreLog_wrapper(chartInstance->S, 9);
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 177U);
    for (c1_i88 = 0; c1_i88 < 864; c1_i88++) {
      c1_set_target(chartInstance, c1_i88, c1_get_Q(chartInstance, c1_i88));
    }

    ssUpdateDataStoreLog_wrapper(chartInstance->S, 17);
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 178U);
    for (c1_i89 = 0; c1_i89 < 864; c1_i89++) {
      c1_Q_prev[c1_i89] = c1_get_Q(chartInstance, c1_i89);
    }

    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 180U);
    c1_set_target(chartInstance, ((int16_T)sf_eml_array_bounds_check
      (sfGlobalDebugInstanceStruct, chartInstance->S, 1U, 5557, 28, 18U,
       (int32_T)c1_get_state_idx(chartInstance, 0), 1, 72) + 72 * (c1__s32_d_
      (chartInstance, c1_action_idx, 1U, 5557U, 28U) - 1)) - 1, c1_get_Q_targ
                  (chartInstance, ((int16_T)sf_eml_array_bounds_check
      (sfGlobalDebugInstanceStruct, chartInstance->S, 1U, 5586, 28, 10U,
       (int32_T)c1_get_state_idx(chartInstance, 0), 1, 72) + 72 * (c1__s32_d_
      (chartInstance, c1_action_idx, 1U, 5586U, 28U) - 1)) - 1));
    ssUpdateDataStoreLog_wrapper(chartInstance->S, 17);
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 202U);
  for (c1_i38 = 0; c1_i38 < 864; c1_i38++) {
    c1_set_Q(chartInstance, c1_i38, c1_get_target(chartInstance, c1_i38));
  }

  ssUpdateDataStoreLog_wrapper(chartInstance->S, 7);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 206U);
  c1_set_ove_count(chartInstance, 0, c1_get_ove_count(chartInstance, 0) + 1.0);
  ssUpdateDataStoreLog_wrapper(chartInstance->S, 14);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 214U);
  c1_k = 1.0;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 239U);
  c1_d_y = 1.0;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 240U);
  c1_set_ove_count(chartInstance, 0, c1_get_ove_count(chartInstance, 0) + 1.0);
  ssUpdateDataStoreLog_wrapper(chartInstance->S, 14);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, -240);
  _SFD_SYMBOL_SCOPE_POP();
  sf_mex_destroy(&c1_net);
  sf_mex_destroy(&c1_tr);
  sf_mex_destroy(&c1_y_prime);
  sf_mex_destroy(&c1_b_net);
  sf_mex_destroy(&c1_b_tr);
  *chartInstance->c1_c_y = c1_d_y;
  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 0U, chartInstance->c1_sfEvent);
  c1_emxFree_real_T(chartInstance, &c1_Input);
  c1_emxFree_real_T(chartInstance, &c1_Out);
  c1_emxFree_real_T(chartInstance, &c1_Output);
  c1_emxFree_real_T(chartInstance, &c1_In);
}

static void initSimStructsc1_RM3_copy(SFc1_RM3_copyInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void init_script_number_translation(uint32_T c1_machineNumber, uint32_T
  c1_chartNumber, uint32_T c1_instanceNumber)
{
  (void)(c1_machineNumber);
  (void)(c1_chartNumber);
  (void)(c1_instanceNumber);
}

static const mxArray *c1_sf_marshallOut(void *chartInstanceVoid, void *c1_inData)
{
  const mxArray *c1_mxArrayOutData;
  real_T c1_b_u;
  const mxArray *c1_d_y = NULL;
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  c1_b_u = *(real_T *)c1_inData;
  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", &c1_b_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_d_y, false);
  return c1_mxArrayOutData;
}

static void c1_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_datasample;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_d_y;
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_datasample = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_d_y = c1_m_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_datasample),
    &c1_thisId);
  sf_mex_destroy(&c1_datasample);
  *(real_T *)c1_outData = c1_d_y;
  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_b_sf_marshallOut(void *chartInstanceVoid, real_T
  c1_inData_data[], int32_T c1_inData_size[2])
{
  const mxArray *c1_mxArrayOutData;
  int32_T c1_u_size[2];
  int32_T c1_b_u;
  int32_T c1_c_u;
  int32_T c1_i90;
  const mxArray *c1_d_y = NULL;
  real_T c1_u_data[4];
  SFc1_RM3_copyInstanceStruct *chartInstance;
  (void)c1_inData_size;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  c1_u_size[0] = 1;
  c1_u_size[1] = 4;
  c1_b_u = c1_u_size[0];
  c1_c_u = c1_u_size[1];
  for (c1_i90 = 0; c1_i90 < 4; c1_i90++) {
    c1_u_data[c1_i90] = c1_inData_data[c1_i90];
  }

  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", (void *)&c1_u_data, 0, 0U, 1U, 0U, 2,
    c1_u_size[0], c1_u_size[1]), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_d_y, false);
  return c1_mxArrayOutData;
}

static void c1_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T
  c1_y_data[], int32_T c1_y_size[2])
{
  int32_T c1_i91;
  int32_T c1_tmp_size[2];
  uint32_T c1_uv0[2];
  int32_T c1_i92;
  real_T c1_b_tmp_data[4];
  boolean_T c1_bv3[2];
  int32_T c1_d_y;
  int32_T c1_e_y;
  int32_T c1_i93;
  (void)chartInstance;
  for (c1_i91 = 0; c1_i91 < 2; c1_i91++) {
    c1_uv0[c1_i91] = 1U + 3U * (uint32_T)c1_i91;
  }

  c1_tmp_size[0] = sf_mex_get_dimension(c1_b_u, 0);
  c1_tmp_size[1] = sf_mex_get_dimension(c1_b_u, 1);
  for (c1_i92 = 0; c1_i92 < 2; c1_i92++) {
    c1_bv3[c1_i92] = false;
  }

  sf_mex_import_vs(c1_parentId, sf_mex_dup(c1_b_u), (void *)&c1_b_tmp_data, 1, 0,
                   0U, 1, 0U, 2, c1_bv3, c1_uv0, c1_tmp_size);
  c1_y_size[0] = 1;
  c1_y_size[1] = 4;
  c1_d_y = c1_y_size[0];
  c1_e_y = c1_y_size[1];
  for (c1_i93 = 0; c1_i93 < 4; c1_i93++) {
    c1_y_data[c1_i93] = c1_b_tmp_data[c1_i93];
  }

  sf_mex_destroy(&c1_b_u);
}

static void c1_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, real_T c1_outData_data[], int32_T
  c1_outData_size[2])
{
  const mxArray *c1_Input;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y_data[4];
  int32_T c1_y_size[2];
  int32_T c1_i94;
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_Input = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_Input), &c1_thisId, c1_y_data,
                      c1_y_size);
  sf_mex_destroy(&c1_Input);
  c1_outData_size[0] = 1;
  c1_outData_size[1] = 4;
  for (c1_i94 = 0; c1_i94 < 4; c1_i94++) {
    c1_outData_data[c1_outData_size[0] * c1_i94] = c1_y_data[c1_y_size[0] *
      c1_i94];
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_c_sf_marshallOut(void *chartInstanceVoid,
  c1_emxArray_real_T *c1_inData)
{
  const mxArray *c1_mxArrayOutData;
  c1_emxArray_real_T *c1_b_u;
  int32_T c1_i95;
  int32_T c1_c_u;
  int32_T c1_d_u;
  int32_T c1_loop_ub;
  int32_T c1_i96;
  const mxArray *c1_d_y = NULL;
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_emxInit_real_T(chartInstance, &c1_b_u, 2, (emlrtRTEInfo *)NULL);
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  c1_i95 = c1_b_u->size[0] * c1_b_u->size[1];
  c1_b_u->size[0] = c1_inData->size[0];
  c1_b_u->size[1] = 4;
  c1_emxEnsureCapacity_real_T(chartInstance, c1_b_u, c1_i95, (emlrtRTEInfo *)
    NULL);
  c1_c_u = c1_b_u->size[0];
  c1_d_u = c1_b_u->size[1];
  c1_loop_ub = c1_inData->size[0] * c1_inData->size[1] - 1;
  for (c1_i96 = 0; c1_i96 <= c1_loop_ub; c1_i96++) {
    c1_b_u->data[c1_i96] = c1_inData->data[c1_i96];
  }

  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", c1_b_u->data, 0, 0U, 1U, 0U, 2,
    c1_b_u->size[0], c1_b_u->size[1]), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_d_y, false);
  c1_emxFree_real_T(chartInstance, &c1_b_u);
  return c1_mxArrayOutData;
}

static void c1_b_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId,
  c1_emxArray_real_T *c1_d_y)
{
  c1_emxArray_real_T *c1_r0;
  int32_T c1_i97;
  int32_T c1_i98;
  uint32_T c1_uv1[2];
  int32_T c1_i99;
  boolean_T c1_bv4[2];
  static boolean_T c1_bv5[2] = { true, false };

  int32_T c1_i100;
  int32_T c1_e_y;
  int32_T c1_f_y;
  int32_T c1_loop_ub;
  int32_T c1_i101;
  c1_emxInit_real_T(chartInstance, &c1_r0, 2, (emlrtRTEInfo *)NULL);
  for (c1_i97 = 0; c1_i97 < 2; c1_i97++) {
    c1_uv1[c1_i97] = 5000U + (uint32_T)(-4996 * c1_i97);
  }

  c1_i98 = c1_r0->size[0] * c1_r0->size[1];
  c1_r0->size[0] = sf_mex_get_dimension(c1_b_u, 0);
  c1_r0->size[1] = sf_mex_get_dimension(c1_b_u, 1);
  c1_emxEnsureCapacity_real_T(chartInstance, c1_r0, c1_i98, (emlrtRTEInfo *)NULL);
  for (c1_i99 = 0; c1_i99 < 2; c1_i99++) {
    c1_bv4[c1_i99] = c1_bv5[c1_i99];
  }

  sf_mex_import_vs(c1_parentId, sf_mex_dup(c1_b_u), c1_r0->data, 1, 0, 0U, 1, 0U,
                   2, c1_bv4, c1_uv1, c1_r0->size);
  c1_i100 = c1_d_y->size[0] * c1_d_y->size[1];
  c1_d_y->size[0] = c1_r0->size[0];
  c1_d_y->size[1] = 4;
  c1_emxEnsureCapacity_real_T(chartInstance, c1_d_y, c1_i100, (emlrtRTEInfo *)
    NULL);
  c1_e_y = c1_d_y->size[0];
  c1_f_y = c1_d_y->size[1];
  c1_loop_ub = c1_r0->size[0] * c1_r0->size[1] - 1;
  for (c1_i101 = 0; c1_i101 <= c1_loop_ub; c1_i101++) {
    c1_d_y->data[c1_i101] = c1_r0->data[c1_i101];
  }

  sf_mex_destroy(&c1_b_u);
  c1_emxFree_real_T(chartInstance, &c1_r0);
}

static void c1_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, c1_emxArray_real_T *c1_outData)
{
  c1_emxArray_real_T *c1_d_y;
  const mxArray *c1_Input;
  emlrtMsgIdentifier c1_thisId;
  int32_T c1_i102;
  int32_T c1_i103;
  int32_T c1_loop_ub;
  int32_T c1_i104;
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_emxInit_real_T(chartInstance, &c1_d_y, 2, (emlrtRTEInfo *)NULL);
  c1_Input = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_Input), &c1_thisId, c1_d_y);
  sf_mex_destroy(&c1_Input);
  c1_i102 = c1_outData->size[0] * c1_outData->size[1];
  c1_outData->size[0] = c1_d_y->size[0];
  c1_outData->size[1] = 4;
  c1_emxEnsureCapacity_real_T(chartInstance, c1_outData, c1_i102, (emlrtRTEInfo *)
    NULL);
  for (c1_i103 = 0; c1_i103 < 4; c1_i103++) {
    c1_loop_ub = c1_d_y->size[0] - 1;
    for (c1_i104 = 0; c1_i104 <= c1_loop_ub; c1_i104++) {
      c1_outData->data[c1_i104 + c1_outData->size[0] * c1_i103] = c1_d_y->
        data[c1_i104 + c1_d_y->size[0] * c1_i103];
    }
  }

  c1_emxFree_real_T(chartInstance, &c1_d_y);
  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_d_sf_marshallOut(void *chartInstanceVoid, real_T
  c1_inData_data[], int32_T c1_inData_size[1])
{
  const mxArray *c1_mxArrayOutData;
  int32_T c1_u_size[1];
  int32_T c1_loop_ub;
  int32_T c1_i105;
  const mxArray *c1_d_y = NULL;
  real_T c1_u_data[5000];
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  c1_u_size[0] = c1_inData_size[0];
  c1_loop_ub = c1_inData_size[0] - 1;
  for (c1_i105 = 0; c1_i105 <= c1_loop_ub; c1_i105++) {
    c1_u_data[c1_i105] = c1_inData_data[c1_i105];
  }

  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", (void *)&c1_u_data, 0, 0U, 1U, 0U, 1,
    c1_u_size[0]), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_d_y, false);
  return c1_mxArrayOutData;
}

static void c1_c_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T
  c1_y_data[], int32_T c1_y_size[1])
{
  uint32_T c1_uv2[1];
  int32_T c1_tmp_size[1];
  boolean_T c1_bv6[1];
  real_T c1_b_tmp_data[5000];
  int32_T c1_loop_ub;
  int32_T c1_i106;
  (void)chartInstance;
  c1_uv2[0] = 5000U;
  c1_tmp_size[0] = sf_mex_get_dimension(c1_b_u, 0);
  c1_bv6[0] = true;
  sf_mex_import_vs(c1_parentId, sf_mex_dup(c1_b_u), (void *)&c1_b_tmp_data, 1, 0,
                   0U, 1, 0U, 1, c1_bv6, c1_uv2, c1_tmp_size);
  c1_y_size[0] = c1_tmp_size[0];
  c1_loop_ub = c1_tmp_size[0] - 1;
  for (c1_i106 = 0; c1_i106 <= c1_loop_ub; c1_i106++) {
    c1_y_data[c1_i106] = c1_b_tmp_data[c1_i106];
  }

  sf_mex_destroy(&c1_b_u);
}

static void c1_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, real_T c1_outData_data[], int32_T
  c1_outData_size[1])
{
  const mxArray *c1_Out;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y_data[5000];
  int32_T c1_y_size[1];
  int32_T c1_loop_ub;
  int32_T c1_i107;
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_Out = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_Out), &c1_thisId, c1_y_data,
                        c1_y_size);
  sf_mex_destroy(&c1_Out);
  c1_outData_size[0] = c1_y_size[0];
  c1_loop_ub = c1_y_size[0] - 1;
  for (c1_i107 = 0; c1_i107 <= c1_loop_ub; c1_i107++) {
    c1_outData_data[c1_i107] = c1_y_data[c1_i107];
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_e_sf_marshallOut(void *chartInstanceVoid,
  c1_emxArray_real_T *c1_inData)
{
  const mxArray *c1_mxArrayOutData;
  c1_emxArray_real_T *c1_b_u;
  int32_T c1_i108;
  int32_T c1_c_u;
  int32_T c1_d_u;
  int32_T c1_loop_ub;
  int32_T c1_i109;
  const mxArray *c1_d_y = NULL;
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_emxInit_real_T(chartInstance, &c1_b_u, 2, (emlrtRTEInfo *)NULL);
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  c1_i108 = c1_b_u->size[0] * c1_b_u->size[1];
  c1_b_u->size[0] = c1_inData->size[0];
  c1_b_u->size[1] = 22;
  c1_emxEnsureCapacity_real_T(chartInstance, c1_b_u, c1_i108, (emlrtRTEInfo *)
    NULL);
  c1_c_u = c1_b_u->size[0];
  c1_d_u = c1_b_u->size[1];
  c1_loop_ub = c1_inData->size[0] * c1_inData->size[1] - 1;
  for (c1_i109 = 0; c1_i109 <= c1_loop_ub; c1_i109++) {
    c1_b_u->data[c1_i109] = c1_inData->data[c1_i109];
  }

  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", c1_b_u->data, 0, 0U, 1U, 0U, 2,
    c1_b_u->size[0], c1_b_u->size[1]), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_d_y, false);
  c1_emxFree_real_T(chartInstance, &c1_b_u);
  return c1_mxArrayOutData;
}

static void c1_d_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId,
  c1_emxArray_real_T *c1_d_y)
{
  c1_emxArray_real_T *c1_r1;
  int32_T c1_i110;
  int32_T c1_i111;
  uint32_T c1_uv3[2];
  int32_T c1_i112;
  boolean_T c1_bv7[2];
  static boolean_T c1_bv8[2] = { true, false };

  int32_T c1_i113;
  int32_T c1_e_y;
  int32_T c1_f_y;
  int32_T c1_loop_ub;
  int32_T c1_i114;
  c1_emxInit_real_T(chartInstance, &c1_r1, 2, (emlrtRTEInfo *)NULL);
  for (c1_i110 = 0; c1_i110 < 2; c1_i110++) {
    c1_uv3[c1_i110] = 5000U + (uint32_T)(-4978 * c1_i110);
  }

  c1_i111 = c1_r1->size[0] * c1_r1->size[1];
  c1_r1->size[0] = sf_mex_get_dimension(c1_b_u, 0);
  c1_r1->size[1] = sf_mex_get_dimension(c1_b_u, 1);
  c1_emxEnsureCapacity_real_T(chartInstance, c1_r1, c1_i111, (emlrtRTEInfo *)
    NULL);
  for (c1_i112 = 0; c1_i112 < 2; c1_i112++) {
    c1_bv7[c1_i112] = c1_bv8[c1_i112];
  }

  sf_mex_import_vs(c1_parentId, sf_mex_dup(c1_b_u), c1_r1->data, 1, 0, 0U, 1, 0U,
                   2, c1_bv7, c1_uv3, c1_r1->size);
  c1_i113 = c1_d_y->size[0] * c1_d_y->size[1];
  c1_d_y->size[0] = c1_r1->size[0];
  c1_d_y->size[1] = 22;
  c1_emxEnsureCapacity_real_T(chartInstance, c1_d_y, c1_i113, (emlrtRTEInfo *)
    NULL);
  c1_e_y = c1_d_y->size[0];
  c1_f_y = c1_d_y->size[1];
  c1_loop_ub = c1_r1->size[0] * c1_r1->size[1] - 1;
  for (c1_i114 = 0; c1_i114 <= c1_loop_ub; c1_i114++) {
    c1_d_y->data[c1_i114] = c1_r1->data[c1_i114];
  }

  sf_mex_destroy(&c1_b_u);
  c1_emxFree_real_T(chartInstance, &c1_r1);
}

static void c1_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, c1_emxArray_real_T *c1_outData)
{
  c1_emxArray_real_T *c1_d_y;
  const mxArray *c1_Out;
  emlrtMsgIdentifier c1_thisId;
  int32_T c1_i115;
  int32_T c1_i116;
  int32_T c1_loop_ub;
  int32_T c1_i117;
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_emxInit_real_T(chartInstance, &c1_d_y, 2, (emlrtRTEInfo *)NULL);
  c1_Out = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_Out), &c1_thisId, c1_d_y);
  sf_mex_destroy(&c1_Out);
  c1_i115 = c1_outData->size[0] * c1_outData->size[1];
  c1_outData->size[0] = c1_d_y->size[0];
  c1_outData->size[1] = 22;
  c1_emxEnsureCapacity_real_T(chartInstance, c1_outData, c1_i115, (emlrtRTEInfo *)
    NULL);
  for (c1_i116 = 0; c1_i116 < 22; c1_i116++) {
    c1_loop_ub = c1_d_y->size[0] - 1;
    for (c1_i117 = 0; c1_i117 <= c1_loop_ub; c1_i117++) {
      c1_outData->data[c1_i117 + c1_outData->size[0] * c1_i116] = c1_d_y->
        data[c1_i117 + c1_d_y->size[0] * c1_i116];
    }
  }

  c1_emxFree_real_T(chartInstance, &c1_d_y);
  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_f_sf_marshallOut(void *chartInstanceVoid, real_T
  c1_inData_data[], int32_T c1_inData_size[1])
{
  const mxArray *c1_mxArrayOutData;
  int32_T c1_u_size[1];
  int32_T c1_loop_ub;
  int32_T c1_i118;
  const mxArray *c1_d_y = NULL;
  real_T c1_u_data[72];
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  c1_u_size[0] = c1_inData_size[0];
  c1_loop_ub = c1_inData_size[0] - 1;
  for (c1_i118 = 0; c1_i118 <= c1_loop_ub; c1_i118++) {
    c1_u_data[c1_i118] = c1_inData_data[c1_i118];
  }

  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", (void *)&c1_u_data, 0, 0U, 1U, 0U, 1,
    c1_u_size[0]), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_d_y, false);
  return c1_mxArrayOutData;
}

static void c1_e_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T
  c1_y_data[], int32_T c1_y_size[1])
{
  uint32_T c1_uv4[1];
  int32_T c1_tmp_size[1];
  boolean_T c1_bv9[1];
  real_T c1_b_tmp_data[72];
  int32_T c1_loop_ub;
  int32_T c1_i119;
  (void)chartInstance;
  c1_uv4[0] = 72U;
  c1_tmp_size[0] = sf_mex_get_dimension(c1_b_u, 0);
  c1_bv9[0] = true;
  sf_mex_import_vs(c1_parentId, sf_mex_dup(c1_b_u), (void *)&c1_b_tmp_data, 1, 0,
                   0U, 1, 0U, 1, c1_bv9, c1_uv4, c1_tmp_size);
  c1_y_size[0] = c1_tmp_size[0];
  c1_loop_ub = c1_tmp_size[0] - 1;
  for (c1_i119 = 0; c1_i119 <= c1_loop_ub; c1_i119++) {
    c1_y_data[c1_i119] = c1_b_tmp_data[c1_i119];
  }

  sf_mex_destroy(&c1_b_u);
}

static void c1_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, real_T c1_outData_data[], int32_T
  c1_outData_size[1])
{
  const mxArray *c1_Height;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y_data[72];
  int32_T c1_y_size[1];
  int32_T c1_loop_ub;
  int32_T c1_i120;
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_Height = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_Height), &c1_thisId,
                        c1_y_data, c1_y_size);
  sf_mex_destroy(&c1_Height);
  c1_outData_size[0] = c1_y_size[0];
  c1_loop_ub = c1_y_size[0] - 1;
  for (c1_i120 = 0; c1_i120 <= c1_loop_ub; c1_i120++) {
    c1_outData_data[c1_i120] = c1_y_data[c1_i120];
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_g_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData;
  int32_T c1_i121;
  int32_T c1_i122;
  const mxArray *c1_d_y = NULL;
  int32_T c1_i123;
  real_T c1_b_u[864];
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  c1_i121 = 0;
  for (c1_i122 = 0; c1_i122 < 12; c1_i122++) {
    for (c1_i123 = 0; c1_i123 < 72; c1_i123++) {
      c1_b_u[c1_i123 + c1_i121] = (*(real_T (*)[864])c1_inData)[c1_i123 +
        c1_i121];
    }

    c1_i121 += 72;
  }

  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", c1_b_u, 0, 0U, 1U, 0U, 2, 72, 12),
                false);
  sf_mex_assign(&c1_mxArrayOutData, c1_d_y, false);
  return c1_mxArrayOutData;
}

static void c1_f_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_d_y
  [864])
{
  real_T c1_dv5[864];
  int32_T c1_i124;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_b_u), c1_dv5, 1, 0, 0U, 1, 0U, 2, 72,
                12);
  for (c1_i124 = 0; c1_i124 < 864; c1_i124++) {
    c1_d_y[c1_i124] = c1_dv5[c1_i124];
  }

  sf_mex_destroy(&c1_b_u);
}

static void c1_g_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_Q_prev;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_d_y[864];
  int32_T c1_i125;
  int32_T c1_i126;
  int32_T c1_i127;
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_Q_prev = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_Q_prev), &c1_thisId, c1_d_y);
  sf_mex_destroy(&c1_Q_prev);
  c1_i125 = 0;
  for (c1_i126 = 0; c1_i126 < 12; c1_i126++) {
    for (c1_i127 = 0; c1_i127 < 72; c1_i127++) {
      (*(real_T (*)[864])c1_outData)[c1_i127 + c1_i125] = c1_d_y[c1_i127 +
        c1_i125];
    }

    c1_i125 += 72;
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_h_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  const mxArray *c1_b_u;
  const mxArray *c1_d_y = NULL;
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_b_u = sf_mex_dup(*(const mxArray **)c1_inData);
  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_duplicatearraysafe(&c1_b_u), false);
  sf_mex_destroy(&c1_b_u);
  sf_mex_assign(&c1_mxArrayOutData, c1_d_y, false);
  return c1_mxArrayOutData;
}

static const mxArray *c1_i_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  const mxArray *c1_d_y = NULL;
  SFc1_RM3_copyInstanceStruct *chartInstance;
  (void)c1_inData;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_d_y, false);
  return c1_mxArrayOutData;
}

static void c1_g_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId)
{
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_b_u), NULL, 1, 0, 0U, 1, 0U, 2, 0, 0);
  sf_mex_destroy(&c1_b_u);
}

static void c1_h_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_Input;
  emlrtMsgIdentifier c1_thisId;
  SFc1_RM3_copyInstanceStruct *chartInstance;
  (void)c1_outData;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_Input = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_g_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_Input), &c1_thisId);
  sf_mex_destroy(&c1_Input);
  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_j_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData;
  int32_T c1_i128;
  const mxArray *c1_d_y = NULL;
  real_T c1_b_u[2];
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  for (c1_i128 = 0; c1_i128 < 2; c1_i128++) {
    c1_b_u[c1_i128] = (*(real_T (*)[2])c1_inData)[c1_i128];
  }

  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", c1_b_u, 0, 0U, 1U, 0U, 2, 1, 2),
                false);
  sf_mex_assign(&c1_mxArrayOutData, c1_d_y, false);
  return c1_mxArrayOutData;
}

static void c1_i_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_datasample;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_d_y[2];
  int32_T c1_i129;
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_datasample = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_k_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_datasample), &c1_thisId,
                        c1_d_y);
  sf_mex_destroy(&c1_datasample);
  for (c1_i129 = 0; c1_i129 < 2; c1_i129++) {
    (*(real_T (*)[2])c1_outData)[c1_i129] = c1_d_y[c1_i129];
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_k_sf_marshallOut(void *chartInstanceVoid, real_T
  c1_inData_data[], int32_T c1_inData_size[1])
{
  const mxArray *c1_mxArrayOutData;
  int32_T c1_u_size[1];
  int32_T c1_loop_ub;
  int32_T c1_i130;
  const mxArray *c1_d_y = NULL;
  real_T c1_u_data[12];
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  c1_u_size[0] = c1_inData_size[0];
  c1_loop_ub = c1_inData_size[0] - 1;
  for (c1_i130 = 0; c1_i130 <= c1_loop_ub; c1_i130++) {
    c1_u_data[c1_i130] = c1_inData_data[c1_i130];
  }

  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", (void *)&c1_u_data, 0, 0U, 1U, 0U, 1,
    c1_u_size[0]), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_d_y, false);
  return c1_mxArrayOutData;
}

static void c1_h_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T
  c1_y_data[], int32_T c1_y_size[1])
{
  uint32_T c1_uv5[1];
  int32_T c1_tmp_size[1];
  boolean_T c1_bv10[1];
  real_T c1_b_tmp_data[12];
  int32_T c1_loop_ub;
  int32_T c1_i131;
  (void)chartInstance;
  c1_uv5[0] = 12U;
  c1_tmp_size[0] = sf_mex_get_dimension(c1_b_u, 0);
  c1_bv10[0] = true;
  sf_mex_import_vs(c1_parentId, sf_mex_dup(c1_b_u), (void *)&c1_b_tmp_data, 1, 0,
                   0U, 1, 0U, 1, c1_bv10, c1_uv5, c1_tmp_size);
  c1_y_size[0] = c1_tmp_size[0];
  c1_loop_ub = c1_tmp_size[0] - 1;
  for (c1_i131 = 0; c1_i131 <= c1_loop_ub; c1_i131++) {
    c1_y_data[c1_i131] = c1_b_tmp_data[c1_i131];
  }

  sf_mex_destroy(&c1_b_u);
}

static void c1_j_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, real_T c1_outData_data[], int32_T
  c1_outData_size[1])
{
  const mxArray *c1_action_idx_pr;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y_data[12];
  int32_T c1_y_size[1];
  int32_T c1_loop_ub;
  int32_T c1_i132;
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_action_idx_pr = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_h_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_action_idx_pr), &c1_thisId,
                        c1_y_data, c1_y_size);
  sf_mex_destroy(&c1_action_idx_pr);
  c1_outData_size[0] = c1_y_size[0];
  c1_loop_ub = c1_y_size[0] - 1;
  for (c1_i132 = 0; c1_i132 <= c1_loop_ub; c1_i132++) {
    c1_outData_data[c1_i132] = c1_y_data[c1_i132];
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_l_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData;
  c1_sI3an6DNyfQAOCkXY61B1nE c1_b_u;
  const mxArray *c1_d_y = NULL;
  real_T c1_c_u;
  const mxArray *c1_e_y = NULL;
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  c1_b_u = *(c1_sI3an6DNyfQAOCkXY61B1nE *)c1_inData;
  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_createstruct("structure", 2, 1, 1), false);
  c1_c_u = c1_b_u.StepRatio;
  c1_e_y = NULL;
  sf_mex_assign(&c1_e_y, sf_mex_create("y", &c1_c_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c1_d_y, c1_e_y, "StepRatio", "StepRatio", 0);
  sf_mex_assign(&c1_mxArrayOutData, c1_d_y, false);
  return c1_mxArrayOutData;
}

static c1_sI3an6DNyfQAOCkXY61B1nE c1_i_emlrt_marshallIn
  (SFc1_RM3_copyInstanceStruct *chartInstance, const mxArray *c1_b_u, const
   emlrtMsgIdentifier *c1_parentId)
{
  c1_sI3an6DNyfQAOCkXY61B1nE c1_d_y;
  emlrtMsgIdentifier c1_thisId;
  static const char * c1_fieldNames[1] = { "StepRatio" };

  c1_thisId.fParent = c1_parentId;
  c1_thisId.bParentIsCell = false;
  sf_mex_check_struct(c1_parentId, c1_b_u, 1, c1_fieldNames, 0U, NULL);
  c1_thisId.fIdentifier = "StepRatio";
  c1_d_y.StepRatio = c1_m_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getfield(c1_b_u, "StepRatio", "StepRatio", 0)), &c1_thisId);
  sf_mex_destroy(&c1_b_u);
  return c1_d_y;
}

static void c1_k_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_opts;
  emlrtMsgIdentifier c1_thisId;
  c1_sI3an6DNyfQAOCkXY61B1nE c1_d_y;
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_opts = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_d_y = c1_i_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_opts), &c1_thisId);
  sf_mex_destroy(&c1_opts);
  *(c1_sI3an6DNyfQAOCkXY61B1nE *)c1_outData = c1_d_y;
  sf_mex_destroy(&c1_mxArrayInData);
}

const mxArray *sf_c1_RM3_copy_get_eml_resolved_functions_info(void)
{
  const mxArray *c1_nameCaptureInfo = NULL;
  c1_nameCaptureInfo = NULL;
  sf_mex_assign(&c1_nameCaptureInfo, sf_mex_create("nameCaptureInfo", NULL, 0,
    0U, 1U, 0U, 2, 0, 1), false);
  return c1_nameCaptureInfo;
}

static real_T c1_sum(SFc1_RM3_copyInstanceStruct *chartInstance, real_T c1_x[50])
{
  real_T c1_d_y;
  int32_T c1_k;
  (void)chartInstance;
  c1_d_y = c1_x[0];
  for (c1_k = 0; c1_k < 49; c1_k++) {
    c1_d_y += c1_x[c1_k + 1];
  }

  return c1_d_y;
}

static int16_T c1_rdivide(SFc1_RM3_copyInstanceStruct *chartInstance, real_T
  c1_x, int16_T c1_d_y)
{
  real_T c1_d0;
  int16_T c1_i133;
  (void)chartInstance;
  c1_d0 = muDoubleScalarRound(c1_x / (real_T)c1_d_y);
  if (c1_d0 < 32768.0) {
    if (c1_d0 >= -32768.0) {
      c1_i133 = (int16_T)c1_d0;
    } else {
      c1_i133 = MIN_int16_T;
      _SFD_OVERFLOW_DETECTION(SFDB_SATURATE, 1U, 0U, 0U);
    }
  } else if (c1_d0 >= 32768.0) {
    c1_i133 = MAX_int16_T;
    _SFD_OVERFLOW_DETECTION(SFDB_SATURATE, 1U, 0U, 0U);
  } else {
    c1_i133 = 0;
  }

  return c1_i133;
}

static void c1_mpower(SFc1_RM3_copyInstanceStruct *chartInstance, real_T
                      c1_a_data[], int32_T c1_a_size[1], real_T c1_c_data[],
                      int32_T c1_c_size[1])
{
  boolean_T c1_b3;
  const mxArray *c1_d_y = NULL;
  static char_T c1_cv13[19] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T', 'L',
    'A', 'B', ':', 's', 'q', 'u', 'a', 'r', 'e' };

  int32_T c1_b_a_size[1];
  int32_T c1_loop_ub;
  int32_T c1_i134;
  real_T c1_b_a_data[72];
  c1_b3 = ((real_T)c1_a_size[0] == 1.0);
  if (c1_b3) {
  } else {
    c1_d_y = NULL;
    sf_mex_assign(&c1_d_y, sf_mex_create("y", c1_cv13, 10, 0U, 1U, 0U, 2, 1, 19),
                  false);
    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                      sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message",
      1U, 1U, 14, c1_d_y));
  }

  c1_b_a_size[0] = c1_a_size[0];
  c1_loop_ub = c1_a_size[0] - 1;
  for (c1_i134 = 0; c1_i134 <= c1_loop_ub; c1_i134++) {
    c1_b_a_data[c1_i134] = c1_a_data[c1_i134];
  }

  c1_matrix_to_scalar_power(chartInstance, c1_b_a_data, c1_b_a_size, c1_c_data,
    c1_c_size);
}

static void c1_matrix_to_scalar_power(SFc1_RM3_copyInstanceStruct *chartInstance,
  real_T c1_a_data[], int32_T c1_a_size[1], real_T c1_c_data[], int32_T
  c1_c_size[1])
{
  real_T c1_a;
  const mxArray *c1_d_y = NULL;
  const mxArray *c1_e_y = NULL;
  int32_T c1_loop_ub;
  static char_T c1_cv14[21] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T', 'L',
    'A', 'B', ':', 'i', 'n', 'n', 'e', 'r', 'd', 'i', 'm' };

  static char_T c1_cv15[45] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'm', 't', 'i', 'm', 'e', 's', '_', 'n', 'o', 'D', 'y',
    'n', 'a', 'm', 'i', 'c', 'S', 'c', 'a', 'l', 'a', 'r', 'E', 'x', 'p', 'a',
    'n', 's', 'i', 'o', 'n' };

  int32_T c1_i135;
  real_T c1_c;
  const mxArray *c1_f_y = NULL;
  const mxArray *c1_g_y = NULL;
  int32_T c1_b_loop_ub;
  static char_T c1_cv16[21] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T', 'L',
    'A', 'B', ':', 'i', 'n', 'n', 'e', 'r', 'd', 'i', 'm' };

  static char_T c1_cv17[45] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'm', 't', 'i', 'm', 'e', 's', '_', 'n', 'o', 'D', 'y',
    'n', 'a', 'm', 'i', 'c', 'S', 'c', 'a', 'l', 'a', 'r', 'E', 'x', 'p', 'a',
    'n', 's', 'i', 'o', 'n' };

  int32_T c1_i136;
  real_T c1_b_c;
  const mxArray *c1_h_y = NULL;
  const mxArray *c1_i_y = NULL;
  int32_T c1_c_loop_ub;
  static char_T c1_cv18[21] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T', 'L',
    'A', 'B', ':', 'i', 'n', 'n', 'e', 'r', 'd', 'i', 'm' };

  static char_T c1_cv19[45] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'm', 't', 'i', 'm', 'e', 's', '_', 'n', 'o', 'D', 'y',
    'n', 'a', 'm', 'i', 'c', 'S', 'c', 'a', 'l', 'a', 'r', 'E', 'x', 'p', 'a',
    'n', 's', 'i', 'o', 'n' };

  int32_T c1_i137;
  (void)chartInstance;
  if (!(1.0 == (real_T)c1_a_size[0])) {
    if ((c1_a_size[0] == 1) || (c1_a_size[0] == 1)) {
      c1_e_y = NULL;
      sf_mex_assign(&c1_e_y, sf_mex_create("y", c1_cv15, 10, 0U, 1U, 0U, 2, 1,
        45), false);
      sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                        sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message",
        1U, 1U, 14, c1_e_y));
    } else {
      c1_d_y = NULL;
      sf_mex_assign(&c1_d_y, sf_mex_create("y", c1_cv14, 10, 0U, 1U, 0U, 2, 1,
        21), false);
      sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                        sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message",
        1U, 1U, 14, c1_d_y));
    }
  }

  c1_a = c1_a_data[0];
  c1_c_size[0] = c1_a_size[0];
  c1_loop_ub = c1_a_size[0] - 1;
  for (c1_i135 = 0; c1_i135 <= c1_loop_ub; c1_i135++) {
    c1_c_data[c1_i135] = c1_a_data[c1_i135] * c1_a;
  }

  if (!(1.0 == (real_T)c1_c_size[0])) {
    if ((c1_c_size[0] == 1) || (c1_c_size[0] == 1)) {
      c1_g_y = NULL;
      sf_mex_assign(&c1_g_y, sf_mex_create("y", c1_cv17, 10, 0U, 1U, 0U, 2, 1,
        45), false);
      sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                        sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message",
        1U, 1U, 14, c1_g_y));
    } else {
      c1_f_y = NULL;
      sf_mex_assign(&c1_f_y, sf_mex_create("y", c1_cv16, 10, 0U, 1U, 0U, 2, 1,
        21), false);
      sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                        sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message",
        1U, 1U, 14, c1_f_y));
    }
  }

  c1_c = c1_c_data[0];
  c1_c_size[0];
  c1_b_loop_ub = c1_c_size[0] - 1;
  for (c1_i136 = 0; c1_i136 <= c1_b_loop_ub; c1_i136++) {
    c1_c_data[c1_i136] *= c1_c;
  }

  if (!(1.0 == (real_T)c1_c_size[0])) {
    if ((c1_a_size[0] == 1) || (c1_c_size[0] == 1)) {
      c1_i_y = NULL;
      sf_mex_assign(&c1_i_y, sf_mex_create("y", c1_cv19, 10, 0U, 1U, 0U, 2, 1,
        45), false);
      sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                        sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message",
        1U, 1U, 14, c1_i_y));
    } else {
      c1_h_y = NULL;
      sf_mex_assign(&c1_h_y, sf_mex_create("y", c1_cv18, 10, 0U, 1U, 0U, 2, 1,
        21), false);
      sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                        sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message",
        1U, 1U, 14, c1_h_y));
    }
  }

  c1_b_c = c1_c_data[0];
  c1_c_size[0] = c1_a_size[0];
  c1_c_loop_ub = c1_a_size[0] - 1;
  for (c1_i137 = 0; c1_i137 <= c1_c_loop_ub; c1_i137++) {
    c1_c_data[c1_i137] = c1_a_data[c1_i137] * c1_b_c;
  }
}

static void c1_b_rdivide(SFc1_RM3_copyInstanceStruct *chartInstance, real_T
  c1_x_data[], int32_T c1_x_size[1], real_T c1_y_data[], int32_T c1_y_size[1],
  real_T c1_z_data[], int32_T c1_z_size[1])
{
  real_T c1_a[2];
  real_T c1_b[2];
  boolean_T c1_p;
  boolean_T c1_b_p;
  int32_T c1_k;
  boolean_T c1_c_p;
  const mxArray *c1_d_y = NULL;
  int32_T c1_loop_ub;
  int32_T c1_i138;
  boolean_T exitg1;
  (void)chartInstance;
  c1_a[0] = (real_T)c1_x_size[0];
  c1_a[1] = 1.0;
  c1_b[0] = (real_T)c1_y_size[0];
  c1_b[1] = 1.0;
  c1_p = false;
  c1_b_p = true;
  c1_k = 0;
  exitg1 = false;
  while ((!exitg1) && (c1_k < 2)) {
    if (!(c1_a[c1_k] == c1_b[c1_k])) {
      c1_b_p = false;
      exitg1 = true;
    } else {
      c1_k++;
    }
  }

  if (!c1_b_p) {
  } else {
    c1_p = true;
  }

  c1_c_p = c1_p;
  if (c1_c_p) {
  } else {
    c1_d_y = NULL;
    sf_mex_assign(&c1_d_y, sf_mex_create("y", c1_cv0, 10, 0U, 1U, 0U, 2, 1, 15),
                  false);
    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                      sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message",
      1U, 1U, 14, c1_d_y));
  }

  c1_z_size[0] = c1_x_size[0];
  c1_loop_ub = c1_x_size[0] - 1;
  for (c1_i138 = 0; c1_i138 <= c1_loop_ub; c1_i138++) {
    c1_z_data[c1_i138] = c1_x_data[c1_i138] / c1_y_data[c1_i138];
  }
}

static void c1_c_rdivide(SFc1_RM3_copyInstanceStruct *chartInstance, real_T
  c1_x_data[], int32_T c1_x_size[2], real_T c1_y_data[], int32_T c1_y_size[2],
  real_T c1_z_data[], int32_T c1_z_size[2])
{
  int32_T c1_i139;
  int32_T c1_i140;
  real_T c1_a[2];
  boolean_T c1_p;
  real_T c1_b[2];
  boolean_T c1_b_p;
  int32_T c1_k;
  boolean_T c1_c_p;
  const mxArray *c1_d_y = NULL;
  int32_T c1_z;
  int32_T c1_b_z;
  int32_T c1_i141;
  boolean_T exitg1;
  (void)chartInstance;
  for (c1_i139 = 0; c1_i139 < 2; c1_i139++) {
    c1_a[c1_i139] = (real_T)c1_x_size[c1_i139];
  }

  for (c1_i140 = 0; c1_i140 < 2; c1_i140++) {
    c1_b[c1_i140] = (real_T)c1_y_size[c1_i140];
  }

  c1_p = false;
  c1_b_p = true;
  c1_k = 0;
  exitg1 = false;
  while ((!exitg1) && (c1_k < 2)) {
    if (!(c1_a[c1_k] == c1_b[c1_k])) {
      c1_b_p = false;
      exitg1 = true;
    } else {
      c1_k++;
    }
  }

  if (!c1_b_p) {
  } else {
    c1_p = true;
  }

  c1_c_p = c1_p;
  if (c1_c_p) {
  } else {
    c1_d_y = NULL;
    sf_mex_assign(&c1_d_y, sf_mex_create("y", c1_cv0, 10, 0U, 1U, 0U, 2, 1, 15),
                  false);
    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                      sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message",
      1U, 1U, 14, c1_d_y));
  }

  c1_z_size[0] = 1;
  c1_z_size[1] = 4;
  c1_z = c1_z_size[0];
  c1_b_z = c1_z_size[1];
  for (c1_i141 = 0; c1_i141 < 4; c1_i141++) {
    c1_z_data[c1_i141] = c1_x_data[c1_i141] / c1_y_data[c1_i141];
  }
}

static const mxArray *c1_emlrt_marshallOut(SFc1_RM3_copyInstanceStruct
  *chartInstance, const real_T c1_b_u[24])
{
  const mxArray *c1_d_y = NULL;
  (void)chartInstance;
  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", c1_b_u, 0, 0U, 1U, 0U, 2, 12, 2),
                false);
  return c1_d_y;
}

static const mxArray *c1_b_emlrt_marshallOut(SFc1_RM3_copyInstanceStruct
  *chartInstance, const real_T c1_b_u)
{
  const mxArray *c1_d_y = NULL;
  (void)chartInstance;
  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", &c1_b_u, 0, 0U, 0U, 0U, 0), false);
  return c1_d_y;
}

static void c1_j_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_datasample, const char_T *c1_identifier, real_T c1_d_y[2])
{
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = (const char *)c1_identifier;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_k_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_datasample), &c1_thisId,
                        c1_d_y);
  sf_mex_destroy(&c1_datasample);
}

static void c1_k_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_d_y[2])
{
  real_T c1_dv6[2];
  int32_T c1_i142;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_b_u), c1_dv6, 1, 0, 0U, 1, 0U, 2, 1,
                2);
  for (c1_i142 = 0; c1_i142 < 2; c1_i142++) {
    c1_d_y[c1_i142] = c1_dv6[c1_i142];
  }

  sf_mex_destroy(&c1_b_u);
}

static const mxArray *c1_c_emlrt_marshallOut(SFc1_RM3_copyInstanceStruct
  *chartInstance, const real_T c1_b_u[72])
{
  const mxArray *c1_d_y = NULL;
  (void)chartInstance;
  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", c1_b_u, 0, 0U, 1U, 0U, 2, 1, 72),
                false);
  return c1_d_y;
}

static real_T c1_l_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_datasample, const char_T *c1_identifier)
{
  real_T c1_d_y;
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = (const char *)c1_identifier;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_d_y = c1_m_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_datasample),
    &c1_thisId);
  sf_mex_destroy(&c1_datasample);
  return c1_d_y;
}

static real_T c1_m_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId)
{
  real_T c1_d_y;
  real_T c1_d1;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_b_u), &c1_d1, 1, 0, 0U, 0, 0U, 0);
  c1_d_y = c1_d1;
  sf_mex_destroy(&c1_b_u);
  return c1_d_y;
}

static const mxArray *c1_d_emlrt_marshallOut(SFc1_RM3_copyInstanceStruct
  *chartInstance, const char_T c1_b_u[16])
{
  const mxArray *c1_d_y = NULL;
  (void)chartInstance;
  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", c1_b_u, 10, 0U, 1U, 0U, 2, 1, 16),
                false);
  return c1_d_y;
}

static const mxArray *c1_e_emlrt_marshallOut(SFc1_RM3_copyInstanceStruct
  *chartInstance, const int16_T c1_b_u)
{
  const mxArray *c1_d_y = NULL;
  (void)chartInstance;
  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", &c1_b_u, 4, 0U, 0U, 0U, 0), false);
  return c1_d_y;
}

static const mxArray *c1_f_emlrt_marshallOut(SFc1_RM3_copyInstanceStruct
  *chartInstance, const char_T c1_b_u[14])
{
  const mxArray *c1_d_y = NULL;
  (void)chartInstance;
  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", c1_b_u, 10, 0U, 1U, 0U, 2, 1, 14),
                false);
  return c1_d_y;
}

static const mxArray *c1_g_emlrt_marshallOut(SFc1_RM3_copyInstanceStruct
  *chartInstance, const real_T c1_u_data[], const int32_T c1_u_size[1])
{
  const mxArray *c1_d_y = NULL;
  (void)chartInstance;
  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", (void *)c1_u_data, 0, 0U, 1U, 0U, 1,
    c1_u_size[0]), false);
  return c1_d_y;
}

static const mxArray *c1_h_emlrt_marshallOut(SFc1_RM3_copyInstanceStruct
  *chartInstance, const real_T c1_b_u[2])
{
  const mxArray *c1_d_y = NULL;
  (void)chartInstance;
  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", c1_b_u, 0, 0U, 1U, 0U, 2, 1, 2),
                false);
  return c1_d_y;
}

static const mxArray *c1_i_emlrt_marshallOut(SFc1_RM3_copyInstanceStruct
  *chartInstance, const char_T c1_b_u[15])
{
  const mxArray *c1_d_y = NULL;
  (void)chartInstance;
  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", c1_b_u, 10, 0U, 1U, 0U, 2, 1, 15),
                false);
  return c1_d_y;
}

static const mxArray *c1_j_emlrt_marshallOut(SFc1_RM3_copyInstanceStruct
  *chartInstance, const real_T c1_u_data[], const int32_T c1_u_size[1])
{
  const mxArray *c1_d_y = NULL;
  (void)chartInstance;
  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", (void *)c1_u_data, 0, 0U, 1U, 0U, 1,
    c1_u_size[0]), false);
  return c1_d_y;
}

static const mxArray *c1_k_emlrt_marshallOut(SFc1_RM3_copyInstanceStruct
  *chartInstance, const real_T c1_u_data[], const int32_T c1_u_size[2])
{
  const mxArray *c1_d_y = NULL;
  (void)chartInstance;
  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", (void *)c1_u_data, 0, 0U, 1U, 0U, 2,
    c1_u_size[0], c1_u_size[1]), false);
  return c1_d_y;
}

static const mxArray *c1_l_emlrt_marshallOut(SFc1_RM3_copyInstanceStruct
  *chartInstance, const real_T c1_u_data[], const int32_T c1_u_size[2])
{
  const mxArray *c1_d_y = NULL;
  (void)chartInstance;
  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", (void *)c1_u_data, 0, 0U, 1U, 0U, 2,
    c1_u_size[0], c1_u_size[1]), false);
  return c1_d_y;
}

static void c1_n_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_transpose, const char_T *c1_identifier, real_T c1_y_data[],
  int32_T c1_y_size[2])
{
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = (const char *)c1_identifier;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_o_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_transpose), &c1_thisId,
                        c1_y_data, c1_y_size);
  sf_mex_destroy(&c1_transpose);
}

static void c1_o_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T
  c1_y_data[], int32_T c1_y_size[2])
{
  int32_T c1_i143;
  int32_T c1_tmp_size[2];
  uint32_T c1_uv6[2];
  int32_T c1_i144;
  real_T c1_b_tmp_data[864];
  boolean_T c1_bv11[2];
  static boolean_T c1_bv12[2] = { true, false };

  int32_T c1_d_y;
  int32_T c1_e_y;
  int32_T c1_loop_ub;
  int32_T c1_i145;
  (void)chartInstance;
  for (c1_i143 = 0; c1_i143 < 2; c1_i143++) {
    c1_uv6[c1_i143] = 72U + (uint32_T)(-60 * c1_i143);
  }

  c1_tmp_size[0] = sf_mex_get_dimension(c1_b_u, 0);
  c1_tmp_size[1] = sf_mex_get_dimension(c1_b_u, 1);
  for (c1_i144 = 0; c1_i144 < 2; c1_i144++) {
    c1_bv11[c1_i144] = c1_bv12[c1_i144];
  }

  sf_mex_import_vs(c1_parentId, sf_mex_dup(c1_b_u), (void *)&c1_b_tmp_data, 0, 0,
                   0U, 1, 0U, 2, c1_bv11, c1_uv6, c1_tmp_size);
  c1_y_size[0] = c1_tmp_size[0];
  c1_y_size[1] = 12;
  c1_d_y = c1_y_size[0];
  c1_e_y = c1_y_size[1];
  c1_loop_ub = c1_tmp_size[0] * c1_tmp_size[1] - 1;
  for (c1_i145 = 0; c1_i145 <= c1_loop_ub; c1_i145++) {
    c1_y_data[c1_i145] = c1_b_tmp_data[c1_i145];
  }

  sf_mex_destroy(&c1_b_u);
}

static const mxArray *c1_m_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData;
  int32_T c1_b_u;
  const mxArray *c1_d_y = NULL;
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  c1_b_u = *(int32_T *)c1_inData;
  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", &c1_b_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_d_y, false);
  return c1_mxArrayOutData;
}

static int32_T c1_p_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId)
{
  int32_T c1_d_y;
  int32_T c1_i146;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_b_u), &c1_i146, 1, 6, 0U, 0, 0U, 0);
  c1_d_y = c1_i146;
  sf_mex_destroy(&c1_b_u);
  return c1_d_y;
}

static void c1_l_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_b_sfEvent;
  emlrtMsgIdentifier c1_thisId;
  int32_T c1_d_y;
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_b_sfEvent = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_d_y = c1_p_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_sfEvent),
    &c1_thisId);
  sf_mex_destroy(&c1_b_sfEvent);
  *(int32_T *)c1_outData = c1_d_y;
  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_n_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData;
  int32_T c1_i147;
  int32_T c1_i148;
  const mxArray *c1_d_y = NULL;
  int32_T c1_i149;
  real_T c1_b_u[24];
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  c1_i147 = 0;
  for (c1_i148 = 0; c1_i148 < 2; c1_i148++) {
    for (c1_i149 = 0; c1_i149 < 12; c1_i149++) {
      c1_b_u[c1_i149 + c1_i147] = (*(real_T (*)[24])c1_inData)[c1_i149 + c1_i147];
    }

    c1_i147 += 12;
  }

  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", c1_b_u, 0, 0U, 1U, 0U, 2, 12, 2),
                false);
  sf_mex_assign(&c1_mxArrayOutData, c1_d_y, false);
  return c1_mxArrayOutData;
}

static void c1_q_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_d_y[24])
{
  real_T c1_dv7[24];
  int32_T c1_i150;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_b_u), c1_dv7, 1, 0, 0U, 1, 0U, 2, 12,
                2);
  for (c1_i150 = 0; c1_i150 < 24; c1_i150++) {
    c1_d_y[c1_i150] = c1_dv7[c1_i150];
  }

  sf_mex_destroy(&c1_b_u);
}

static void c1_m_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_action;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_d_y[24];
  int32_T c1_i151;
  int32_T c1_i152;
  int32_T c1_i153;
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_action = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_q_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_action), &c1_thisId, c1_d_y);
  sf_mex_destroy(&c1_action);
  c1_i151 = 0;
  for (c1_i152 = 0; c1_i152 < 2; c1_i152++) {
    for (c1_i153 = 0; c1_i153 < 12; c1_i153++) {
      (*(real_T (*)[24])c1_outData)[c1_i153 + c1_i151] = c1_d_y[c1_i153 +
        c1_i151];
    }

    c1_i151 += 12;
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_o_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData;
  int32_T c1_i154;
  int32_T c1_i155;
  const mxArray *c1_d_y = NULL;
  int32_T c1_i156;
  real_T c1_b_u[15000];
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  c1_i154 = 0;
  for (c1_i155 = 0; c1_i155 < 5; c1_i155++) {
    for (c1_i156 = 0; c1_i156 < 3000; c1_i156++) {
      c1_b_u[c1_i156 + c1_i154] = (*(real_T (*)[15000])c1_inData)[c1_i156 +
        c1_i154];
    }

    c1_i154 += 3000;
  }

  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", c1_b_u, 0, 0U, 1U, 0U, 2, 3000, 5),
                false);
  sf_mex_assign(&c1_mxArrayOutData, c1_d_y, false);
  return c1_mxArrayOutData;
}

static void c1_r_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_d_y
  [15000])
{
  real_T c1_dv8[15000];
  int32_T c1_i157;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_b_u), c1_dv8, 1, 0, 0U, 1, 0U, 2,
                3000, 5);
  for (c1_i157 = 0; c1_i157 < 15000; c1_i157++) {
    c1_d_y[c1_i157] = c1_dv8[c1_i157];
  }

  sf_mex_destroy(&c1_b_u);
}

static void c1_n_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_Stor_eta;
  emlrtMsgIdentifier c1_thisId;
  int32_T c1_i158;
  int32_T c1_i159;
  int32_T c1_i160;
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_Stor_eta = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_r_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_Stor_eta), &c1_thisId,
                        chartInstance->c1_b_y);
  sf_mex_destroy(&c1_Stor_eta);
  c1_i158 = 0;
  for (c1_i159 = 0; c1_i159 < 5; c1_i159++) {
    for (c1_i160 = 0; c1_i160 < 3000; c1_i160++) {
      (*(real_T (*)[15000])c1_outData)[c1_i160 + c1_i158] =
        chartInstance->c1_b_y[c1_i160 + c1_i158];
    }

    c1_i158 += 3000;
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_p_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData;
  int16_T c1_b_u;
  const mxArray *c1_d_y = NULL;
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  c1_b_u = *(int16_T *)c1_inData;
  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", &c1_b_u, 4, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_d_y, false);
  return c1_mxArrayOutData;
}

static int16_T c1_s_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId)
{
  int16_T c1_d_y;
  int16_T c1_i161;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_b_u), &c1_i161, 1, 4, 0U, 0, 0U, 0);
  c1_d_y = c1_i161;
  sf_mex_destroy(&c1_b_u);
  return c1_d_y;
}

static void c1_o_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_state_idx;
  emlrtMsgIdentifier c1_thisId;
  int16_T c1_d_y;
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_state_idx = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_d_y = c1_s_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_state_idx),
    &c1_thisId);
  sf_mex_destroy(&c1_state_idx);
  *(int16_T *)c1_outData = c1_d_y;
  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_q_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData;
  int32_T c1_i162;
  int32_T c1_i163;
  const mxArray *c1_d_y = NULL;
  int32_T c1_i164;
  real_T c1_b_u[288];
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  c1_i162 = 0;
  for (c1_i163 = 0; c1_i163 < 4; c1_i163++) {
    for (c1_i164 = 0; c1_i164 < 72; c1_i164++) {
      c1_b_u[c1_i164 + c1_i162] = (*(real_T (*)[288])c1_inData)[c1_i164 +
        c1_i162];
    }

    c1_i162 += 72;
  }

  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", c1_b_u, 0, 0U, 1U, 0U, 2, 72, 4),
                false);
  sf_mex_assign(&c1_mxArrayOutData, c1_d_y, false);
  return c1_mxArrayOutData;
}

static void c1_t_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_d_y
  [288])
{
  real_T c1_dv9[288];
  int32_T c1_i165;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_b_u), c1_dv9, 1, 0, 0U, 1, 0U, 2, 72,
                4);
  for (c1_i165 = 0; c1_i165 < 288; c1_i165++) {
    c1_d_y[c1_i165] = c1_dv9[c1_i165];
  }

  sf_mex_destroy(&c1_b_u);
}

static void c1_p_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_states_space;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_d_y[288];
  int32_T c1_i166;
  int32_T c1_i167;
  int32_T c1_i168;
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_states_space = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_t_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_states_space), &c1_thisId,
                        c1_d_y);
  sf_mex_destroy(&c1_states_space);
  c1_i166 = 0;
  for (c1_i167 = 0; c1_i167 < 4; c1_i167++) {
    for (c1_i168 = 0; c1_i168 < 72; c1_i168++) {
      (*(real_T (*)[288])c1_outData)[c1_i168 + c1_i166] = c1_d_y[c1_i168 +
        c1_i166];
    }

    c1_i166 += 72;
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_r_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData;
  int32_T c1_i169;
  int32_T c1_i170;
  const mxArray *c1_d_y = NULL;
  int32_T c1_i171;
  int16_T c1_b_u[10000];
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  c1_i169 = 0;
  for (c1_i170 = 0; c1_i170 < 100; c1_i170++) {
    for (c1_i171 = 0; c1_i171 < 100; c1_i171++) {
      c1_b_u[c1_i171 + c1_i169] = (*(int16_T (*)[10000])c1_inData)[c1_i171 +
        c1_i169];
    }

    c1_i169 += 100;
  }

  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", c1_b_u, 4, 0U, 1U, 0U, 2, 100, 100),
                false);
  sf_mex_assign(&c1_mxArrayOutData, c1_d_y, false);
  return c1_mxArrayOutData;
}

static void c1_u_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, int16_T c1_d_y
  [10000])
{
  int16_T c1_iv0[10000];
  int32_T c1_i172;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_b_u), c1_iv0, 1, 4, 0U, 1, 0U, 2, 100,
                100);
  for (c1_i172 = 0; c1_i172 < 10000; c1_i172++) {
    c1_d_y[c1_i172] = c1_iv0[c1_i172];
  }

  sf_mex_destroy(&c1_b_u);
}

static void c1_q_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_num_R;
  emlrtMsgIdentifier c1_thisId;
  int16_T c1_d_y[10000];
  int32_T c1_i173;
  int32_T c1_i174;
  int32_T c1_i175;
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_num_R = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_u_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_num_R), &c1_thisId, c1_d_y);
  sf_mex_destroy(&c1_num_R);
  c1_i173 = 0;
  for (c1_i174 = 0; c1_i174 < 100; c1_i174++) {
    for (c1_i175 = 0; c1_i175 < 100; c1_i175++) {
      (*(int16_T (*)[10000])c1_outData)[c1_i175 + c1_i173] = c1_d_y[c1_i175 +
        c1_i173];
    }

    c1_i173 += 100;
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_s_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData;
  int32_T c1_i176;
  int32_T c1_i177;
  const mxArray *c1_d_y = NULL;
  int32_T c1_i178;
  real_T c1_b_u[3600];
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  c1_i176 = 0;
  for (c1_i177 = 0; c1_i177 < 50; c1_i177++) {
    for (c1_i178 = 0; c1_i178 < 72; c1_i178++) {
      c1_b_u[c1_i178 + c1_i176] = (*(real_T (*)[3600])c1_inData)[c1_i178 +
        c1_i176];
    }

    c1_i176 += 72;
  }

  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", c1_b_u, 0, 0U, 1U, 0U, 2, 72, 50),
                false);
  sf_mex_assign(&c1_mxArrayOutData, c1_d_y, false);
  return c1_mxArrayOutData;
}

static void c1_v_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_d_y
  [3600])
{
  real_T c1_dv10[3600];
  int32_T c1_i179;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_b_u), c1_dv10, 1, 0, 0U, 1, 0U, 2, 72,
                50);
  for (c1_i179 = 0; c1_i179 < 3600; c1_i179++) {
    c1_d_y[c1_i179] = c1_dv10[c1_i179];
  }

  sf_mex_destroy(&c1_b_u);
}

static void c1_r_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_R;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_d_y[3600];
  int32_T c1_i180;
  int32_T c1_i181;
  int32_T c1_i182;
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_R = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_v_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_R), &c1_thisId, c1_d_y);
  sf_mex_destroy(&c1_R);
  c1_i180 = 0;
  for (c1_i181 = 0; c1_i181 < 50; c1_i181++) {
    for (c1_i182 = 0; c1_i182 < 72; c1_i182++) {
      (*(real_T (*)[3600])c1_outData)[c1_i182 + c1_i180] = c1_d_y[c1_i182 +
        c1_i180];
    }

    c1_i180 += 72;
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_t_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData;
  int32_T c1_i183;
  int32_T c1_i184;
  const mxArray *c1_d_y = NULL;
  int32_T c1_i185;
  real_T c1_b_u[10000];
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  c1_i183 = 0;
  for (c1_i184 = 0; c1_i184 < 100; c1_i184++) {
    for (c1_i185 = 0; c1_i185 < 100; c1_i185++) {
      c1_b_u[c1_i185 + c1_i183] = (*(real_T (*)[10000])c1_inData)[c1_i185 +
        c1_i183];
    }

    c1_i183 += 100;
  }

  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", c1_b_u, 0, 0U, 1U, 0U, 2, 100, 100),
                false);
  sf_mex_assign(&c1_mxArrayOutData, c1_d_y, false);
  return c1_mxArrayOutData;
}

static void c1_w_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_d_y
  [10000])
{
  real_T c1_dv11[10000];
  int32_T c1_i186;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_b_u), c1_dv11, 1, 0, 0U, 1, 0U, 2,
                100, 100);
  for (c1_i186 = 0; c1_i186 < 10000; c1_i186++) {
    c1_d_y[c1_i186] = c1_dv11[c1_i186];
  }

  sf_mex_destroy(&c1_b_u);
}

static void c1_s_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_MeanR;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_d_y[10000];
  int32_T c1_i187;
  int32_T c1_i188;
  int32_T c1_i189;
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_MeanR = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_w_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_MeanR), &c1_thisId, c1_d_y);
  sf_mex_destroy(&c1_MeanR);
  c1_i187 = 0;
  for (c1_i188 = 0; c1_i188 < 100; c1_i188++) {
    for (c1_i189 = 0; c1_i189 < 100; c1_i189++) {
      (*(real_T (*)[10000])c1_outData)[c1_i189 + c1_i187] = c1_d_y[c1_i189 +
        c1_i187];
    }

    c1_i187 += 100;
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_u_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData;
  int32_T c1_i190;
  int32_T c1_i191;
  const mxArray *c1_d_y = NULL;
  int32_T c1_i192;
  real_T c1_b_u[1000];
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  c1_i190 = 0;
  for (c1_i191 = 0; c1_i191 < 10; c1_i191++) {
    for (c1_i192 = 0; c1_i192 < 100; c1_i192++) {
      c1_b_u[c1_i192 + c1_i190] = (*(real_T (*)[1000])c1_inData)[c1_i192 +
        c1_i190];
    }

    c1_i190 += 100;
  }

  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", c1_b_u, 0, 0U, 1U, 0U, 2, 100, 10),
                false);
  sf_mex_assign(&c1_mxArrayOutData, c1_d_y, false);
  return c1_mxArrayOutData;
}

static void c1_x_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_d_y
  [1000])
{
  real_T c1_dv12[1000];
  int32_T c1_i193;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_b_u), c1_dv12, 1, 0, 0U, 1, 0U, 2,
                100, 10);
  for (c1_i193 = 0; c1_i193 < 1000; c1_i193++) {
    c1_d_y[c1_i193] = c1_dv12[c1_i193];
  }

  sf_mex_destroy(&c1_b_u);
}

static void c1_t_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_Final_Rew_matt;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_d_y[1000];
  int32_T c1_i194;
  int32_T c1_i195;
  int32_T c1_i196;
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_Final_Rew_matt = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_x_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_Final_Rew_matt), &c1_thisId,
                        c1_d_y);
  sf_mex_destroy(&c1_Final_Rew_matt);
  c1_i194 = 0;
  for (c1_i195 = 0; c1_i195 < 10; c1_i195++) {
    for (c1_i196 = 0; c1_i196 < 100; c1_i196++) {
      (*(real_T (*)[1000])c1_outData)[c1_i196 + c1_i194] = c1_d_y[c1_i196 +
        c1_i194];
    }

    c1_i194 += 100;
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_v_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData;
  int32_T c1_i197;
  int32_T c1_i198;
  const mxArray *c1_d_y = NULL;
  int32_T c1_i199;
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  c1_i197 = 0;
  for (c1_i198 = 0; c1_i198 < 30; c1_i198++) {
    for (c1_i199 = 0; c1_i199 < 5000; c1_i199++) {
      chartInstance->c1_u[c1_i199 + c1_i197] = (*(real_T (*)[150000])c1_inData)
        [c1_i199 + c1_i197];
    }

    c1_i197 += 5000;
  }

  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", chartInstance->c1_u, 0, 0U, 1U, 0U,
    2, 5000, 30), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_d_y, false);
  return c1_mxArrayOutData;
}

static void c1_y_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_d_y
  [150000])
{
  int32_T c1_i200;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_b_u), chartInstance->c1_dv13, 1, 0,
                0U, 1, 0U, 2, 5000, 30);
  for (c1_i200 = 0; c1_i200 < 150000; c1_i200++) {
    c1_d_y[c1_i200] = chartInstance->c1_dv13[c1_i200];
  }

  sf_mex_destroy(&c1_b_u);
}

static void c1_u_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_Memory_Stor;
  emlrtMsgIdentifier c1_thisId;
  int32_T c1_i201;
  int32_T c1_i202;
  int32_T c1_i203;
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)chartInstanceVoid;
  c1_Memory_Stor = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_y_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_Memory_Stor), &c1_thisId,
                        chartInstance->c1_y);
  sf_mex_destroy(&c1_Memory_Stor);
  c1_i201 = 0;
  for (c1_i202 = 0; c1_i202 < 30; c1_i202++) {
    for (c1_i203 = 0; c1_i203 < 5000; c1_i203++) {
      (*(real_T (*)[150000])c1_outData)[c1_i203 + c1_i201] = chartInstance->
        c1_y[c1_i203 + c1_i201];
    }

    c1_i201 += 5000;
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static uint8_T c1_ab_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_is_active_c1_RM3_copy, const char_T *c1_identifier)
{
  uint8_T c1_d_y;
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = (const char *)c1_identifier;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_d_y = c1_bb_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c1_b_is_active_c1_RM3_copy), &c1_thisId);
  sf_mex_destroy(&c1_b_is_active_c1_RM3_copy);
  return c1_d_y;
}

static uint8_T c1_bb_emlrt_marshallIn(SFc1_RM3_copyInstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId)
{
  uint8_T c1_d_y;
  uint8_T c1_u0;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_b_u), &c1_u0, 1, 3, 0U, 0, 0U, 0);
  c1_d_y = c1_u0;
  sf_mex_destroy(&c1_b_u);
  return c1_d_y;
}

static void c1_emxEnsureCapacity_real_T(SFc1_RM3_copyInstanceStruct
  *chartInstance, c1_emxArray_real_T *c1_emxArray, int32_T c1_oldNumel, const
  emlrtRTEInfo *c1_srcLocation)
{
  int32_T c1_newNumel;
  int32_T c1_i;
  int32_T c1_newCapacity;
  void *c1_newData;
  if (c1_oldNumel < 0) {
    c1_oldNumel = 0;
  }

  c1_newNumel = 1;
  for (c1_i = 0; c1_i < c1_emxArray->numDimensions; c1_i++) {
    c1_newNumel = (int32_T)emlrtSizeMulR2012b((uint32_T)c1_newNumel, (uint32_T)
      c1_emxArray->size[c1_i], c1_srcLocation, chartInstance->c1_fEmlrtCtx);
  }

  if (c1_newNumel > c1_emxArray->allocatedSize) {
    c1_newCapacity = c1_emxArray->allocatedSize;
    if (c1_newCapacity < 16) {
      c1_newCapacity = 16;
    }

    while (c1_newCapacity < c1_newNumel) {
      if (c1_newCapacity > 1073741823) {
        c1_newCapacity = MAX_int32_T;
      } else {
        c1_newCapacity <<= 1;
      }
    }

    c1_newData = emlrtCallocMex((uint32_T)c1_newCapacity, sizeof(real_T));
    if (c1_newData == NULL) {
      emlrtHeapAllocationErrorR2012b(c1_srcLocation, chartInstance->c1_fEmlrtCtx);
    }

    if (c1_emxArray->data != NULL) {
      memcpy(c1_newData, (void *)c1_emxArray->data, sizeof(real_T) * (uint32_T)
             c1_oldNumel);
      if (c1_emxArray->canFreeData) {
        emlrtFreeMex((void *)c1_emxArray->data);
      }
    }

    c1_emxArray->data = (real_T *)c1_newData;
    c1_emxArray->allocatedSize = c1_newCapacity;
    c1_emxArray->canFreeData = true;
  }
}

static void c1_emxInit_real_T(SFc1_RM3_copyInstanceStruct *chartInstance,
  c1_emxArray_real_T **c1_pEmxArray, int32_T c1_numDimensions, const
  emlrtRTEInfo *c1_srcLocation)
{
  c1_emxArray_real_T *c1_emxArray;
  int32_T c1_i;
  *c1_pEmxArray = (c1_emxArray_real_T *)emlrtMallocMex(sizeof(c1_emxArray_real_T));
  if ((void *)*c1_pEmxArray == NULL) {
    emlrtHeapAllocationErrorR2012b(c1_srcLocation, chartInstance->c1_fEmlrtCtx);
  }

  c1_emxArray = *c1_pEmxArray;
  c1_emxArray->data = (real_T *)NULL;
  c1_emxArray->numDimensions = c1_numDimensions;
  c1_emxArray->size = (int32_T *)emlrtMallocMex((uint32_T)(sizeof(int32_T)
    * c1_numDimensions));
  if ((void *)c1_emxArray->size == NULL) {
    emlrtHeapAllocationErrorR2012b(c1_srcLocation, chartInstance->c1_fEmlrtCtx);
  }

  c1_emxArray->allocatedSize = 0;
  c1_emxArray->canFreeData = true;
  for (c1_i = 0; c1_i < c1_numDimensions; c1_i++) {
    c1_emxArray->size[c1_i] = 0;
  }
}

static void c1_emxFree_real_T(SFc1_RM3_copyInstanceStruct *chartInstance,
  c1_emxArray_real_T **c1_pEmxArray)
{
  (void)chartInstance;
  if (*c1_pEmxArray != (c1_emxArray_real_T *)NULL) {
    if (((*c1_pEmxArray)->data != (real_T *)NULL) && (*c1_pEmxArray)
        ->canFreeData) {
      emlrtFreeMex((void *)(*c1_pEmxArray)->data);
    }

    emlrtFreeMex((void *)(*c1_pEmxArray)->size);
    emlrtFreeMex((void *)*c1_pEmxArray);
    *c1_pEmxArray = (c1_emxArray_real_T *)NULL;
  }
}

static int32_T c1__s32_d_(SFc1_RM3_copyInstanceStruct *chartInstance, real_T
  c1_b, uint32_T c1_ssid_src_loc, int32_T c1_offset_src_loc, int32_T
  c1_length_src_loc)
{
  int32_T c1_a;
  real_T c1_b_b;
  (void)chartInstance;
  c1_a = (int32_T)c1_b;
  if (c1_b < 0.0) {
    c1_b_b = muDoubleScalarCeil(c1_b);
  } else {
    c1_b_b = muDoubleScalarFloor(c1_b);
  }

  if ((real_T)c1_a != c1_b_b) {
    _SFD_OVERFLOW_DETECTION(SFDB_OVERFLOW, c1_ssid_src_loc, c1_offset_src_loc,
      c1_length_src_loc);
  }

  return c1_a;
}

static real_T c1_get_B2(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_B2) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 0, "B2",
      c1_elementIndex);
  }

  return *chartInstance->c1_B2_address;
}

static void c1_set_B2(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
                      c1_elementIndex, real_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_B2) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 0, "B2", c1_elementIndex);
  }

  *chartInstance->c1_B2_address = c1_elementValue;
}

static real_T *c1_access_B2(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_B2) {
    ssAccessDataStore_wrapper(chartInstance->S, 0, "B2", c1_rdOnly);
  }

  return chartInstance->c1_B2_address;
}

static real_T c1_get_Final_Rew_matt(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_Final_Rew_matt) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 1, "Final_Rew_matt",
      c1_elementIndex);
  }

  return (*chartInstance->c1__indexFinal_Rew_matt_address)[c1_elementIndex];
}

static void c1_set_Final_Rew_matt(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_elementIndex, real_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_Final_Rew_matt) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 1, "Final_Rew_matt",
      c1_elementIndex);
  }

  (*chartInstance->c1__indexFinal_Rew_matt_address)[c1_elementIndex] =
    c1_elementValue;
}

static real_T *c1_access_Final_Rew_matt(SFc1_RM3_copyInstanceStruct
  *chartInstance, uint32_T c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_Final_Rew_matt) {
    ssAccessDataStore_wrapper(chartInstance->S, 1, "Final_Rew_matt", c1_rdOnly);
  }

  return &(*chartInstance->c1__indexFinal_Rew_matt_address)[0U];
}

static real_T c1_get_K2(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_K2) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 2, "K2",
      c1_elementIndex);
  }

  return *chartInstance->c1__indexK2_address;
}

static void c1_set_K2(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
                      c1_elementIndex, real_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_K2) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 2, "K2", c1_elementIndex);
  }

  *chartInstance->c1__indexK2_address = c1_elementValue;
}

static real_T *c1_access_K2(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_K2) {
    ssAccessDataStore_wrapper(chartInstance->S, 2, "K2", c1_rdOnly);
  }

  return chartInstance->c1__indexK2_address;
}

static real_T c1_get_MeanR(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_MeanR) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 3, "MeanR",
      c1_elementIndex);
  }

  return (*chartInstance->c1__indexMeanR_address)[c1_elementIndex];
}

static void c1_set_MeanR(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex, real_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_MeanR) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 3, "MeanR",
      c1_elementIndex);
  }

  (*chartInstance->c1__indexMeanR_address)[c1_elementIndex] = c1_elementValue;
}

static real_T *c1_access_MeanR(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_MeanR) {
    ssAccessDataStore_wrapper(chartInstance->S, 3, "MeanR", c1_rdOnly);
  }

  return &(*chartInstance->c1__indexMeanR_address)[0U];
}

static real_T c1_get_Memory_D_size(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_Memory_D_size) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 4, "Memory_D_size",
      c1_elementIndex);
  }

  return *chartInstance->c1__indexMemory_D_size_address;
}

static void c1_set_Memory_D_size(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_elementIndex, real_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_Memory_D_size) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 4, "Memory_D_size",
      c1_elementIndex);
  }

  *chartInstance->c1__indexMemory_D_size_address = c1_elementValue;
}

static real_T *c1_access_Memory_D_size(SFc1_RM3_copyInstanceStruct
  *chartInstance, uint32_T c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_Memory_D_size) {
    ssAccessDataStore_wrapper(chartInstance->S, 4, "Memory_D_size", c1_rdOnly);
  }

  return chartInstance->c1__indexMemory_D_size_address;
}

static real_T c1_get_Memory_Stor(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_Memory_Stor) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 5, "Memory_Stor",
      c1_elementIndex);
  }

  return (*chartInstance->c1__indexMemory_Stor_address)[c1_elementIndex];
}

static void c1_set_Memory_Stor(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_elementIndex, real_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_Memory_Stor) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 5, "Memory_Stor",
      c1_elementIndex);
  }

  (*chartInstance->c1__indexMemory_Stor_address)[c1_elementIndex] =
    c1_elementValue;
}

static real_T *c1_access_Memory_Stor(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_Memory_Stor) {
    ssAccessDataStore_wrapper(chartInstance->S, 5, "Memory_Stor", c1_rdOnly);
  }

  return &(*chartInstance->c1__indexMemory_Stor_address)[0U];
}

static real_T c1_get_NN_dd(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_NN_dd) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 6, "NN_dd",
      c1_elementIndex);
  }

  return *chartInstance->c1__indexNN_dd_address;
}

static void c1_set_NN_dd(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex, real_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_NN_dd) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 6, "NN_dd",
      c1_elementIndex);
  }

  *chartInstance->c1__indexNN_dd_address = c1_elementValue;
}

static real_T *c1_access_NN_dd(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_NN_dd) {
    ssAccessDataStore_wrapper(chartInstance->S, 6, "NN_dd", c1_rdOnly);
  }

  return chartInstance->c1__indexNN_dd_address;
}

static real_T c1_get_Q(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_Q) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 7, "Q", c1_elementIndex);
  }

  return (*chartInstance->c1__indexQ_address)[c1_elementIndex];
}

static void c1_set_Q(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
                     c1_elementIndex, real_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_Q) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 7, "Q", c1_elementIndex);
  }

  (*chartInstance->c1__indexQ_address)[c1_elementIndex] = c1_elementValue;
}

static real_T *c1_access_Q(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_Q) {
    ssAccessDataStore_wrapper(chartInstance->S, 7, "Q", c1_rdOnly);
  }

  return &(*chartInstance->c1__indexQ_address)[0U];
}

static real_T c1_get_Q_prime(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_Q_prime) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 8, "Q_prime",
      c1_elementIndex);
  }

  return (*chartInstance->c1__indexQ_prime_address)[c1_elementIndex];
}

static void c1_set_Q_prime(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex, real_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_Q_prime) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 8, "Q_prime",
      c1_elementIndex);
  }

  (*chartInstance->c1__indexQ_prime_address)[c1_elementIndex] = c1_elementValue;
}

static real_T *c1_access_Q_prime(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_Q_prime) {
    ssAccessDataStore_wrapper(chartInstance->S, 8, "Q_prime", c1_rdOnly);
  }

  return &(*chartInstance->c1__indexQ_prime_address)[0U];
}

static real_T c1_get_Q_targ(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_Q_targ) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 9, "Q_targ",
      c1_elementIndex);
  }

  return (*chartInstance->c1__indexQ_targ_address)[c1_elementIndex];
}

static void c1_set_Q_targ(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex, real_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_Q_targ) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 9, "Q_targ",
      c1_elementIndex);
  }

  (*chartInstance->c1__indexQ_targ_address)[c1_elementIndex] = c1_elementValue;
}

static real_T *c1_access_Q_targ(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_Q_targ) {
    ssAccessDataStore_wrapper(chartInstance->S, 9, "Q_targ", c1_rdOnly);
  }

  return &(*chartInstance->c1__indexQ_targ_address)[0U];
}

static real_T c1_get_R(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_R) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 10, "R",
      c1_elementIndex);
  }

  return (*chartInstance->c1__indexR_address)[c1_elementIndex];
}

static void c1_set_R(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
                     c1_elementIndex, real_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_R) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 10, "R", c1_elementIndex);
  }

  (*chartInstance->c1__indexR_address)[c1_elementIndex] = c1_elementValue;
}

static real_T *c1_access_R(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_R) {
    ssAccessDataStore_wrapper(chartInstance->S, 10, "R", c1_rdOnly);
  }

  return &(*chartInstance->c1__indexR_address)[0U];
}

static real_T c1_get_Stor_eta(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_Stor_eta) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 11, "Stor_eta",
      c1_elementIndex);
  }

  return (*chartInstance->c1__indexStor_eta_address)[c1_elementIndex];
}

static void c1_set_Stor_eta(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex, real_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_Stor_eta) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 11, "Stor_eta",
      c1_elementIndex);
  }

  (*chartInstance->c1__indexStor_eta_address)[c1_elementIndex] = c1_elementValue;
}

static real_T *c1_access_Stor_eta(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_Stor_eta) {
    ssAccessDataStore_wrapper(chartInstance->S, 11, "Stor_eta", c1_rdOnly);
  }

  return &(*chartInstance->c1__indexStor_eta_address)[0U];
}

static real_T c1_get_action(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_action) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 12, "action",
      c1_elementIndex);
  }

  return (*chartInstance->c1__indexaction_address)[c1_elementIndex];
}

static void c1_set_action(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex, real_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_action) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 12, "action",
      c1_elementIndex);
  }

  (*chartInstance->c1__indexaction_address)[c1_elementIndex] = c1_elementValue;
}

static real_T *c1_access_action(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_action) {
    ssAccessDataStore_wrapper(chartInstance->S, 12, "action", c1_rdOnly);
  }

  return &(*chartInstance->c1__indexaction_address)[0U];
}

static int16_T c1_get_num_R(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_num_R) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 13, "num_R",
      c1_elementIndex);
  }

  return (*chartInstance->c1__indexnum_R_address)[c1_elementIndex];
}

static void c1_set_num_R(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex, int16_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_num_R) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 13, "num_R",
      c1_elementIndex);
  }

  (*chartInstance->c1__indexnum_R_address)[c1_elementIndex] = c1_elementValue;
}

static int16_T *c1_access_num_R(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_num_R) {
    ssAccessDataStore_wrapper(chartInstance->S, 13, "num_R", c1_rdOnly);
  }

  return &(*chartInstance->c1__indexnum_R_address)[0U];
}

static real_T c1_get_ove_count(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_ove_count) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 14, "ove_count",
      c1_elementIndex);
  }

  return *chartInstance->c1__indexove_count_address;
}

static void c1_set_ove_count(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_elementIndex, real_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_ove_count) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 14, "ove_count",
      c1_elementIndex);
  }

  *chartInstance->c1__indexove_count_address = c1_elementValue;
}

static real_T *c1_access_ove_count(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_ove_count) {
    ssAccessDataStore_wrapper(chartInstance->S, 14, "ove_count", c1_rdOnly);
  }

  return chartInstance->c1__indexove_count_address;
}

static int16_T c1_get_state_idx(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_state_idx) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 15, "state_idx",
      c1_elementIndex);
  }

  return *chartInstance->c1__indexstate_idx_address;
}

static void c1_set_state_idx(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_elementIndex, int16_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_state_idx) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 15, "state_idx",
      c1_elementIndex);
  }

  *chartInstance->c1__indexstate_idx_address = c1_elementValue;
}

static int16_T *c1_access_state_idx(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_state_idx) {
    ssAccessDataStore_wrapper(chartInstance->S, 15, "state_idx", c1_rdOnly);
  }

  return chartInstance->c1__indexstate_idx_address;
}

static real_T c1_get_states_space(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_states_space) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 16, "states_space",
      c1_elementIndex);
  }

  return (*chartInstance->c1__indexstates_space_address)[c1_elementIndex];
}

static void c1_set_states_space(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_elementIndex, real_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_states_space) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 16, "states_space",
      c1_elementIndex);
  }

  (*chartInstance->c1__indexstates_space_address)[c1_elementIndex] =
    c1_elementValue;
}

static real_T *c1_access_states_space(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_states_space) {
    ssAccessDataStore_wrapper(chartInstance->S, 16, "states_space", c1_rdOnly);
  }

  return &(*chartInstance->c1__indexstates_space_address)[0U];
}

static real_T c1_get_target(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_target) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 17, "target",
      c1_elementIndex);
  }

  return (*chartInstance->c1__indextarget_address)[c1_elementIndex];
}

static void c1_set_target(SFc1_RM3_copyInstanceStruct *chartInstance, uint32_T
  c1_elementIndex, real_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_target) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 17, "target",
      c1_elementIndex);
  }

  (*chartInstance->c1__indextarget_address)[c1_elementIndex] = c1_elementValue;
}

static real_T *c1_access_target(SFc1_RM3_copyInstanceStruct *chartInstance,
  uint32_T c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_target) {
    ssAccessDataStore_wrapper(chartInstance->S, 17, "target", c1_rdOnly);
  }

  return &(*chartInstance->c1__indextarget_address)[0U];
}

static void init_dsm_address_info(SFc1_RM3_copyInstanceStruct *chartInstance)
{
  ssGetSFcnDataStoreNameAddrIdx_wrapper(chartInstance->S, "B2", (void **)
    &chartInstance->c1_B2_address, &chartInstance->c1__index);
  ssGetSFcnDataStoreNameAddrIdx_wrapper(chartInstance->S, "Final_Rew_matt",
    (void **)&chartInstance->c1__indexFinal_Rew_matt_address,
    &chartInstance->c1__indexRew_matt);
  ssGetSFcnDataStoreNameAddrIdx_wrapper(chartInstance->S, "K2", (void **)
    &chartInstance->c1__indexK2_address, &chartInstance->c1_b__index);
  ssGetSFcnDataStoreNameAddrIdx_wrapper(chartInstance->S, "MeanR", (void **)
    &chartInstance->c1__indexMeanR_address, &chartInstance->c1_c__index);
  ssGetSFcnDataStoreNameAddrIdx_wrapper(chartInstance->S, "Memory_D_size", (void
    **)&chartInstance->c1__indexMemory_D_size_address,
    &chartInstance->c1__index_D_size);
  ssGetSFcnDataStoreNameAddrIdx_wrapper(chartInstance->S, "Memory_Stor", (void **)
    &chartInstance->c1__indexMemory_Stor_address, &chartInstance->c1__index_Stor);
  ssGetSFcnDataStoreNameAddrIdx_wrapper(chartInstance->S, "NN_dd", (void **)
    &chartInstance->c1__indexNN_dd_address, &chartInstance->c1_d__index);
  ssGetSFcnDataStoreNameAddrIdx_wrapper(chartInstance->S, "Q", (void **)
    &chartInstance->c1__indexQ_address, &chartInstance->c1_e__index);
  ssGetSFcnDataStoreNameAddrIdx_wrapper(chartInstance->S, "Q_prime", (void **)
    &chartInstance->c1__indexQ_prime_address, &chartInstance->c1__indexe);
  ssGetSFcnDataStoreNameAddrIdx_wrapper(chartInstance->S, "Q_targ", (void **)
    &chartInstance->c1__indexQ_targ_address, &chartInstance->c1_f__index);
  ssGetSFcnDataStoreNameAddrIdx_wrapper(chartInstance->S, "R", (void **)
    &chartInstance->c1__indexR_address, &chartInstance->c1_g__index);
  ssGetSFcnDataStoreNameAddrIdx_wrapper(chartInstance->S, "Stor_eta", (void **)
    &chartInstance->c1__indexStor_eta_address, &chartInstance->c1__indexta);
  ssGetSFcnDataStoreNameAddrIdx_wrapper(chartInstance->S, "action", (void **)
    &chartInstance->c1__indexaction_address, &chartInstance->c1_h__index);
  ssGetSFcnDataStoreNameAddrIdx_wrapper(chartInstance->S, "num_R", (void **)
    &chartInstance->c1__indexnum_R_address, &chartInstance->c1_i__index);
  ssGetSFcnDataStoreNameAddrIdx_wrapper(chartInstance->S, "ove_count", (void **)
    &chartInstance->c1__indexove_count_address, &chartInstance->c1__indexunt);
  ssGetSFcnDataStoreNameAddrIdx_wrapper(chartInstance->S, "state_idx", (void **)
    &chartInstance->c1__indexstate_idx_address, &chartInstance->c1__indexidx);
  ssGetSFcnDataStoreNameAddrIdx_wrapper(chartInstance->S, "states_space", (void **)
    &chartInstance->c1__indexstates_space_address,
    &chartInstance->c1__index_space);
  ssGetSFcnDataStoreNameAddrIdx_wrapper(chartInstance->S, "target", (void **)
    &chartInstance->c1__indextarget_address, &chartInstance->c1_j__index);
}

static void init_simulink_io_address(SFc1_RM3_copyInstanceStruct *chartInstance)
{
  chartInstance->c1_fEmlrtCtx = (void *)sfrtGetEmlrtCtx(chartInstance->S);
  chartInstance->c1_c_y = (real_T *)ssGetOutputPortSignal_wrapper
    (chartInstance->S, 1);
  chartInstance->c1_dsmdiag_B2 = (boolean_T)
    ssGetDSMBlockDiagnosticsEnabled_wrapper(chartInstance->S, 0, "B2");
  chartInstance->c1_dsmdiag_Final_Rew_matt = (boolean_T)
    ssGetDSMBlockDiagnosticsEnabled_wrapper(chartInstance->S, 1,
    "Final_Rew_matt");
  chartInstance->c1_dsmdiag_K2 = (boolean_T)
    ssGetDSMBlockDiagnosticsEnabled_wrapper(chartInstance->S, 2, "K2");
  chartInstance->c1_dsmdiag_MeanR = (boolean_T)
    ssGetDSMBlockDiagnosticsEnabled_wrapper(chartInstance->S, 3, "MeanR");
  chartInstance->c1_dsmdiag_Memory_D_size = (boolean_T)
    ssGetDSMBlockDiagnosticsEnabled_wrapper(chartInstance->S, 4, "Memory_D_size");
  chartInstance->c1_dsmdiag_Memory_Stor = (boolean_T)
    ssGetDSMBlockDiagnosticsEnabled_wrapper(chartInstance->S, 5, "Memory_Stor");
  chartInstance->c1_dsmdiag_NN_dd = (boolean_T)
    ssGetDSMBlockDiagnosticsEnabled_wrapper(chartInstance->S, 6, "NN_dd");
  chartInstance->c1_dsmdiag_Q = (boolean_T)
    ssGetDSMBlockDiagnosticsEnabled_wrapper(chartInstance->S, 7, "Q");
  chartInstance->c1_dsmdiag_Q_prime = (boolean_T)
    ssGetDSMBlockDiagnosticsEnabled_wrapper(chartInstance->S, 8, "Q_prime");
  chartInstance->c1_dsmdiag_Q_targ = (boolean_T)
    ssGetDSMBlockDiagnosticsEnabled_wrapper(chartInstance->S, 9, "Q_targ");
  chartInstance->c1_dsmdiag_R = (boolean_T)
    ssGetDSMBlockDiagnosticsEnabled_wrapper(chartInstance->S, 10, "R");
  chartInstance->c1_dsmdiag_Stor_eta = (boolean_T)
    ssGetDSMBlockDiagnosticsEnabled_wrapper(chartInstance->S, 11, "Stor_eta");
  chartInstance->c1_dsmdiag_action = (boolean_T)
    ssGetDSMBlockDiagnosticsEnabled_wrapper(chartInstance->S, 12, "action");
  chartInstance->c1_dsmdiag_num_R = (boolean_T)
    ssGetDSMBlockDiagnosticsEnabled_wrapper(chartInstance->S, 13, "num_R");
  chartInstance->c1_dsmdiag_ove_count = (boolean_T)
    ssGetDSMBlockDiagnosticsEnabled_wrapper(chartInstance->S, 14, "ove_count");
  chartInstance->c1_dsmdiag_state_idx = (boolean_T)
    ssGetDSMBlockDiagnosticsEnabled_wrapper(chartInstance->S, 15, "state_idx");
  chartInstance->c1_dsmdiag_states_space = (boolean_T)
    ssGetDSMBlockDiagnosticsEnabled_wrapper(chartInstance->S, 16, "states_space");
  chartInstance->c1_dsmdiag_target = (boolean_T)
    ssGetDSMBlockDiagnosticsEnabled_wrapper(chartInstance->S, 17, "target");
}

/* SFunction Glue Code */
#ifdef utFree
#undef utFree
#endif

#ifdef utMalloc
#undef utMalloc
#endif

#ifdef __cplusplus

extern "C" void *utMalloc(size_t size);
extern "C" void utFree(void*);

#else

extern void *utMalloc(size_t size);
extern void utFree(void*);

#endif

void sf_c1_RM3_copy_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(1259117718U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(811005451U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3184310428U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(2705518353U);
}

mxArray* sf_c1_RM3_copy_get_post_codegen_info(void);
mxArray *sf_c1_RM3_copy_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals", "postCodegenInfo" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1, 1, sizeof
    (autoinheritanceFields)/sizeof(autoinheritanceFields[0]),
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("ujv2UnreK0sQ9vfDTLD0TG");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxCreateDoubleMatrix(0,0,
                mxREAL));
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,1,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,0,mxREAL);
      double *pr = mxGetPr(mxSize);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt", "isFixedPointType" };

      mxArray *mxType = mxCreateStructMatrix(1,1,sizeof(typeFields)/sizeof
        (typeFields[0]),typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxType,0,"isFixedPointType",mxCreateDoubleScalar(0));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  {
    mxArray* mxPostCodegenInfo = sf_c1_RM3_copy_get_post_codegen_info();
    mxSetField(mxAutoinheritanceInfo,0,"postCodegenInfo",mxPostCodegenInfo);
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c1_RM3_copy_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c1_RM3_copy_jit_fallback_info(void)
{
  const char *infoFields[] = { "fallbackType", "fallbackReason",
    "hiddenFallbackType", "hiddenFallbackReason", "incompatibleSymbol" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 5, infoFields);
  mxArray *fallbackType = mxCreateString("pre");
  mxArray *fallbackReason = mxCreateString("hasBreakpoints");
  mxArray *hiddenFallbackType = mxCreateString("none");
  mxArray *hiddenFallbackReason = mxCreateString("");
  mxArray *incompatibleSymbol = mxCreateString("");
  mxSetField(mxInfo, 0, infoFields[0], fallbackType);
  mxSetField(mxInfo, 0, infoFields[1], fallbackReason);
  mxSetField(mxInfo, 0, infoFields[2], hiddenFallbackType);
  mxSetField(mxInfo, 0, infoFields[3], hiddenFallbackReason);
  mxSetField(mxInfo, 0, infoFields[4], incompatibleSymbol);
  return mxInfo;
}

mxArray *sf_c1_RM3_copy_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

mxArray* sf_c1_RM3_copy_get_post_codegen_info(void)
{
  const char* fieldNames[] = { "exportedFunctionsUsedByThisChart",
    "exportedFunctionsChecksum" };

  mwSize dims[2] = { 1, 1 };

  mxArray* mxPostCodegenInfo = mxCreateStructArray(2, dims, sizeof(fieldNames)/
    sizeof(fieldNames[0]), fieldNames);

  {
    mxArray* mxExportedFunctionsChecksum = mxCreateString("");
    mwSize exp_dims[2] = { 0, 1 };

    mxArray* mxExportedFunctionsUsedByThisChart = mxCreateCellArray(2, exp_dims);
    mxSetField(mxPostCodegenInfo, 0, "exportedFunctionsUsedByThisChart",
               mxExportedFunctionsUsedByThisChart);
    mxSetField(mxPostCodegenInfo, 0, "exportedFunctionsChecksum",
               mxExportedFunctionsChecksum);
  }

  return mxPostCodegenInfo;
}

static const mxArray *sf_get_sim_state_info_c1_RM3_copy(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[5],T\"y\",},{M[8],M[0],T\"is_active_c1_RM3_copy\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c1_RM3_copy_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc1_RM3_copyInstanceStruct *chartInstance = (SFc1_RM3_copyInstanceStruct *)
      sf_get_chart_instance_ptr(S);
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _RM3_copyMachineNumber_,
           1,
           1,
           1,
           0,
           19,
           0,
           0,
           0,
           0,
           0,
           &chartInstance->chartNumber,
           &chartInstance->instanceNumber,
           (void *)S);

        /* Each instance must initialize its own list of scripts */
        init_script_number_translation(_RM3_copyMachineNumber_,
          chartInstance->chartNumber,chartInstance->instanceNumber);
        if (chartAlreadyPresent==0) {
          /* this is the first instance */
          sf_debug_set_chart_disable_implicit_casting
            (sfGlobalDebugInstanceStruct,_RM3_copyMachineNumber_,
             chartInstance->chartNumber,1);
          sf_debug_set_chart_event_thresholds(sfGlobalDebugInstanceStruct,
            _RM3_copyMachineNumber_,
            chartInstance->chartNumber,
            0,
            0,
            0);
          _SFD_SET_DATA_PROPS(0,2,0,1,"y");
          _SFD_SET_DATA_PROPS(1,11,0,0,"B2");
          _SFD_SET_DATA_PROPS(2,11,0,0,"Final_Rew_matt");
          _SFD_SET_DATA_PROPS(3,11,0,0,"K2");
          _SFD_SET_DATA_PROPS(4,11,0,0,"MeanR");
          _SFD_SET_DATA_PROPS(5,11,0,0,"Memory_D_size");
          _SFD_SET_DATA_PROPS(6,11,0,0,"Memory_Stor");
          _SFD_SET_DATA_PROPS(7,11,0,0,"NN_dd");
          _SFD_SET_DATA_PROPS(8,11,0,0,"Q");
          _SFD_SET_DATA_PROPS(9,11,0,0,"Q_prime");
          _SFD_SET_DATA_PROPS(10,11,0,0,"Q_targ");
          _SFD_SET_DATA_PROPS(11,11,0,0,"R");
          _SFD_SET_DATA_PROPS(12,11,0,0,"Stor_eta");
          _SFD_SET_DATA_PROPS(13,11,0,0,"action");
          _SFD_SET_DATA_PROPS(14,11,0,0,"num_R");
          _SFD_SET_DATA_PROPS(15,11,0,0,"ove_count");
          _SFD_SET_DATA_PROPS(16,11,0,0,"state_idx");
          _SFD_SET_DATA_PROPS(17,11,0,0,"states_space");
          _SFD_SET_DATA_PROPS(18,11,0,0,"target");
          _SFD_STATE_INFO(0,0,2);
          _SFD_CH_SUBSTATE_COUNT(0);
          _SFD_CH_SUBSTATE_DECOMP(0);
        }

        _SFD_CV_INIT_CHART(0,0,0,0);

        {
          _SFD_CV_INIT_STATE(0,0,0,0,0,0,NULL,NULL);
        }

        _SFD_CV_INIT_TRANS(0,0,NULL,NULL,0,NULL);

        /* Initialization of MATLAB Function Model Coverage */
        _SFD_CV_INIT_EML(0,1,1,0,4,0,1,0,0,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,6993);
        _SFD_CV_INIT_EML_SATURATION(0,1,0,2111,-1,2137);
        _SFD_CV_INIT_EML_IF(0,1,0,299,316,396,495);
        _SFD_CV_INIT_EML_IF(0,1,1,531,548,606,671);
        _SFD_CV_INIT_EML_IF(0,1,2,3420,3439,5965,6415);
        _SFD_CV_INIT_EML_IF(0,1,3,3556,3582,3855,4080);
        _SFD_CV_INIT_EML_RELATIONAL(0,1,0,302,316,-1,5);
        _SFD_CV_INIT_EML_RELATIONAL(0,1,1,534,548,-1,5);
        _SFD_CV_INIT_EML_RELATIONAL(0,1,2,3423,3439,-1,4);
        _SFD_CV_INIT_EML_RELATIONAL(0,1,3,3559,3582,-1,4);
        _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)c1_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)c1_sf_marshallIn);

        {
          unsigned int dimVector[2];
          dimVector[0]= 100U;
          dimVector[1]= 10U;
          _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_u_sf_marshallOut,(MexInFcnForType)
            c1_t_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)c1_sf_marshallIn);

        {
          unsigned int dimVector[2];
          dimVector[0]= 100U;
          dimVector[1]= 100U;
          _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_t_sf_marshallOut,(MexInFcnForType)
            c1_s_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)c1_sf_marshallIn);

        {
          unsigned int dimVector[2];
          dimVector[0]= 5000U;
          dimVector[1]= 30U;
          _SFD_SET_DATA_COMPILED_PROPS(6,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_v_sf_marshallOut,(MexInFcnForType)
            c1_u_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(7,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)c1_sf_marshallIn);

        {
          unsigned int dimVector[2];
          dimVector[0]= 72U;
          dimVector[1]= 12U;
          _SFD_SET_DATA_COMPILED_PROPS(8,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_g_sf_marshallOut,(MexInFcnForType)
            c1_g_sf_marshallIn);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 72U;
          dimVector[1]= 12U;
          _SFD_SET_DATA_COMPILED_PROPS(9,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_g_sf_marshallOut,(MexInFcnForType)
            c1_g_sf_marshallIn);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 72U;
          dimVector[1]= 12U;
          _SFD_SET_DATA_COMPILED_PROPS(10,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_g_sf_marshallOut,(MexInFcnForType)
            c1_g_sf_marshallIn);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 72U;
          dimVector[1]= 50U;
          _SFD_SET_DATA_COMPILED_PROPS(11,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_s_sf_marshallOut,(MexInFcnForType)
            c1_r_sf_marshallIn);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 3000U;
          dimVector[1]= 5U;
          _SFD_SET_DATA_COMPILED_PROPS(12,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_o_sf_marshallOut,(MexInFcnForType)
            c1_n_sf_marshallIn);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 12U;
          dimVector[1]= 2U;
          _SFD_SET_DATA_COMPILED_PROPS(13,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_n_sf_marshallOut,(MexInFcnForType)
            c1_m_sf_marshallIn);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 100U;
          dimVector[1]= 100U;
          _SFD_SET_DATA_COMPILED_PROPS(14,SF_INT16,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_r_sf_marshallOut,(MexInFcnForType)
            c1_q_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(15,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)c1_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(16,SF_INT16,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_p_sf_marshallOut,(MexInFcnForType)c1_o_sf_marshallIn);

        {
          unsigned int dimVector[2];
          dimVector[0]= 72U;
          dimVector[1]= 4U;
          _SFD_SET_DATA_COMPILED_PROPS(17,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_q_sf_marshallOut,(MexInFcnForType)
            c1_p_sf_marshallIn);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 72U;
          dimVector[1]= 12U;
          _SFD_SET_DATA_COMPILED_PROPS(18,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_g_sf_marshallOut,(MexInFcnForType)
            c1_g_sf_marshallIn);
        }
      }
    } else {
      sf_debug_reset_current_state_configuration(sfGlobalDebugInstanceStruct,
        _RM3_copyMachineNumber_,chartInstance->chartNumber,
        chartInstance->instanceNumber);
    }
  }
}

static void chart_debug_initialize_data_addresses(SimStruct *S)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc1_RM3_copyInstanceStruct *chartInstance = (SFc1_RM3_copyInstanceStruct *)
      sf_get_chart_instance_ptr(S);
    if (ssIsFirstInitCond(S)) {
      /* do this only if simulation is starting and after we know the addresses of all data */
      {
        _SFD_SET_DATA_VALUE_PTR(0U, (void *)chartInstance->c1_c_y);
        _SFD_SET_DATA_VALUE_PTR(13U, (void *)
          chartInstance->c1__indexaction_address);
        _SFD_SET_DATA_VALUE_PTR(15U, (void *)
          chartInstance->c1__indexove_count_address);
        _SFD_SET_DATA_VALUE_PTR(1U, (void *)chartInstance->c1_B2_address);
        _SFD_SET_DATA_VALUE_PTR(3U, (void *)chartInstance->c1__indexK2_address);
        _SFD_SET_DATA_VALUE_PTR(12U, (void *)
          chartInstance->c1__indexStor_eta_address);
        _SFD_SET_DATA_VALUE_PTR(16U, (void *)
          chartInstance->c1__indexstate_idx_address);
        _SFD_SET_DATA_VALUE_PTR(17U, (void *)
          chartInstance->c1__indexstates_space_address);
        _SFD_SET_DATA_VALUE_PTR(14U, (void *)
          chartInstance->c1__indexnum_R_address);
        _SFD_SET_DATA_VALUE_PTR(11U, (void *)chartInstance->c1__indexR_address);
        _SFD_SET_DATA_VALUE_PTR(4U, (void *)
          chartInstance->c1__indexMeanR_address);
        _SFD_SET_DATA_VALUE_PTR(2U, (void *)
          chartInstance->c1__indexFinal_Rew_matt_address);
        _SFD_SET_DATA_VALUE_PTR(7U, (void *)
          chartInstance->c1__indexNN_dd_address);
        _SFD_SET_DATA_VALUE_PTR(5U, (void *)
          chartInstance->c1__indexMemory_D_size_address);
        _SFD_SET_DATA_VALUE_PTR(6U, (void *)
          chartInstance->c1__indexMemory_Stor_address);
        _SFD_SET_DATA_VALUE_PTR(9U, (void *)
          chartInstance->c1__indexQ_prime_address);
        _SFD_SET_DATA_VALUE_PTR(10U, (void *)
          chartInstance->c1__indexQ_targ_address);
        _SFD_SET_DATA_VALUE_PTR(8U, (void *)chartInstance->c1__indexQ_address);
        _SFD_SET_DATA_VALUE_PTR(18U, (void *)
          chartInstance->c1__indextarget_address);
      }
    }
  }
}

static const char* sf_get_instance_specialization(void)
{
  return "sCk8R4g2tLvMpXpmpGG9evF";
}

static void sf_opaque_initialize_c1_RM3_copy(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc1_RM3_copyInstanceStruct*) chartInstanceVar)
    ->S,0);
  initialize_params_c1_RM3_copy((SFc1_RM3_copyInstanceStruct*) chartInstanceVar);
  initialize_c1_RM3_copy((SFc1_RM3_copyInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c1_RM3_copy(void *chartInstanceVar)
{
  enable_c1_RM3_copy((SFc1_RM3_copyInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c1_RM3_copy(void *chartInstanceVar)
{
  disable_c1_RM3_copy((SFc1_RM3_copyInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c1_RM3_copy(void *chartInstanceVar)
{
  sf_gateway_c1_RM3_copy((SFc1_RM3_copyInstanceStruct*) chartInstanceVar);
}

static const mxArray* sf_opaque_get_sim_state_c1_RM3_copy(SimStruct* S)
{
  return get_sim_state_c1_RM3_copy((SFc1_RM3_copyInstanceStruct *)
    sf_get_chart_instance_ptr(S));     /* raw sim ctx */
}

static void sf_opaque_set_sim_state_c1_RM3_copy(SimStruct* S, const mxArray *st)
{
  set_sim_state_c1_RM3_copy((SFc1_RM3_copyInstanceStruct*)
    sf_get_chart_instance_ptr(S), st);
}

static void sf_opaque_terminate_c1_RM3_copy(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc1_RM3_copyInstanceStruct*) chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_RM3_copy_optimization_info();
    }

    finalize_c1_RM3_copy((SFc1_RM3_copyInstanceStruct*) chartInstanceVar);
    utFree(chartInstanceVar);
    if (ssGetUserData(S)!= NULL) {
      sf_free_ChartRunTimeInfo(S);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc1_RM3_copy((SFc1_RM3_copyInstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c1_RM3_copy(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c1_RM3_copy((SFc1_RM3_copyInstanceStruct*)
      sf_get_chart_instance_ptr(S));
  }
}

static void mdlSetWorkWidths_c1_RM3_copy(SimStruct *S)
{
  /* Set overwritable ports for inplace optimization */
  ssSetInputPortDirectFeedThrough(S, 0, 1);
  ssSetStatesModifiedOnlyInUpdate(S, 0);
  ssSetBlockIsPurelyCombinatorial_wrapper(S, 0);
  ssMdlUpdateIsEmpty(S, 1);
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_RM3_copy_optimization_info(sim_mode_is_rtw_gen(S),
      sim_mode_is_modelref_sim(S), sim_mode_is_external(S));
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,1);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,1);
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,1,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_set_chart_accesses_machine_info(S, sf_get_instance_specialization(),
      infoStruct, 1);
    sf_update_buildInfo(S, sf_get_instance_specialization(),infoStruct,1);
    if (chartIsInlinable) {
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,1,1);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=1; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,1);
    sf_register_codegen_names_for_scoped_functions_defined_by_chart(S);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(3065903991U));
  ssSetChecksum1(S,(2733176318U));
  ssSetChecksum2(S,(2694920348U));
  ssSetChecksum3(S,(3674888359U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSetStateSemanticsClassicAndSynchronous(S, true);
  ssSupportsMultipleExecInstances(S,0);
}

static void mdlRTW_c1_RM3_copy(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c1_RM3_copy(SimStruct *S)
{
  SFc1_RM3_copyInstanceStruct *chartInstance;
  chartInstance = (SFc1_RM3_copyInstanceStruct *)utMalloc(sizeof
    (SFc1_RM3_copyInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  memset(chartInstance, 0, sizeof(SFc1_RM3_copyInstanceStruct));
  chartInstance->chartInfo.chartInstance = chartInstance;
  if (ssGetSampleTime(S, 0) == CONTINUOUS_SAMPLE_TIME && ssGetOffsetTime(S, 0) ==
      0 && ssGetNumContStates(ssGetRootSS(S)) > 0) {
    sf_error_out_about_continuous_sample_time_with_persistent_vars(S);
  }

  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c1_RM3_copy;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c1_RM3_copy;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c1_RM3_copy;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c1_RM3_copy;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c1_RM3_copy;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c1_RM3_copy;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c1_RM3_copy;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c1_RM3_copy;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c1_RM3_copy;
  chartInstance->chartInfo.mdlStart = mdlStart_c1_RM3_copy;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c1_RM3_copy;
  chartInstance->chartInfo.callGetHoverDataForMsg = NULL;
  chartInstance->chartInfo.extModeExec = NULL;
  chartInstance->chartInfo.restoreLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.restoreBeforeLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.storeCurrentConfiguration = NULL;
  chartInstance->chartInfo.callAtomicSubchartUserFcn = NULL;
  chartInstance->chartInfo.callAtomicSubchartAutoFcn = NULL;
  chartInstance->chartInfo.debugInstance = sfGlobalDebugInstanceStruct;
  chartInstance->S = S;
  sf_init_ChartRunTimeInfo(S, &(chartInstance->chartInfo), false, 0);
  init_dsm_address_info(chartInstance);
  init_simulink_io_address(chartInstance);
  if (!sim_mode_is_rtw_gen(S)) {
  }

  chart_debug_initialization(S,1);
  mdl_start_c1_RM3_copy(chartInstance);
}

void c1_RM3_copy_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c1_RM3_copy(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c1_RM3_copy(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c1_RM3_copy(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c1_RM3_copy_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
