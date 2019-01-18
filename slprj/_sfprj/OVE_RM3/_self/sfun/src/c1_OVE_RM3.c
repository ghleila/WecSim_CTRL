/* Include files */

#include "OVE_RM3_sfun.h"
#include "c1_OVE_RM3.h"
#include "mwmathutil.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "OVE_RM3_sfun_debug_macros.h"
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
static const char * c1_debug_family_names[20] = { "PptoAvg", "gamma", "opts",
  "alpha", "current_action", "action_idx", "Height", "next_stPrimee_idx",
  "st_chosen", "Height_prime", "Ts_prime", "next_state_idx", "h", "nargin",
  "nargout", "power", "H", "Ts", "B3", "K3" };

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

static emlrtRTEInfo c1_c_emlrtRTEI = { 110,/* lineNo */
  10,                                  /* colNo */
  "MATLAB Function1",                  /* fName */
  "#OVE_RM3:3"                         /* pName */
};

static emlrtRTEInfo c1_d_emlrtRTEI = { 137,/* lineNo */
  9,                                   /* colNo */
  "MATLAB Function1",                  /* fName */
  "#OVE_RM3:3"                         /* pName */
};

static emlrtRTEInfo c1_e_emlrtRTEI = { 165,/* lineNo */
  9,                                   /* colNo */
  "MATLAB Function1",                  /* fName */
  "#OVE_RM3:3"                         /* pName */
};

/* Function Declarations */
static void initialize_c1_OVE_RM3(SFc1_OVE_RM3InstanceStruct *chartInstance);
static void initialize_params_c1_OVE_RM3(SFc1_OVE_RM3InstanceStruct
  *chartInstance);
static void enable_c1_OVE_RM3(SFc1_OVE_RM3InstanceStruct *chartInstance);
static void disable_c1_OVE_RM3(SFc1_OVE_RM3InstanceStruct *chartInstance);
static void c1_update_debugger_state_c1_OVE_RM3(SFc1_OVE_RM3InstanceStruct
  *chartInstance);
static const mxArray *get_sim_state_c1_OVE_RM3(SFc1_OVE_RM3InstanceStruct
  *chartInstance);
static void set_sim_state_c1_OVE_RM3(SFc1_OVE_RM3InstanceStruct *chartInstance,
  const mxArray *c1_st);
static void finalize_c1_OVE_RM3(SFc1_OVE_RM3InstanceStruct *chartInstance);
static void sf_gateway_c1_OVE_RM3(SFc1_OVE_RM3InstanceStruct *chartInstance);
static void mdl_start_c1_OVE_RM3(SFc1_OVE_RM3InstanceStruct *chartInstance);
static void c1_chartstep_c1_OVE_RM3(SFc1_OVE_RM3InstanceStruct *chartInstance);
static void initSimStructsc1_OVE_RM3(SFc1_OVE_RM3InstanceStruct *chartInstance);
static void init_script_number_translation(uint32_T c1_machineNumber, uint32_T
  c1_chartNumber, uint32_T c1_instanceNumber);
static const mxArray *c1_sf_marshallOut(void *chartInstanceVoid, void *c1_inData);
static void c1_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_b_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static const mxArray *c1_c_sf_marshallOut(void *chartInstanceVoid, real_T
  c1_inData_data[], int32_T c1_inData_size[1]);
static void c1_emlrt_marshallIn(SFc1_OVE_RM3InstanceStruct *chartInstance, const
  mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y_data[],
  int32_T c1_y_size[1]);
static void c1_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, real_T c1_outData_data[], int32_T
  c1_outData_size[1]);
static const mxArray *c1_d_sf_marshallOut(void *chartInstanceVoid, real_T
  c1_inData_data[], int32_T c1_inData_size[1]);
static void c1_b_emlrt_marshallIn(SFc1_OVE_RM3InstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T
  c1_y_data[], int32_T c1_y_size[1]);
static void c1_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, real_T c1_outData_data[], int32_T
  c1_outData_size[1]);
static const mxArray *c1_e_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static const mxArray *c1_f_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static c1_sI3an6DNyfQAOCkXY61B1nE c1_c_emlrt_marshallIn
  (SFc1_OVE_RM3InstanceStruct *chartInstance, const mxArray *c1_b_u, const
   emlrtMsgIdentifier *c1_parentId);
static void c1_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_g_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_d_emlrt_marshallIn(SFc1_OVE_RM3InstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_c_y[6]);
static void c1_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static real_T c1_e_emlrt_marshallIn(SFc1_OVE_RM3InstanceStruct *chartInstance,
  const mxArray *c1_coder_internal_mxSubscript, const char_T *c1_identifier);
static real_T c1_f_emlrt_marshallIn(SFc1_OVE_RM3InstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId);
static const mxArray *c1_h_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static int32_T c1_g_emlrt_marshallIn(SFc1_OVE_RM3InstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_i_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_h_emlrt_marshallIn(SFc1_OVE_RM3InstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_c_y[8]);
static void c1_g_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_j_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_i_emlrt_marshallIn(SFc1_OVE_RM3InstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_c_y[2]);
static void c1_h_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_k_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_j_emlrt_marshallIn(SFc1_OVE_RM3InstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_c_y
  [20000]);
static void c1_i_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_l_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_k_emlrt_marshallIn(SFc1_OVE_RM3InstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_c_y
  [288]);
static void c1_j_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_m_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_l_emlrt_marshallIn(SFc1_OVE_RM3InstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_c_y[72]);
static void c1_k_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_n_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_m_emlrt_marshallIn(SFc1_OVE_RM3InstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_c_y
  [7200]);
static void c1_l_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_o_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_n_emlrt_marshallIn(SFc1_OVE_RM3InstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_c_y[72]);
static void c1_m_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_p_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_o_emlrt_marshallIn(SFc1_OVE_RM3InstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_c_y
  [48000]);
static void c1_n_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static uint8_T c1_p_emlrt_marshallIn(SFc1_OVE_RM3InstanceStruct *chartInstance,
  const mxArray *c1_b_is_active_c1_OVE_RM3, const char_T *c1_identifier);
static uint8_T c1_q_emlrt_marshallIn(SFc1_OVE_RM3InstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId);
static int32_T c1__s32_d_(SFc1_OVE_RM3InstanceStruct *chartInstance, real_T c1_b,
  uint32_T c1_ssid_src_loc, int32_T c1_offset_src_loc, int32_T c1_length_src_loc);
static real_T c1_get_B2(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex);
static void c1_set_B2(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
                      c1_elementIndex, real_T c1_elementValue);
static real_T *c1_access_B2(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_rdOnly);
static real_T c1_get_Final_Rew_matt(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_elementIndex);
static void c1_set_Final_Rew_matt(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_elementIndex, real_T c1_elementValue);
static real_T *c1_access_Final_Rew_matt(SFc1_OVE_RM3InstanceStruct
  *chartInstance, uint32_T c1_rdOnly);
static real_T c1_get_K2(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex);
static void c1_set_K2(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
                      c1_elementIndex, real_T c1_elementValue);
static real_T *c1_access_K2(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_rdOnly);
static real_T c1_get_MeanR(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex);
static void c1_set_MeanR(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex, real_T c1_elementValue);
static real_T *c1_access_MeanR(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_rdOnly);
static real_T c1_get_Memory_D_size(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_elementIndex);
static void c1_set_Memory_D_size(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_elementIndex, real_T c1_elementValue);
static real_T *c1_access_Memory_D_size(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_rdOnly);
static real_T c1_get_Memory_Stor(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_elementIndex);
static void c1_set_Memory_Stor(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_elementIndex, real_T c1_elementValue);
static real_T *c1_access_Memory_Stor(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_rdOnly);
static real_T c1_get_NN_dd(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex);
static void c1_set_NN_dd(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex, real_T c1_elementValue);
static real_T *c1_access_NN_dd(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_rdOnly);
static real_T c1_get_Q(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex);
static void c1_set_Q(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
                     c1_elementIndex, real_T c1_elementValue);
static real_T *c1_access_Q(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_rdOnly);
static real_T c1_get_Q_prime(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex);
static void c1_set_Q_prime(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex, real_T c1_elementValue);
static real_T *c1_access_Q_prime(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_rdOnly);
static real_T c1_get_Q_targ(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex);
static void c1_set_Q_targ(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex, real_T c1_elementValue);
static real_T *c1_access_Q_targ(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_rdOnly);
static real_T c1_get_R(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex);
static void c1_set_R(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
                     c1_elementIndex, real_T c1_elementValue);
static real_T *c1_access_R(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_rdOnly);
static real_T c1_get_Stor_eta(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_elementIndex);
static void c1_set_Stor_eta(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex, real_T c1_elementValue);
static real_T *c1_access_Stor_eta(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_rdOnly);
static real_T c1_get_action(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex);
static void c1_set_action(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex, real_T c1_elementValue);
static real_T *c1_access_action(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_rdOnly);
static real_T c1_get_num_R(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex);
static void c1_set_num_R(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex, real_T c1_elementValue);
static real_T *c1_access_num_R(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_rdOnly);
static real_T c1_get_ove_count(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_elementIndex);
static void c1_set_ove_count(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex, real_T c1_elementValue);
static real_T *c1_access_ove_count(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_rdOnly);
static real_T c1_get_state_idx(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_elementIndex);
static void c1_set_state_idx(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex, real_T c1_elementValue);
static real_T *c1_access_state_idx(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_rdOnly);
static real_T c1_get_states_space(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_elementIndex);
static void c1_set_states_space(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_elementIndex, real_T c1_elementValue);
static real_T *c1_access_states_space(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_rdOnly);
static real_T c1_get_xp(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex);
static void c1_set_xp(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
                      c1_elementIndex, real_T c1_elementValue);
static real_T *c1_access_xp(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_rdOnly);
static void init_dsm_address_info(SFc1_OVE_RM3InstanceStruct *chartInstance);
static void init_simulink_io_address(SFc1_OVE_RM3InstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c1_OVE_RM3(SFc1_OVE_RM3InstanceStruct *chartInstance)
{
  if (sf_is_first_init_cond(chartInstance->S)) {
    initSimStructsc1_OVE_RM3(chartInstance);
    chart_debug_initialize_data_addresses(chartInstance->S);
  }

  chartInstance->c1_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c1_is_active_c1_OVE_RM3 = 0U;
}

static void initialize_params_c1_OVE_RM3(SFc1_OVE_RM3InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void enable_c1_OVE_RM3(SFc1_OVE_RM3InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c1_OVE_RM3(SFc1_OVE_RM3InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c1_update_debugger_state_c1_OVE_RM3(SFc1_OVE_RM3InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c1_OVE_RM3(SFc1_OVE_RM3InstanceStruct
  *chartInstance)
{
  const mxArray *c1_st;
  const mxArray *c1_c_y = NULL;
  real_T c1_hoistedGlobal;
  const mxArray *c1_d_y = NULL;
  real_T c1_b_hoistedGlobal;
  const mxArray *c1_e_y = NULL;
  real_T c1_c_hoistedGlobal;
  const mxArray *c1_f_y = NULL;
  real_T c1_d_hoistedGlobal;
  const mxArray *c1_g_y = NULL;
  uint8_T c1_e_hoistedGlobal;
  const mxArray *c1_h_y = NULL;
  c1_st = NULL;
  c1_st = NULL;
  c1_c_y = NULL;
  sf_mex_assign(&c1_c_y, sf_mex_createcellmatrix(5, 1), false);
  c1_hoistedGlobal = *chartInstance->c1_B3;
  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", &c1_hoistedGlobal, 0, 0U, 0U, 0U, 0),
                false);
  sf_mex_setcell(c1_c_y, 0, c1_d_y);
  c1_b_hoistedGlobal = *chartInstance->c1_H;
  c1_e_y = NULL;
  sf_mex_assign(&c1_e_y, sf_mex_create("y", &c1_b_hoistedGlobal, 0, 0U, 0U, 0U,
    0), false);
  sf_mex_setcell(c1_c_y, 1, c1_e_y);
  c1_c_hoistedGlobal = *chartInstance->c1_K3;
  c1_f_y = NULL;
  sf_mex_assign(&c1_f_y, sf_mex_create("y", &c1_c_hoistedGlobal, 0, 0U, 0U, 0U,
    0), false);
  sf_mex_setcell(c1_c_y, 2, c1_f_y);
  c1_d_hoistedGlobal = *chartInstance->c1_Ts;
  c1_g_y = NULL;
  sf_mex_assign(&c1_g_y, sf_mex_create("y", &c1_d_hoistedGlobal, 0, 0U, 0U, 0U,
    0), false);
  sf_mex_setcell(c1_c_y, 3, c1_g_y);
  c1_e_hoistedGlobal = chartInstance->c1_is_active_c1_OVE_RM3;
  c1_h_y = NULL;
  sf_mex_assign(&c1_h_y, sf_mex_create("y", &c1_e_hoistedGlobal, 3, 0U, 0U, 0U,
    0), false);
  sf_mex_setcell(c1_c_y, 4, c1_h_y);
  sf_mex_assign(&c1_st, c1_c_y, false);
  return c1_st;
}

static void set_sim_state_c1_OVE_RM3(SFc1_OVE_RM3InstanceStruct *chartInstance,
  const mxArray *c1_st)
{
  const mxArray *c1_b_u;
  chartInstance->c1_doneDoubleBufferReInit = true;
  c1_b_u = sf_mex_dup(c1_st);
  *chartInstance->c1_B3 = c1_e_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getcell(c1_b_u, 0)), "B3");
  *chartInstance->c1_H = c1_e_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getcell(c1_b_u, 1)), "H");
  *chartInstance->c1_K3 = c1_e_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getcell(c1_b_u, 2)), "K3");
  *chartInstance->c1_Ts = c1_e_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getcell(c1_b_u, 3)), "Ts");
  chartInstance->c1_is_active_c1_OVE_RM3 = c1_p_emlrt_marshallIn(chartInstance,
    sf_mex_dup(sf_mex_getcell(c1_b_u, 4)), "is_active_c1_OVE_RM3");
  sf_mex_destroy(&c1_b_u);
  c1_update_debugger_state_c1_OVE_RM3(chartInstance);
  sf_mex_destroy(&c1_st);
}

static void finalize_c1_OVE_RM3(SFc1_OVE_RM3InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c1_OVE_RM3(SFc1_OVE_RM3InstanceStruct *chartInstance)
{
  int32_T c1_i0;
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 0U, chartInstance->c1_sfEvent);
  for (c1_i0 = 0; c1_i0 < 6; c1_i0++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_power)[c1_i0], 0U);
  }

  chartInstance->c1_sfEvent = CALL_EVENT;
  c1_chartstep_c1_OVE_RM3(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_OVE_RM3MachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  _SFD_DATA_RANGE_CHECK(*chartInstance->c1_H, 1U);
  _SFD_DATA_RANGE_CHECK(*chartInstance->c1_Ts, 2U);
  _SFD_DATA_RANGE_CHECK(*chartInstance->c1_B3, 3U);
  _SFD_DATA_RANGE_CHECK(*chartInstance->c1_K3, 4U);
}

static void mdl_start_c1_OVE_RM3(SFc1_OVE_RM3InstanceStruct *chartInstance)
{
  sim_mode_is_external(chartInstance->S);
}

static void c1_chartstep_c1_OVE_RM3(SFc1_OVE_RM3InstanceStruct *chartInstance)
{
  int32_T c1_i1;
  uint32_T c1_debug_family_var_map[20];
  real_T c1_b_power[6];
  real_T c1_PptoAvg[6];
  real_T c1_gamma;
  c1_sI3an6DNyfQAOCkXY61B1nE c1_opts;
  real_T c1_alpha;
  const mxArray *c1_current_action = NULL;
  real_T c1_action_idx_data[4];
  int32_T c1_action_idx_size[1];
  real_T c1_Height;
  real_T c1_next_stPrimee_idx_data[72];
  int32_T c1_next_stPrimee_idx_size[1];
  real_T c1_st_chosen;
  real_T c1_Height_prime;
  real_T c1_Ts_prime;
  real_T c1_next_state_idx_data[72];
  int32_T c1_next_state_idx_size[1];
  real_T c1_h;
  real_T c1_action_idx;
  real_T c1_nargin = 1.0;
  real_T c1_nargout = 4.0;
  real_T c1_b_H;
  real_T c1_b_Ts;
  real_T c1_b_B3;
  real_T c1_b_K3;
  int32_T c1_i2;
  const mxArray *c1_c_y = NULL;
  static char_T c1_cv0[11] = { 'i', 't', 'e', 'r', 'a', 't', 'i', 'o', 'n', ':',
    ' ' };

  real_T c1_hoistedGlobal;
  const mxArray *c1_d_y = NULL;
  const mxArray *c1_e_y = NULL;
  real_T c1_b_u;
  const mxArray *c1_f_y = NULL;
  real_T c1_c_u;
  const mxArray *c1_g_y = NULL;
  real_T c1_d_u;
  const mxArray *c1_h_y = NULL;
  real_T c1_e_u;
  const mxArray *c1_i_y = NULL;
  real_T c1_f_u;
  const mxArray *c1_j_y = NULL;
  int32_T c1_i3;
  int32_T c1_idx;
  boolean_T c1_x[4];
  int32_T c1_ii_size[1];
  int32_T c1_ii;
  int32_T c1_i4;
  int32_T c1_ii_data[4];
  int32_T c1_loop_ub;
  int32_T c1_i5;
  int32_T c1_b_action_idx[2];
  int32_T c1_i6;
  int32_T c1_i7;
  boolean_T c1_b_x[72];
  int32_T c1_i8;
  boolean_T c1_bv0[72];
  int32_T c1_b_idx;
  int32_T c1_b_ii_size[1];
  int32_T c1_b_ii;
  int32_T c1_i9;
  int32_T c1_b_ii_data[72];
  int32_T c1_b_loop_ub;
  int32_T c1_i10;
  const mxArray *c1_k_y = NULL;
  static real_T c1_dv0[72] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
    11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0,
    24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0,
    37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0,
    50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0, 60.0, 61.0, 62.0,
    63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70.0, 71.0, 72.0 };

  real_T c1_g_u;
  const mxArray *c1_l_y = NULL;
  int32_T c1_i11;
  int32_T c1_i12;
  int32_T c1_i13;
  int32_T c1_c_idx;
  int32_T c1_c_ii;
  int32_T c1_i14;
  int32_T c1_c_loop_ub;
  int32_T c1_i15;
  boolean_T exitg1;
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 0U, chartInstance->c1_sfEvent);
  for (c1_i1 = 0; c1_i1 < 6; c1_i1++) {
    c1_b_power[c1_i1] = (*chartInstance->c1_power)[c1_i1];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 20U, 21U, c1_debug_family_names,
    c1_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_PptoAvg, 0U, c1_g_sf_marshallOut,
    c1_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_gamma, 1U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_opts, 2U, c1_f_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_alpha, 3U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c1_current_action, 4U, c1_e_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_DYN_EMX_IMPORTABLE((void *)&c1_action_idx_data, (
    const int32_T *)&c1_action_idx_size, NULL, 0, -1, (void *)
    c1_d_sf_marshallOut, (void *)c1_c_sf_marshallIn, (void *)&c1_action_idx_data,
    false);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_Height, 6U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_DYN_EMX_IMPORTABLE((void *)
    &c1_next_stPrimee_idx_data, (const int32_T *)&c1_next_stPrimee_idx_size,
    NULL, 0, 7, (void *)c1_c_sf_marshallOut, (void *)c1_b_sf_marshallIn, (void *)
    &c1_next_stPrimee_idx_data, false);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_st_chosen, 8U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_Height_prime, 9U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_Ts_prime, 10U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_DYN_EMX_IMPORTABLE((void *)&c1_next_state_idx_data,
    (const int32_T *)&c1_next_state_idx_size, NULL, 0, 11, (void *)
    c1_c_sf_marshallOut, (void *)c1_b_sf_marshallIn, (void *)
    &c1_next_state_idx_data, false);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_h, 12U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_action_idx, MAX_uint32_T,
    c1_sf_marshallOut, c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargin, 13U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargout, 14U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_b_power, 15U, c1_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_b_H, 16U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_b_Ts, 17U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_b_B3, 18U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_b_K3, 19U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 3);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 12);
  for (c1_i2 = 0; c1_i2 < 6; c1_i2++) {
    c1_PptoAvg[c1_i2] = c1_b_power[c1_i2];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 29);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 33);
  c1_c_y = NULL;
  sf_mex_assign(&c1_c_y, sf_mex_create("y", c1_cv0, 10, 0U, 1U, 0U, 2, 1, 11),
                false);
  c1_hoistedGlobal = c1_get_ove_count(chartInstance, 0);
  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", &c1_hoistedGlobal, 0, 0U, 0U, 0U, 0),
                false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "disp", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "horzcat", 1U,
    2U, 14, c1_c_y, 14, sf_mex_call_debug(sfGlobalDebugInstanceStruct, "num2str",
    1U, 1U, 14, c1_d_y)));
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 44);
  if (CV_EML_IF(0, 1, 0, CV_RELATIONAL_EVAL(4U, 0U, 0, c1_get_ove_count
        (chartInstance, 0), 400.0, -1, 5U, c1_get_ove_count(chartInstance, 0) >=
        400.0))) {
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 46);
    c1_gamma = 0.6;
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 47);
    c1_opts.StepRatio = 0.01;
  } else {
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 50);
    c1_gamma = 0.4;
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 51);
    c1_opts.StepRatio = 0.1;
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 57);
  if (CV_EML_IF(0, 1, 1, CV_RELATIONAL_EVAL(4U, 0U, 1, c1_get_ove_count
        (chartInstance, 0), 200.0, -1, 5U, c1_get_ove_count(chartInstance, 0) >=
        200.0))) {
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 59);
    c1_alpha = 0.6;
  } else {
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 63);
    c1_alpha = 0.4;
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 104);
  c1_e_y = NULL;
  sf_mex_assign(&c1_e_y, sf_mex_create("y", c1_access_action(chartInstance, true),
    0, 0U, 1U, 0U, 2, 4, 2), false);
  c1_b_u = 1.0;
  c1_f_y = NULL;
  sf_mex_assign(&c1_f_y, sf_mex_create("y", &c1_b_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c1_current_action, sf_mex_call_debug
                (sfGlobalDebugInstanceStruct, "datasample", 1U, 2U, 14, c1_e_y,
                 14, c1_f_y), false);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 106);
  c1_c_u = 1.0;
  c1_g_y = NULL;
  sf_mex_assign(&c1_g_y, sf_mex_create("y", &c1_c_u, 0, 0U, 0U, 0U, 0), false);
  c1_d_u = 1.0;
  c1_h_y = NULL;
  sf_mex_assign(&c1_h_y, sf_mex_create("y", &c1_d_u, 0, 0U, 0U, 0U, 0), false);
  c1_set_B2(chartInstance, 0, c1_e_emlrt_marshallIn(chartInstance,
             sf_mex_call_debug(sfGlobalDebugInstanceStruct,
              "coder.internal.mxSubscript", 1U, 3U, 14, sf_mex_dup
              (c1_current_action), 14, c1_g_y, 14, c1_h_y),
             "coder.internal.mxSubscript"));
  ssUpdateDataStoreLog_wrapper(chartInstance->S, 0);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 108);
  c1_e_u = 1.0;
  c1_i_y = NULL;
  sf_mex_assign(&c1_i_y, sf_mex_create("y", &c1_e_u, 0, 0U, 0U, 0U, 0), false);
  c1_f_u = 2.0;
  c1_j_y = NULL;
  sf_mex_assign(&c1_j_y, sf_mex_create("y", &c1_f_u, 0, 0U, 0U, 0U, 0), false);
  c1_set_K2(chartInstance, 0, c1_e_emlrt_marshallIn(chartInstance,
             sf_mex_call_debug(sfGlobalDebugInstanceStruct,
              "coder.internal.mxSubscript", 1U, 3U, 14, sf_mex_dup
              (c1_current_action), 14, c1_i_y, 14, c1_j_y),
             "coder.internal.mxSubscript"));
  ssUpdateDataStoreLog_wrapper(chartInstance->S, 2);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 110);
  for (c1_i3 = 0; c1_i3 < 4; c1_i3++) {
    c1_x[c1_i3] = ((c1_get_action(chartInstance, c1_i3) == c1_get_B2
                    (chartInstance, 0)) && (c1_get_action(chartInstance, c1_i3 +
      4) == c1_get_K2(chartInstance, 0)));
  }

  c1_idx = 0;
  c1_ii_size[0] = 4;
  c1_ii = 1;
  exitg1 = false;
  while ((!exitg1) && (c1_ii < 5)) {
    if (c1_x[c1_ii - 1]) {
      c1_idx++;
      c1_ii_data[c1_idx - 1] = c1_ii;
      if (c1_idx >= 4) {
        exitg1 = true;
      } else {
        c1_ii++;
      }
    } else {
      c1_ii++;
    }
  }

  if (1 > c1_idx) {
    c1_i4 = 0;
  } else {
    c1_i4 = c1_idx;
  }

  c1_ii_size[0] = c1_i4;
  c1_action_idx_size[0] = c1_ii_size[0];
  c1_loop_ub = c1_ii_size[0] - 1;
  for (c1_i5 = 0; c1_i5 <= c1_loop_ub; c1_i5++) {
    c1_action_idx_data[c1_i5] = (real_T)c1_ii_data[c1_i5];
  }

  _SFD_SYMBOL_SWITCH(5U, 5U);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 115);
  c1_set_Stor_eta(chartInstance, sf_eml_array_bounds_check
                  (sfGlobalDebugInstanceStruct, chartInstance->S, 1U, 2721, 9,
                   16U, (int32_T)sf_integer_check(chartInstance->S, 1U, 2721U,
    9U, c1_get_ove_count(chartInstance, 0)), 1, 4000) - 1, c1_get_state_idx
                  (chartInstance, 0));
  ssUpdateDataStoreLog_wrapper(chartInstance->S, 11);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 119);
  (real_T)sf_eml_array_bounds_check(sfGlobalDebugInstanceStruct,
    chartInstance->S, 1U, 2844, 1, MAX_uint32_T, 1, 1, c1_action_idx_size[0]);
  c1_b_action_idx[0] = c1_action_idx_size[0];
  c1_b_action_idx[1] = 1;
  c1_action_idx = c1_action_idx_data[0];
  _SFD_SYMBOL_SWITCH(5U, 13U);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 121);
  c1_set_Stor_eta(chartInstance, sf_eml_array_bounds_check
                  (sfGlobalDebugInstanceStruct, chartInstance->S, 1U, 2874, 9,
                   16U, (int32_T)sf_integer_check(chartInstance->S, 1U, 2874U,
    9U, c1_get_ove_count(chartInstance, 0)), 1, 4000) + 3999, c1_action_idx);
  ssUpdateDataStoreLog_wrapper(chartInstance->S, 11);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 131U);
  c1_Height = c1_get_states_space(chartInstance, sf_eml_array_bounds_check
    (sfGlobalDebugInstanceStruct, chartInstance->S, 1U, 3098, 26, 21U, (int32_T)
     sf_integer_check(chartInstance->S, 1U, 3098U, 26U, c1_get_state_idx
                      (chartInstance, 0)), 1, 72) - 1);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 133U);
  c1_b_Ts = c1_get_states_space(chartInstance, sf_eml_array_bounds_check
    (sfGlobalDebugInstanceStruct, chartInstance->S, 1U, 3192, 26, 21U, (int32_T)
     sf_integer_check(chartInstance->S, 1U, 3192U, 26U, c1_get_state_idx
                      (chartInstance, 0)), 1, 72) + 71);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 137U);
  for (c1_i6 = 0; c1_i6 < 72; c1_i6++) {
    c1_b_x[c1_i6] = (c1_get_states_space(chartInstance, c1_i6) == c1_Height);
  }

  for (c1_i7 = 0; c1_i7 < 72; c1_i7++) {
    c1_bv0[c1_i7] = (c1_get_states_space(chartInstance, c1_i7 + 72) == c1_b_Ts);
  }

  for (c1_i8 = 0; c1_i8 < 72; c1_i8++) {
    c1_b_x[c1_i8] = (c1_b_x[c1_i8] && c1_bv0[c1_i8] && (c1_get_states_space
      (chartInstance, c1_i8 + 144) == c1_get_B2(chartInstance, 0)) &&
                     (c1_get_states_space(chartInstance, c1_i8 + 216) ==
                      c1_get_K2(chartInstance, 0)));
  }

  c1_b_idx = 0;
  c1_b_ii_size[0] = 72;
  c1_b_ii = 1;
  exitg1 = false;
  while ((!exitg1) && (c1_b_ii < 73)) {
    if (c1_b_x[c1_b_ii - 1]) {
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
    c1_i9 = 0;
  } else {
    c1_i9 = c1_b_idx;
  }

  c1_b_ii_size[0] = c1_i9;
  c1_next_stPrimee_idx_size[0] = c1_b_ii_size[0];
  c1_b_loop_ub = c1_b_ii_size[0] - 1;
  for (c1_i10 = 0; c1_i10 <= c1_b_loop_ub; c1_i10++) {
    c1_next_stPrimee_idx_data[c1_i10] = (real_T)c1_b_ii_data[c1_i10];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 158U);
  c1_st_chosen = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 159U);
  c1_k_y = NULL;
  sf_mex_assign(&c1_k_y, sf_mex_create("y", c1_dv0, 0, 0U, 1U, 0U, 2, 1, 72),
                false);
  c1_g_u = 1.0;
  c1_l_y = NULL;
  sf_mex_assign(&c1_l_y, sf_mex_create("y", &c1_g_u, 0, 0U, 0U, 0U, 0), false);
  c1_st_chosen = c1_e_emlrt_marshallIn(chartInstance, sf_mex_call_debug
    (sfGlobalDebugInstanceStruct, "datasample", 1U, 2U, 14, c1_k_y, 14, c1_l_y),
    "datasample");
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 162U);
  c1_Height_prime = c1_get_states_space(chartInstance, sf_eml_array_bounds_check
    (sfGlobalDebugInstanceStruct, chartInstance->S, 1U, 4112, 26, 21U, (int32_T)
     sf_integer_check(chartInstance->S, 1U, 4112U, 26U, c1_st_chosen), 1, 72) -
    1);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 163U);
  c1_Ts_prime = c1_get_states_space(chartInstance, c1__s32_d_(chartInstance,
    c1_st_chosen, 1U, 4203U, 26U) + 71);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 165U);
  for (c1_i11 = 0; c1_i11 < 72; c1_i11++) {
    c1_b_x[c1_i11] = (c1_get_states_space(chartInstance, c1_i11) ==
                      c1_Height_prime);
  }

  for (c1_i12 = 0; c1_i12 < 72; c1_i12++) {
    c1_bv0[c1_i12] = (c1_get_states_space(chartInstance, c1_i12 + 72) ==
                      c1_Ts_prime);
  }

  for (c1_i13 = 0; c1_i13 < 72; c1_i13++) {
    c1_b_x[c1_i13] = (c1_b_x[c1_i13] && c1_bv0[c1_i13] && (c1_get_states_space
      (chartInstance, c1_i13 + 144) == c1_get_B2(chartInstance, 0)) &&
                      (c1_get_states_space(chartInstance, c1_i13 + 216) ==
                       c1_get_K2(chartInstance, 0)));
  }

  c1_c_idx = 0;
  c1_b_ii_size[0] = 72;
  c1_c_ii = 1;
  exitg1 = false;
  while ((!exitg1) && (c1_c_ii < 73)) {
    if (c1_b_x[c1_c_ii - 1]) {
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
    c1_i14 = 0;
  } else {
    c1_i14 = c1_c_idx;
  }

  c1_b_ii_size[0] = c1_i14;
  c1_next_state_idx_size[0] = c1_b_ii_size[0];
  c1_c_loop_ub = c1_b_ii_size[0] - 1;
  for (c1_i15 = 0; c1_i15 <= c1_c_loop_ub; c1_i15++) {
    c1_next_state_idx_data[c1_i15] = (real_T)c1_b_ii_data[c1_i15];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 334);
  c1_b_B3 = 11.0;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 335);
  c1_b_K3 = 3.0;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 338);
  c1_b_H = 5.0;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 339);
  c1_b_Ts = 6.0;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 340);
  c1_h = 1.0;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, -340);
  _SFD_SYMBOL_SCOPE_POP();
  sf_mex_destroy(&c1_current_action);
  *chartInstance->c1_H = c1_b_H;
  *chartInstance->c1_Ts = c1_b_Ts;
  *chartInstance->c1_B3 = c1_b_B3;
  *chartInstance->c1_K3 = c1_b_K3;
  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 0U, chartInstance->c1_sfEvent);
}

static void initSimStructsc1_OVE_RM3(SFc1_OVE_RM3InstanceStruct *chartInstance)
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
  const mxArray *c1_c_y = NULL;
  SFc1_OVE_RM3InstanceStruct *chartInstance;
  chartInstance = (SFc1_OVE_RM3InstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  c1_b_u = *(real_T *)c1_inData;
  c1_c_y = NULL;
  sf_mex_assign(&c1_c_y, sf_mex_create("y", &c1_b_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_c_y, false);
  return c1_mxArrayOutData;
}

static void c1_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_coder_internal_mxSubscript;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_c_y;
  SFc1_OVE_RM3InstanceStruct *chartInstance;
  chartInstance = (SFc1_OVE_RM3InstanceStruct *)chartInstanceVoid;
  c1_coder_internal_mxSubscript = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_c_y = c1_f_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c1_coder_internal_mxSubscript), &c1_thisId);
  sf_mex_destroy(&c1_coder_internal_mxSubscript);
  *(real_T *)c1_outData = c1_c_y;
  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_b_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData;
  int32_T c1_i16;
  const mxArray *c1_c_y = NULL;
  real_T c1_b_u[6];
  SFc1_OVE_RM3InstanceStruct *chartInstance;
  chartInstance = (SFc1_OVE_RM3InstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  for (c1_i16 = 0; c1_i16 < 6; c1_i16++) {
    c1_b_u[c1_i16] = (*(real_T (*)[6])c1_inData)[c1_i16];
  }

  c1_c_y = NULL;
  sf_mex_assign(&c1_c_y, sf_mex_create("y", c1_b_u, 0, 0U, 1U, 0U, 2, 6, 1),
                false);
  sf_mex_assign(&c1_mxArrayOutData, c1_c_y, false);
  return c1_mxArrayOutData;
}

static const mxArray *c1_c_sf_marshallOut(void *chartInstanceVoid, real_T
  c1_inData_data[], int32_T c1_inData_size[1])
{
  const mxArray *c1_mxArrayOutData;
  int32_T c1_u_size[1];
  int32_T c1_loop_ub;
  int32_T c1_i17;
  const mxArray *c1_c_y = NULL;
  real_T c1_u_data[72];
  SFc1_OVE_RM3InstanceStruct *chartInstance;
  chartInstance = (SFc1_OVE_RM3InstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  c1_u_size[0] = c1_inData_size[0];
  c1_loop_ub = c1_inData_size[0] - 1;
  for (c1_i17 = 0; c1_i17 <= c1_loop_ub; c1_i17++) {
    c1_u_data[c1_i17] = c1_inData_data[c1_i17];
  }

  c1_c_y = NULL;
  sf_mex_assign(&c1_c_y, sf_mex_create("y", (void *)&c1_u_data, 0, 0U, 1U, 0U, 1,
    c1_u_size[0]), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_c_y, false);
  return c1_mxArrayOutData;
}

static void c1_emlrt_marshallIn(SFc1_OVE_RM3InstanceStruct *chartInstance, const
  mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y_data[],
  int32_T c1_y_size[1])
{
  uint32_T c1_uv0[1];
  int32_T c1_tmp_size[1];
  boolean_T c1_bv1[1];
  real_T c1_tmp_data[72];
  int32_T c1_loop_ub;
  int32_T c1_i18;
  (void)chartInstance;
  c1_uv0[0] = 72U;
  c1_tmp_size[0] = sf_mex_get_dimension(c1_b_u, 0);
  c1_bv1[0] = true;
  sf_mex_import_vs(c1_parentId, sf_mex_dup(c1_b_u), (void *)&c1_tmp_data, 1, 0,
                   0U, 1, 0U, 1, c1_bv1, c1_uv0, c1_tmp_size);
  c1_y_size[0] = c1_tmp_size[0];
  c1_loop_ub = c1_tmp_size[0] - 1;
  for (c1_i18 = 0; c1_i18 <= c1_loop_ub; c1_i18++) {
    c1_y_data[c1_i18] = c1_tmp_data[c1_i18];
  }

  sf_mex_destroy(&c1_b_u);
}

static void c1_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, real_T c1_outData_data[], int32_T
  c1_outData_size[1])
{
  const mxArray *c1_next_state_idx;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y_data[72];
  int32_T c1_y_size[1];
  int32_T c1_loop_ub;
  int32_T c1_i19;
  SFc1_OVE_RM3InstanceStruct *chartInstance;
  chartInstance = (SFc1_OVE_RM3InstanceStruct *)chartInstanceVoid;
  c1_next_state_idx = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_next_state_idx), &c1_thisId,
                      c1_y_data, c1_y_size);
  sf_mex_destroy(&c1_next_state_idx);
  c1_outData_size[0] = c1_y_size[0];
  c1_loop_ub = c1_y_size[0] - 1;
  for (c1_i19 = 0; c1_i19 <= c1_loop_ub; c1_i19++) {
    c1_outData_data[c1_i19] = c1_y_data[c1_i19];
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_d_sf_marshallOut(void *chartInstanceVoid, real_T
  c1_inData_data[], int32_T c1_inData_size[1])
{
  const mxArray *c1_mxArrayOutData;
  int32_T c1_u_size[1];
  int32_T c1_loop_ub;
  int32_T c1_i20;
  const mxArray *c1_c_y = NULL;
  real_T c1_u_data[4];
  SFc1_OVE_RM3InstanceStruct *chartInstance;
  chartInstance = (SFc1_OVE_RM3InstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  c1_u_size[0] = c1_inData_size[0];
  c1_loop_ub = c1_inData_size[0] - 1;
  for (c1_i20 = 0; c1_i20 <= c1_loop_ub; c1_i20++) {
    c1_u_data[c1_i20] = c1_inData_data[c1_i20];
  }

  c1_c_y = NULL;
  sf_mex_assign(&c1_c_y, sf_mex_create("y", (void *)&c1_u_data, 0, 0U, 1U, 0U, 1,
    c1_u_size[0]), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_c_y, false);
  return c1_mxArrayOutData;
}

static void c1_b_emlrt_marshallIn(SFc1_OVE_RM3InstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T
  c1_y_data[], int32_T c1_y_size[1])
{
  uint32_T c1_uv1[1];
  int32_T c1_tmp_size[1];
  boolean_T c1_bv2[1];
  real_T c1_tmp_data[4];
  int32_T c1_loop_ub;
  int32_T c1_i21;
  (void)chartInstance;
  c1_uv1[0] = 4U;
  c1_tmp_size[0] = sf_mex_get_dimension(c1_b_u, 0);
  c1_bv2[0] = true;
  sf_mex_import_vs(c1_parentId, sf_mex_dup(c1_b_u), (void *)&c1_tmp_data, 1, 0,
                   0U, 1, 0U, 1, c1_bv2, c1_uv1, c1_tmp_size);
  c1_y_size[0] = c1_tmp_size[0];
  c1_loop_ub = c1_tmp_size[0] - 1;
  for (c1_i21 = 0; c1_i21 <= c1_loop_ub; c1_i21++) {
    c1_y_data[c1_i21] = c1_tmp_data[c1_i21];
  }

  sf_mex_destroy(&c1_b_u);
}

static void c1_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, real_T c1_outData_data[], int32_T
  c1_outData_size[1])
{
  const mxArray *c1_action_idx;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y_data[4];
  int32_T c1_y_size[1];
  int32_T c1_loop_ub;
  int32_T c1_i22;
  SFc1_OVE_RM3InstanceStruct *chartInstance;
  chartInstance = (SFc1_OVE_RM3InstanceStruct *)chartInstanceVoid;
  c1_action_idx = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_action_idx), &c1_thisId,
                        c1_y_data, c1_y_size);
  sf_mex_destroy(&c1_action_idx);
  c1_outData_size[0] = c1_y_size[0];
  c1_loop_ub = c1_y_size[0] - 1;
  for (c1_i22 = 0; c1_i22 <= c1_loop_ub; c1_i22++) {
    c1_outData_data[c1_i22] = c1_y_data[c1_i22];
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_e_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  const mxArray *c1_b_u;
  const mxArray *c1_c_y = NULL;
  SFc1_OVE_RM3InstanceStruct *chartInstance;
  chartInstance = (SFc1_OVE_RM3InstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_b_u = sf_mex_dup(*(const mxArray **)c1_inData);
  c1_c_y = NULL;
  sf_mex_assign(&c1_c_y, sf_mex_duplicatearraysafe(&c1_b_u), false);
  sf_mex_destroy(&c1_b_u);
  sf_mex_assign(&c1_mxArrayOutData, c1_c_y, false);
  return c1_mxArrayOutData;
}

static const mxArray *c1_f_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData;
  c1_sI3an6DNyfQAOCkXY61B1nE c1_b_u;
  const mxArray *c1_c_y = NULL;
  real_T c1_c_u;
  const mxArray *c1_d_y = NULL;
  SFc1_OVE_RM3InstanceStruct *chartInstance;
  chartInstance = (SFc1_OVE_RM3InstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  c1_b_u = *(c1_sI3an6DNyfQAOCkXY61B1nE *)c1_inData;
  c1_c_y = NULL;
  sf_mex_assign(&c1_c_y, sf_mex_createstruct("structure", 2, 1, 1), false);
  c1_c_u = c1_b_u.StepRatio;
  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", &c1_c_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c1_c_y, c1_d_y, "StepRatio", "StepRatio", 0);
  sf_mex_assign(&c1_mxArrayOutData, c1_c_y, false);
  return c1_mxArrayOutData;
}

static c1_sI3an6DNyfQAOCkXY61B1nE c1_c_emlrt_marshallIn
  (SFc1_OVE_RM3InstanceStruct *chartInstance, const mxArray *c1_b_u, const
   emlrtMsgIdentifier *c1_parentId)
{
  c1_sI3an6DNyfQAOCkXY61B1nE c1_c_y;
  emlrtMsgIdentifier c1_thisId;
  static const char * c1_fieldNames[1] = { "StepRatio" };

  c1_thisId.fParent = c1_parentId;
  c1_thisId.bParentIsCell = false;
  sf_mex_check_struct(c1_parentId, c1_b_u, 1, c1_fieldNames, 0U, NULL);
  c1_thisId.fIdentifier = "StepRatio";
  c1_c_y.StepRatio = c1_f_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getfield(c1_b_u, "StepRatio", "StepRatio", 0)), &c1_thisId);
  sf_mex_destroy(&c1_b_u);
  return c1_c_y;
}

static void c1_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_opts;
  emlrtMsgIdentifier c1_thisId;
  c1_sI3an6DNyfQAOCkXY61B1nE c1_c_y;
  SFc1_OVE_RM3InstanceStruct *chartInstance;
  chartInstance = (SFc1_OVE_RM3InstanceStruct *)chartInstanceVoid;
  c1_opts = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_c_y = c1_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_opts), &c1_thisId);
  sf_mex_destroy(&c1_opts);
  *(c1_sI3an6DNyfQAOCkXY61B1nE *)c1_outData = c1_c_y;
  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_g_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData;
  int32_T c1_i23;
  const mxArray *c1_c_y = NULL;
  real_T c1_b_u[6];
  SFc1_OVE_RM3InstanceStruct *chartInstance;
  chartInstance = (SFc1_OVE_RM3InstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  for (c1_i23 = 0; c1_i23 < 6; c1_i23++) {
    c1_b_u[c1_i23] = (*(real_T (*)[6])c1_inData)[c1_i23];
  }

  c1_c_y = NULL;
  sf_mex_assign(&c1_c_y, sf_mex_create("y", c1_b_u, 0, 0U, 1U, 0U, 1, 6), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_c_y, false);
  return c1_mxArrayOutData;
}

static void c1_d_emlrt_marshallIn(SFc1_OVE_RM3InstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_c_y[6])
{
  real_T c1_dv1[6];
  int32_T c1_i24;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_b_u), c1_dv1, 1, 0, 0U, 1, 0U, 1, 6);
  for (c1_i24 = 0; c1_i24 < 6; c1_i24++) {
    c1_c_y[c1_i24] = c1_dv1[c1_i24];
  }

  sf_mex_destroy(&c1_b_u);
}

static void c1_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_PptoAvg;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_c_y[6];
  int32_T c1_i25;
  SFc1_OVE_RM3InstanceStruct *chartInstance;
  chartInstance = (SFc1_OVE_RM3InstanceStruct *)chartInstanceVoid;
  c1_PptoAvg = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_PptoAvg), &c1_thisId,
                        c1_c_y);
  sf_mex_destroy(&c1_PptoAvg);
  for (c1_i25 = 0; c1_i25 < 6; c1_i25++) {
    (*(real_T (*)[6])c1_outData)[c1_i25] = c1_c_y[c1_i25];
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

const mxArray *sf_c1_OVE_RM3_get_eml_resolved_functions_info(void)
{
  const mxArray *c1_nameCaptureInfo = NULL;
  c1_nameCaptureInfo = NULL;
  sf_mex_assign(&c1_nameCaptureInfo, sf_mex_create("nameCaptureInfo", NULL, 0,
    0U, 1U, 0U, 2, 0, 1), false);
  return c1_nameCaptureInfo;
}

static real_T c1_e_emlrt_marshallIn(SFc1_OVE_RM3InstanceStruct *chartInstance,
  const mxArray *c1_coder_internal_mxSubscript, const char_T *c1_identifier)
{
  real_T c1_c_y;
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = (const char *)c1_identifier;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_c_y = c1_f_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c1_coder_internal_mxSubscript), &c1_thisId);
  sf_mex_destroy(&c1_coder_internal_mxSubscript);
  return c1_c_y;
}

static real_T c1_f_emlrt_marshallIn(SFc1_OVE_RM3InstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId)
{
  real_T c1_c_y;
  real_T c1_d0;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_b_u), &c1_d0, 1, 0, 0U, 0, 0U, 0);
  c1_c_y = c1_d0;
  sf_mex_destroy(&c1_b_u);
  return c1_c_y;
}

static const mxArray *c1_h_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData;
  int32_T c1_b_u;
  const mxArray *c1_c_y = NULL;
  SFc1_OVE_RM3InstanceStruct *chartInstance;
  chartInstance = (SFc1_OVE_RM3InstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  c1_b_u = *(int32_T *)c1_inData;
  c1_c_y = NULL;
  sf_mex_assign(&c1_c_y, sf_mex_create("y", &c1_b_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_c_y, false);
  return c1_mxArrayOutData;
}

static int32_T c1_g_emlrt_marshallIn(SFc1_OVE_RM3InstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId)
{
  int32_T c1_c_y;
  int32_T c1_i26;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_b_u), &c1_i26, 1, 6, 0U, 0, 0U, 0);
  c1_c_y = c1_i26;
  sf_mex_destroy(&c1_b_u);
  return c1_c_y;
}

static void c1_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_b_sfEvent;
  emlrtMsgIdentifier c1_thisId;
  int32_T c1_c_y;
  SFc1_OVE_RM3InstanceStruct *chartInstance;
  chartInstance = (SFc1_OVE_RM3InstanceStruct *)chartInstanceVoid;
  c1_b_sfEvent = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_c_y = c1_g_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_sfEvent),
    &c1_thisId);
  sf_mex_destroy(&c1_b_sfEvent);
  *(int32_T *)c1_outData = c1_c_y;
  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_i_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData;
  int32_T c1_i27;
  int32_T c1_i28;
  const mxArray *c1_c_y = NULL;
  int32_T c1_i29;
  real_T c1_b_u[8];
  SFc1_OVE_RM3InstanceStruct *chartInstance;
  chartInstance = (SFc1_OVE_RM3InstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  c1_i27 = 0;
  for (c1_i28 = 0; c1_i28 < 2; c1_i28++) {
    for (c1_i29 = 0; c1_i29 < 4; c1_i29++) {
      c1_b_u[c1_i29 + c1_i27] = (*(real_T (*)[8])c1_inData)[c1_i29 + c1_i27];
    }

    c1_i27 += 4;
  }

  c1_c_y = NULL;
  sf_mex_assign(&c1_c_y, sf_mex_create("y", c1_b_u, 0, 0U, 1U, 0U, 2, 4, 2),
                false);
  sf_mex_assign(&c1_mxArrayOutData, c1_c_y, false);
  return c1_mxArrayOutData;
}

static void c1_h_emlrt_marshallIn(SFc1_OVE_RM3InstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_c_y[8])
{
  real_T c1_dv2[8];
  int32_T c1_i30;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_b_u), c1_dv2, 1, 0, 0U, 1, 0U, 2, 4,
                2);
  for (c1_i30 = 0; c1_i30 < 8; c1_i30++) {
    c1_c_y[c1_i30] = c1_dv2[c1_i30];
  }

  sf_mex_destroy(&c1_b_u);
}

static void c1_g_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_action;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_c_y[8];
  int32_T c1_i31;
  int32_T c1_i32;
  int32_T c1_i33;
  SFc1_OVE_RM3InstanceStruct *chartInstance;
  chartInstance = (SFc1_OVE_RM3InstanceStruct *)chartInstanceVoid;
  c1_action = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_h_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_action), &c1_thisId, c1_c_y);
  sf_mex_destroy(&c1_action);
  c1_i31 = 0;
  for (c1_i32 = 0; c1_i32 < 2; c1_i32++) {
    for (c1_i33 = 0; c1_i33 < 4; c1_i33++) {
      (*(real_T (*)[8])c1_outData)[c1_i33 + c1_i31] = c1_c_y[c1_i33 + c1_i31];
    }

    c1_i31 += 4;
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_j_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData;
  int32_T c1_i34;
  const mxArray *c1_c_y = NULL;
  real_T c1_b_u[2];
  SFc1_OVE_RM3InstanceStruct *chartInstance;
  chartInstance = (SFc1_OVE_RM3InstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  for (c1_i34 = 0; c1_i34 < 2; c1_i34++) {
    c1_b_u[c1_i34] = (*(real_T (*)[2])c1_inData)[c1_i34];
  }

  c1_c_y = NULL;
  sf_mex_assign(&c1_c_y, sf_mex_create("y", c1_b_u, 0, 0U, 1U, 0U, 2, 1, 2),
                false);
  sf_mex_assign(&c1_mxArrayOutData, c1_c_y, false);
  return c1_mxArrayOutData;
}

static void c1_i_emlrt_marshallIn(SFc1_OVE_RM3InstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_c_y[2])
{
  real_T c1_dv3[2];
  int32_T c1_i35;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_b_u), c1_dv3, 1, 0, 0U, 1, 0U, 2, 1,
                2);
  for (c1_i35 = 0; c1_i35 < 2; c1_i35++) {
    c1_c_y[c1_i35] = c1_dv3[c1_i35];
  }

  sf_mex_destroy(&c1_b_u);
}

static void c1_h_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_xp;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_c_y[2];
  int32_T c1_i36;
  SFc1_OVE_RM3InstanceStruct *chartInstance;
  chartInstance = (SFc1_OVE_RM3InstanceStruct *)chartInstanceVoid;
  c1_xp = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_i_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_xp), &c1_thisId, c1_c_y);
  sf_mex_destroy(&c1_xp);
  for (c1_i36 = 0; c1_i36 < 2; c1_i36++) {
    (*(real_T (*)[2])c1_outData)[c1_i36] = c1_c_y[c1_i36];
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_k_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData;
  int32_T c1_i37;
  int32_T c1_i38;
  const mxArray *c1_c_y = NULL;
  int32_T c1_i39;
  real_T c1_b_u[20000];
  SFc1_OVE_RM3InstanceStruct *chartInstance;
  chartInstance = (SFc1_OVE_RM3InstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  c1_i37 = 0;
  for (c1_i38 = 0; c1_i38 < 5; c1_i38++) {
    for (c1_i39 = 0; c1_i39 < 4000; c1_i39++) {
      c1_b_u[c1_i39 + c1_i37] = (*(real_T (*)[20000])c1_inData)[c1_i39 + c1_i37];
    }

    c1_i37 += 4000;
  }

  c1_c_y = NULL;
  sf_mex_assign(&c1_c_y, sf_mex_create("y", c1_b_u, 0, 0U, 1U, 0U, 2, 4000, 5),
                false);
  sf_mex_assign(&c1_mxArrayOutData, c1_c_y, false);
  return c1_mxArrayOutData;
}

static void c1_j_emlrt_marshallIn(SFc1_OVE_RM3InstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_c_y
  [20000])
{
  real_T c1_dv4[20000];
  int32_T c1_i40;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_b_u), c1_dv4, 1, 0, 0U, 1, 0U, 2,
                4000, 5);
  for (c1_i40 = 0; c1_i40 < 20000; c1_i40++) {
    c1_c_y[c1_i40] = c1_dv4[c1_i40];
  }

  sf_mex_destroy(&c1_b_u);
}

static void c1_i_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_Stor_eta;
  emlrtMsgIdentifier c1_thisId;
  int32_T c1_i41;
  int32_T c1_i42;
  int32_T c1_i43;
  SFc1_OVE_RM3InstanceStruct *chartInstance;
  chartInstance = (SFc1_OVE_RM3InstanceStruct *)chartInstanceVoid;
  c1_Stor_eta = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_j_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_Stor_eta), &c1_thisId,
                        chartInstance->c1_b_y);
  sf_mex_destroy(&c1_Stor_eta);
  c1_i41 = 0;
  for (c1_i42 = 0; c1_i42 < 5; c1_i42++) {
    for (c1_i43 = 0; c1_i43 < 4000; c1_i43++) {
      (*(real_T (*)[20000])c1_outData)[c1_i43 + c1_i41] = chartInstance->
        c1_b_y[c1_i43 + c1_i41];
    }

    c1_i41 += 4000;
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_l_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData;
  int32_T c1_i44;
  int32_T c1_i45;
  const mxArray *c1_c_y = NULL;
  int32_T c1_i46;
  real_T c1_b_u[288];
  SFc1_OVE_RM3InstanceStruct *chartInstance;
  chartInstance = (SFc1_OVE_RM3InstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  c1_i44 = 0;
  for (c1_i45 = 0; c1_i45 < 4; c1_i45++) {
    for (c1_i46 = 0; c1_i46 < 72; c1_i46++) {
      c1_b_u[c1_i46 + c1_i44] = (*(real_T (*)[288])c1_inData)[c1_i46 + c1_i44];
    }

    c1_i44 += 72;
  }

  c1_c_y = NULL;
  sf_mex_assign(&c1_c_y, sf_mex_create("y", c1_b_u, 0, 0U, 1U, 0U, 2, 72, 4),
                false);
  sf_mex_assign(&c1_mxArrayOutData, c1_c_y, false);
  return c1_mxArrayOutData;
}

static void c1_k_emlrt_marshallIn(SFc1_OVE_RM3InstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_c_y
  [288])
{
  real_T c1_dv5[288];
  int32_T c1_i47;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_b_u), c1_dv5, 1, 0, 0U, 1, 0U, 2, 72,
                4);
  for (c1_i47 = 0; c1_i47 < 288; c1_i47++) {
    c1_c_y[c1_i47] = c1_dv5[c1_i47];
  }

  sf_mex_destroy(&c1_b_u);
}

static void c1_j_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_states_space;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_c_y[288];
  int32_T c1_i48;
  int32_T c1_i49;
  int32_T c1_i50;
  SFc1_OVE_RM3InstanceStruct *chartInstance;
  chartInstance = (SFc1_OVE_RM3InstanceStruct *)chartInstanceVoid;
  c1_states_space = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_k_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_states_space), &c1_thisId,
                        c1_c_y);
  sf_mex_destroy(&c1_states_space);
  c1_i48 = 0;
  for (c1_i49 = 0; c1_i49 < 4; c1_i49++) {
    for (c1_i50 = 0; c1_i50 < 72; c1_i50++) {
      (*(real_T (*)[288])c1_outData)[c1_i50 + c1_i48] = c1_c_y[c1_i50 + c1_i48];
    }

    c1_i48 += 72;
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_m_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData;
  int32_T c1_i51;
  const mxArray *c1_c_y = NULL;
  real_T c1_b_u[72];
  SFc1_OVE_RM3InstanceStruct *chartInstance;
  chartInstance = (SFc1_OVE_RM3InstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  for (c1_i51 = 0; c1_i51 < 72; c1_i51++) {
    c1_b_u[c1_i51] = (*(real_T (*)[72])c1_inData)[c1_i51];
  }

  c1_c_y = NULL;
  sf_mex_assign(&c1_c_y, sf_mex_create("y", c1_b_u, 0, 0U, 1U, 0U, 2, 72, 1),
                false);
  sf_mex_assign(&c1_mxArrayOutData, c1_c_y, false);
  return c1_mxArrayOutData;
}

static void c1_l_emlrt_marshallIn(SFc1_OVE_RM3InstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_c_y[72])
{
  real_T c1_dv6[72];
  int32_T c1_i52;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_b_u), c1_dv6, 1, 0, 0U, 1, 0U, 2, 72,
                1);
  for (c1_i52 = 0; c1_i52 < 72; c1_i52++) {
    c1_c_y[c1_i52] = c1_dv6[c1_i52];
  }

  sf_mex_destroy(&c1_b_u);
}

static void c1_k_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_num_R;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_c_y[72];
  int32_T c1_i53;
  SFc1_OVE_RM3InstanceStruct *chartInstance;
  chartInstance = (SFc1_OVE_RM3InstanceStruct *)chartInstanceVoid;
  c1_num_R = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_l_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_num_R), &c1_thisId, c1_c_y);
  sf_mex_destroy(&c1_num_R);
  for (c1_i53 = 0; c1_i53 < 72; c1_i53++) {
    (*(real_T (*)[72])c1_outData)[c1_i53] = c1_c_y[c1_i53];
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_n_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData;
  int32_T c1_i54;
  int32_T c1_i55;
  const mxArray *c1_c_y = NULL;
  int32_T c1_i56;
  real_T c1_b_u[7200];
  SFc1_OVE_RM3InstanceStruct *chartInstance;
  chartInstance = (SFc1_OVE_RM3InstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  c1_i54 = 0;
  for (c1_i55 = 0; c1_i55 < 100; c1_i55++) {
    for (c1_i56 = 0; c1_i56 < 72; c1_i56++) {
      c1_b_u[c1_i56 + c1_i54] = (*(real_T (*)[7200])c1_inData)[c1_i56 + c1_i54];
    }

    c1_i54 += 72;
  }

  c1_c_y = NULL;
  sf_mex_assign(&c1_c_y, sf_mex_create("y", c1_b_u, 0, 0U, 1U, 0U, 2, 72, 100),
                false);
  sf_mex_assign(&c1_mxArrayOutData, c1_c_y, false);
  return c1_mxArrayOutData;
}

static void c1_m_emlrt_marshallIn(SFc1_OVE_RM3InstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_c_y
  [7200])
{
  real_T c1_dv7[7200];
  int32_T c1_i57;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_b_u), c1_dv7, 1, 0, 0U, 1, 0U, 2, 72,
                100);
  for (c1_i57 = 0; c1_i57 < 7200; c1_i57++) {
    c1_c_y[c1_i57] = c1_dv7[c1_i57];
  }

  sf_mex_destroy(&c1_b_u);
}

static void c1_l_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_R;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_c_y[7200];
  int32_T c1_i58;
  int32_T c1_i59;
  int32_T c1_i60;
  SFc1_OVE_RM3InstanceStruct *chartInstance;
  chartInstance = (SFc1_OVE_RM3InstanceStruct *)chartInstanceVoid;
  c1_R = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_m_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_R), &c1_thisId, c1_c_y);
  sf_mex_destroy(&c1_R);
  c1_i58 = 0;
  for (c1_i59 = 0; c1_i59 < 100; c1_i59++) {
    for (c1_i60 = 0; c1_i60 < 72; c1_i60++) {
      (*(real_T (*)[7200])c1_outData)[c1_i60 + c1_i58] = c1_c_y[c1_i60 + c1_i58];
    }

    c1_i58 += 72;
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_o_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData;
  int32_T c1_i61;
  const mxArray *c1_c_y = NULL;
  real_T c1_b_u[72];
  SFc1_OVE_RM3InstanceStruct *chartInstance;
  chartInstance = (SFc1_OVE_RM3InstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  for (c1_i61 = 0; c1_i61 < 72; c1_i61++) {
    c1_b_u[c1_i61] = (*(real_T (*)[72])c1_inData)[c1_i61];
  }

  c1_c_y = NULL;
  sf_mex_assign(&c1_c_y, sf_mex_create("y", c1_b_u, 0, 0U, 1U, 0U, 1, 72), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_c_y, false);
  return c1_mxArrayOutData;
}

static void c1_n_emlrt_marshallIn(SFc1_OVE_RM3InstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_c_y[72])
{
  real_T c1_dv8[72];
  int32_T c1_i62;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_b_u), c1_dv8, 1, 0, 0U, 1, 0U, 1, 72);
  for (c1_i62 = 0; c1_i62 < 72; c1_i62++) {
    c1_c_y[c1_i62] = c1_dv8[c1_i62];
  }

  sf_mex_destroy(&c1_b_u);
}

static void c1_m_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_Final_Rew_matt;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_c_y[72];
  int32_T c1_i63;
  SFc1_OVE_RM3InstanceStruct *chartInstance;
  chartInstance = (SFc1_OVE_RM3InstanceStruct *)chartInstanceVoid;
  c1_Final_Rew_matt = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_n_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_Final_Rew_matt), &c1_thisId,
                        c1_c_y);
  sf_mex_destroy(&c1_Final_Rew_matt);
  for (c1_i63 = 0; c1_i63 < 72; c1_i63++) {
    (*(real_T (*)[72])c1_outData)[c1_i63] = c1_c_y[c1_i63];
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_p_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData;
  int32_T c1_i64;
  int32_T c1_i65;
  const mxArray *c1_c_y = NULL;
  int32_T c1_i66;
  SFc1_OVE_RM3InstanceStruct *chartInstance;
  chartInstance = (SFc1_OVE_RM3InstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_mxArrayOutData = NULL;
  c1_i64 = 0;
  for (c1_i65 = 0; c1_i65 < 12; c1_i65++) {
    for (c1_i66 = 0; c1_i66 < 4000; c1_i66++) {
      chartInstance->c1_u[c1_i66 + c1_i64] = (*(real_T (*)[48000])c1_inData)
        [c1_i66 + c1_i64];
    }

    c1_i64 += 4000;
  }

  c1_c_y = NULL;
  sf_mex_assign(&c1_c_y, sf_mex_create("y", chartInstance->c1_u, 0, 0U, 1U, 0U,
    2, 4000, 12), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_c_y, false);
  return c1_mxArrayOutData;
}

static void c1_o_emlrt_marshallIn(SFc1_OVE_RM3InstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_c_y
  [48000])
{
  int32_T c1_i67;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_b_u), chartInstance->c1_dv9, 1, 0, 0U,
                1, 0U, 2, 4000, 12);
  for (c1_i67 = 0; c1_i67 < 48000; c1_i67++) {
    c1_c_y[c1_i67] = chartInstance->c1_dv9[c1_i67];
  }

  sf_mex_destroy(&c1_b_u);
}

static void c1_n_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_Memory_Stor;
  emlrtMsgIdentifier c1_thisId;
  int32_T c1_i68;
  int32_T c1_i69;
  int32_T c1_i70;
  SFc1_OVE_RM3InstanceStruct *chartInstance;
  chartInstance = (SFc1_OVE_RM3InstanceStruct *)chartInstanceVoid;
  c1_Memory_Stor = sf_mex_dup(c1_mxArrayInData);
  c1_thisId.fIdentifier = (const char *)c1_varName;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_o_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_Memory_Stor), &c1_thisId,
                        chartInstance->c1_y);
  sf_mex_destroy(&c1_Memory_Stor);
  c1_i68 = 0;
  for (c1_i69 = 0; c1_i69 < 12; c1_i69++) {
    for (c1_i70 = 0; c1_i70 < 4000; c1_i70++) {
      (*(real_T (*)[48000])c1_outData)[c1_i70 + c1_i68] = chartInstance->
        c1_y[c1_i70 + c1_i68];
    }

    c1_i68 += 4000;
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static uint8_T c1_p_emlrt_marshallIn(SFc1_OVE_RM3InstanceStruct *chartInstance,
  const mxArray *c1_b_is_active_c1_OVE_RM3, const char_T *c1_identifier)
{
  uint8_T c1_c_y;
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = (const char *)c1_identifier;
  c1_thisId.fParent = NULL;
  c1_thisId.bParentIsCell = false;
  c1_c_y = c1_q_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c1_b_is_active_c1_OVE_RM3), &c1_thisId);
  sf_mex_destroy(&c1_b_is_active_c1_OVE_RM3);
  return c1_c_y;
}

static uint8_T c1_q_emlrt_marshallIn(SFc1_OVE_RM3InstanceStruct *chartInstance,
  const mxArray *c1_b_u, const emlrtMsgIdentifier *c1_parentId)
{
  uint8_T c1_c_y;
  uint8_T c1_u0;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_b_u), &c1_u0, 1, 3, 0U, 0, 0U, 0);
  c1_c_y = c1_u0;
  sf_mex_destroy(&c1_b_u);
  return c1_c_y;
}

static int32_T c1__s32_d_(SFc1_OVE_RM3InstanceStruct *chartInstance, real_T c1_b,
  uint32_T c1_ssid_src_loc, int32_T c1_offset_src_loc, int32_T c1_length_src_loc)
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

static real_T c1_get_B2(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_B2) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 0, "B2",
      c1_elementIndex);
  }

  return *chartInstance->c1_B2_address;
}

static void c1_set_B2(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
                      c1_elementIndex, real_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_B2) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 0, "B2", c1_elementIndex);
  }

  *chartInstance->c1_B2_address = c1_elementValue;
}

static real_T *c1_access_B2(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_B2) {
    ssAccessDataStore_wrapper(chartInstance->S, 0, "B2", c1_rdOnly);
  }

  return chartInstance->c1_B2_address;
}

static real_T c1_get_Final_Rew_matt(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_Final_Rew_matt) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 1, "Final_Rew_matt",
      c1_elementIndex);
  }

  return (*chartInstance->c1__indexFinal_Rew_matt_address)[c1_elementIndex];
}

static void c1_set_Final_Rew_matt(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_elementIndex, real_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_Final_Rew_matt) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 1, "Final_Rew_matt",
      c1_elementIndex);
  }

  (*chartInstance->c1__indexFinal_Rew_matt_address)[c1_elementIndex] =
    c1_elementValue;
}

static real_T *c1_access_Final_Rew_matt(SFc1_OVE_RM3InstanceStruct
  *chartInstance, uint32_T c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_Final_Rew_matt) {
    ssAccessDataStore_wrapper(chartInstance->S, 1, "Final_Rew_matt", c1_rdOnly);
  }

  return &(*chartInstance->c1__indexFinal_Rew_matt_address)[0U];
}

static real_T c1_get_K2(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_K2) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 2, "K2",
      c1_elementIndex);
  }

  return *chartInstance->c1__indexK2_address;
}

static void c1_set_K2(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
                      c1_elementIndex, real_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_K2) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 2, "K2", c1_elementIndex);
  }

  *chartInstance->c1__indexK2_address = c1_elementValue;
}

static real_T *c1_access_K2(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_K2) {
    ssAccessDataStore_wrapper(chartInstance->S, 2, "K2", c1_rdOnly);
  }

  return chartInstance->c1__indexK2_address;
}

static real_T c1_get_MeanR(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_MeanR) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 3, "MeanR",
      c1_elementIndex);
  }

  return (*chartInstance->c1__indexMeanR_address)[c1_elementIndex];
}

static void c1_set_MeanR(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex, real_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_MeanR) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 3, "MeanR",
      c1_elementIndex);
  }

  (*chartInstance->c1__indexMeanR_address)[c1_elementIndex] = c1_elementValue;
}

static real_T *c1_access_MeanR(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_MeanR) {
    ssAccessDataStore_wrapper(chartInstance->S, 3, "MeanR", c1_rdOnly);
  }

  return &(*chartInstance->c1__indexMeanR_address)[0U];
}

static real_T c1_get_Memory_D_size(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_Memory_D_size) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 4, "Memory_D_size",
      c1_elementIndex);
  }

  return *chartInstance->c1__indexMemory_D_size_address;
}

static void c1_set_Memory_D_size(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_elementIndex, real_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_Memory_D_size) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 4, "Memory_D_size",
      c1_elementIndex);
  }

  *chartInstance->c1__indexMemory_D_size_address = c1_elementValue;
}

static real_T *c1_access_Memory_D_size(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_Memory_D_size) {
    ssAccessDataStore_wrapper(chartInstance->S, 4, "Memory_D_size", c1_rdOnly);
  }

  return chartInstance->c1__indexMemory_D_size_address;
}

static real_T c1_get_Memory_Stor(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_Memory_Stor) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 5, "Memory_Stor",
      c1_elementIndex);
  }

  return (*chartInstance->c1__indexMemory_Stor_address)[c1_elementIndex];
}

static void c1_set_Memory_Stor(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_elementIndex, real_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_Memory_Stor) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 5, "Memory_Stor",
      c1_elementIndex);
  }

  (*chartInstance->c1__indexMemory_Stor_address)[c1_elementIndex] =
    c1_elementValue;
}

static real_T *c1_access_Memory_Stor(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_Memory_Stor) {
    ssAccessDataStore_wrapper(chartInstance->S, 5, "Memory_Stor", c1_rdOnly);
  }

  return &(*chartInstance->c1__indexMemory_Stor_address)[0U];
}

static real_T c1_get_NN_dd(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_NN_dd) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 6, "NN_dd",
      c1_elementIndex);
  }

  return *chartInstance->c1__indexNN_dd_address;
}

static void c1_set_NN_dd(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex, real_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_NN_dd) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 6, "NN_dd",
      c1_elementIndex);
  }

  *chartInstance->c1__indexNN_dd_address = c1_elementValue;
}

static real_T *c1_access_NN_dd(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_NN_dd) {
    ssAccessDataStore_wrapper(chartInstance->S, 6, "NN_dd", c1_rdOnly);
  }

  return chartInstance->c1__indexNN_dd_address;
}

static real_T c1_get_Q(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_Q) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 7, "Q", c1_elementIndex);
  }

  return (*chartInstance->c1__indexQ_address)[c1_elementIndex];
}

static void c1_set_Q(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
                     c1_elementIndex, real_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_Q) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 7, "Q", c1_elementIndex);
  }

  (*chartInstance->c1__indexQ_address)[c1_elementIndex] = c1_elementValue;
}

static real_T *c1_access_Q(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_Q) {
    ssAccessDataStore_wrapper(chartInstance->S, 7, "Q", c1_rdOnly);
  }

  return &(*chartInstance->c1__indexQ_address)[0U];
}

static real_T c1_get_Q_prime(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_Q_prime) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 8, "Q_prime",
      c1_elementIndex);
  }

  return (*chartInstance->c1__indexQ_prime_address)[c1_elementIndex];
}

static void c1_set_Q_prime(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex, real_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_Q_prime) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 8, "Q_prime",
      c1_elementIndex);
  }

  (*chartInstance->c1__indexQ_prime_address)[c1_elementIndex] = c1_elementValue;
}

static real_T *c1_access_Q_prime(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_Q_prime) {
    ssAccessDataStore_wrapper(chartInstance->S, 8, "Q_prime", c1_rdOnly);
  }

  return &(*chartInstance->c1__indexQ_prime_address)[0U];
}

static real_T c1_get_Q_targ(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_Q_targ) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 9, "Q_targ",
      c1_elementIndex);
  }

  return (*chartInstance->c1__indexQ_targ_address)[c1_elementIndex];
}

static void c1_set_Q_targ(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex, real_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_Q_targ) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 9, "Q_targ",
      c1_elementIndex);
  }

  (*chartInstance->c1__indexQ_targ_address)[c1_elementIndex] = c1_elementValue;
}

static real_T *c1_access_Q_targ(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_Q_targ) {
    ssAccessDataStore_wrapper(chartInstance->S, 9, "Q_targ", c1_rdOnly);
  }

  return &(*chartInstance->c1__indexQ_targ_address)[0U];
}

static real_T c1_get_R(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_R) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 10, "R",
      c1_elementIndex);
  }

  return (*chartInstance->c1__indexR_address)[c1_elementIndex];
}

static void c1_set_R(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
                     c1_elementIndex, real_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_R) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 10, "R", c1_elementIndex);
  }

  (*chartInstance->c1__indexR_address)[c1_elementIndex] = c1_elementValue;
}

static real_T *c1_access_R(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_R) {
    ssAccessDataStore_wrapper(chartInstance->S, 10, "R", c1_rdOnly);
  }

  return &(*chartInstance->c1__indexR_address)[0U];
}

static real_T c1_get_Stor_eta(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_Stor_eta) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 11, "Stor_eta",
      c1_elementIndex);
  }

  return (*chartInstance->c1__indexStor_eta_address)[c1_elementIndex];
}

static void c1_set_Stor_eta(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex, real_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_Stor_eta) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 11, "Stor_eta",
      c1_elementIndex);
  }

  (*chartInstance->c1__indexStor_eta_address)[c1_elementIndex] = c1_elementValue;
}

static real_T *c1_access_Stor_eta(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_Stor_eta) {
    ssAccessDataStore_wrapper(chartInstance->S, 11, "Stor_eta", c1_rdOnly);
  }

  return &(*chartInstance->c1__indexStor_eta_address)[0U];
}

static real_T c1_get_action(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_action) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 12, "action",
      c1_elementIndex);
  }

  return (*chartInstance->c1__indexaction_address)[c1_elementIndex];
}

static void c1_set_action(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex, real_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_action) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 12, "action",
      c1_elementIndex);
  }

  (*chartInstance->c1__indexaction_address)[c1_elementIndex] = c1_elementValue;
}

static real_T *c1_access_action(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_action) {
    ssAccessDataStore_wrapper(chartInstance->S, 12, "action", c1_rdOnly);
  }

  return &(*chartInstance->c1__indexaction_address)[0U];
}

static real_T c1_get_num_R(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_num_R) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 13, "num_R",
      c1_elementIndex);
  }

  return (*chartInstance->c1__indexnum_R_address)[c1_elementIndex];
}

static void c1_set_num_R(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex, real_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_num_R) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 13, "num_R",
      c1_elementIndex);
  }

  (*chartInstance->c1__indexnum_R_address)[c1_elementIndex] = c1_elementValue;
}

static real_T *c1_access_num_R(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_num_R) {
    ssAccessDataStore_wrapper(chartInstance->S, 13, "num_R", c1_rdOnly);
  }

  return &(*chartInstance->c1__indexnum_R_address)[0U];
}

static real_T c1_get_ove_count(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_ove_count) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 14, "ove_count",
      c1_elementIndex);
  }

  return *chartInstance->c1__indexove_count_address;
}

static void c1_set_ove_count(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex, real_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_ove_count) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 14, "ove_count",
      c1_elementIndex);
  }

  *chartInstance->c1__indexove_count_address = c1_elementValue;
}

static real_T *c1_access_ove_count(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_ove_count) {
    ssAccessDataStore_wrapper(chartInstance->S, 14, "ove_count", c1_rdOnly);
  }

  return chartInstance->c1__indexove_count_address;
}

static real_T c1_get_state_idx(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_state_idx) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 15, "state_idx",
      c1_elementIndex);
  }

  return *chartInstance->c1__indexstate_idx_address;
}

static void c1_set_state_idx(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex, real_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_state_idx) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 15, "state_idx",
      c1_elementIndex);
  }

  *chartInstance->c1__indexstate_idx_address = c1_elementValue;
}

static real_T *c1_access_state_idx(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_state_idx) {
    ssAccessDataStore_wrapper(chartInstance->S, 15, "state_idx", c1_rdOnly);
  }

  return chartInstance->c1__indexstate_idx_address;
}

static real_T c1_get_states_space(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_states_space) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 16, "states_space",
      c1_elementIndex);
  }

  return (*chartInstance->c1__indexstates_space_address)[c1_elementIndex];
}

static void c1_set_states_space(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_elementIndex, real_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_states_space) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 16, "states_space",
      c1_elementIndex);
  }

  (*chartInstance->c1__indexstates_space_address)[c1_elementIndex] =
    c1_elementValue;
}

static real_T *c1_access_states_space(SFc1_OVE_RM3InstanceStruct *chartInstance,
  uint32_T c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_states_space) {
    ssAccessDataStore_wrapper(chartInstance->S, 16, "states_space", c1_rdOnly);
  }

  return &(*chartInstance->c1__indexstates_space_address)[0U];
}

static real_T c1_get_xp(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_elementIndex)
{
  if (chartInstance->c1_dsmdiag_xp) {
    ssReadFromDataStoreElement_wrapper(chartInstance->S, 17, "xp",
      c1_elementIndex);
  }

  return (*chartInstance->c1__indexxp_address)[c1_elementIndex];
}

static void c1_set_xp(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
                      c1_elementIndex, real_T c1_elementValue)
{
  if (chartInstance->c1_dsmdiag_xp) {
    ssWriteToDataStoreElement_wrapper(chartInstance->S, 17, "xp",
      c1_elementIndex);
  }

  (*chartInstance->c1__indexxp_address)[c1_elementIndex] = c1_elementValue;
}

static real_T *c1_access_xp(SFc1_OVE_RM3InstanceStruct *chartInstance, uint32_T
  c1_rdOnly)
{
  if (chartInstance->c1_dsmdiag_xp) {
    ssAccessDataStore_wrapper(chartInstance->S, 17, "xp", c1_rdOnly);
  }

  return &(*chartInstance->c1__indexxp_address)[0U];
}

static void init_dsm_address_info(SFc1_OVE_RM3InstanceStruct *chartInstance)
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
  ssGetSFcnDataStoreNameAddrIdx_wrapper(chartInstance->S, "xp", (void **)
    &chartInstance->c1__indexxp_address, &chartInstance->c1_j__index);
}

static void init_simulink_io_address(SFc1_OVE_RM3InstanceStruct *chartInstance)
{
  chartInstance->c1_fEmlrtCtx = (void *)sfrtGetEmlrtCtx(chartInstance->S);
  chartInstance->c1_H = (real_T *)ssGetOutputPortSignal_wrapper(chartInstance->S,
    1);
  chartInstance->c1_Ts = (real_T *)ssGetOutputPortSignal_wrapper
    (chartInstance->S, 2);
  chartInstance->c1_power = (real_T (*)[6])ssGetInputPortSignal_wrapper
    (chartInstance->S, 0);
  chartInstance->c1_B3 = (real_T *)ssGetOutputPortSignal_wrapper
    (chartInstance->S, 3);
  chartInstance->c1_K3 = (real_T *)ssGetOutputPortSignal_wrapper
    (chartInstance->S, 4);
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
  chartInstance->c1_dsmdiag_xp = (boolean_T)
    ssGetDSMBlockDiagnosticsEnabled_wrapper(chartInstance->S, 17, "xp");
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

void sf_c1_OVE_RM3_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(2515993147U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(1426346805U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(150936956U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(4157173700U);
}

mxArray* sf_c1_OVE_RM3_get_post_codegen_info(void);
mxArray *sf_c1_OVE_RM3_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals", "postCodegenInfo" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1, 1, sizeof
    (autoinheritanceFields)/sizeof(autoinheritanceFields[0]),
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("txsZtjRHQJ0VXiDKFnWtrF");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,1,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(6);
      pr[1] = (double)(1);
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
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxCreateDoubleMatrix(0,0,
                mxREAL));
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,4,3,dataFields);

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

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,0,mxREAL);
      double *pr = mxGetPr(mxSize);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt", "isFixedPointType" };

      mxArray *mxType = mxCreateStructMatrix(1,1,sizeof(typeFields)/sizeof
        (typeFields[0]),typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxType,0,"isFixedPointType",mxCreateDoubleScalar(0));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,0,mxREAL);
      double *pr = mxGetPr(mxSize);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt", "isFixedPointType" };

      mxArray *mxType = mxCreateStructMatrix(1,1,sizeof(typeFields)/sizeof
        (typeFields[0]),typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxType,0,"isFixedPointType",mxCreateDoubleScalar(0));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,0,mxREAL);
      double *pr = mxGetPr(mxSize);
      mxSetField(mxData,3,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt", "isFixedPointType" };

      mxArray *mxType = mxCreateStructMatrix(1,1,sizeof(typeFields)/sizeof
        (typeFields[0]),typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxType,0,"isFixedPointType",mxCreateDoubleScalar(0));
      mxSetField(mxData,3,"type",mxType);
    }

    mxSetField(mxData,3,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  {
    mxArray* mxPostCodegenInfo = sf_c1_OVE_RM3_get_post_codegen_info();
    mxSetField(mxAutoinheritanceInfo,0,"postCodegenInfo",mxPostCodegenInfo);
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c1_OVE_RM3_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c1_OVE_RM3_jit_fallback_info(void)
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

mxArray *sf_c1_OVE_RM3_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

mxArray* sf_c1_OVE_RM3_get_post_codegen_info(void)
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

static const mxArray *sf_get_sim_state_info_c1_OVE_RM3(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x5'type','srcId','name','auxInfo'{{M[1],M[27],T\"B3\",},{M[1],M[29],T\"H\",},{M[1],M[28],T\"K3\",},{M[1],M[30],T\"Ts\",},{M[8],M[0],T\"is_active_c1_OVE_RM3\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 5, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c1_OVE_RM3_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc1_OVE_RM3InstanceStruct *chartInstance = (SFc1_OVE_RM3InstanceStruct *)
      sf_get_chart_instance_ptr(S);
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _OVE_RM3MachineNumber_,
           1,
           1,
           1,
           0,
           23,
           0,
           0,
           0,
           0,
           0,
           &chartInstance->chartNumber,
           &chartInstance->instanceNumber,
           (void *)S);

        /* Each instance must initialize its own list of scripts */
        init_script_number_translation(_OVE_RM3MachineNumber_,
          chartInstance->chartNumber,chartInstance->instanceNumber);
        if (chartAlreadyPresent==0) {
          /* this is the first instance */
          sf_debug_set_chart_disable_implicit_casting
            (sfGlobalDebugInstanceStruct,_OVE_RM3MachineNumber_,
             chartInstance->chartNumber,1);
          sf_debug_set_chart_event_thresholds(sfGlobalDebugInstanceStruct,
            _OVE_RM3MachineNumber_,
            chartInstance->chartNumber,
            0,
            0,
            0);
          _SFD_SET_DATA_PROPS(0,1,1,0,"power");
          _SFD_SET_DATA_PROPS(1,2,0,1,"H");
          _SFD_SET_DATA_PROPS(2,2,0,1,"Ts");
          _SFD_SET_DATA_PROPS(3,2,0,1,"B3");
          _SFD_SET_DATA_PROPS(4,2,0,1,"K3");
          _SFD_SET_DATA_PROPS(5,11,0,0,"B2");
          _SFD_SET_DATA_PROPS(6,11,0,0,"Final_Rew_matt");
          _SFD_SET_DATA_PROPS(7,11,0,0,"K2");
          _SFD_SET_DATA_PROPS(8,11,0,0,"MeanR");
          _SFD_SET_DATA_PROPS(9,11,0,0,"Memory_D_size");
          _SFD_SET_DATA_PROPS(10,11,0,0,"Memory_Stor");
          _SFD_SET_DATA_PROPS(11,11,0,0,"NN_dd");
          _SFD_SET_DATA_PROPS(12,11,0,0,"Q");
          _SFD_SET_DATA_PROPS(13,11,0,0,"Q_prime");
          _SFD_SET_DATA_PROPS(14,11,0,0,"Q_targ");
          _SFD_SET_DATA_PROPS(15,11,0,0,"R");
          _SFD_SET_DATA_PROPS(16,11,0,0,"Stor_eta");
          _SFD_SET_DATA_PROPS(17,11,0,0,"action");
          _SFD_SET_DATA_PROPS(18,11,0,0,"num_R");
          _SFD_SET_DATA_PROPS(19,11,0,0,"ove_count");
          _SFD_SET_DATA_PROPS(20,11,0,0,"state_idx");
          _SFD_SET_DATA_PROPS(21,11,0,0,"states_space");
          _SFD_SET_DATA_PROPS(22,11,0,0,"xp");
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
        _SFD_CV_INIT_EML(0,1,1,0,2,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,10440);
        _SFD_CV_INIT_EML_IF(0,1,0,924,941,1021,1120);
        _SFD_CV_INIT_EML_IF(0,1,1,1156,1173,1231,1296);
        _SFD_CV_INIT_EML_RELATIONAL(0,1,0,927,941,-1,5);
        _SFD_CV_INIT_EML_RELATIONAL(0,1,1,1159,1173,-1,5);

        {
          unsigned int dimVector[2];
          dimVector[0]= 6U;
          dimVector[1]= 1U;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)c1_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)c1_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)c1_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)c1_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)c1_sf_marshallIn);

        {
          unsigned int dimVector[1];
          dimVector[0]= 72U;
          _SFD_SET_DATA_COMPILED_PROPS(6,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_o_sf_marshallOut,(MexInFcnForType)
            c1_m_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(7,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)c1_sf_marshallIn);

        {
          unsigned int dimVector[2];
          dimVector[0]= 72U;
          dimVector[1]= 1U;
          _SFD_SET_DATA_COMPILED_PROPS(8,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_m_sf_marshallOut,(MexInFcnForType)
            c1_k_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(9,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)c1_sf_marshallIn);

        {
          unsigned int dimVector[2];
          dimVector[0]= 4000U;
          dimVector[1]= 12U;
          _SFD_SET_DATA_COMPILED_PROPS(10,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_p_sf_marshallOut,(MexInFcnForType)
            c1_n_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(11,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)c1_sf_marshallIn);

        {
          unsigned int dimVector[2];
          dimVector[0]= 72U;
          dimVector[1]= 4U;
          _SFD_SET_DATA_COMPILED_PROPS(12,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_l_sf_marshallOut,(MexInFcnForType)
            c1_j_sf_marshallIn);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 72U;
          dimVector[1]= 4U;
          _SFD_SET_DATA_COMPILED_PROPS(13,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_l_sf_marshallOut,(MexInFcnForType)
            c1_j_sf_marshallIn);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 72U;
          dimVector[1]= 4U;
          _SFD_SET_DATA_COMPILED_PROPS(14,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_l_sf_marshallOut,(MexInFcnForType)
            c1_j_sf_marshallIn);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 72U;
          dimVector[1]= 100U;
          _SFD_SET_DATA_COMPILED_PROPS(15,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_n_sf_marshallOut,(MexInFcnForType)
            c1_l_sf_marshallIn);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 4000U;
          dimVector[1]= 5U;
          _SFD_SET_DATA_COMPILED_PROPS(16,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_k_sf_marshallOut,(MexInFcnForType)
            c1_i_sf_marshallIn);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 4U;
          dimVector[1]= 2U;
          _SFD_SET_DATA_COMPILED_PROPS(17,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_i_sf_marshallOut,(MexInFcnForType)
            c1_g_sf_marshallIn);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 72U;
          dimVector[1]= 1U;
          _SFD_SET_DATA_COMPILED_PROPS(18,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_m_sf_marshallOut,(MexInFcnForType)
            c1_k_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(19,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)c1_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(20,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)c1_sf_marshallIn);

        {
          unsigned int dimVector[2];
          dimVector[0]= 72U;
          dimVector[1]= 4U;
          _SFD_SET_DATA_COMPILED_PROPS(21,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_l_sf_marshallOut,(MexInFcnForType)
            c1_j_sf_marshallIn);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 1U;
          dimVector[1]= 2U;
          _SFD_SET_DATA_COMPILED_PROPS(22,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_j_sf_marshallOut,(MexInFcnForType)
            c1_h_sf_marshallIn);
        }
      }
    } else {
      sf_debug_reset_current_state_configuration(sfGlobalDebugInstanceStruct,
        _OVE_RM3MachineNumber_,chartInstance->chartNumber,
        chartInstance->instanceNumber);
    }
  }
}

static void chart_debug_initialize_data_addresses(SimStruct *S)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc1_OVE_RM3InstanceStruct *chartInstance = (SFc1_OVE_RM3InstanceStruct *)
      sf_get_chart_instance_ptr(S);
    if (ssIsFirstInitCond(S)) {
      /* do this only if simulation is starting and after we know the addresses of all data */
      {
        _SFD_SET_DATA_VALUE_PTR(1U, (void *)chartInstance->c1_H);
        _SFD_SET_DATA_VALUE_PTR(2U, (void *)chartInstance->c1_Ts);
        _SFD_SET_DATA_VALUE_PTR(19U, (void *)
          chartInstance->c1__indexove_count_address);
        _SFD_SET_DATA_VALUE_PTR(17U, (void *)
          chartInstance->c1__indexaction_address);
        _SFD_SET_DATA_VALUE_PTR(22U, (void *)chartInstance->c1__indexxp_address);
        _SFD_SET_DATA_VALUE_PTR(5U, (void *)chartInstance->c1_B2_address);
        _SFD_SET_DATA_VALUE_PTR(7U, (void *)chartInstance->c1__indexK2_address);
        _SFD_SET_DATA_VALUE_PTR(16U, (void *)
          chartInstance->c1__indexStor_eta_address);
        _SFD_SET_DATA_VALUE_PTR(20U, (void *)
          chartInstance->c1__indexstate_idx_address);
        _SFD_SET_DATA_VALUE_PTR(21U, (void *)
          chartInstance->c1__indexstates_space_address);
        _SFD_SET_DATA_VALUE_PTR(18U, (void *)
          chartInstance->c1__indexnum_R_address);
        _SFD_SET_DATA_VALUE_PTR(8U, (void *)
          chartInstance->c1__indexMeanR_address);
        _SFD_SET_DATA_VALUE_PTR(15U, (void *)chartInstance->c1__indexR_address);
        _SFD_SET_DATA_VALUE_PTR(6U, (void *)
          chartInstance->c1__indexFinal_Rew_matt_address);
        _SFD_SET_DATA_VALUE_PTR(11U, (void *)
          chartInstance->c1__indexNN_dd_address);
        _SFD_SET_DATA_VALUE_PTR(9U, (void *)
          chartInstance->c1__indexMemory_D_size_address);
        _SFD_SET_DATA_VALUE_PTR(10U, (void *)
          chartInstance->c1__indexMemory_Stor_address);
        _SFD_SET_DATA_VALUE_PTR(13U, (void *)
          chartInstance->c1__indexQ_prime_address);
        _SFD_SET_DATA_VALUE_PTR(14U, (void *)
          chartInstance->c1__indexQ_targ_address);
        _SFD_SET_DATA_VALUE_PTR(12U, (void *)chartInstance->c1__indexQ_address);
        _SFD_SET_DATA_VALUE_PTR(0U, (void *)chartInstance->c1_power);
        _SFD_SET_DATA_VALUE_PTR(3U, (void *)chartInstance->c1_B3);
        _SFD_SET_DATA_VALUE_PTR(4U, (void *)chartInstance->c1_K3);
      }
    }
  }
}

static const char* sf_get_instance_specialization(void)
{
  return "s0hyMAMnLtbdHXY2v4G8nH";
}

static void sf_opaque_initialize_c1_OVE_RM3(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc1_OVE_RM3InstanceStruct*) chartInstanceVar)->S,
    0);
  initialize_params_c1_OVE_RM3((SFc1_OVE_RM3InstanceStruct*) chartInstanceVar);
  initialize_c1_OVE_RM3((SFc1_OVE_RM3InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c1_OVE_RM3(void *chartInstanceVar)
{
  enable_c1_OVE_RM3((SFc1_OVE_RM3InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c1_OVE_RM3(void *chartInstanceVar)
{
  disable_c1_OVE_RM3((SFc1_OVE_RM3InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c1_OVE_RM3(void *chartInstanceVar)
{
  sf_gateway_c1_OVE_RM3((SFc1_OVE_RM3InstanceStruct*) chartInstanceVar);
}

static const mxArray* sf_opaque_get_sim_state_c1_OVE_RM3(SimStruct* S)
{
  return get_sim_state_c1_OVE_RM3((SFc1_OVE_RM3InstanceStruct *)
    sf_get_chart_instance_ptr(S));     /* raw sim ctx */
}

static void sf_opaque_set_sim_state_c1_OVE_RM3(SimStruct* S, const mxArray *st)
{
  set_sim_state_c1_OVE_RM3((SFc1_OVE_RM3InstanceStruct*)
    sf_get_chart_instance_ptr(S), st);
}

static void sf_opaque_terminate_c1_OVE_RM3(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc1_OVE_RM3InstanceStruct*) chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_OVE_RM3_optimization_info();
    }

    finalize_c1_OVE_RM3((SFc1_OVE_RM3InstanceStruct*) chartInstanceVar);
    utFree(chartInstanceVar);
    if (ssGetUserData(S)!= NULL) {
      sf_free_ChartRunTimeInfo(S);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc1_OVE_RM3((SFc1_OVE_RM3InstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c1_OVE_RM3(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c1_OVE_RM3((SFc1_OVE_RM3InstanceStruct*)
      sf_get_chart_instance_ptr(S));
  }
}

static void mdlSetWorkWidths_c1_OVE_RM3(SimStruct *S)
{
  /* Set overwritable ports for inplace optimization */
  ssSetInputPortDirectFeedThrough(S, 0, 1);
  ssSetStatesModifiedOnlyInUpdate(S, 0);
  ssSetBlockIsPurelyCombinatorial_wrapper(S, 0);
  ssMdlUpdateIsEmpty(S, 1);
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_OVE_RM3_optimization_info(sim_mode_is_rtw_gen(S),
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
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,1,1);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,1,4);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=4; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 1; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,1);
    sf_register_codegen_names_for_scoped_functions_defined_by_chart(S);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(3161466152U));
  ssSetChecksum1(S,(4019633090U));
  ssSetChecksum2(S,(2426859392U));
  ssSetChecksum3(S,(2110554876U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSetStateSemanticsClassicAndSynchronous(S, true);
  ssSupportsMultipleExecInstances(S,0);
}

static void mdlRTW_c1_OVE_RM3(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c1_OVE_RM3(SimStruct *S)
{
  SFc1_OVE_RM3InstanceStruct *chartInstance;
  chartInstance = (SFc1_OVE_RM3InstanceStruct *)utMalloc(sizeof
    (SFc1_OVE_RM3InstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  memset(chartInstance, 0, sizeof(SFc1_OVE_RM3InstanceStruct));
  chartInstance->chartInfo.chartInstance = chartInstance;
  if (ssGetSampleTime(S, 0) == CONTINUOUS_SAMPLE_TIME && ssGetOffsetTime(S, 0) ==
      0 && ssGetNumContStates(ssGetRootSS(S)) > 0) {
    sf_error_out_about_continuous_sample_time_with_persistent_vars(S);
  }

  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c1_OVE_RM3;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c1_OVE_RM3;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c1_OVE_RM3;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c1_OVE_RM3;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c1_OVE_RM3;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c1_OVE_RM3;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c1_OVE_RM3;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c1_OVE_RM3;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c1_OVE_RM3;
  chartInstance->chartInfo.mdlStart = mdlStart_c1_OVE_RM3;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c1_OVE_RM3;
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
  mdl_start_c1_OVE_RM3(chartInstance);
}

void c1_OVE_RM3_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c1_OVE_RM3(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c1_OVE_RM3(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c1_OVE_RM3(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c1_OVE_RM3_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
