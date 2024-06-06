/*
 * reaction_pendulum.h
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "reaction_pendulum".
 *
 * Model version              : 7.4
 * Simulink Coder version : 9.6 (R2021b) 14-May-2021
 * C source code generated on : Thu Jun  6 17:09:01 2024
 *
 * Target selection: rtcon_rpend_usb2.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_reaction_pendulum_h_
#define RTW_HEADER_reaction_pendulum_h_
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <float.h>
#ifndef reaction_pendulum_COMMON_INCLUDES_
#define reaction_pendulum_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "simstruc.h"
#include "fixedpoint.h"
#include "rt_logging.h"
#include "dt_info.h"
#include "ext_work.h"
#endif                                 /* reaction_pendulum_COMMON_INCLUDES_ */

#include "reaction_pendulum_types.h"

/* Shared type includes */
#include "multiword_types.h"
#include "rtGetInf.h"
#include "rt_nonfinite.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm) ((rtm)->CTOutputIncnstWithState)
#endif

#ifndef rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm, val) ((rtm)->CTOutputIncnstWithState = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
#define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
#define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetFinalTime
#define rtmGetFinalTime(rtm)           ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetRTWExtModeInfo
#define rtmGetRTWExtModeInfo(rtm)      ((rtm)->extModeInfo)
#endif

#ifndef rtmGetRTWLogInfo
#define rtmGetRTWLogInfo(rtm)          ((rtm)->rtwLogInfo)
#endif

#ifndef rtmGetSampleHitArray
#define rtmGetSampleHitArray(rtm)      ((rtm)->Timing.sampleHitArray)
#endif

#ifndef rtmGetStepSize
#define rtmGetStepSize(rtm)            ((rtm)->Timing.stepSize)
#endif

#ifndef rtmGetZCCacheNeedsReset
#define rtmGetZCCacheNeedsReset(rtm)   ((rtm)->zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
#define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->zCCacheNeedsReset = (val))
#endif

#ifndef rtmGet_TimeOfLastOutput
#define rtmGet_TimeOfLastOutput(rtm)   ((rtm)->Timing.timeOfLastOutput)
#endif

#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
#define rtmGetStopRequested(rtm)       ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
#define rtmSetStopRequested(rtm, val)  ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
#define rtmGetStopRequestedPtr(rtm)    (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
#define rtmGetT(rtm)                   (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTFinal
#define rtmGetTFinal(rtm)              ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetTPtr
#define rtmGetTPtr(rtm)                ((rtm)->Timing.t)
#endif

#ifndef rtmGetTStart
#define rtmGetTStart(rtm)              ((rtm)->Timing.tStart)
#endif

#ifndef rtmGetTimeOfLastOutput
#define rtmGetTimeOfLastOutput(rtm)    ((rtm)->Timing.timeOfLastOutput)
#endif

/* Block signals (default storage) */
typedef struct {
  real_T DelayOneStep[3];              /* '<Root>/Delay One Step' */
  real_T Sum;                          /* '<Root>/Sum' */
  real_T Saturation;                   /* '<S2>/Saturation' */
  real_T Gain[2];                      /* '<S2>/Gain' */
  real_T Prescaler;                    /* '<S2>/Prescaler' */
  real_T ThermFlag;                    /* '<S2>/ThermFlag' */
  real_T SFunction_o1;                 /* '<S2>/S-Function' */
  real_T SFunction_o2;                 /* '<S2>/S-Function' */
  real_T SFunction_o3;                 /* '<S2>/S-Function' */
  real_T SFunction_o4;                 /* '<S2>/S-Function' */
  real_T SFunction_o5;                 /* '<S2>/S-Function' */
  real_T SFunction_o6;                 /* '<S2>/S-Function' */
  real_T SFunction_o7[2];              /* '<S2>/S-Function' */
  real_T PendulumAnglerad;             /* '<S2>/Pendulum Convert to rad' */
  real_T Periodms;                     /* '<S2>/Gain1' */
  real_T DCVelrads;                    /* '<S2>/Divide' */
  real_T Reshapexhat[3];               /* '<S1>/Reshapexhat' */
  real_T DCConverttoA1;                /* '<S2>/DC Convert to [A]1' */
  real_T Gain_i;                       /* '<Root>/Gain' */
  real_T Sum1;                         /* '<Root>/Sum1' */
  real_T thetadotdiff;                 /* '<Root>/Derivative' */
  real_T Product2[3];                  /* '<S56>/Product2' */
  real_T Product3[3];                  /* '<S54>/Product3' */
} B_reaction_pendulum_T;

/* Block states (default storage) for system '<Root>' */
typedef struct {
  real_T DelayOneStep_DSTATE[3];       /* '<Root>/Delay One Step' */
  real_T MemoryX_DSTATE[3];            /* '<S1>/MemoryX' */
  real_T Memory1_PreviousInput;        /* '<S2>/Memory1' */
  real_T Memory_PreviousInput;         /* '<S2>/Memory' */
  real_T TimeStampA;                   /* '<Root>/Derivative' */
  real_T LastUAtTimeA;                 /* '<Root>/Derivative' */
  real_T TimeStampB;                   /* '<Root>/Derivative' */
  real_T LastUAtTimeB;                 /* '<Root>/Derivative' */
  struct {
    void *LoggedData[2];
  } Phidot_PWORK;                      /* '<Root>/Phi dot' */

  struct {
    void *LoggedData[5];
  } PlotState_PWORK;                   /* '<Root>/PlotState' */

  struct {
    void *LoggedData;
  } Scope1_PWORK;                      /* '<Root>/Scope1' */

  struct {
    void *LoggedData;
  } Scope2_PWORK;                      /* '<Root>/Scope2' */

  struct {
    void *LoggedData;
  } Scope3_PWORK;                      /* '<Root>/Scope3' */

  struct {
    void *LoggedData;
  } Scope4_PWORK;                      /* '<Root>/Scope4' */

  struct {
    void *LoggedData;
  } Scope5_PWORK;                      /* '<Root>/Scope5' */

  struct {
    void *LoggedData;
  } Scope6_PWORK;                      /* '<Root>/Scope6' */

  struct {
    void *LoggedData;
  } Scope7_PWORK;                      /* '<Root>/Scope7' */

  struct {
    void *LoggedData;
  } theta_PWORK;                       /* '<Root>/theta' */

  struct {
    void *LoggedData;
  } thetadot_PWORK;                    /* '<Root>/theta dot' */

  struct {
    void *LoggedData;
  } theta1_PWORK;                      /* '<Root>/theta1' */

  int8_T EnabledSubsystem_SubsysRanBC; /* '<S30>/Enabled Subsystem' */
  int8_T MeasurementUpdate_SubsysRanBC;/* '<S23>/MeasurementUpdate' */
  boolean_T icLoad;                    /* '<S1>/MemoryX' */
  boolean_T EnabledSubsystem_MODE;     /* '<S30>/Enabled Subsystem' */
  boolean_T MeasurementUpdate_MODE;    /* '<S23>/MeasurementUpdate' */
} DW_reaction_pendulum_T;

/* Parameters (default storage) */
struct P_reaction_pendulum_T_ {
  real_T K[3];                         /* Variable: K
                                        * Referenced by: '<Root>/Gain'
                                        */
  real_T SFunction_P2_Size[2];         /* Computed Parameter: SFunction_P2_Size
                                        * Referenced by: '<S2>/S-Function'
                                        */
  real_T h;                            /* Variable: h
                                        * Referenced by: '<S2>/S-Function'
                                        */
  real_T Lykyhatkk1_Y0;                /* Expression: 0
                                        * Referenced by: '<S54>/L*(y[k]-yhat[k|k-1])'
                                        */
  real_T deltax_Y0;                    /* Expression: 0
                                        * Referenced by: '<S56>/deltax'
                                        */
  real_T Gain1_Gain;                   /* Expression: 0.3
                                        * Referenced by: '<Root>/Gain1'
                                        */
  real_T Reset_Value;                  /* Expression: 1
                                        * Referenced by: '<Root>/Reset'
                                        */
  real_T Normal_Value;                 /* Expression: 0
                                        * Referenced by: '<Root>/Normal'
                                        */
  real_T DelayOneStep_InitialCondition;/* Expression: 0.0
                                        * Referenced by: '<Root>/Delay One Step'
                                        */
  real_T KalmanGainM_Value[6];         /* Expression: pInitialization.M
                                        * Referenced by: '<S3>/KalmanGainM'
                                        */
  real_T C_Value[6];                   /* Expression: pInitialization.C
                                        * Referenced by: '<S1>/C'
                                        */
  real_T DC_Ctrl1_Value;               /* Expression: 0
                                        * Referenced by: '<Root>/DC_Ctrl1'
                                        */
  real_T SignalGenerator_Amplitude;    /* Expression: 1
                                        * Referenced by: '<Root>/Signal Generator'
                                        */
  real_T SignalGenerator_Frequency;    /* Expression: 1/2
                                        * Referenced by: '<Root>/Signal Generator'
                                        */
  real_T Gain2_Gain;                   /* Expression: 1
                                        * Referenced by: '<S2>/Gain2'
                                        */
  real_T Saturation_UpperSat;          /* Expression: 1
                                        * Referenced by: '<S2>/Saturation'
                                        */
  real_T Saturation_LowerSat;          /* Expression: -1
                                        * Referenced by: '<S2>/Saturation'
                                        */
  real_T Gain_Gain[2];                 /* Expression: [1;1]
                                        * Referenced by: '<S2>/Gain'
                                        */
  real_T Prescaler_Value;              /* Expression: 0
                                        * Referenced by: '<S2>/Prescaler'
                                        */
  real_T ThermFlag_Value;              /* Expression: 0
                                        * Referenced by: '<S2>/ThermFlag'
                                        */
  real_T SFunction_P1_Size[2];         /* Computed Parameter: SFunction_P1_Size
                                        * Referenced by: '<S2>/S-Function'
                                        */
  real_T SFunction_P1;                 /* Expression: 1
                                        * Referenced by: '<S2>/S-Function'
                                        */
  real_T PendulumConverttorad_Gain;    /* Expression: 2*pi/20000
                                        * Referenced by: '<S2>/Pendulum Convert to rad'
                                        */
  real_T DCConverttorad_Gain;          /* Expression: 2*pi/4096
                                        * Referenced by: '<S2>/DC Convert to rad'
                                        */
  real_T Memory1_InitialCondition;     /* Expression: 0
                                        * Referenced by: '<S2>/Memory1'
                                        */
  real_T rad2RPM_Gain;                 /* Expression: 1
                                        * Referenced by: '<S2>/rad2RPM'
                                        */
  real_T Memory_InitialCondition;      /* Expression: 0
                                        * Referenced by: '<S2>/Memory'
                                        */
  real_T Gain1_Gain_h;                 /* Expression: 1/20000000
                                        * Referenced by: '<S2>/Gain1'
                                        */
  real_T X0_Value[3];                  /* Expression: pInitialization.X0
                                        * Referenced by: '<S1>/X0'
                                        */
  real_T DCConverttoA1_Gain;           /* Expression: 0.8333
                                        * Referenced by: '<S2>/DC Convert to [A]1'
                                        */
  real_T Constant_Value;               /* Expression: 0
                                        * Referenced by: '<Root>/Constant'
                                        */
  real_T A_Value[9];                   /* Expression: pInitialization.A
                                        * Referenced by: '<S1>/A'
                                        */
  real_T B_Value[3];                   /* Expression: pInitialization.B
                                        * Referenced by: '<S1>/B'
                                        */
  real_T CovarianceZ_Value[9];         /* Expression: pInitialization.Z
                                        * Referenced by: '<S3>/CovarianceZ'
                                        */
  real_T KalmanGainL_Value[6];         /* Expression: pInitialization.L
                                        * Referenced by: '<S3>/KalmanGainL'
                                        */
  real_T D_Value[2];                   /* Expression: pInitialization.D
                                        * Referenced by: '<S1>/D'
                                        */
  boolean_T Enable_Value;              /* Expression: true()
                                        * Referenced by: '<S1>/Enable'
                                        */
  boolean_T isSqrtUsed_Value;          /* Expression: pInitialization.isSqrtUsed
                                        * Referenced by: '<S52>/isSqrtUsed'
                                        */
  uint8_T ResetEncoders2_CurrentSetting;
                            /* Computed Parameter: ResetEncoders2_CurrentSetting
                             * Referenced by: '<Root>/Reset Encoders2'
                             */
  uint8_T ResetEncoders_CurrentSetting;
                             /* Computed Parameter: ResetEncoders_CurrentSetting
                              * Referenced by: '<Root>/Reset Encoders'
                              */
};

/* Real-time Model Data Structure */
struct tag_RTM_reaction_pendulum_T {
  struct SimStruct_tag * *childSfunctions;
  const char_T *errorStatus;
  SS_SimMode simMode;
  RTWLogInfo *rtwLogInfo;
  RTWExtModeInfo *extModeInfo;
  RTWSolverInfo solverInfo;
  RTWSolverInfo *solverInfoPtr;
  void *sfcnInfo;

  /*
   * NonInlinedSFcns:
   * The following substructure contains information regarding
   * non-inlined s-functions used in the model.
   */
  struct {
    RTWSfcnInfo sfcnInfo;
    time_T *taskTimePtrs[2];
    SimStruct childSFunctions[1];
    SimStruct *childSFunctionPtrs[1];
    struct _ssBlkInfo2 blkInfo2[1];
    struct _ssSFcnModelMethods2 methods2[1];
    struct _ssSFcnModelMethods3 methods3[1];
    struct _ssSFcnModelMethods4 methods4[1];
    struct _ssStatesInfo2 statesInfo2[1];
    ssPeriodicStatesInfo periodicStatesInfo[1];
    struct _ssPortInfo2 inputOutputPortInfo2[1];
    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[4];
      struct _ssInPortUnit inputPortUnits[4];
      struct _ssInPortCoSimAttribute inputPortCoSimAttribute[4];
      real_T const *UPtrs0[1];
      real_T const *UPtrs1[2];
      real_T const *UPtrs2[1];
      real_T const *UPtrs3[1];
      struct _ssPortOutputs outputPortInfo[7];
      struct _ssOutPortUnit outputPortUnits[7];
      struct _ssOutPortCoSimAttribute outputPortCoSimAttribute[7];
      uint_T attribs[2];
      mxArray *params[2];
    } Sfcn0;
  } NonInlinedSFcns;

  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T CTOutputIncnstWithState;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    uint32_T checksums[4];
    uint32_T options;
    int_T numContStates;
    int_T numU;
    int_T numY;
    int_T numSampTimes;
    int_T numBlocks;
    int_T numBlockIO;
    int_T numBlockPrms;
    int_T numDwork;
    int_T numSFcnPrms;
    int_T numSFcns;
    int_T numIports;
    int_T numOports;
    int_T numNonSampZCs;
    int_T sysDirFeedThru;
    int_T rtwGenSfcn;
  } Sizes;

  /*
   * SpecialInfo:
   * The following substructure contains special information
   * related to other components that are dependent on RTW.
   */
  struct {
    const void *mappingInfo;
  } SpecialInfo;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    time_T stepSize;
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    uint32_T clockTick1;
    uint32_T clockTickH1;
    time_T stepSize1;
    time_T tStart;
    time_T tFinal;
    time_T timeOfLastOutput;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *sampleTimes;
    time_T *offsetTimes;
    int_T *sampleTimeTaskIDPtr;
    int_T *sampleHits;
    int_T *perTaskSampleHits;
    time_T *t;
    time_T sampleTimesArray[2];
    time_T offsetTimesArray[2];
    int_T sampleTimeTaskIDArray[2];
    int_T sampleHitArray[2];
    int_T perTaskSampleHitsArray[4];
    time_T tArray[2];
  } Timing;
};

/* Block parameters (default storage) */
extern P_reaction_pendulum_T reaction_pendulum_P;

/* Block signals (default storage) */
extern B_reaction_pendulum_T reaction_pendulum_B;

/* Block states (default storage) */
extern DW_reaction_pendulum_T reaction_pendulum_DW;

/* Model entry point functions */
extern void reaction_pendulum_initialize(void);
extern void reaction_pendulum_step(void);
extern void reaction_pendulum_terminate(void);

/* Real-time Model object */
extern RT_MODEL_reaction_pendulum_T *const reaction_pendulum_M;

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'reaction_pendulum'
 * '<S1>'   : 'reaction_pendulum/Kalman Filter'
 * '<S2>'   : 'reaction_pendulum/RP Driver'
 * '<S3>'   : 'reaction_pendulum/Kalman Filter/CalculatePL'
 * '<S4>'   : 'reaction_pendulum/Kalman Filter/CalculateYhat'
 * '<S5>'   : 'reaction_pendulum/Kalman Filter/CovarianceOutputConfigurator'
 * '<S6>'   : 'reaction_pendulum/Kalman Filter/DataTypeConversionA'
 * '<S7>'   : 'reaction_pendulum/Kalman Filter/DataTypeConversionB'
 * '<S8>'   : 'reaction_pendulum/Kalman Filter/DataTypeConversionC'
 * '<S9>'   : 'reaction_pendulum/Kalman Filter/DataTypeConversionD'
 * '<S10>'  : 'reaction_pendulum/Kalman Filter/DataTypeConversionEnable'
 * '<S11>'  : 'reaction_pendulum/Kalman Filter/DataTypeConversionG'
 * '<S12>'  : 'reaction_pendulum/Kalman Filter/DataTypeConversionH'
 * '<S13>'  : 'reaction_pendulum/Kalman Filter/DataTypeConversionN'
 * '<S14>'  : 'reaction_pendulum/Kalman Filter/DataTypeConversionP'
 * '<S15>'  : 'reaction_pendulum/Kalman Filter/DataTypeConversionP0'
 * '<S16>'  : 'reaction_pendulum/Kalman Filter/DataTypeConversionQ'
 * '<S17>'  : 'reaction_pendulum/Kalman Filter/DataTypeConversionR'
 * '<S18>'  : 'reaction_pendulum/Kalman Filter/DataTypeConversionReset'
 * '<S19>'  : 'reaction_pendulum/Kalman Filter/DataTypeConversionX'
 * '<S20>'  : 'reaction_pendulum/Kalman Filter/DataTypeConversionX0'
 * '<S21>'  : 'reaction_pendulum/Kalman Filter/DataTypeConversionu'
 * '<S22>'  : 'reaction_pendulum/Kalman Filter/MemoryP'
 * '<S23>'  : 'reaction_pendulum/Kalman Filter/Observer'
 * '<S24>'  : 'reaction_pendulum/Kalman Filter/ReducedQRN'
 * '<S25>'  : 'reaction_pendulum/Kalman Filter/Reset'
 * '<S26>'  : 'reaction_pendulum/Kalman Filter/Reshapeyhat'
 * '<S27>'  : 'reaction_pendulum/Kalman Filter/ScalarExpansionP0'
 * '<S28>'  : 'reaction_pendulum/Kalman Filter/ScalarExpansionQ'
 * '<S29>'  : 'reaction_pendulum/Kalman Filter/ScalarExpansionR'
 * '<S30>'  : 'reaction_pendulum/Kalman Filter/UseCurrentEstimator'
 * '<S31>'  : 'reaction_pendulum/Kalman Filter/checkA'
 * '<S32>'  : 'reaction_pendulum/Kalman Filter/checkB'
 * '<S33>'  : 'reaction_pendulum/Kalman Filter/checkC'
 * '<S34>'  : 'reaction_pendulum/Kalman Filter/checkD'
 * '<S35>'  : 'reaction_pendulum/Kalman Filter/checkEnable'
 * '<S36>'  : 'reaction_pendulum/Kalman Filter/checkG'
 * '<S37>'  : 'reaction_pendulum/Kalman Filter/checkH'
 * '<S38>'  : 'reaction_pendulum/Kalman Filter/checkN'
 * '<S39>'  : 'reaction_pendulum/Kalman Filter/checkP0'
 * '<S40>'  : 'reaction_pendulum/Kalman Filter/checkQ'
 * '<S41>'  : 'reaction_pendulum/Kalman Filter/checkR'
 * '<S42>'  : 'reaction_pendulum/Kalman Filter/checkReset'
 * '<S43>'  : 'reaction_pendulum/Kalman Filter/checkX0'
 * '<S44>'  : 'reaction_pendulum/Kalman Filter/checku'
 * '<S45>'  : 'reaction_pendulum/Kalman Filter/checky'
 * '<S46>'  : 'reaction_pendulum/Kalman Filter/CalculatePL/DataTypeConversionL'
 * '<S47>'  : 'reaction_pendulum/Kalman Filter/CalculatePL/DataTypeConversionM'
 * '<S48>'  : 'reaction_pendulum/Kalman Filter/CalculatePL/DataTypeConversionP'
 * '<S49>'  : 'reaction_pendulum/Kalman Filter/CalculatePL/DataTypeConversionZ'
 * '<S50>'  : 'reaction_pendulum/Kalman Filter/CalculatePL/Ground'
 * '<S51>'  : 'reaction_pendulum/Kalman Filter/CalculateYhat/Ground'
 * '<S52>'  : 'reaction_pendulum/Kalman Filter/CovarianceOutputConfigurator/decideOutput'
 * '<S53>'  : 'reaction_pendulum/Kalman Filter/CovarianceOutputConfigurator/decideOutput/SqrtUsedFcn'
 * '<S54>'  : 'reaction_pendulum/Kalman Filter/Observer/MeasurementUpdate'
 * '<S55>'  : 'reaction_pendulum/Kalman Filter/ReducedQRN/Ground'
 * '<S56>'  : 'reaction_pendulum/Kalman Filter/UseCurrentEstimator/Enabled Subsystem'
 */
#endif                                 /* RTW_HEADER_reaction_pendulum_h_ */
