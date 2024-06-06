/*
 * reaction_pendulum.c
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

#include "reaction_pendulum.h"
#include "reaction_pendulum_private.h"
#include "reaction_pendulum_dt.h"

/* Block signals (default storage) */
B_reaction_pendulum_T reaction_pendulum_B;

/* Block states (default storage) */
DW_reaction_pendulum_T reaction_pendulum_DW;

/* Real-time model */
static RT_MODEL_reaction_pendulum_T reaction_pendulum_M_;
RT_MODEL_reaction_pendulum_T *const reaction_pendulum_M = &reaction_pendulum_M_;

/* Model step function */
void reaction_pendulum_step(void)
{
  /* local block i/o variables */
  real_T rtb_DCAnglerad;
  real_T rtb_Add_j[3];
  real_T rtb_TmpSignalConversionAtDelayO[3];

  {
    real_T rtb_Sum[2];
    real_T tmp[2];
    real_T rtb_Gain2;
    real_T tmp_0;
    real_T *lastU;
    int32_T i;

    /* Reset subsysRan breadcrumbs */
    srClearBC(reaction_pendulum_DW.MeasurementUpdate_SubsysRanBC);

    /* Reset subsysRan breadcrumbs */
    srClearBC(reaction_pendulum_DW.EnabledSubsystem_SubsysRanBC);

    /* Delay: '<Root>/Delay One Step' */
    reaction_pendulum_B.DelayOneStep[0] =
      reaction_pendulum_DW.DelayOneStep_DSTATE[0];
    reaction_pendulum_B.DelayOneStep[1] =
      reaction_pendulum_DW.DelayOneStep_DSTATE[1];
    reaction_pendulum_B.DelayOneStep[2] =
      reaction_pendulum_DW.DelayOneStep_DSTATE[2];

    /* ManualSwitch: '<Root>/Reset Encoders2' */
    if (reaction_pendulum_P.ResetEncoders2_CurrentSetting == 1) {
      /* Sum: '<Root>/Sum' incorporates:
       *  Constant: '<Root>/DC_Ctrl1'
       */
      reaction_pendulum_B.Sum = reaction_pendulum_P.DC_Ctrl1_Value;
    } else {
      /* Sum: '<Root>/Sum' incorporates:
       *  Gain: '<Root>/Gain1'
       *  SignalGenerator: '<Root>/Signal Generator'
       */
      reaction_pendulum_B.Sum = sin(6.2831853071795862 *
        reaction_pendulum_M->Timing.t[0] *
        reaction_pendulum_P.SignalGenerator_Frequency) *
        reaction_pendulum_P.SignalGenerator_Amplitude *
        reaction_pendulum_P.Gain1_Gain;
    }

    /* End of ManualSwitch: '<Root>/Reset Encoders2' */

    /* Gain: '<S2>/Gain2' */
    rtb_Gain2 = reaction_pendulum_P.Gain2_Gain * reaction_pendulum_B.Sum;

    /* Saturate: '<S2>/Saturation' */
    if (rtb_Gain2 > reaction_pendulum_P.Saturation_UpperSat) {
      /* Saturate: '<S2>/Saturation' */
      reaction_pendulum_B.Saturation = reaction_pendulum_P.Saturation_UpperSat;
    } else if (rtb_Gain2 < reaction_pendulum_P.Saturation_LowerSat) {
      /* Saturate: '<S2>/Saturation' */
      reaction_pendulum_B.Saturation = reaction_pendulum_P.Saturation_LowerSat;
    } else {
      /* Saturate: '<S2>/Saturation' */
      reaction_pendulum_B.Saturation = rtb_Gain2;
    }

    /* End of Saturate: '<S2>/Saturation' */

    /* ManualSwitch: '<Root>/Reset Encoders' incorporates:
     *  Constant: '<Root>/Normal'
     *  Constant: '<Root>/Reset'
     */
    if (reaction_pendulum_P.ResetEncoders_CurrentSetting == 1) {
      rtb_Gain2 = reaction_pendulum_P.Reset_Value;
    } else {
      rtb_Gain2 = reaction_pendulum_P.Normal_Value;
    }

    /* End of ManualSwitch: '<Root>/Reset Encoders' */

    /* Gain: '<S2>/Gain' */
    reaction_pendulum_B.Gain[0] = reaction_pendulum_P.Gain_Gain[0] * rtb_Gain2;
    reaction_pendulum_B.Gain[1] = reaction_pendulum_P.Gain_Gain[1] * rtb_Gain2;

    /* Constant: '<S2>/Prescaler' */
    reaction_pendulum_B.Prescaler = reaction_pendulum_P.Prescaler_Value;

    /* Constant: '<S2>/ThermFlag' */
    reaction_pendulum_B.ThermFlag = reaction_pendulum_P.ThermFlag_Value;

    /* S-Function (rtdacusb2_rpend_dd): '<S2>/S-Function' */

    /* Level2 S-Function Block: '<S2>/S-Function' (rtdacusb2_rpend_dd) */
    {
      SimStruct *rts = reaction_pendulum_M->childSfunctions[0];
      sfcnOutputs(rts,0);
    }

    /* Gain: '<S2>/Pendulum Convert to rad' */
    reaction_pendulum_B.PendulumAnglerad =
      reaction_pendulum_P.PendulumConverttorad_Gain *
      reaction_pendulum_B.SFunction_o2;

    /* Gain: '<S2>/DC Convert to rad' */
    rtb_DCAnglerad = reaction_pendulum_P.DCConverttorad_Gain *
      reaction_pendulum_B.SFunction_o3;

    /* Gain: '<S2>/Gain1' incorporates:
     *  Memory: '<S2>/Memory'
     *  Sum: '<S2>/Add'
     */
    reaction_pendulum_B.Periodms = (reaction_pendulum_B.SFunction_o6 -
      reaction_pendulum_DW.Memory_PreviousInput) *
      reaction_pendulum_P.Gain1_Gain_h;

    /* Product: '<S2>/Divide' incorporates:
     *  Gain: '<S2>/rad2RPM'
     *  Memory: '<S2>/Memory1'
     *  Sum: '<S2>/Add1'
     */
    reaction_pendulum_B.DCVelrads = (rtb_DCAnglerad -
      reaction_pendulum_DW.Memory1_PreviousInput) *
      reaction_pendulum_P.rad2RPM_Gain / reaction_pendulum_B.Periodms;

    /* Delay: '<S1>/MemoryX' incorporates:
     *  Constant: '<S1>/X0'
     */
    if (reaction_pendulum_DW.icLoad) {
      reaction_pendulum_DW.MemoryX_DSTATE[0] = reaction_pendulum_P.X0_Value[0];
      reaction_pendulum_DW.MemoryX_DSTATE[1] = reaction_pendulum_P.X0_Value[1];
      reaction_pendulum_DW.MemoryX_DSTATE[2] = reaction_pendulum_P.X0_Value[2];
    }

    /* Outputs for Enabled SubSystem: '<S30>/Enabled Subsystem' incorporates:
     *  EnablePort: '<S56>/Enable'
     */
    /* Constant: '<S1>/Enable' */
    if (reaction_pendulum_P.Enable_Value) {
      reaction_pendulum_DW.EnabledSubsystem_MODE = true;

      /* Product: '<S56>/Product' incorporates:
       *  Constant: '<S1>/C'
       *  Delay: '<S1>/MemoryX'
       */
      for (i = 0; i < 2; i++) {
        rtb_Sum[i] = (reaction_pendulum_P.C_Value[i + 2] *
                      reaction_pendulum_DW.MemoryX_DSTATE[1] +
                      reaction_pendulum_P.C_Value[i] *
                      reaction_pendulum_DW.MemoryX_DSTATE[0]) +
          reaction_pendulum_P.C_Value[i + 4] *
          reaction_pendulum_DW.MemoryX_DSTATE[2];
      }

      /* End of Product: '<S56>/Product' */

      /* Sum: '<S56>/Add1' */
      rtb_Gain2 = reaction_pendulum_B.PendulumAnglerad - rtb_Sum[0];
      tmp_0 = reaction_pendulum_B.DCVelrads - rtb_Sum[1];
      for (i = 0; i < 3; i++) {
        /* Product: '<S56>/Product2' incorporates:
         *  Constant: '<S3>/KalmanGainM'
         */
        reaction_pendulum_B.Product2[i] = 0.0;
        reaction_pendulum_B.Product2[i] +=
          reaction_pendulum_P.KalmanGainM_Value[i] * rtb_Gain2;
        reaction_pendulum_B.Product2[i] +=
          reaction_pendulum_P.KalmanGainM_Value[i + 3] * tmp_0;
      }

      srUpdateBC(reaction_pendulum_DW.EnabledSubsystem_SubsysRanBC);
    } else if (reaction_pendulum_DW.EnabledSubsystem_MODE) {
      /* Disable for Product: '<S56>/Product2' incorporates:
       *  Outport: '<S56>/deltax'
       */
      reaction_pendulum_B.Product2[0] = reaction_pendulum_P.deltax_Y0;
      reaction_pendulum_B.Product2[1] = reaction_pendulum_P.deltax_Y0;
      reaction_pendulum_B.Product2[2] = reaction_pendulum_P.deltax_Y0;
      reaction_pendulum_DW.EnabledSubsystem_MODE = false;
    }

    /* End of Outputs for SubSystem: '<S30>/Enabled Subsystem' */

    /* Reshape: '<S1>/Reshapexhat' incorporates:
     *  Delay: '<S1>/MemoryX'
     *  Sum: '<S30>/Add'
     */
    reaction_pendulum_B.Reshapexhat[0] = reaction_pendulum_B.Product2[0] +
      reaction_pendulum_DW.MemoryX_DSTATE[0];
    reaction_pendulum_B.Reshapexhat[1] = reaction_pendulum_B.Product2[1] +
      reaction_pendulum_DW.MemoryX_DSTATE[1];
    reaction_pendulum_B.Reshapexhat[2] = reaction_pendulum_B.Product2[2] +
      reaction_pendulum_DW.MemoryX_DSTATE[2];

    /* Gain: '<S2>/DC Convert to [A]1' */
    reaction_pendulum_B.DCConverttoA1 = reaction_pendulum_P.DCConverttoA1_Gain *
      reaction_pendulum_B.SFunction_o4;

    /* Gain: '<Root>/Gain' */
    reaction_pendulum_B.Gain_i = (-reaction_pendulum_P.K[0] *
      reaction_pendulum_B.DelayOneStep[0] + -reaction_pendulum_P.K[1] *
      reaction_pendulum_B.DelayOneStep[1]) + -reaction_pendulum_P.K[2] *
      reaction_pendulum_B.DelayOneStep[2];

    /* Sum: '<Root>/Sum1' incorporates:
     *  Constant: '<Root>/Constant'
     */
    reaction_pendulum_B.Sum1 = reaction_pendulum_P.Constant_Value -
      reaction_pendulum_B.PendulumAnglerad;

    /* Derivative: '<Root>/Derivative' */
    rtb_Gain2 = reaction_pendulum_M->Timing.t[0];
    if ((reaction_pendulum_DW.TimeStampA >= rtb_Gain2) &&
        (reaction_pendulum_DW.TimeStampB >= rtb_Gain2)) {
      /* Derivative: '<Root>/Derivative' */
      reaction_pendulum_B.thetadotdiff = 0.0;
    } else {
      tmp_0 = reaction_pendulum_DW.TimeStampA;
      lastU = &reaction_pendulum_DW.LastUAtTimeA;
      if (reaction_pendulum_DW.TimeStampA < reaction_pendulum_DW.TimeStampB) {
        if (reaction_pendulum_DW.TimeStampB < rtb_Gain2) {
          tmp_0 = reaction_pendulum_DW.TimeStampB;
          lastU = &reaction_pendulum_DW.LastUAtTimeB;
        }
      } else if (reaction_pendulum_DW.TimeStampA >= rtb_Gain2) {
        tmp_0 = reaction_pendulum_DW.TimeStampB;
        lastU = &reaction_pendulum_DW.LastUAtTimeB;
      }

      /* Derivative: '<Root>/Derivative' */
      reaction_pendulum_B.thetadotdiff = (reaction_pendulum_B.PendulumAnglerad -
        *lastU) / (rtb_Gain2 - tmp_0);
    }

    /* End of Derivative: '<Root>/Derivative' */
    /* Outputs for Enabled SubSystem: '<S23>/MeasurementUpdate' incorporates:
     *  EnablePort: '<S54>/Enable'
     */
    /* Constant: '<S1>/Enable' */
    /* MATLAB Function 'Kalman Filter/CovarianceOutputConfigurator/decideOutput/SqrtUsedFcn': '<S53>:1' */
    if (reaction_pendulum_P.Enable_Value) {
      reaction_pendulum_DW.MeasurementUpdate_MODE = true;

      /* Sum: '<S54>/Sum' */
      rtb_Sum[0] = reaction_pendulum_B.PendulumAnglerad;
      rtb_Sum[1] = reaction_pendulum_B.DCVelrads;
      for (i = 0; i < 2; i++) {
        /* Sum: '<S54>/Sum' incorporates:
         *  Constant: '<S1>/C'
         *  Constant: '<S1>/D'
         *  Delay: '<S1>/MemoryX'
         *  Product: '<S54>/C[k]*xhat[k|k-1]'
         *  Product: '<S54>/D[k]*u[k]'
         *  Sum: '<S54>/Add1'
         */
        tmp[i] = rtb_Sum[i] - (((reaction_pendulum_P.C_Value[i + 2] *
          reaction_pendulum_DW.MemoryX_DSTATE[1] + reaction_pendulum_P.C_Value[i]
          * reaction_pendulum_DW.MemoryX_DSTATE[0]) +
          reaction_pendulum_P.C_Value[i + 4] *
          reaction_pendulum_DW.MemoryX_DSTATE[2]) +
          reaction_pendulum_P.D_Value[i] * reaction_pendulum_B.Sum);
      }

      for (i = 0; i < 3; i++) {
        /* Product: '<S54>/Product3' incorporates:
         *  Constant: '<S3>/KalmanGainL'
         */
        reaction_pendulum_B.Product3[i] = 0.0;
        reaction_pendulum_B.Product3[i] +=
          reaction_pendulum_P.KalmanGainL_Value[i] * tmp[0];
        reaction_pendulum_B.Product3[i] +=
          reaction_pendulum_P.KalmanGainL_Value[i + 3] * tmp[1];
      }

      srUpdateBC(reaction_pendulum_DW.MeasurementUpdate_SubsysRanBC);
    } else if (reaction_pendulum_DW.MeasurementUpdate_MODE) {
      /* Disable for Product: '<S54>/Product3' incorporates:
       *  Outport: '<S54>/L*(y[k]-yhat[k|k-1])'
       */
      reaction_pendulum_B.Product3[0] = reaction_pendulum_P.Lykyhatkk1_Y0;
      reaction_pendulum_B.Product3[1] = reaction_pendulum_P.Lykyhatkk1_Y0;
      reaction_pendulum_B.Product3[2] = reaction_pendulum_P.Lykyhatkk1_Y0;
      reaction_pendulum_DW.MeasurementUpdate_MODE = false;
    }

    /* End of Outputs for SubSystem: '<S23>/MeasurementUpdate' */
    for (i = 0; i < 3; i++) {
      /* Sum: '<S23>/Add' incorporates:
       *  Constant: '<S1>/A'
       *  Constant: '<S1>/B'
       *  Delay: '<S1>/MemoryX'
       *  Product: '<S23>/A[k]*xhat[k|k-1]'
       *  Product: '<S23>/B[k]*u[k]'
       *  Product: '<S54>/Product3'
       */
      rtb_Add_j[i] = (((reaction_pendulum_P.A_Value[i + 3] *
                        reaction_pendulum_DW.MemoryX_DSTATE[1] +
                        reaction_pendulum_P.A_Value[i] *
                        reaction_pendulum_DW.MemoryX_DSTATE[0]) +
                       reaction_pendulum_P.A_Value[i + 6] *
                       reaction_pendulum_DW.MemoryX_DSTATE[2]) +
                      reaction_pendulum_P.B_Value[i] * reaction_pendulum_B.Sum)
        + reaction_pendulum_B.Product3[i];
    }

    /* SignalConversion generated from: '<Root>/Delay One Step' incorporates:
     *  ZeroOrderHold: '<Root>/Zero-Order Hold'
     */
    rtb_TmpSignalConversionAtDelayO[0] = reaction_pendulum_B.PendulumAnglerad;
    rtb_TmpSignalConversionAtDelayO[1] = reaction_pendulum_B.thetadotdiff;
    rtb_TmpSignalConversionAtDelayO[2] = reaction_pendulum_B.DCVelrads;
  }

  /* Matfile logging */
  rt_UpdateTXYLogVars(reaction_pendulum_M->rtwLogInfo,
                      (reaction_pendulum_M->Timing.t));

  {
    real_T *lastU;

    /* Update for Memory: '<S2>/Memory1' */
    reaction_pendulum_DW.Memory1_PreviousInput = rtb_DCAnglerad;

    /* Update for Memory: '<S2>/Memory' */
    reaction_pendulum_DW.Memory_PreviousInput = reaction_pendulum_B.SFunction_o6;

    /* Update for Delay: '<S1>/MemoryX' */
    reaction_pendulum_DW.icLoad = false;

    /* Update for Delay: '<Root>/Delay One Step' */
    reaction_pendulum_DW.DelayOneStep_DSTATE[0] =
      rtb_TmpSignalConversionAtDelayO[0];

    /* Update for Delay: '<S1>/MemoryX' */
    reaction_pendulum_DW.MemoryX_DSTATE[0] = rtb_Add_j[0];

    /* Update for Delay: '<Root>/Delay One Step' */
    reaction_pendulum_DW.DelayOneStep_DSTATE[1] =
      rtb_TmpSignalConversionAtDelayO[1];

    /* Update for Delay: '<S1>/MemoryX' */
    reaction_pendulum_DW.MemoryX_DSTATE[1] = rtb_Add_j[1];

    /* Update for Delay: '<Root>/Delay One Step' */
    reaction_pendulum_DW.DelayOneStep_DSTATE[2] =
      rtb_TmpSignalConversionAtDelayO[2];

    /* Update for Delay: '<S1>/MemoryX' */
    reaction_pendulum_DW.MemoryX_DSTATE[2] = rtb_Add_j[2];

    /* Update for Derivative: '<Root>/Derivative' */
    if (reaction_pendulum_DW.TimeStampA == (rtInf)) {
      reaction_pendulum_DW.TimeStampA = reaction_pendulum_M->Timing.t[0];
      lastU = &reaction_pendulum_DW.LastUAtTimeA;
    } else if (reaction_pendulum_DW.TimeStampB == (rtInf)) {
      reaction_pendulum_DW.TimeStampB = reaction_pendulum_M->Timing.t[0];
      lastU = &reaction_pendulum_DW.LastUAtTimeB;
    } else if (reaction_pendulum_DW.TimeStampA < reaction_pendulum_DW.TimeStampB)
    {
      reaction_pendulum_DW.TimeStampA = reaction_pendulum_M->Timing.t[0];
      lastU = &reaction_pendulum_DW.LastUAtTimeA;
    } else {
      reaction_pendulum_DW.TimeStampB = reaction_pendulum_M->Timing.t[0];
      lastU = &reaction_pendulum_DW.LastUAtTimeB;
    }

    *lastU = reaction_pendulum_B.PendulumAnglerad;

    /* End of Update for Derivative: '<Root>/Derivative' */
  }

  /* External mode */
  rtExtModeUploadCheckTrigger(2);

  {                                    /* Sample time: [0.0s, 0.0s] */
    rtExtModeUpload(0, (real_T)reaction_pendulum_M->Timing.t[0]);
  }

  {                                    /* Sample time: [0.05s, 0.0s] */
    rtExtModeUpload(1, (real_T)reaction_pendulum_M->Timing.t[1]);
  }

  /* signal main to stop simulation */
  {                                    /* Sample time: [0.0s, 0.0s] */
    if ((rtmGetTFinal(reaction_pendulum_M)!=-1) &&
        !((rtmGetTFinal(reaction_pendulum_M)-reaction_pendulum_M->Timing.t[0]) >
          reaction_pendulum_M->Timing.t[0] * (DBL_EPSILON))) {
      rtmSetErrorStatus(reaction_pendulum_M, "Simulation finished");
    }

    if (rtmGetStopRequested(reaction_pendulum_M)) {
      rtmSetErrorStatus(reaction_pendulum_M, "Simulation finished");
    }
  }

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++reaction_pendulum_M->Timing.clockTick0)) {
    ++reaction_pendulum_M->Timing.clockTickH0;
  }

  reaction_pendulum_M->Timing.t[0] = reaction_pendulum_M->Timing.clockTick0 *
    reaction_pendulum_M->Timing.stepSize0 +
    reaction_pendulum_M->Timing.clockTickH0 *
    reaction_pendulum_M->Timing.stepSize0 * 4294967296.0;

  {
    /* Update absolute timer for sample time: [0.05s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick1"
     * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++reaction_pendulum_M->Timing.clockTick1)) {
      ++reaction_pendulum_M->Timing.clockTickH1;
    }

    reaction_pendulum_M->Timing.t[1] = reaction_pendulum_M->Timing.clockTick1 *
      reaction_pendulum_M->Timing.stepSize1 +
      reaction_pendulum_M->Timing.clockTickH1 *
      reaction_pendulum_M->Timing.stepSize1 * 4294967296.0;
  }
}

/* Model initialize function */
void reaction_pendulum_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* initialize real-time model */
  (void) memset((void *)reaction_pendulum_M, 0,
                sizeof(RT_MODEL_reaction_pendulum_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&reaction_pendulum_M->solverInfo,
                          &reaction_pendulum_M->Timing.simTimeStep);
    rtsiSetTPtr(&reaction_pendulum_M->solverInfo, &rtmGetTPtr
                (reaction_pendulum_M));
    rtsiSetStepSizePtr(&reaction_pendulum_M->solverInfo,
                       &reaction_pendulum_M->Timing.stepSize0);
    rtsiSetErrorStatusPtr(&reaction_pendulum_M->solverInfo, (&rtmGetErrorStatus
      (reaction_pendulum_M)));
    rtsiSetRTModelPtr(&reaction_pendulum_M->solverInfo, reaction_pendulum_M);
  }

  rtsiSetSimTimeStep(&reaction_pendulum_M->solverInfo, MAJOR_TIME_STEP);
  rtsiSetSolverName(&reaction_pendulum_M->solverInfo,"FixedStepDiscrete");
  reaction_pendulum_M->solverInfoPtr = (&reaction_pendulum_M->solverInfo);

  /* Initialize timing info */
  {
    int_T *mdlTsMap = reaction_pendulum_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;

    /* polyspace +2 MISRA2012:D4.1 [Justified:Low] "reaction_pendulum_M points to
       static memory which is guaranteed to be non-NULL" */
    reaction_pendulum_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    reaction_pendulum_M->Timing.sampleTimes =
      (&reaction_pendulum_M->Timing.sampleTimesArray[0]);
    reaction_pendulum_M->Timing.offsetTimes =
      (&reaction_pendulum_M->Timing.offsetTimesArray[0]);

    /* task periods */
    reaction_pendulum_M->Timing.sampleTimes[0] = (0.0);
    reaction_pendulum_M->Timing.sampleTimes[1] = (0.05);

    /* task offsets */
    reaction_pendulum_M->Timing.offsetTimes[0] = (0.0);
    reaction_pendulum_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(reaction_pendulum_M, &reaction_pendulum_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = reaction_pendulum_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    reaction_pendulum_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(reaction_pendulum_M, 100.0);
  reaction_pendulum_M->Timing.stepSize0 = 0.05;
  reaction_pendulum_M->Timing.stepSize1 = 0.05;

  /* Setup for data logging */
  {
    static RTWLogInfo rt_DataLoggingInfo;
    rt_DataLoggingInfo.loggingInterval = (NULL);
    reaction_pendulum_M->rtwLogInfo = &rt_DataLoggingInfo;
  }

  /* Setup for data logging */
  {
    rtliSetLogXSignalInfo(reaction_pendulum_M->rtwLogInfo, (NULL));
    rtliSetLogXSignalPtrs(reaction_pendulum_M->rtwLogInfo, (NULL));
    rtliSetLogT(reaction_pendulum_M->rtwLogInfo, "");
    rtliSetLogX(reaction_pendulum_M->rtwLogInfo, "");
    rtliSetLogXFinal(reaction_pendulum_M->rtwLogInfo, "");
    rtliSetLogVarNameModifier(reaction_pendulum_M->rtwLogInfo, "rt_");
    rtliSetLogFormat(reaction_pendulum_M->rtwLogInfo, 0);
    rtliSetLogMaxRows(reaction_pendulum_M->rtwLogInfo, 0);
    rtliSetLogDecimation(reaction_pendulum_M->rtwLogInfo, 1);
    rtliSetLogY(reaction_pendulum_M->rtwLogInfo, "");
    rtliSetLogYSignalInfo(reaction_pendulum_M->rtwLogInfo, (NULL));
    rtliSetLogYSignalPtrs(reaction_pendulum_M->rtwLogInfo, (NULL));
  }

  /* External mode info */
  reaction_pendulum_M->Sizes.checksums[0] = (4043391257U);
  reaction_pendulum_M->Sizes.checksums[1] = (3157747759U);
  reaction_pendulum_M->Sizes.checksums[2] = (4109200965U);
  reaction_pendulum_M->Sizes.checksums[3] = (4046214105U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[7];
    reaction_pendulum_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    systemRan[1] = &rtAlwaysEnabled;
    systemRan[2] = (sysRanDType *)
      &reaction_pendulum_DW.MeasurementUpdate_SubsysRanBC;
    systemRan[3] = (sysRanDType *)
      &reaction_pendulum_DW.EnabledSubsystem_SubsysRanBC;
    systemRan[4] = &rtAlwaysEnabled;
    systemRan[5] = &rtAlwaysEnabled;
    systemRan[6] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(reaction_pendulum_M->extModeInfo,
      &reaction_pendulum_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(reaction_pendulum_M->extModeInfo,
                        reaction_pendulum_M->Sizes.checksums);
    rteiSetTPtr(reaction_pendulum_M->extModeInfo, rtmGetTPtr(reaction_pendulum_M));
  }

  reaction_pendulum_M->solverInfoPtr = (&reaction_pendulum_M->solverInfo);
  reaction_pendulum_M->Timing.stepSize = (0.05);
  rtsiSetFixedStepSize(&reaction_pendulum_M->solverInfo, 0.05);
  rtsiSetSolverMode(&reaction_pendulum_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  (void) memset(((void *) &reaction_pendulum_B), 0,
                sizeof(B_reaction_pendulum_T));

  /* states (dwork) */
  (void) memset((void *)&reaction_pendulum_DW, 0,
                sizeof(DW_reaction_pendulum_T));

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    reaction_pendulum_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 18;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.BTransTable = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.PTransTable = &rtPTransTable;
  }

  /* child S-Function registration */
  {
    RTWSfcnInfo *sfcnInfo = &reaction_pendulum_M->NonInlinedSFcns.sfcnInfo;
    reaction_pendulum_M->sfcnInfo = (sfcnInfo);
    rtssSetErrorStatusPtr(sfcnInfo, (&rtmGetErrorStatus(reaction_pendulum_M)));
    reaction_pendulum_M->Sizes.numSampTimes = (2);
    rtssSetNumRootSampTimesPtr(sfcnInfo,
      &reaction_pendulum_M->Sizes.numSampTimes);
    reaction_pendulum_M->NonInlinedSFcns.taskTimePtrs[0] = &(rtmGetTPtr
      (reaction_pendulum_M)[0]);
    reaction_pendulum_M->NonInlinedSFcns.taskTimePtrs[1] = &(rtmGetTPtr
      (reaction_pendulum_M)[1]);
    rtssSetTPtrPtr(sfcnInfo,reaction_pendulum_M->NonInlinedSFcns.taskTimePtrs);
    rtssSetTStartPtr(sfcnInfo, &rtmGetTStart(reaction_pendulum_M));
    rtssSetTFinalPtr(sfcnInfo, &rtmGetTFinal(reaction_pendulum_M));
    rtssSetTimeOfLastOutputPtr(sfcnInfo, &rtmGetTimeOfLastOutput
      (reaction_pendulum_M));
    rtssSetStepSizePtr(sfcnInfo, &reaction_pendulum_M->Timing.stepSize);
    rtssSetStopRequestedPtr(sfcnInfo, &rtmGetStopRequested(reaction_pendulum_M));
    rtssSetDerivCacheNeedsResetPtr(sfcnInfo,
      &reaction_pendulum_M->derivCacheNeedsReset);
    rtssSetZCCacheNeedsResetPtr(sfcnInfo,
      &reaction_pendulum_M->zCCacheNeedsReset);
    rtssSetContTimeOutputInconsistentWithStateAtMajorStepPtr(sfcnInfo,
      &reaction_pendulum_M->CTOutputIncnstWithState);
    rtssSetSampleHitsPtr(sfcnInfo, &reaction_pendulum_M->Timing.sampleHits);
    rtssSetPerTaskSampleHitsPtr(sfcnInfo,
      &reaction_pendulum_M->Timing.perTaskSampleHits);
    rtssSetSimModePtr(sfcnInfo, &reaction_pendulum_M->simMode);
    rtssSetSolverInfoPtr(sfcnInfo, &reaction_pendulum_M->solverInfoPtr);
  }

  reaction_pendulum_M->Sizes.numSFcns = (1);

  /* register each child */
  {
    (void) memset((void *)&reaction_pendulum_M->NonInlinedSFcns.childSFunctions
                  [0], 0,
                  1*sizeof(SimStruct));
    reaction_pendulum_M->childSfunctions =
      (&reaction_pendulum_M->NonInlinedSFcns.childSFunctionPtrs[0]);
    reaction_pendulum_M->childSfunctions[0] =
      (&reaction_pendulum_M->NonInlinedSFcns.childSFunctions[0]);

    /* Level2 S-Function Block: reaction_pendulum/<S2>/S-Function (rtdacusb2_rpend_dd) */
    {
      SimStruct *rts = reaction_pendulum_M->childSfunctions[0];

      /* timing info */
      time_T *sfcnPeriod = reaction_pendulum_M->NonInlinedSFcns.Sfcn0.sfcnPeriod;
      time_T *sfcnOffset = reaction_pendulum_M->NonInlinedSFcns.Sfcn0.sfcnOffset;
      int_T *sfcnTsMap = reaction_pendulum_M->NonInlinedSFcns.Sfcn0.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      {
        ssSetBlkInfo2Ptr(rts, &reaction_pendulum_M->NonInlinedSFcns.blkInfo2[0]);
      }

      _ssSetBlkInfo2PortInfo2Ptr(rts,
        &reaction_pendulum_M->NonInlinedSFcns.inputOutputPortInfo2[0]);

      /* Set up the mdlInfo pointer */
      ssSetRTWSfcnInfo(rts, reaction_pendulum_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts, &reaction_pendulum_M->NonInlinedSFcns.methods2[0]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts, &reaction_pendulum_M->NonInlinedSFcns.methods3[0]);
      }

      /* Allocate memory of model methods 4 */
      {
        ssSetModelMethods4(rts, &reaction_pendulum_M->NonInlinedSFcns.methods4[0]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts, &reaction_pendulum_M->NonInlinedSFcns.statesInfo2
                         [0]);
        ssSetPeriodicStatesInfo(rts,
          &reaction_pendulum_M->NonInlinedSFcns.periodicStatesInfo[0]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 4);
        ssSetPortInfoForInputs(rts,
          &reaction_pendulum_M->NonInlinedSFcns.Sfcn0.inputPortInfo[0]);
        _ssSetPortInfo2ForInputUnits(rts,
          &reaction_pendulum_M->NonInlinedSFcns.Sfcn0.inputPortUnits[0]);
        ssSetInputPortUnit(rts, 0, 0);
        ssSetInputPortUnit(rts, 1, 0);
        ssSetInputPortUnit(rts, 2, 0);
        ssSetInputPortUnit(rts, 3, 0);
        _ssSetPortInfo2ForInputCoSimAttribute(rts,
          &reaction_pendulum_M->NonInlinedSFcns.Sfcn0.inputPortCoSimAttribute[0]);
        ssSetInputPortIsContinuousQuantity(rts, 0, 0);
        ssSetInputPortIsContinuousQuantity(rts, 1, 0);
        ssSetInputPortIsContinuousQuantity(rts, 2, 0);
        ssSetInputPortIsContinuousQuantity(rts, 3, 0);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &reaction_pendulum_M->NonInlinedSFcns.Sfcn0.UPtrs0;
          sfcnUPtrs[0] = &reaction_pendulum_B.Saturation;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }

        /* port 1 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &reaction_pendulum_M->NonInlinedSFcns.Sfcn0.UPtrs1;
          sfcnUPtrs[0] = reaction_pendulum_B.Gain;
          sfcnUPtrs[1] = &reaction_pendulum_B.Gain[1];
          ssSetInputPortSignalPtrs(rts, 1, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 2);
        }

        /* port 2 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &reaction_pendulum_M->NonInlinedSFcns.Sfcn0.UPtrs2;
          sfcnUPtrs[0] = &reaction_pendulum_B.Prescaler;
          ssSetInputPortSignalPtrs(rts, 2, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 2, 1);
          ssSetInputPortWidth(rts, 2, 1);
        }

        /* port 3 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &reaction_pendulum_M->NonInlinedSFcns.Sfcn0.UPtrs3;
          sfcnUPtrs[0] = &reaction_pendulum_B.ThermFlag;
          ssSetInputPortSignalPtrs(rts, 3, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 3, 1);
          ssSetInputPortWidth(rts, 3, 1);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &reaction_pendulum_M->NonInlinedSFcns.Sfcn0.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 7);
        _ssSetPortInfo2ForOutputUnits(rts,
          &reaction_pendulum_M->NonInlinedSFcns.Sfcn0.outputPortUnits[0]);
        ssSetOutputPortUnit(rts, 0, 0);
        ssSetOutputPortUnit(rts, 1, 0);
        ssSetOutputPortUnit(rts, 2, 0);
        ssSetOutputPortUnit(rts, 3, 0);
        ssSetOutputPortUnit(rts, 4, 0);
        ssSetOutputPortUnit(rts, 5, 0);
        ssSetOutputPortUnit(rts, 6, 0);
        _ssSetPortInfo2ForOutputCoSimAttribute(rts,
          &reaction_pendulum_M->NonInlinedSFcns.Sfcn0.outputPortCoSimAttribute[0]);
        ssSetOutputPortIsContinuousQuantity(rts, 0, 0);
        ssSetOutputPortIsContinuousQuantity(rts, 1, 0);
        ssSetOutputPortIsContinuousQuantity(rts, 2, 0);
        ssSetOutputPortIsContinuousQuantity(rts, 3, 0);
        ssSetOutputPortIsContinuousQuantity(rts, 4, 0);
        ssSetOutputPortIsContinuousQuantity(rts, 5, 0);
        ssSetOutputPortIsContinuousQuantity(rts, 6, 0);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            &reaction_pendulum_B.SFunction_o1));
        }

        /* port 1 */
        {
          _ssSetOutputPortNumDimensions(rts, 1, 1);
          ssSetOutputPortWidth(rts, 1, 1);
          ssSetOutputPortSignal(rts, 1, ((real_T *)
            &reaction_pendulum_B.SFunction_o2));
        }

        /* port 2 */
        {
          _ssSetOutputPortNumDimensions(rts, 2, 1);
          ssSetOutputPortWidth(rts, 2, 1);
          ssSetOutputPortSignal(rts, 2, ((real_T *)
            &reaction_pendulum_B.SFunction_o3));
        }

        /* port 3 */
        {
          _ssSetOutputPortNumDimensions(rts, 3, 1);
          ssSetOutputPortWidth(rts, 3, 1);
          ssSetOutputPortSignal(rts, 3, ((real_T *)
            &reaction_pendulum_B.SFunction_o4));
        }

        /* port 4 */
        {
          _ssSetOutputPortNumDimensions(rts, 4, 1);
          ssSetOutputPortWidth(rts, 4, 1);
          ssSetOutputPortSignal(rts, 4, ((real_T *)
            &reaction_pendulum_B.SFunction_o5));
        }

        /* port 5 */
        {
          _ssSetOutputPortNumDimensions(rts, 5, 1);
          ssSetOutputPortWidth(rts, 5, 1);
          ssSetOutputPortSignal(rts, 5, ((real_T *)
            &reaction_pendulum_B.SFunction_o6));
        }

        /* port 6 */
        {
          _ssSetOutputPortNumDimensions(rts, 6, 1);
          ssSetOutputPortWidth(rts, 6, 2);
          ssSetOutputPortSignal(rts, 6, ((real_T *)
            reaction_pendulum_B.SFunction_o7));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts, "reaction_pendulum/RP Driver/S-Function");
      ssSetRTModel(rts,reaction_pendulum_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &reaction_pendulum_M->NonInlinedSFcns.Sfcn0.params;
        ssSetSFcnParamsCount(rts, 2);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)reaction_pendulum_P.SFunction_P1_Size);
        ssSetSFcnParam(rts, 1, (mxArray*)reaction_pendulum_P.SFunction_P2_Size);
      }

      /* registration */
      rtdacusb2_rpend_dd(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 0.05);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetInputPortConnected(rts, 1, 1);
      _ssSetInputPortConnected(rts, 2, 1);
      _ssSetInputPortConnected(rts, 3, 1);
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortConnected(rts, 1, 1);
      _ssSetOutputPortConnected(rts, 2, 1);
      _ssSetOutputPortConnected(rts, 3, 1);
      _ssSetOutputPortConnected(rts, 4, 1);
      _ssSetOutputPortConnected(rts, 5, 1);
      _ssSetOutputPortConnected(rts, 6, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);
      _ssSetOutputPortBeingMerged(rts, 1, 0);
      _ssSetOutputPortBeingMerged(rts, 2, 0);
      _ssSetOutputPortBeingMerged(rts, 3, 0);
      _ssSetOutputPortBeingMerged(rts, 4, 0);
      _ssSetOutputPortBeingMerged(rts, 5, 0);
      _ssSetOutputPortBeingMerged(rts, 6, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
      ssSetInputPortBufferDstPort(rts, 1, -1);
      ssSetInputPortBufferDstPort(rts, 2, -1);
      ssSetInputPortBufferDstPort(rts, 3, -1);
    }
  }

  /* Matfile logging */
  rt_StartDataLoggingWithStartTime(reaction_pendulum_M->rtwLogInfo, 0.0,
    rtmGetTFinal(reaction_pendulum_M), reaction_pendulum_M->Timing.stepSize0,
    (&rtmGetErrorStatus(reaction_pendulum_M)));

  /* Start for Constant: '<S2>/Prescaler' */
  reaction_pendulum_B.Prescaler = reaction_pendulum_P.Prescaler_Value;

  /* Start for Constant: '<S2>/ThermFlag' */
  reaction_pendulum_B.ThermFlag = reaction_pendulum_P.ThermFlag_Value;

  /* Start for S-Function (rtdacusb2_rpend_dd): '<S2>/S-Function' */
  /* Level2 S-Function Block: '<S2>/S-Function' (rtdacusb2_rpend_dd) */
  {
    SimStruct *rts = reaction_pendulum_M->childSfunctions[0];
    sfcnStart(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* Start for Enabled SubSystem: '<S30>/Enabled Subsystem' */
  reaction_pendulum_DW.EnabledSubsystem_MODE = false;

  /* End of Start for SubSystem: '<S30>/Enabled Subsystem' */

  /* Start for Enabled SubSystem: '<S23>/MeasurementUpdate' */
  reaction_pendulum_DW.MeasurementUpdate_MODE = false;

  /* End of Start for SubSystem: '<S23>/MeasurementUpdate' */

  /* InitializeConditions for Memory: '<S2>/Memory1' */
  reaction_pendulum_DW.Memory1_PreviousInput =
    reaction_pendulum_P.Memory1_InitialCondition;

  /* InitializeConditions for Memory: '<S2>/Memory' */
  reaction_pendulum_DW.Memory_PreviousInput =
    reaction_pendulum_P.Memory_InitialCondition;

  /* InitializeConditions for Delay: '<S1>/MemoryX' */
  reaction_pendulum_DW.icLoad = true;

  /* InitializeConditions for Derivative: '<Root>/Derivative' */
  reaction_pendulum_DW.TimeStampA = (rtInf);
  reaction_pendulum_DW.TimeStampB = (rtInf);

  /* InitializeConditions for Delay: '<Root>/Delay One Step' */
  reaction_pendulum_DW.DelayOneStep_DSTATE[0] =
    reaction_pendulum_P.DelayOneStep_InitialCondition;

  /* SystemInitialize for Enabled SubSystem: '<S30>/Enabled Subsystem' */
  /* SystemInitialize for Product: '<S56>/Product2' incorporates:
   *  Outport: '<S56>/deltax'
   */
  reaction_pendulum_B.Product2[0] = reaction_pendulum_P.deltax_Y0;

  /* End of SystemInitialize for SubSystem: '<S30>/Enabled Subsystem' */

  /* SystemInitialize for Enabled SubSystem: '<S23>/MeasurementUpdate' */
  /* SystemInitialize for Product: '<S54>/Product3' incorporates:
   *  Outport: '<S54>/L*(y[k]-yhat[k|k-1])'
   */
  reaction_pendulum_B.Product3[0] = reaction_pendulum_P.Lykyhatkk1_Y0;

  /* End of SystemInitialize for SubSystem: '<S23>/MeasurementUpdate' */

  /* InitializeConditions for Delay: '<Root>/Delay One Step' */
  reaction_pendulum_DW.DelayOneStep_DSTATE[1] =
    reaction_pendulum_P.DelayOneStep_InitialCondition;

  /* SystemInitialize for Enabled SubSystem: '<S30>/Enabled Subsystem' */
  /* SystemInitialize for Product: '<S56>/Product2' incorporates:
   *  Outport: '<S56>/deltax'
   */
  reaction_pendulum_B.Product2[1] = reaction_pendulum_P.deltax_Y0;

  /* End of SystemInitialize for SubSystem: '<S30>/Enabled Subsystem' */

  /* SystemInitialize for Enabled SubSystem: '<S23>/MeasurementUpdate' */
  /* SystemInitialize for Product: '<S54>/Product3' incorporates:
   *  Outport: '<S54>/L*(y[k]-yhat[k|k-1])'
   */
  reaction_pendulum_B.Product3[1] = reaction_pendulum_P.Lykyhatkk1_Y0;

  /* End of SystemInitialize for SubSystem: '<S23>/MeasurementUpdate' */

  /* InitializeConditions for Delay: '<Root>/Delay One Step' */
  reaction_pendulum_DW.DelayOneStep_DSTATE[2] =
    reaction_pendulum_P.DelayOneStep_InitialCondition;

  /* SystemInitialize for Enabled SubSystem: '<S30>/Enabled Subsystem' */
  /* SystemInitialize for Product: '<S56>/Product2' incorporates:
   *  Outport: '<S56>/deltax'
   */
  reaction_pendulum_B.Product2[2] = reaction_pendulum_P.deltax_Y0;

  /* End of SystemInitialize for SubSystem: '<S30>/Enabled Subsystem' */

  /* SystemInitialize for Enabled SubSystem: '<S23>/MeasurementUpdate' */
  /* SystemInitialize for Product: '<S54>/Product3' incorporates:
   *  Outport: '<S54>/L*(y[k]-yhat[k|k-1])'
   */
  reaction_pendulum_B.Product3[2] = reaction_pendulum_P.Lykyhatkk1_Y0;

  /* End of SystemInitialize for SubSystem: '<S23>/MeasurementUpdate' */
}

/* Model terminate function */
void reaction_pendulum_terminate(void)
{
  /* Terminate for S-Function (rtdacusb2_rpend_dd): '<S2>/S-Function' */
  /* Level2 S-Function Block: '<S2>/S-Function' (rtdacusb2_rpend_dd) */
  {
    SimStruct *rts = reaction_pendulum_M->childSfunctions[0];
    sfcnTerminate(rts);
  }
}

#include <stdio.h>

/* Final time from "Simulation Parameters"   */
real_T finaltime = 100.0;
real_T StepSize = 0.05;

////////////////////////////////////////////////
//
//  Return compilation date and time
//
char *GetDateAndTime( void )
{
  static char AuxStr[ 128 ];
  sprintf( AuxStr, "%s %s", __DATE__, __TIME__ );
  return( AuxStr );
}
