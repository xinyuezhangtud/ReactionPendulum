/*
 * reaction_pendulum_data.c
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "reaction_pendulum".
 *
 * Model version              : 5.1
 * Simulink Coder version : 9.6 (R2021b) 14-May-2021
 * C source code generated on : Tue Apr 30 15:09:29 2024
 *
 * Target selection: rtcon_rpend_usb2.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "reaction_pendulum.h"
#include "reaction_pendulum_private.h"

/* Block parameters (default storage) */
P_reaction_pendulum_T reaction_pendulum_P = {
  /* Computed Parameter: SFunction_P2_Size
   * Referenced by: '<S1>/S-Function'
   */
  { 1.0, 1.0 },

  /* Variable: h
   * Referenced by: '<S1>/S-Function'
   */
  0.05,

  /* Expression: 0.3
   * Referenced by: '<Root>/Gain1'
   */
  0.3,

  /* Expression: 1
   * Referenced by: '<Root>/Reset'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<Root>/Normal'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<Root>/DC_Ctrl1'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<Root>/Signal Generator'
   */
  1.0,

  /* Expression: 1/2
   * Referenced by: '<Root>/Signal Generator'
   */
  0.5,

  /* Expression: 1
   * Referenced by: '<S1>/Gain2'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S1>/Saturation'
   */
  1.0,

  /* Expression: -1
   * Referenced by: '<S1>/Saturation'
   */
  -1.0,

  /* Expression: [1;1]
   * Referenced by: '<S1>/Gain'
   */
  { 1.0, 1.0 },

  /* Expression: 0
   * Referenced by: '<S1>/Prescaler'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S1>/ThermFlag'
   */
  0.0,

  /* Computed Parameter: SFunction_P1_Size
   * Referenced by: '<S1>/S-Function'
   */
  { 1.0, 1.0 },

  /* Expression: 1
   * Referenced by: '<S1>/S-Function'
   */
  1.0,

  /* Expression: 2*pi/20000
   * Referenced by: '<S1>/Pendulum Convert to rad'
   */
  0.00031415926535897931,

  /* Expression: 2*pi/4096
   * Referenced by: '<S1>/DC Convert to rad'
   */
  0.0015339807878856412,

  /* Expression: 0
   * Referenced by: '<S1>/Memory1'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S1>/rad2RPM'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S1>/Memory'
   */
  0.0,

  /* Expression: 1/20000000
   * Referenced by: '<S1>/Gain1'
   */
  5.0E-8,

  /* Expression: 0.8333
   * Referenced by: '<S1>/DC Convert to [A]1'
   */
  0.8333,

  /* Computed Parameter: ResetEncoders2_CurrentSetting
   * Referenced by: '<Root>/Reset Encoders2'
   */
  0U,

  /* Computed Parameter: ResetEncoders_CurrentSetting
   * Referenced by: '<Root>/Reset Encoders'
   */
  0U
};
