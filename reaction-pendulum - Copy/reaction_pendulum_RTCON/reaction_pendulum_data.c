/*
 * reaction_pendulum_data.c
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

/* Block parameters (default storage) */
P_reaction_pendulum_T reaction_pendulum_P = {
  /* Variable: K
   * Referenced by: '<Root>/Gain'
   */
  { 0.028103205953057014, -0.0016146829609239347, 0.055788562358802807 },

  /* Computed Parameter: SFunction_P2_Size
   * Referenced by: '<S2>/S-Function'
   */
  { 1.0, 1.0 },

  /* Variable: h
   * Referenced by: '<S2>/S-Function'
   */
  0.05,

  /* Expression: 0
   * Referenced by: '<S54>/L*(y[k]-yhat[k|k-1])'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S56>/deltax'
   */
  0.0,

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

  /* Expression: 0.0
   * Referenced by: '<Root>/Delay One Step'
   */
  0.0,

  /* Expression: pInitialization.M
   * Referenced by: '<S3>/KalmanGainM'
   */
  { 0.95448452084725366, 0.01708018358799358, 0.0011132279607907576,
    0.0011132279607907504, -0.00066593889619775363, 0.9548737864654484 },

  /* Expression: pInitialization.C
   * Referenced by: '<S1>/C'
   */
  { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0 },

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
   * Referenced by: '<S2>/Gain2'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S2>/Saturation'
   */
  1.0,

  /* Expression: -1
   * Referenced by: '<S2>/Saturation'
   */
  -1.0,

  /* Expression: [1;1]
   * Referenced by: '<S2>/Gain'
   */
  { 1.0, 1.0 },

  /* Expression: 0
   * Referenced by: '<S2>/Prescaler'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S2>/ThermFlag'
   */
  0.0,

  /* Computed Parameter: SFunction_P1_Size
   * Referenced by: '<S2>/S-Function'
   */
  { 1.0, 1.0 },

  /* Expression: 1
   * Referenced by: '<S2>/S-Function'
   */
  1.0,

  /* Expression: 2*pi/20000
   * Referenced by: '<S2>/Pendulum Convert to rad'
   */
  0.00031415926535897931,

  /* Expression: 2*pi/4096
   * Referenced by: '<S2>/DC Convert to rad'
   */
  0.0015339807878856412,

  /* Expression: 0
   * Referenced by: '<S2>/Memory1'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S2>/rad2RPM'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S2>/Memory'
   */
  0.0,

  /* Expression: 1/20000000
   * Referenced by: '<S2>/Gain1'
   */
  5.0E-8,

  /* Expression: pInitialization.X0
   * Referenced by: '<S1>/X0'
   */
  { 0.0, 0.0, 0.0 },

  /* Expression: 0.8333
   * Referenced by: '<S2>/DC Convert to [A]1'
   */
  0.8333,

  /* Expression: 0
   * Referenced by: '<Root>/Constant'
   */
  0.0,

  /* Expression: pInitialization.A
   * Referenced by: '<S1>/A'
   */
  { 0.98938551434434385, -0.39121037208628806, 0.5588570044781821,
    0.039199074802239262, 0.5980450437143634, 0.010561489196473029,
    1.1302891669337189E-5, 0.00041288918162421236, 0.95478117398525364 },

  /* Expression: pInitialization.B
   * Referenced by: '<S1>/B'
   */
  { -0.0031881710089001576, -0.11610905714770475, 17.023365225908535 },

  /* Expression: pInitialization.Z
   * Referenced by: '<S3>/CovarianceZ'
   */
  { 0.047724226042362684, 0.0008540091793996789, 5.5661398039537645E-5,
    0.00085400917939967911, 1.5670524462048108, -3.3296944809887683E-5,
    5.5661398039537645E-5, -3.3296944809887683E-5, 0.047743689323272435 },

  /* Expression: pInitialization.L
   * Referenced by: '<S3>/KalmanGainL'
   */
  { 0.94502269856897136, -0.3631890657709444, 0.53466364141524447,
    0.0010861002649298461, -0.00043951072479786533, 0.91231061678623637 },

  /* Expression: pInitialization.D
   * Referenced by: '<S1>/D'
   */
  { 0.0, 0.0 },

  /* Expression: true()
   * Referenced by: '<S1>/Enable'
   */
  true,

  /* Expression: pInitialization.isSqrtUsed
   * Referenced by: '<S52>/isSqrtUsed'
   */
  false,

  /* Computed Parameter: ResetEncoders2_CurrentSetting
   * Referenced by: '<Root>/Reset Encoders2'
   */
  1U,

  /* Computed Parameter: ResetEncoders_CurrentSetting
   * Referenced by: '<Root>/Reset Encoders'
   */
  0U
};
