    function targMap = targDataMap(),

    ;%***********************
    ;% Create Parameter Map *
    ;%***********************
    
        nTotData      = 0; %add to this count as we go
        nTotSects     = 3;
        sectIdxOffset = 0;

        ;%
        ;% Define dummy sections & preallocate arrays
        ;%
        dumSection.nData = -1;
        dumSection.data  = [];

        dumData.logicalSrcIdx = -1;
        dumData.dtTransOffset = -1;

        ;%
        ;% Init/prealloc paramMap
        ;%
        paramMap.nSections           = nTotSects;
        paramMap.sectIdxOffset       = sectIdxOffset;
            paramMap.sections(nTotSects) = dumSection; %prealloc
        paramMap.nTotData            = -1;

        ;%
        ;% Auto data (reaction_pendulum_P)
        ;%
            section.nData     = 36;
            section.data(36)  = dumData; %prealloc

                    ;% reaction_pendulum_P.K
                    section.data(1).logicalSrcIdx = 0;
                    section.data(1).dtTransOffset = 0;

                    ;% reaction_pendulum_P.SFunction_P2_Size
                    section.data(2).logicalSrcIdx = 1;
                    section.data(2).dtTransOffset = 3;

                    ;% reaction_pendulum_P.h
                    section.data(3).logicalSrcIdx = 2;
                    section.data(3).dtTransOffset = 5;

                    ;% reaction_pendulum_P.Lykyhatkk1_Y0
                    section.data(4).logicalSrcIdx = 3;
                    section.data(4).dtTransOffset = 6;

                    ;% reaction_pendulum_P.deltax_Y0
                    section.data(5).logicalSrcIdx = 4;
                    section.data(5).dtTransOffset = 7;

                    ;% reaction_pendulum_P.Gain1_Gain
                    section.data(6).logicalSrcIdx = 5;
                    section.data(6).dtTransOffset = 8;

                    ;% reaction_pendulum_P.Reset_Value
                    section.data(7).logicalSrcIdx = 6;
                    section.data(7).dtTransOffset = 9;

                    ;% reaction_pendulum_P.Normal_Value
                    section.data(8).logicalSrcIdx = 7;
                    section.data(8).dtTransOffset = 10;

                    ;% reaction_pendulum_P.DelayOneStep_InitialCondition
                    section.data(9).logicalSrcIdx = 8;
                    section.data(9).dtTransOffset = 11;

                    ;% reaction_pendulum_P.KalmanGainM_Value
                    section.data(10).logicalSrcIdx = 9;
                    section.data(10).dtTransOffset = 12;

                    ;% reaction_pendulum_P.C_Value
                    section.data(11).logicalSrcIdx = 10;
                    section.data(11).dtTransOffset = 18;

                    ;% reaction_pendulum_P.DC_Ctrl1_Value
                    section.data(12).logicalSrcIdx = 11;
                    section.data(12).dtTransOffset = 24;

                    ;% reaction_pendulum_P.SignalGenerator_Amplitude
                    section.data(13).logicalSrcIdx = 12;
                    section.data(13).dtTransOffset = 25;

                    ;% reaction_pendulum_P.SignalGenerator_Frequency
                    section.data(14).logicalSrcIdx = 13;
                    section.data(14).dtTransOffset = 26;

                    ;% reaction_pendulum_P.Gain2_Gain
                    section.data(15).logicalSrcIdx = 14;
                    section.data(15).dtTransOffset = 27;

                    ;% reaction_pendulum_P.Saturation_UpperSat
                    section.data(16).logicalSrcIdx = 15;
                    section.data(16).dtTransOffset = 28;

                    ;% reaction_pendulum_P.Saturation_LowerSat
                    section.data(17).logicalSrcIdx = 16;
                    section.data(17).dtTransOffset = 29;

                    ;% reaction_pendulum_P.Gain_Gain
                    section.data(18).logicalSrcIdx = 17;
                    section.data(18).dtTransOffset = 30;

                    ;% reaction_pendulum_P.Prescaler_Value
                    section.data(19).logicalSrcIdx = 18;
                    section.data(19).dtTransOffset = 32;

                    ;% reaction_pendulum_P.ThermFlag_Value
                    section.data(20).logicalSrcIdx = 19;
                    section.data(20).dtTransOffset = 33;

                    ;% reaction_pendulum_P.SFunction_P1_Size
                    section.data(21).logicalSrcIdx = 20;
                    section.data(21).dtTransOffset = 34;

                    ;% reaction_pendulum_P.SFunction_P1
                    section.data(22).logicalSrcIdx = 21;
                    section.data(22).dtTransOffset = 36;

                    ;% reaction_pendulum_P.PendulumConverttorad_Gain
                    section.data(23).logicalSrcIdx = 22;
                    section.data(23).dtTransOffset = 37;

                    ;% reaction_pendulum_P.DCConverttorad_Gain
                    section.data(24).logicalSrcIdx = 23;
                    section.data(24).dtTransOffset = 38;

                    ;% reaction_pendulum_P.Memory1_InitialCondition
                    section.data(25).logicalSrcIdx = 24;
                    section.data(25).dtTransOffset = 39;

                    ;% reaction_pendulum_P.rad2RPM_Gain
                    section.data(26).logicalSrcIdx = 25;
                    section.data(26).dtTransOffset = 40;

                    ;% reaction_pendulum_P.Memory_InitialCondition
                    section.data(27).logicalSrcIdx = 26;
                    section.data(27).dtTransOffset = 41;

                    ;% reaction_pendulum_P.Gain1_Gain_h
                    section.data(28).logicalSrcIdx = 27;
                    section.data(28).dtTransOffset = 42;

                    ;% reaction_pendulum_P.X0_Value
                    section.data(29).logicalSrcIdx = 28;
                    section.data(29).dtTransOffset = 43;

                    ;% reaction_pendulum_P.DCConverttoA1_Gain
                    section.data(30).logicalSrcIdx = 29;
                    section.data(30).dtTransOffset = 46;

                    ;% reaction_pendulum_P.Constant_Value
                    section.data(31).logicalSrcIdx = 30;
                    section.data(31).dtTransOffset = 47;

                    ;% reaction_pendulum_P.A_Value
                    section.data(32).logicalSrcIdx = 31;
                    section.data(32).dtTransOffset = 48;

                    ;% reaction_pendulum_P.B_Value
                    section.data(33).logicalSrcIdx = 32;
                    section.data(33).dtTransOffset = 57;

                    ;% reaction_pendulum_P.CovarianceZ_Value
                    section.data(34).logicalSrcIdx = 33;
                    section.data(34).dtTransOffset = 60;

                    ;% reaction_pendulum_P.KalmanGainL_Value
                    section.data(35).logicalSrcIdx = 34;
                    section.data(35).dtTransOffset = 69;

                    ;% reaction_pendulum_P.D_Value
                    section.data(36).logicalSrcIdx = 35;
                    section.data(36).dtTransOffset = 75;

            nTotData = nTotData + section.nData;
            paramMap.sections(1) = section;
            clear section

            section.nData     = 2;
            section.data(2)  = dumData; %prealloc

                    ;% reaction_pendulum_P.Enable_Value
                    section.data(1).logicalSrcIdx = 36;
                    section.data(1).dtTransOffset = 0;

                    ;% reaction_pendulum_P.isSqrtUsed_Value
                    section.data(2).logicalSrcIdx = 37;
                    section.data(2).dtTransOffset = 1;

            nTotData = nTotData + section.nData;
            paramMap.sections(2) = section;
            clear section

            section.nData     = 2;
            section.data(2)  = dumData; %prealloc

                    ;% reaction_pendulum_P.ResetEncoders2_CurrentSetting
                    section.data(1).logicalSrcIdx = 38;
                    section.data(1).dtTransOffset = 0;

                    ;% reaction_pendulum_P.ResetEncoders_CurrentSetting
                    section.data(2).logicalSrcIdx = 39;
                    section.data(2).dtTransOffset = 1;

            nTotData = nTotData + section.nData;
            paramMap.sections(3) = section;
            clear section


            ;%
            ;% Non-auto Data (parameter)
            ;%


        ;%
        ;% Add final counts to struct.
        ;%
        paramMap.nTotData = nTotData;



    ;%**************************
    ;% Create Block Output Map *
    ;%**************************
    
        nTotData      = 0; %add to this count as we go
        nTotSects     = 1;
        sectIdxOffset = 0;

        ;%
        ;% Define dummy sections & preallocate arrays
        ;%
        dumSection.nData = -1;
        dumSection.data  = [];

        dumData.logicalSrcIdx = -1;
        dumData.dtTransOffset = -1;

        ;%
        ;% Init/prealloc sigMap
        ;%
        sigMap.nSections           = nTotSects;
        sigMap.sectIdxOffset       = sectIdxOffset;
            sigMap.sections(nTotSects) = dumSection; %prealloc
        sigMap.nTotData            = -1;

        ;%
        ;% Auto data (reaction_pendulum_B)
        ;%
            section.nData     = 23;
            section.data(23)  = dumData; %prealloc

                    ;% reaction_pendulum_B.DelayOneStep
                    section.data(1).logicalSrcIdx = 0;
                    section.data(1).dtTransOffset = 0;

                    ;% reaction_pendulum_B.Sum
                    section.data(2).logicalSrcIdx = 1;
                    section.data(2).dtTransOffset = 3;

                    ;% reaction_pendulum_B.Saturation
                    section.data(3).logicalSrcIdx = 2;
                    section.data(3).dtTransOffset = 4;

                    ;% reaction_pendulum_B.Gain
                    section.data(4).logicalSrcIdx = 3;
                    section.data(4).dtTransOffset = 5;

                    ;% reaction_pendulum_B.Prescaler
                    section.data(5).logicalSrcIdx = 4;
                    section.data(5).dtTransOffset = 7;

                    ;% reaction_pendulum_B.ThermFlag
                    section.data(6).logicalSrcIdx = 5;
                    section.data(6).dtTransOffset = 8;

                    ;% reaction_pendulum_B.SFunction_o1
                    section.data(7).logicalSrcIdx = 6;
                    section.data(7).dtTransOffset = 9;

                    ;% reaction_pendulum_B.SFunction_o2
                    section.data(8).logicalSrcIdx = 7;
                    section.data(8).dtTransOffset = 10;

                    ;% reaction_pendulum_B.SFunction_o3
                    section.data(9).logicalSrcIdx = 8;
                    section.data(9).dtTransOffset = 11;

                    ;% reaction_pendulum_B.SFunction_o4
                    section.data(10).logicalSrcIdx = 9;
                    section.data(10).dtTransOffset = 12;

                    ;% reaction_pendulum_B.SFunction_o5
                    section.data(11).logicalSrcIdx = 10;
                    section.data(11).dtTransOffset = 13;

                    ;% reaction_pendulum_B.SFunction_o6
                    section.data(12).logicalSrcIdx = 11;
                    section.data(12).dtTransOffset = 14;

                    ;% reaction_pendulum_B.SFunction_o7
                    section.data(13).logicalSrcIdx = 12;
                    section.data(13).dtTransOffset = 15;

                    ;% reaction_pendulum_B.PendulumAnglerad
                    section.data(14).logicalSrcIdx = 13;
                    section.data(14).dtTransOffset = 17;

                    ;% reaction_pendulum_B.Periodms
                    section.data(15).logicalSrcIdx = 14;
                    section.data(15).dtTransOffset = 18;

                    ;% reaction_pendulum_B.DCVelrads
                    section.data(16).logicalSrcIdx = 15;
                    section.data(16).dtTransOffset = 19;

                    ;% reaction_pendulum_B.Reshapexhat
                    section.data(17).logicalSrcIdx = 16;
                    section.data(17).dtTransOffset = 20;

                    ;% reaction_pendulum_B.DCConverttoA1
                    section.data(18).logicalSrcIdx = 17;
                    section.data(18).dtTransOffset = 23;

                    ;% reaction_pendulum_B.Gain_i
                    section.data(19).logicalSrcIdx = 18;
                    section.data(19).dtTransOffset = 24;

                    ;% reaction_pendulum_B.Sum1
                    section.data(20).logicalSrcIdx = 19;
                    section.data(20).dtTransOffset = 25;

                    ;% reaction_pendulum_B.thetadotdiff
                    section.data(21).logicalSrcIdx = 20;
                    section.data(21).dtTransOffset = 26;

                    ;% reaction_pendulum_B.Product2
                    section.data(22).logicalSrcIdx = 21;
                    section.data(22).dtTransOffset = 27;

                    ;% reaction_pendulum_B.Product3
                    section.data(23).logicalSrcIdx = 22;
                    section.data(23).dtTransOffset = 30;

            nTotData = nTotData + section.nData;
            sigMap.sections(1) = section;
            clear section


            ;%
            ;% Non-auto Data (signal)
            ;%


        ;%
        ;% Add final counts to struct.
        ;%
        sigMap.nTotData = nTotData;



    ;%*******************
    ;% Create DWork Map *
    ;%*******************
    
        nTotData      = 0; %add to this count as we go
        nTotSects     = 4;
        sectIdxOffset = 1;

        ;%
        ;% Define dummy sections & preallocate arrays
        ;%
        dumSection.nData = -1;
        dumSection.data  = [];

        dumData.logicalSrcIdx = -1;
        dumData.dtTransOffset = -1;

        ;%
        ;% Init/prealloc dworkMap
        ;%
        dworkMap.nSections           = nTotSects;
        dworkMap.sectIdxOffset       = sectIdxOffset;
            dworkMap.sections(nTotSects) = dumSection; %prealloc
        dworkMap.nTotData            = -1;

        ;%
        ;% Auto data (reaction_pendulum_DW)
        ;%
            section.nData     = 8;
            section.data(8)  = dumData; %prealloc

                    ;% reaction_pendulum_DW.DelayOneStep_DSTATE
                    section.data(1).logicalSrcIdx = 0;
                    section.data(1).dtTransOffset = 0;

                    ;% reaction_pendulum_DW.MemoryX_DSTATE
                    section.data(2).logicalSrcIdx = 1;
                    section.data(2).dtTransOffset = 3;

                    ;% reaction_pendulum_DW.Memory1_PreviousInput
                    section.data(3).logicalSrcIdx = 2;
                    section.data(3).dtTransOffset = 6;

                    ;% reaction_pendulum_DW.Memory_PreviousInput
                    section.data(4).logicalSrcIdx = 3;
                    section.data(4).dtTransOffset = 7;

                    ;% reaction_pendulum_DW.TimeStampA
                    section.data(5).logicalSrcIdx = 4;
                    section.data(5).dtTransOffset = 8;

                    ;% reaction_pendulum_DW.LastUAtTimeA
                    section.data(6).logicalSrcIdx = 5;
                    section.data(6).dtTransOffset = 9;

                    ;% reaction_pendulum_DW.TimeStampB
                    section.data(7).logicalSrcIdx = 6;
                    section.data(7).dtTransOffset = 10;

                    ;% reaction_pendulum_DW.LastUAtTimeB
                    section.data(8).logicalSrcIdx = 7;
                    section.data(8).dtTransOffset = 11;

            nTotData = nTotData + section.nData;
            dworkMap.sections(1) = section;
            clear section

            section.nData     = 12;
            section.data(12)  = dumData; %prealloc

                    ;% reaction_pendulum_DW.Phidot_PWORK.LoggedData
                    section.data(1).logicalSrcIdx = 8;
                    section.data(1).dtTransOffset = 0;

                    ;% reaction_pendulum_DW.PlotState_PWORK.LoggedData
                    section.data(2).logicalSrcIdx = 9;
                    section.data(2).dtTransOffset = 2;

                    ;% reaction_pendulum_DW.Scope1_PWORK.LoggedData
                    section.data(3).logicalSrcIdx = 10;
                    section.data(3).dtTransOffset = 7;

                    ;% reaction_pendulum_DW.Scope2_PWORK.LoggedData
                    section.data(4).logicalSrcIdx = 11;
                    section.data(4).dtTransOffset = 8;

                    ;% reaction_pendulum_DW.Scope3_PWORK.LoggedData
                    section.data(5).logicalSrcIdx = 12;
                    section.data(5).dtTransOffset = 9;

                    ;% reaction_pendulum_DW.Scope4_PWORK.LoggedData
                    section.data(6).logicalSrcIdx = 13;
                    section.data(6).dtTransOffset = 10;

                    ;% reaction_pendulum_DW.Scope5_PWORK.LoggedData
                    section.data(7).logicalSrcIdx = 14;
                    section.data(7).dtTransOffset = 11;

                    ;% reaction_pendulum_DW.Scope6_PWORK.LoggedData
                    section.data(8).logicalSrcIdx = 15;
                    section.data(8).dtTransOffset = 12;

                    ;% reaction_pendulum_DW.Scope7_PWORK.LoggedData
                    section.data(9).logicalSrcIdx = 16;
                    section.data(9).dtTransOffset = 13;

                    ;% reaction_pendulum_DW.theta_PWORK.LoggedData
                    section.data(10).logicalSrcIdx = 17;
                    section.data(10).dtTransOffset = 14;

                    ;% reaction_pendulum_DW.thetadot_PWORK.LoggedData
                    section.data(11).logicalSrcIdx = 18;
                    section.data(11).dtTransOffset = 15;

                    ;% reaction_pendulum_DW.theta1_PWORK.LoggedData
                    section.data(12).logicalSrcIdx = 19;
                    section.data(12).dtTransOffset = 16;

            nTotData = nTotData + section.nData;
            dworkMap.sections(2) = section;
            clear section

            section.nData     = 2;
            section.data(2)  = dumData; %prealloc

                    ;% reaction_pendulum_DW.EnabledSubsystem_SubsysRanBC
                    section.data(1).logicalSrcIdx = 20;
                    section.data(1).dtTransOffset = 0;

                    ;% reaction_pendulum_DW.MeasurementUpdate_SubsysRanBC
                    section.data(2).logicalSrcIdx = 21;
                    section.data(2).dtTransOffset = 1;

            nTotData = nTotData + section.nData;
            dworkMap.sections(3) = section;
            clear section

            section.nData     = 3;
            section.data(3)  = dumData; %prealloc

                    ;% reaction_pendulum_DW.icLoad
                    section.data(1).logicalSrcIdx = 22;
                    section.data(1).dtTransOffset = 0;

                    ;% reaction_pendulum_DW.EnabledSubsystem_MODE
                    section.data(2).logicalSrcIdx = 23;
                    section.data(2).dtTransOffset = 1;

                    ;% reaction_pendulum_DW.MeasurementUpdate_MODE
                    section.data(3).logicalSrcIdx = 24;
                    section.data(3).dtTransOffset = 2;

            nTotData = nTotData + section.nData;
            dworkMap.sections(4) = section;
            clear section


            ;%
            ;% Non-auto Data (dwork)
            ;%


        ;%
        ;% Add final counts to struct.
        ;%
        dworkMap.nTotData = nTotData;



    ;%
    ;% Add individual maps to base struct.
    ;%

    targMap.paramMap  = paramMap;
    targMap.signalMap = sigMap;
    targMap.dworkMap  = dworkMap;

    ;%
    ;% Add checksums to base struct.
    ;%


    targMap.checksum0 = 4043391257;
    targMap.checksum1 = 3157747759;
    targMap.checksum2 = 4109200965;
    targMap.checksum3 = 4046214105;

