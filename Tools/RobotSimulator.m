classdef RobotSimulator < handle
    %ROBOTMCT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        headOrientation = 0
        audio                     % Source audio [nSamples x nSources]
        fsHz                      % Sampling rate
        azRef                     % Source reference positions
        room                      % Label for BRIRs
        nTotalSamples             % Number of audio samples
        isFinished = false;       % Return true if end of audio reached
        sampleIndex = 0;          % Current processing index
        bUseSurrey;               % Whether to use the Surrey azimuth convention
    end
    
    methods
        function obj = RobotSimulator(audio, fsHz, azRef, room)
            obj.audio = audio;
            obj.fsHz = fsHz;
            obj.azRef = azRef;
            obj.room = room;
            obj.nTotalSamples = size(audio,1);
            if isempty(strfind(room,'SURREY'))
                obj.bUseSurrey = false;
            else
                obj.bUseSurrey = true;
            end  
        end
        
        function rotateHead(obj, angleDeg)
            if angleDeg ~= 0
                obj.headOrientation = mod(obj.headOrientation+angleDeg, 360);
            end
        end
        
        % Returns the duration of signal in seconds if there is less than
        % durSec seconds left
        function durSec = canGetSignal(obj, durSec)
            if obj.isFinished
                durSec = 0;
            else
                durSamples = floor(durSec * obj.fsHz);
                leftSamples = obj.nTotalSamples - obj.fsHz * 10E-3 - obj.sampleIndex;
                durSec = min(durSamples, leftSamples) / obj.fsHz;
            end
        end

        function durSec = remainingDuration(obj)
            if obj.isFinished
                durSec = 0;
            else
                durSec = (obj.nTotalSamples - obj.sampleIndex) / obj.fsHz;
            end
        end
        
        function [sig,fsHz] = getSignal(obj, durSec)
            fsHz = obj.fsHz;
            if obj.isFinished
                sig = [];
            else
                durSamples = floor(durSec * obj.fsHz);
                blockSamples = (1:durSamples) + obj.sampleIndex;
                
                % Make sure the last block does not exceed nTotalSamples
                blockSamples = blockSamples(blockSamples<=obj.nTotalSamples);
                
                obj.sampleIndex = obj.sampleIndex + durSamples;
                
                % If there is less than 10 ms we stop processing
                if obj.sampleIndex >= obj.nTotalSamples - obj.fsHz * 10E-3
                    obj.isFinished = true;
                end
                
                relAzRef = mod(obj.azRef - obj.headOrientation, 360);
                if obj.bUseSurrey
                    % Transform azimuth to SURREY representation
                    relAzRef = convertAzimuthsWP1ToSurrey(relAzRef);
                end
                sig = spatializeAudio(obj.audio(blockSamples,:), obj.fsHz, relAzRef, obj.room);
            end
        end
        
    end
    
end


