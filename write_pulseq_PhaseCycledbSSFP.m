%%
% write_pulseq_PhaseCycledbSSFP.m
% Written by Kelly Bonzon, Yasaman Safarkhanlo
% Upload date: 19/08/2024 Last modified: 19/08/2024

%% Gradient Modes
% This section sets different gradient modes that determine the maximum
% gradient strength and slew rate, which are important for controlling the
% speed and smoothness of the gradient transitions.

grad_mode = 'Normal'; % default gradient mode

% Set maximum gradient strength and slew rate based on the chosen gradient mode
switch grad_mode
    case 'Fast'
        max_grad = 26;      % Max gradient strength [mT/m]
        max_slew = 180;     % Maximum slew rate [mT/m/ms]
    case 'Normal'
        max_grad = 22;      % Max gradient strength [mT/m]
        max_slew = 100;     % Maximum slew rate [mT/m/ms]
    case 'Whisper'
        max_grad = 22;      % Max gradient strength [mT/m]
        max_slew = 50;      % Maximum slew rate [mT/m/ms]
    case 'Performance'
        max_grad = 37;      % Max gradient strength [mT/m]
        max_slew = 188;     % Maximum slew rate [mT/m/ms]
end

%% System Limits
% Define the system hardware limits using the Pulseq options.
sys = mr.opts('MaxGrad', max_grad, 'GradUnit', 'mT/m', ...
              'MaxSlew', max_slew, 'SlewUnit', 'T/m/s', ...
              'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
              'adcDeadTime', 20e-6);

%% Start a Sequence
% Initialize a new Pulseq sequence object with the defined system limits.
seq = mr.Sequence(sys);

% Define the resolution and field of view (FOV).
res = [1 1 3];                % Resolution in x, y, z directions
fov = [384e-3 192e-3 192e-3]; % Field of view in meters

% Calculate the delta k-space values
deltak = 1 ./ fov;

% Matrix sizes in x, y, z directions.
Nx = 384;
Ny = 192;
Nz = 64;

% Calculate the area covered per kernel in Z and Y directions.
areaZ = ((0:Nz-1) - Nz/2) * deltak(3);
areaY = ((0:Ny-1) - Ny/2) * deltak(2);

% RF excitation parameters
alpha = 20;         % RF excitation flip angle in degrees
TR = 5.02;          % Repetition time in milliseconds
adc_dur = 384 * 3;  % ADC (Analog-to-Digital Converter) duration in microseconds

% Number of dummy pulses to reach steady state
Ndummy = 1000;

% Phase cycle increments in degrees
pc = 0:20:359;

% RF pulse parameters
rf_dur = 500;       % RF pulse duration in microseconds
rf_apo = 1;         % RF pulse apodization factor
rf_bwt = 3;         % RF pulse bandwidth-time product                             

% Create the RF pulse with the given parameters
[rf, rfDelay] = mr.makeBlockPulse(alpha*pi/180,sys,'Duration',rf_dur*1e-6);

% Create the gradient along the x-axis for the readout
gx = mr.makeTrapezoid('x','FlatArea',...
    Nx*deltak(1),...
    'FlatTime',adc_dur*1e-6,...
    'system',sys);        

% Create the ADC object to sample the signal during the readout gradient
adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime,'system',sys);

% Create the pre-phasing gradient (negative lobe) for the readout
RO_pre = mr.makeTrapezoid('x','Area',-gx.area/2,'system',sys);
RO_pre.delay=0;

% Splitting R.O. gradient to rise_and_flat + fall to be able to start
% rephasers of PE and SE gradients right after flat time of R.O. grad
gx_parts=mr.splitGradientAt(gx,ceil(mr.calcDuration(adc)/sys.gradRasterTime)*sys.gradRasterTime);
RO_rise_and_flat = gx_parts(1);
RO_rise_and_flat.delay = mr.calcDuration(RO_pre);
RO_fall = gx_parts(2);

% adding RO_pre, RO_rise_and_flat to make them consistently concat 
% (pulseq related)
RO_pre_rise_and_flat = mr.addGradients({RO_pre, RO_rise_and_flat});

% giving ADC additional delay to wait for the prephaser or R.O. gradient
adc.delay=adc.delay+mr.calcDuration(RO_pre); 

% Removing any delay of fall of RO 
% to be able start right after flat of RO
RO_fall.delay=0;

% adding RO_fall and RO_pre to make them consistently concat 
% (pulseq related)
RO_fall_and_pre = mr.addGradients({RO_fall, RO_pre});

%% Phase Encoding Gradients
% Phase encoding gradients for Y and Z directions

% Calculating the maximum required PE pre and rephasers to use their
% duration later (for consistent TR)
PE_pre = mr.makeTrapezoid('y','Area',-1*max(abs(areaY)),'system',sys);
PE_reph = mr.makeTrapezoid('y','Area',max(abs(areaY)),'system',sys);

% Calculating the maximum required SE pre and rephasers to use their
% duration later (for consistent TR)
SE_pre = mr.makeTrapezoid('z','Area',-1*max(abs(areaZ)),'system',sys);
SE_reph = mr.makeTrapezoid('z','Area',max(abs(areaZ)),'system',sys);

% Using duration of the maximum duration required by the 
% largest SE or PE prephaser gradients
max_T_pre = max([mr.calcDuration(SE_pre), mr.calcDuration(PE_pre)]);


% Creating one kernel and checking TR
seq.addBlock(rf)
seq.addBlock(PE_pre, SE_pre,RO_pre_rise_and_flat, adc)
seq.addBlock(RO_fall_and_pre, PE_reph,SE_reph)

TR_kernel = (seq.duration);
fprintf("Minimum TR obtainable: %2.2f ms\n", 1e3*TR_kernel);

% Create a delay object to fix to the desired TR
if TR>TR_kernel
    delay = mr.makeDelay((TR-1e3*TR_kernel)*1e-3);
    delay_flag = 1;
    seq.addBlock(delay);
    TR_kernel = (seq.duration);
    fprintf("TR: %2.2f ms\n", 1e3*TR_kernel);
end

clear seq
%% 
seq=mr.Sequence(sys);           
kernel_counter = 0;             % Keep the kernel counter check the TR 
seq_duration_each_kernel = [];  % Save the seq duration after each kernel 

%% Loop Over Phase Cycles
% This loop goes through each phase cycle and creates the sequence blocks.
for ipc = 1:length(pc)          
    pc_increment = pc(ipc);                 
    rf_phase_off = pc_increment * pi/180;    % Set the RF phase offset
    
    % Dummy Scans for Steady State
    % No ADC during dummy scans
    % This section adds dummy pulses to reach the steady state before the actual acquisition.
    for i=1:Ndummy
        rf_phase_off = mod(rf_phase_off + pc_increment * pi/180, 2*pi); 
        rf.phaseOffset = rf_phase_off;

        % Define the sequence blocks for each TR
        seq.addBlock(rf);  
        seq.addBlock(RO_pre_rise_and_flat)
        seq.addBlock(RO_fall_and_pre)
        
        % Optional delay for TR timing adjustment
        if delay_flag
            seq.addBlock(delay);
        end
    end

    % Loop over slice encoding direction     
    for iZ=1:Nz
        SE_pre = mr.makeTrapezoid('z','Area',areaZ(iZ), ...
                                    "duration",max_T_pre,'system',sys);
        SE_reph = mr.makeTrapezoid('z','Area',-areaZ(iZ), ...
                                    "duration",max_T_pre,'system',sys);
        
        % Loop over phase encoding direction     
        for i=1:Ny
            rf_phase_off = mod(rf_phase_off + pc_increment * pi/180, 2*pi);
            rf.phaseOffset=rf_phase_off;
            adc.phaseOffset=mod(rf_phase_off + pi/2, 2*pi);
            seq.addBlock(rf);
  
            PE_pre = mr.makeTrapezoid('y','Area',-1*areaY(i),"duration", max_T_pre,'system',sys);
            PE_reph = mr.makeTrapezoid('y','Area',areaY(i),"duration", max_T_pre,'system',sys);
    
            % Prephase and rephase and sample (ADC on)
            seq.addBlock(RO_pre_rise_and_flat,PE_pre, SE_pre, adc)  
            seq.addBlock(RO_fall_and_pre, PE_reph,SE_reph)

            % If we want a TR larger then min. obtainable TR, add delay
            if delay_flag
                seq.addBlock(delay);
            end
    
        end

    end

    fprintf("Current PC increment: %i\n", pc_increment);    

end

% Check and output the sequence timing
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

calc_TR = 1e3*seq.duration/(length(pc)*(Ny*Nz+Ndummy));
receive_BW = (1/adc.dwell)/Nx;

%% Export the Sequence
% Export the created sequence to a .seq file which can be used on an MRI scanner.
seq.setDefinition('FOV', fov);
seq.setDefinition('Name', 'trufi');

TR_str = strrep(sprintf('%.3f', TR),'.','p');
alpha_str = sprintf('%d',alpha);
adc_dur_str = sprintf('%d',adc_dur);
receive_BW_str = sprintf('%.2f', receive_BW);
fov_str = sprintf('%dx%dx%d', fov(1)*1000,fov(2)*1000,fov(3)*1000);

date = datetime('now');

date_str = datestr(date,'ddmmyyyy_HHMM');
seq_name = sprintf('pulseq_bca_TR%s_FA%s_adc%s_Nz%d_receive_BW%s_FOV%s_RE1x1x3_18PC_%s',...
                    TR_str,alpha_str,adc_dur_str,Nz,receive_BW_str,fov_str,date_str);

sprintf('sequence name: %s',seq_name)

seq.write( strcat(seq_name,'.seq') )    

%% Export the sequence parameters
params.ADC_duration = adc_dur;
params.MaxGrad = max_grad;
params.MaxSlew = max_slew;
params.FA = alpha;
params.RF_duration = rf_dur;
params.FOV = fov;
params.RES = res;
save(strcat(seq_name, '_PARAMS.mat'),'params')

%% Plot the sequence
seq.plot('TimeRange',[Ndummy-5 Ndummy+10]*1e-3*TR)
