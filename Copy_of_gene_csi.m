function csi_data=gene_csi(filepath)
    csi_trace=read_bf_file(filepath);
    for ii=1:length(csi_trace)
        csi_entry = csi_trace{ii};
        csi = get_scaled_csi(csi_entry);
        % Only consider measurements for transmitting on one antenna
        csi = csi(1, :, :);
        % Remove the single element dimension
        csi = squeeze(csi);
    %     [csi, phase_matrix] = remove_sto(csi, FI);
        csi=csi';
        csi_data(ii,:)=csi(:); 
    end
end

function ret = read_bf_file(filename)
%% Input check
narginchk(1,1);

%% Open file
f = fopen(filename, 'rb');
if (f < 0)
    error('Couldn''t open file %s', filename);
    return;
end

status = fseek(f, 0, 'eof');
if status ~= 0
    [msg, errno] = ferror(f);
    error('Error %d seeking: %s', errno, msg);
    fclose(f);
    return;
end
len = ftell(f);

status = fseek(f, 0, 'bof');
if status ~= 0
    [msg, errno] = ferror(f);
    error('Error %d seeking: %s', errno, msg);
    fclose(f);
    return;
end

%% Initialize variables
ret = cell(ceil(len/95),1);     % Holds the return values - 1x1 CSI is 95 bytes big, so this should be upper bound
cur = 0;                        % Current offset into file
count = 0;                      % Number of records output
broken_perm = 0;                % Flag marking whether we've encountered a broken CSI yet
triangle = [1 3 6];             % What perm should sum to for 1,2,3 antennas

%% Process all entries in file
% Need 3 bytes -- 2 byte size field and 1 byte code
while cur < (len - 3)
    % Read size and code
    field_len = fread(f, 1, 'uint16', 0, 'ieee-be');
    code = fread(f,1);
    cur = cur+3;
    
    % If unhandled code, skip (seek over) the record and continue
    if (code == 187) % get beamforming or phy data
        bytes = fread(f, field_len-1, 'uint8=>uint8');
        cur = cur + field_len - 1;
        if (length(bytes) ~= field_len-1)
            fclose(f);
            return;
        end
    else % skip all other info
        fseek(f, field_len - 1, 'cof');
        cur = cur + field_len - 1;
        continue;
    end
    
    if (code == 187) %hex2dec('bb')) Beamforming matrix -- output a record
        count = count + 1;
        ret{count} = read_bfee(bytes);
        
        perm = ret{count}.perm;
        Nrx = ret{count}.Nrx;
        if Nrx == 1 % No permuting needed for only 1 antenna
            continue;
        end
        if sum(perm) ~= triangle(Nrx) % matrix does not contain default values
            if broken_perm == 0
                broken_perm = 1;
                fprintf('WARN ONCE: Found CSI (%s) with Nrx=%d and invalid perm=[%s]\n', filename, Nrx, int2str(perm));
            end
        else
            ret{count}.csi(:,perm(1:Nrx),:) = ret{count}.csi(:,1:Nrx,:);
        end
    end
end
ret = ret(1:count);

%% Close file
fclose(f);
end


function ret = get_scaled_csi(csi_st)
    % Pull out CSI
    csi = csi_st.csi;

    % Calculate the scale factor between normalized CSI and RSSI (mW)
    csi_sq = csi .* conj(csi);
    csi_pwr = sum(csi_sq(:));
    rssi_pwr = dbinv(get_total_rss(csi_st));
%     rssi_pwr = 10;
    %   Scale CSI -> Signal power : rssi_pwr / (mean of csi_pwr)
    scale = rssi_pwr / (csi_pwr / 30);

    % Thermal noise might be undefined if the trace was
    % captured in monitor mode.
    % ... If so, set it to -92
    if (csi_st.noise == -127)
        noise_db = -92;
    else
        noise_db = csi_st.noise;
    end
     thermal_noise_pwr = dbinv(noise_db);
%      thermal_noise_pwr=0.1;
    
    % Quantization error: the coefficients in the matrices are
    % 8-bit signed numbers, max 127/-128 to min 0/1. Given that Intel
    % only uses a 6-bit ADC, I expect every entry to be off by about
    % +/- 1 (total across real & complex parts) per entry.
    %
    % The total power is then 1^2 = 1 per entry, and there are
    % Nrx*Ntx entries per carrier. We only want one carrier's worth of
    % error, since we only computed one carrier's worth of signal above.
    quant_error_pwr = scale * (csi_st.Nrx * csi_st.Ntx);

    % Total noise and error power
    total_noise_pwr = thermal_noise_pwr + quant_error_pwr;

    % Ret now has units of sqrt(SNR) just like H in textbooks
    ret = csi * sqrt(scale / total_noise_pwr);
    if csi_st.Ntx == 2
        ret = ret * sqrt(2);
    elseif csi_st.Ntx == 3
        % Note: this should be sqrt(3)~ 4.77 dB. But, 4.5 dB is how
        % Intel (and some other chip makers) approximate a factor of 3
        %
        % You may need to change this if your card does the right thing.
        ret = ret * sqrt(dbinv(4.5));
    end
end

function ret = get_total_rss(csi_st)
    narginchk(1,1);

    % Careful here: rssis could be zero
    rssi_mag = 0;
    if csi_st.rssi_a ~= 0
        rssi_mag = rssi_mag + dbinv(csi_st.rssi_a);
    end
    if csi_st.rssi_b ~= 0
        rssi_mag = rssi_mag + dbinv(csi_st.rssi_b);
    end
    if csi_st.rssi_c ~= 0
        rssi_mag = rssi_mag + dbinv(csi_st.rssi_c);
    end
    
    ret = db(rssi_mag, 'pow') - 44 - csi_st.agc;
end

function csi_db = dbinv(csi_st)
       csi_db = 10^(csi_st/10);
  
end