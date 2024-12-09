% 设定采样频率
Fs = 500;

% 设计滤波器
Wn = [5 15] / (Fs/2);
[b, a] = butter(3, Wn, 'bandpass');

% 文件夹路径数组
folder_paths = {'D:\Sensing\CW1\ECG_AF', 'D:\Sensing\CW1\ECG_normal', 'D:\Sensing\CW1\ECG_STD'};

% 定义ST段的开始和结束阈值
ST_start_threshold = round(0.06 * Fs); % 假设ST段在S波结束后80ms开始
ST_end_threshold = round(0.06 * Fs);   % 假设ST段在T波开始前80ms结束

% 遍历每个文件夹路径下的文件进行处理
for folder_idx = 1:length(folder_paths)
    folder_path = folder_paths{folder_idx};
    files = dir(fullfile(folder_path, '*.mat'));
    
    % 遍历每个文件
    for i = 1:length(files)
        file_name = files(i).name;
        variable_name = file_name(1:end-4);  % 去除.mat扩展名
        file_data = load(fullfile(folder_path, file_name));  % 加载.mat文件
        
        % 假设数据在ECG.data结构体中
        ECG_data = file_data.ECG.data;
        
        % 对所有导联进行滤波
        filtered_data_all_leads = zeros(size(ECG_data));
        for lead = 1:size(ECG_data, 1)
            filtered_data_all_leads(lead, :) = filtfilt(b, a, ECG_data(lead, :));
        end
        
        % 选择第二导联来进行QRS检测
        lead2_data = filtered_data_all_leads(2, :);
        
        % R波检测阈值计算
        [pks, locs] = findpeaks(lead2_data, 'MinPeakDistance', 0.2*Fs);
        thr = mean(pks) * 0.5; % 50%的平均峰值作为阈值
        r_locs = locs(pks > thr);
        
        % 初始化Q波和S波位置数组
        q_locs = zeros(size(r_locs));
        s_locs = zeros(size(r_locs));
        
        % 搜索Q波和S波
        for j = 1:length(r_locs)
            % Q波通常位于R波之前
            q_window = max(1, r_locs(j) - round(0.12*Fs)):r_locs(j);
            [~, q_index] = min(lead2_data(q_window));
            q_locs(j) = q_window(q_index);
            
            % S波通常位于R波之后
            s_window = r_locs(j):min(length(lead2_data), r_locs(j) + round(0.12*Fs));
            [~, s_index] = min(lead2_data(s_window));
            s_locs(j) = s_window(s_index);
        end
        
        % 搜索T波，作为S波之后的第四个峰值
        t_locs = zeros(size(r_locs));
        for j = 1:length(s_locs)
            if j < length(s_locs)
                % T波通常在S波之后，但在下一个R波之前
                t_window = s_locs(j):r_locs(j+1);
            else
                % 对于最后一个S波，我们可能没有下一个R波的位置
                t_window = s_locs(j):min(length(lead2_data), s_locs(j) + round(0.4*Fs));
            end
            
            % 确保窗口内有足够的数据点来查找峰值
            if length(t_window) >= 3
                [pks_window, locs_window] = findpeaks(lead2_data(t_window));
                if length(locs_window) >= 4
                    t_locs(j) = t_window(locs_window(4)); % 第四个峰值作为T波
                elseif length(locs_window) > 0
                    t_locs(j) = t_window(locs_window(end)); % 如果没有四个峰值，则取最后一个峰值
                end
            end
            % 分析ST段
        for j = 1:length(s_locs)-1
            st_start = s_locs(j) + ST_start_threshold;
            st_end = t_locs(j) - ST_end_threshold;
            st_start = max(st_start, 1);
            st_end = min(st_end, length(lead2_data));
            
           % Calculate ST segment mean
            st_mean = mean(lead2_data(st_start:st_end));
            
            % 计算基线水平
            baseline_start = max(1, q_locs(j) - ST_start_threshold);
            baseline_end = s_locs(j);
            baseline_level = mean(lead2_data(baseline_start:baseline_end));
            
%             % 计算ST均值与基线的差值
%             st_difference = st_mean - baseline_level;
%             st_diff_abs = (st_difference);
            % Save the difference into the ST feature array
            st_features(j) = st_mean;
        end
        
        % Save ST segment characteristics to ST folder
        ST_folder = fullfile('D:\Sensing\CW1\Feature\ST'); % 指定保存ST段特征的文件夹路径
        save(fullfile(ST_folder, [variable_name '_ST.mat']), 'st_features');
        end
        
        % 对于特定文件进行ST段和基线的可视化
        if strcmp(variable_name, 'A0008') || strcmp(variable_name, 'A0002')
            figure;
            plot(lead2_data);
            hold on;
            plot(r_locs, lead2_data(r_locs), 'ro', 'MarkerFaceColor', 'r'); % R波
            plot(s_locs, lead2_data(s_locs), 'gs', 'MarkerFaceColor', 'g'); % S波
            plot(t_locs, lead2_data(t_locs), 'md', 'MarkerFaceColor', 'm'); % T波
            
            % 计算和显示ST段和基线
            for j = 1:length(s_locs)-1
                % 定义ST段的实际开始和结束位置
                st_start = s_locs(j) + ST_start_threshold;
                st_end = t_locs(j) - ST_end_threshold;
                
                % 确保ST段不超过数组的界限
                st_start = max(st_start, 1);
                st_end = min(st_end, length(lead2_data));
                
                % 提取ST段
                st_segment = st_start:st_end;
                
                % 计算基线水平
                baseline = mean(lead2_data(max(1, q_locs(j) - ST_start_threshold):s_locs(j)));
                
                % 画出ST段和基线
                plot(st_segment, lead2_data(st_segment), 'y', 'LineWidth', 2); % ST段
                line([st_start st_end], [baseline baseline], 'Color', 'k', 'LineWidth', 2); % 基线
            end
            
            hold off;
            title(['ECG Lead II with R, S, and T Waves for ' variable_name]);
            xlabel('Sample Number');
            ylabel('Voltage (mV)');
            legend('ECG Signal', 'R Peaks', 'S Waves', 'T Waves', 'ST Segment', 'Baseline');
            grid on;
            axis tight;
        end
    end
end





