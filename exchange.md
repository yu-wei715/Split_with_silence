points2S_sim = [0 500;250 700;550 900; 1000 1000; 1500 1100; 1750 1000; 2100 800; 2400 600; 2700 490; 2880 440];
%points2S_sim = [0 500;500 500;1000 500; 2000 500; 2880 500];
%points2S_sim = [0 500;1000 1200;2000 2400; 2600 2880];

imgW = 2880; imgH = 2880;
margin = 128; 
maxW = 1024; 
maxH = 1024;

%% 2. 呼叫核心演算法
tic;
rois = generate_strict_exclusive_rois(points2S_sim, margin, maxW, maxH);
dt = toc;
fprintf('ROI 生成完畢，耗時: %.4f 秒，共生成 %d 個 ROI\n', dt, size(rois,1));

%% 3. 繪圖驗證
figure('Color', 'w', 'Position', [100, 100, 1000, 800]);
hold on; box on;

% 繪製拼縫趨勢
plot(points2S_sim(:,1), points2S_sim(:,2), 'r.', 'MarkerSize', 10, 'DisplayName', 'Seam Points');

% 繪製各個 ROI
colors = jet(size(rois,1)); % 使用連續色盤
for i = 1:size(rois,1)
    r = rois(i,:);
    w = r(3) - r(1);
    h = r(4) - r(2);
    
    % 畫矩形框
    rectangle('Position', [r(1), r(2), w, h], 'EdgeColor', colors(i,:), 'LineWidth', 1.5);
    
    % 填充透明色以利觀察重疊 (理應完全對接無重疊)
    patch([r(1) r(3) r(3) r(1)], [r(2) r(2) r(4) r(4)], colors(i,:), ...
          'FaceAlpha', 0.15, 'EdgeColor', 'none');
      
    % 標註編號與尺寸
    text(r(1)+5, r(2)+15, sprintf('#%d [%dx%d]', i, w, h), ...
        'Color', colors(i,:), 'FontSize', 8, 'FontWeight', 'bold');
end

% 圖表修飾
set(gca, 'YDir', 'reverse'); % 影像座標
axis equal; 
xlim([0 imgW]); ylim([0 imgH]);
title(['Adaptive Path-based ROI Alignment (Count: ' num2str(size(rois,1)) ')']);
grid on;
legend('Location', 'northeastoutside');


% function finalROIs = generate_grid_rois_optimized(pts2D, margin, maxW, maxH, minSize)
%     % A. 點雲加密 (確保跨網格時不漏點)
%     densePts = interpolate_points(pts2D, margin/2);
% 
%     % B. 建立網格範圍
%     minX = min(densePts(:,1)) - margin;
%     maxX = max(densePts(:,1)) + margin;
%     minY = min(densePts(:,2)) - margin;
%     maxY = max(densePts(:,2)) + margin;
% 
%     xEdges = minX : maxW : (maxX + maxW);
%     yEdges = minY : maxH : (maxY + maxH);
% 
%     % C. Binning 優化: 將點預先分配到網格桶中 (複雜度 O(M))
%     numX = length(xEdges) - 1;
%     numY = length(yEdges) - 1;
%     bins = cell(numX, numY);
% 
%     % 計算每個點屬於哪個網格索引
%     idxX = floor((densePts(:,1) - minX) / maxW) + 1;
%     idxY = floor((densePts(:,2) - minY) / maxH) + 1;
% 
%     for k = 1:size(densePts,1)
%         if idxX(k) >= 1 && idxX(k) <= numX && idxY(k) >= 1 && idxY(k) <= numY
%             bins{idxX(k), idxY(k)} = [bins{idxX(k), idxY(k)}; densePts(k,:)];
%         end
%     end
% 
%     % D. 掃描網格並收縮
%     tempROIs = [];
%     for i = 1:numX
%         for j = 1:numY
%             ptsIn = bins{i,j};
% 
%             % 如果這塊網格內有拼縫點
%             if ~isempty(ptsIn)
%                 % 基礎邊界 (網格物理邊界)
%                 tx1 = xEdges(i); tx2 = xEdges(i+1);
%                 ty1 = yEdges(j); ty2 = yEdges(j+1);
% 
%                 % 收縮邊界至包含點雲 + margin，但不能超出網格
%                 rx1 = max(tx1, min(ptsIn(:,1)) - margin);
%                 ry1 = max(ty1, min(ptsIn(:,2)) - margin);
%                 rx2 = min(tx2, max(ptsIn(:,1)) + margin);
%                 ry2 = min(ty2, max(ptsIn(:,2)) + margin);
% 
%                 % --- 碎片處理邏輯 ---
%                 % 如果收縮後太小，則強行擴張回 minSize (除非會撞到隔壁，但網格內不會)
%                 if (rx2 - rx1) < minSize, rx2 = min(tx2, rx1 + minSize); end
%                 if (ry2 - ry1) < minSize, ry2 = min(ty2, ry1 + minSize); end
% 
%                 tempROIs = [tempROIs; round([rx1, ry1, rx2, ry2])];
%             end
%         end
%     end
% 
%     % E. 最後檢查與修正 (確保完全不重疊)
%     finalROIs = tempROIs;
% end
% 

function finalROIs = generate_strict_exclusive_rois(pts2D, margin, maxW, maxH)
    % --- Step 1: 高密度插值 ---
    densePts = interpolate_points(pts2D, margin/2);

    % --- Step 2: 嚴格互斥生長 ---
    finalROIs = [];
    % 初始化
    curr = [densePts(1,1)-margin, densePts(1,2)-margin, ...
            densePts(1,1)+margin, densePts(1,2)+margin];

    for i = 2:size(densePts,1)
        p = densePts(i,:);
        
        % 預測擴張
        nextX1 = min(curr(1), p(1)-margin);
        nextY1 = min(curr(2), p(2)-margin);
        nextX2 = max(curr(3), p(1)+margin);
        nextY2 = max(curr(4), p(2)+margin);

        if (nextX2 - nextX1) <= maxW && (nextY2 - nextY1) <= maxH
            curr = [nextX1, nextY1, nextX2, nextY2];
        else
            % 【關鍵：強迫切割，不准向回擴張】
            % 1. 將當前 ROI 存入 (已經取過整數)
            curr = round(curr);
            finalROIs = [finalROIs; curr];
            
            % 2. 決定切割方向與邊界
            % 如果拼縫主要水平移動，我們在 X 軸切開
            if abs(p(1) - curr(3)) > abs(p(2) - curr(4))
                x_cut = curr(3); % 切割線就在舊 ROI 的右邊界
                % 新 ROI 的 x1 嚴格從 x_cut 開始，不准往左伸
                newX1 = x_cut;
                newX2 = p(1) + margin;
                newY1 = p(2) - margin;
                newY2 = p(2) + margin;
            else
                y_cut = curr(4); % 切割線就在舊 ROI 的下邊界
                newX1 = p(1) - margin;
                newX2 = p(1) + margin;
                newY1 = y_cut;
                newY2 = p(2) + margin;
            end
            curr = [newX1, newY1, newX2, newY2];
        end
    end
    finalROIs = [finalROIs; round(curr)];
end

function dense = interpolate_points(pts, step)
    dense = [];
    for i = 1:size(pts,1)-1
        p1 = pts(i,:); p2 = pts(i+1,:);
        d = norm(p2-p1);
        n = max(1, ceil(d/step));
        t = linspace(0,1,n+1)';
        seg = (1-t)*p1 + t*p2;
        dense = [dense; seg(1:end-1,:)];
    end
    dense = [dense; pts(end,:)];
end
