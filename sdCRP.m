function z = sdCRP(data, alpha, precision, num, maxIter)
%% Init
if nargin < 4
    num = 100;
end
if nargin < 5
    maxIter = 100;
end

% m is the number of observations, and d is the dimension
[m,d] = size(data);

% table assignment. Init all the observations to be in the first cluster
z = ones(1,m);

% init the centers
mu = mean(data);
mu = repmat(mu, num, 1);

% init the precision of each cluster
Lambda = eye(d);

% init links
c = ones(m,1);
c(2:end) = 1:(m-1);

% get similarity for every pair of observations
S = get_similarity(data);
S(1:(m+1):end) = alpha;

%% Gibbs sampling 
for iter = 1:maxIter
    tic
    %% update links c, sitBehind and table assignment z
    for pos = 1:m
        
        %%%%%%%%%%%%%%%% get the customers that behind pos %%%%%%%%%%%%%%%%
        sitBehind = get_sitBehind(pos, c);
        
        %%%%%%%%%% in the first case, no new table can be generated %%%%%%%
        if ismember(c(pos), sitBehind)
            % compute the likelihood of pos to sit with each table
            t = unique(z);
            p2 = zeros(1, length(t));
            for k = 1:length(t)
                % fetch those sit behind pos
                X = data(sitBehind,:);
                X = X - repmat(mu(t(k),:), length(sitBehind), 1);
                
                % get likelihood
                p2(k) = - 1/2 * sum(sum(X' .* (Lambda * X')));
            end
            
            % compute the probability of pos to sit with j
            p = zeros(1,m);
            for j = 1:m
                % find the table of j
                ix = t == z(j);
                p(j) = p2(ix) + log(S(pos,j));
            end
            
            % normalize p
            p = p - max(p);
            p = exp(p);
            p = p / sum(p);
            
            % update c(pos)
            [~,~,c(pos)] = histcounts(rand(1), [0, cumsum(p)]);
            
            % update z(pos)
            for j = 1:length(sitBehind)
                z(sitBehind(j)) = z(c(pos));
            end
            
        %%%%%%%%%% in the second case, there can be new tables %%%%%%%%%%%%   
        else
            t = [unique(z), 0];
            % the likelihood of each table
            p2 = zeros(1, length(t));
            
            % likelihood for the existing tables
            X = data(sitBehind,:);
            for k = 1:(length(p2)-1)   
                X = X - repmat(mu(t(k),:), length(sitBehind), 1);
                
                % get likelihood
                p2(k) = -1/2 * sum(sum(X' .* (Lambda * X')));
            end
            
            % likelihood for the new table
            nk = length(sitBehind);
            Lambda_1 = nk * Lambda + precision * eye(d);
            if size(X, 1) > 1
                p2(end) = nk/2 * log(det(Lambda)) - d * nk/2 * log(2*pi) ...
                    - 1/2 * log(det(Lambda_1)) ...
                    - 1/2 * sum(sum(X' .* (Lambda * X'))) ...
                    - 1/2 * sum(X) * Lambda' / Lambda_1 * Lambda * (sum(X))';
            else
                p2(end) = nk/2 * log(det(Lambda)) - d * nk/2 * log(2*pi) ...
                    - 1/2 * log(det(Lambda_1)) ...
                    - 1/2 * sum(sum(X' .* (Lambda * X'))) ...
                    - 1/2 * X * Lambda' / Lambda_1 * Lambda * (X)';
            end
             
            % compute the prob of pos to sit with each j
            p = zeros(1,m);
            for j = 1:m
                if ismember(j, sitBehind)
                    p(j) = log(S(pos,j)) + p2(end);
                else
                    ix = t == z(j);
                    p(j) = log(S(pos,j)) + p2(ix);
                end
            end
            
            % normalize p
            p = p - max(p);
            p = exp(p);
            p = p / sum(p);
            
            % compute c(pos)
            [~, ~, c(pos)] = histcounts(rand(1), [0, cumsum(p)]);
            
            % update z(pos)
            null_table = min(setdiff(1:num, z));
            if ismember(c(pos), sitBehind)
                for j = 1:length(sitBehind)
                    z(sitBehind(j)) = null_table;
                end
            else
                for j = 1:length(sitBehind)
                    z(sitBehind(j)) = z(c(pos));
                end
            end
            
        end
    end
    
    %% update parameters
    t = unique(z);
    for k = 1:length(t)
        if size(data(z == t(k),:), 1) > 1
            mu(t(k),:) = mean(data(z == t(k),:));
%             mu(t(k),:) = mvnrnd(zeros(d,1), eye(d) * 5);
        else
            mu(t(k),:) = data(z == t(k),:);
        end
    end
    
    fprintf(['iter ', num2str(iter), ' done\n'])
    toc
end