function point = PAM_GRAY(b)
    % PAM_GRAY maps a binary vector to a PAM constellation point using Gray labeling.
    %
    %   point = PAM_GRAY(b) takes a vector of bits b (each element is 0 or 1)
    %   and returns the corresponding PAM constellation point, which will be one
    %   of the values in {±1, ±3, ..., ±(2^n - 1)} for an n-bit input.
    %
    % The algorithm uses recursion:
    %   - Base case: When the vector has one element, it returns 1 - 2*b(1).
    %   - Recursive case: For b with more than one element, the function uses the
    %     first bit to determine the sign and recursively computes the value for
    %     the remaining bits.
    
    if length(b) > 1
        % MATLAB uses 1-indexing; b(1) is the first element and b(2:end) is the rest.
        point = (1 - 2 * b(1)) * (2^(length(b)-1) - PAM_GRAY(b(2:end)));
    else
        point = 1 - 2 * b(1);
    end
end
% 
% % 1-bit tests:
% result1 = pam_gray([0]);   % expected output: 1
% result2 = pam_gray([1]);   % expected output: -1
% 
% % For a 2-bit vector, the expected outputs are computed as:
% % pam_gray([0,0]): (1 - 2*0) * (2^(1) - pam_gray([0])) = 1 * (2 - 1) = 1
% % pam_gray([0,1]): (1 - 2*0) * (2^(1) - pam_gray([1])) = 1 * (2 - (-1)) = 3
% % pam_gray([1,0]): (1 - 2*1) * (2^(1) - pam_gray([0])) = (-1) * (2 - 1) = -1
% % pam_gray([1,1]): (1 - 2*1) * (2^(1) - pam_gray([1])) = (-1) * (2 - (-1)) = -3
% result3 = pam_gray([0,0]);
% result4 = pam_gray([0,1]);
% result5 = pam_gray([1,0]);
% result6 = pam_gray([1,1])