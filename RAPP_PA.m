function g = RAPP_PA(A,A0,v,p)
   % RAPP PA model implementation
   % Inputs:
   %   A  - Input signal amplitude
   %   A0 - Limiting output amplitude
   %   v  - Small signal gain
   %   p  - Smoothness parameter
    
   % Output:
   %   g  - Output g(A) function 

   % Ensure parameters are valid
   if A0 <= 0 || v < 0 || p <= 0
       error('Invalid parameters: Ensure A0 > 0, v >= 0, and p > 0.');
   end
   % Compute g(A) based on the given RAPP formula
   g = v .* A ./ ((1 + (v.*A ./ A0).^(2 * p)).^(1 / (2 * p)));

end

