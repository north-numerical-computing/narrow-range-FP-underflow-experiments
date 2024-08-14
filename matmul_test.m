% matmul_test.m   Simulation of custom variants of Model 1 in [1].
%
% Requirements:
%   CPFloat (https://github.com/north-numerical-computing/cpfloat/).
%
% References:
%   [1] T. Mary and M. Mikaitis.
%       Error Analysis of Matrix Multiplication with Narrow Range
%       Floating-Point Arithmetic. hal-04671474. Aug. 2024.

% Set up the input format.
options_input.format = input_format;
[~, options] = cpfloat([], options_input);
% Grab various parameters of the format.
t    = options.params(1);
emin = options.params(2);
emax = options.params(3);
% Set up some useful quantities used in the paper.
u = 2^-t;
fmin = 2^emin;
% Special treatment for the maximum value of fp8-e4m3, as per the
% OFP8 format specification.
if (strcmp(input_format, 'fp8-e4m3'))
  fmax = 2^emax*(2-4*u);
else
  fmax = 2^emax*(2-2*u);
end

% Set up the accumulation format.
options_accum.format = accum_format;
[~, options] = cpfloat([], options_accum);
% Grab various parameters of the format.
T    = options.params(1);
Emin = options.params(2);
Emax = options.params(3);
% Set up some useful quantities used in the paper.
U = 2^(-T);
Fmin = 2^Emin;
Fmax = 2^Emax*(2-2*U);

% Matrix dimensions: A is m x n, B is n x q
m = 10;
q = 10;
nlist = floor(logspace(1,6,40));
% Matrix elements: uniformly distributed logarithms in [-l, l].
l = 10;

for subnormals_on = 0:1

    options_input.subnormal = subnormals_on;
    options_accum.subnormal = subnormals_on;

    i = 0;
    for n = nlist
        i=i+1;
        
        % Matrix elements: uniformly distributed logarithms in [-l, l].
        A = (10.^(rand(m,n)*2*l-l));
        B = (10.^(rand(n,q)*2*l-l));
        % Random sign + or - with equal probability
        sA = randi(2,m,n)*2-3;
        sB = randi(2,n,q)*2-3;
        A = A.*sA;
        B = B.*sB;
        
        % Compute a reference result in binary64.
        Ctrue = A*B;
        
        % Compute diagonal scaling matrices L and M such that
        % the elements of L*A and B*M are at most theta.
        theta = min(fmax, sqrt(Fmax/n));
        L = previous_pow2(theta./max(abs(A),[],2));
        M = previous_pow2(theta./max(abs(B)));
        Linv = 1./L;
        Minv = 1./M;
        
        % Round L*A and B*M to the input format.
        temp1 = L.*A;
        temp2 = B.*M;
        LA{1} = cpfloat(temp1, options_input);
        BM{1} = cpfloat(temp2, options_input);
        for j = 2:p
            temp1 = temp1 - LA{j-1}*u^(j-2);
            temp2 = temp2 - BM{j-1}*u^(j-2);
            LA{j} = cpfloat(temp1/u^(j-1), options_input); 
            BM{j} = cpfloat(temp2/u^(j-1), options_input);
        end
        
        % Compute LA*BM in the accumulation format.
        LABM = matmul(LA, BM, options_accum, p, u);
        
        % Scale LABM back to obtain C.
        C = Linv.*LABM.*Minv;
        
        % Compute the error
        err(i) = norm(C-Ctrue,'inf')/norm(A,'inf')/norm(B,'inf');
        
        % Compute the product as if we had no range limitations.
        options_input.explim = 0;
        options_accum.explim = 0;
        temp1 = L.*A;
        temp2 = B.*M;
        LA_nrl{1} = cpfloat(temp1, options_input);
        BM_nrl{1} = cpfloat(temp2, options_input);
        for j = 2:p
            temp1 = temp1 - LA_nrl{j-1}*u^(j-2);
            temp2 = temp2 - BM_nrl{j-1}*u^(j-2);
            LA_nrl{j} = cpfloat(temp1/u^(j-1), options_input); 
            BM_nrl{j} = cpfloat(temp2/u^(j-1), options_input);
        end
        LABM_nrl = matmul(LA_nrl, BM_nrl, options_accum, p, u);
        C_nrl = Linv.*LABM_nrl.*Minv;
        err_nrl(i) =...
            norm(C_nrl-Ctrue,'inf')/norm(A,'inf')/norm(B,'inf');
        options_input.explim = 1;
        options_accum.explim = 1;
        
        % Bound (3.26)
        bound(i) = 2*u + n*U + 4*n^2*fmin/theta + 4*n^2*Fmin/theta^2;
        bound_nrl(i) = 2*u + n*U;
        % Same bound but without dependency on n, which is quite
        % pessimistic.
        %bound(i) = 2*u + U + 4*fmin/theta + 4*Fmin/theta^2;
        %bound_nrl(i) = 2*u + U;
    end

    % Plot errors.
    figure;
    loglog(nlist,err,'-o');
    hold on
    loglog(nlist,bound,'--o');
    loglog(nlist,err_nrl,'-*');
    loglog(nlist,bound_nrl,'--*');
    legend('Error', 'Bound', 'Error (infinite range)',...
        'Bound (infinite range)','location','southwest');
    title(strcat(input_format,...
        '-', accum_format,...
        ' subnormals: ', num2str(subnormals_on), 'words: ', num2str(p)));
    hold off

    % Output various results to .dat files.
    filename = strcat('./data/matmul_test_', input_format,...
        '_', accum_format, '_subnormals',...
        num2str(subnormals_on), '_words_', num2str(p), '.dat');
    fileID = fopen(filename, 'w');
    fprintf(fileID, ...
        ['n error bound error-nrl bound-nrl \n']);
    for j=1:length(nlist)
        fprintf(fileID,'%d %e %e %e %e \n', ...
            nlist(j), err(j), bound(j), err_nrl(j), bound_nrl(j));
    end
end

function y = previous_pow2(x)
% Replace elements of x by the immediately inferior power of two.
y = 2.^floor(log2(x));
end

function C = matmul(A, B, options_accum, p, u)
  C = zeros(size(A{1},1), size(B{1},2));
  for j=1:p
      for k=1:p
          if (j+k-2 < p)
              for i=1:size(A{j},2)
                C = cpfloat(C + u^(j+k-2)*cpfloat(A{j}(:,i)*B{k}(i,:),...
                    options_accum), options_accum);
              end
          end
      end
  end
end
