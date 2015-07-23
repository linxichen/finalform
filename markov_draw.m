function [i_tmr] = markov_draw(i_now,CDF_mat,varargin)
% Given the square transition CDF, index of state today, draw an index of tomorrow.
% You can supply your own random number in [0,1]:
%   markov_draw(i_now,CDF_mat,u), where u = rand.
% Or ask the function to draw one for you if lazy:
%   markov_draw(i_now,CDF_mat)
% Assuming CDF_mat(i,j) = Prob(x=j|x=i)!!!

% Some error checking
n = size(CDF_mat,1);
m = size(CDF_mat,2);
if ( n ~= m || ismatrix(CDF_mat) ~= 1)
	error('CDF not a square matrix!')
end

cdf = [0,CDF_mat(i_now,:)];
switch nargin
	case 2 % user didn't supply a random draw between [0,1]
		u = rand;
		if u == 1
			i_tmr = n;
		else
			i_tmr = find(u<cdf,1,'first')-1;
		end
	case 3 % user supplied something
		u = varargin{1};
		if (u < 0 || u > 1)
			error('random number is not in [0,1]');
		end
		if u == 1
			i_tmr = n;
		else
			i_tmr = find(u<cdf,1,'first')-1;
		end
	otherwise
		error('number of arguments is not right!!!')
end

end
