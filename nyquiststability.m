function [fr_,w_]=nyquiststability(b,a,offset,samples)
%NYQUISTSTABILITY Nyquist diagram with limit closures at infinity for
%closed-loop stability evaluation with negative feedback.
%   
%   NYQUISTSTABILITY(B,A) draws the Nyquist diagram with limit
%   closures of the dynamic system from the transfer function parameters
%   B and A. The function assumes a polynomial transfer function
%   representation.
%   
%             B(s)     s^m*b_m + s^(m-1)*b_(m-1) + ... + s*b_1 + b_0
%     F(s) =  ----  =  ---------------------------------------------
%             A(s)     s^n*a_n + s^(n-1)+a_(n-1) + ... + s*a_1 + a_0
%   
%   NYQUISTSTABILITY(B,A,OFFSET) OFFSET determines a radius of a circle
%   going around the imaginary poles. The default value is set in
%   proportion to the smallest imaginary pole frequency.
%   
%   NYQUISTSTABILITY(B,A,OFFSET,SAMPLES) draws the Nyquist diagram with
%   limit transitions. SAMPLES specify a number of samples for each
%   evaluated frequency segment. The default value is SAMPLES = 1000.
%   
%   [FR,WOUT] = NYQUISTSTABILITY(B,A)
%   [FR,WOUT] = NYQUISTSTABILITY(B,A,OFFSET)
%   [FR,WOUT] = NYQUISTSTABILITY(B,A,OFFSET,N) returns a complex frquency
%   response FR evaluated for frequencies WOUT.

% Author: Jakub Gaier

if isempty(a) | isempty(b)
    error("Either polynom a or b is empty")
end

% poles
p = roots(a);
ad = 3;

eps = 1e-12;
% poles on imaginary axis
imp = p(real(p) < eps & real(p) > -eps);
% use only imaginary part
imp = imag(imp);

% has integrators flag
hi = sum(imp < eps & imp > -eps) > 0;
% remove integrators and complex conjugates
impp = sort(imp(imp > eps), "descend");

if nargin < 4
    % set n to default when not provided
    samples=1000;
end
n = ceil(samples/2);

if nargin < 3
    m = impp;
    % avoid empty array
    if isempty(m)
        m = 1;
    end
    % offset set by smallest imaginary pole
    offset = 10^(round(log10( min(abs(m))/10 )));
end

m = p;
% avoid logarithm of zero
if max(abs(m)) == 0
    m = 1;
end
% maximal frequency set by biggiest pole
wmax = 10^(floor(log10( max(abs(m))*10^ad )));

% frequency bounds to avoid frequencies of imaginary poles
wbounds = zeros(2 + 2*length(impp),1);

wbounds(1) = wmax;

for i=1:length(impp)
    wbounds(i*2) = impp(i) + offset;
    wbounds((i*2)+1) = impp(i) - offset;
end

if hi
    wbounds(end) = offset;
else
    wbounds(end) = 0;
end

% frequencies to evaluate
w = complex(0, zeros(samples*(1+2*length(impp))+hi*n,1));

w(1) = complex(0, -Inf);
% transition - logarithmically spaced values
w(2:samples) = complex(0, -logspace( ...
    log10(wbounds(1)), max(log10(wbounds(2)),-ad), samples-1 ));
    % max to avoid logarithm of zero

% avoid pole and transition
if ~isempty(impp)
    % create half-circle with radius equal to offset
    y = linspace(0,pi,samples+1);
    circ = complex(sin(y), -cos(y))*offset;
    
    for i=1:length(impp)
        % avoid pole
        w((((i*2)-1)*samples):i*2*samples) = circ - complex(0, impp(i));
        % transition
        w((i*2*samples):((i*2)+1)*samples) = complex(0, -logspace( ...
            log10(wbounds((i*2)+1)), ...
            max(log10(wbounds((i+1)*2)),-ad), samples+1) );
    end
end

% avoid integrator
if hi
    % create quarter-circle with radius equal to offset
    y = linspace(0,pi/2,n+1);
    circ = complex(sin(y), -cos(y))*offset;
    w(end-n:end) = circ;
else
    % logarithm can't produce zero
    w(end) = 0;
end

% evaluate frequency response for given (complex) frequencies
fr = polyval(b,w)./polyval(a,w);

% system has relative order r=0
if length(a) == length(b)
    fr(1) = b(1)/a(1);
else
    fr(1) = 0;
end

% return data or plot
if nargout > 0
    % mirrored response for positive frequencies
    fr = flip(complex(real(fr), imag(fr)*-1));
    w = flip(complex(real(w), imag(w)*-1));

    fr_ = fr;
    w_ = w;
else
    clf
    hold on
    style = ["--","-"];
    
    % indicies for arrows so they are in the middle of segments
    idx = zeros(1+hi+2*length(impp),1);
    for i=1:2*length(impp)
        idx(i+1) = samples*i + n;
    end

    if hi
        idx(end) = length(w) - floor(n/2);
    else
        idx(end) = length(w) - samples + floor(n/10);
    end
    if isempty(imp)
        idx(1) = n;
    else
        idx(1) = samples - floor(n/10);
    end
    
    for i=1:2
        % plot frequency response
        pl = plot(fr, "LineStyle",style(i), "Color","#0072BD");
        
        % custom data tips
        wtip = split(sprintf("%.2f ", real(w))) + ...
            split(sprintf("%-+.2fj ", imag(w)));
        pl.DataTipTemplate.DataTipRows(1) = ...
            dataTipTextRow("Re\{F(\omega)\}", real(fr), "%.2f");
        pl.DataTipTemplate.DataTipRows(2) = ...
            dataTipTextRow("Im\{F(\omega)\}", imag(fr), "%.2f");
        pl.DataTipTemplate.DataTipRows(3) = ...
            dataTipTextRow("\omega", wtip(1:end-1));
        
        % plot arrows
        origin = [real(fr(idx)), imag(fr(idx))];
        direction = ([real(fr(idx+1)), imag(fr(idx+1))]-origin)/2;
        for j=1:length(idx)
            annotation("arrow", "Position",[origin(j,:) direction(j,:)], ...
                "Parent",gca, "Color","#0072BD", "HeadStyle","vback1", ...
                "HeadLength",5, "HeadWidth",8);
        end
        
        % mirrored response for positive frequencies
        fr = flip(complex(real(fr), imag(fr)*-1));
        w = flip(complex(real(w), imag(w)*-1));
        idx = length(w)-idx;
    end
    
    xline(0, "LineStyle","-.", "Color","#BBBBBB")
    yline(0, "LineStyle","-.", "Color","#BBBBBB")
    pl = scatter(-1,0, "Marker","+", "MarkerEdgeColor","red", ...
        "LineWidth",1);
    pl.DataTipTemplate.DataTipRows(3:end) = "";
    
    title("Nyquist Stability Diagram")
    xlabel("Real Axis")
    ylabel("Imaginary Axis")
end

end

