%%
% Central moment
% Author: Thanh Phuong NGUYEN, U2IS-ENSTA Paristech
%

%
function [M1, M2, M3, M4] = centralmoment(varargin) 
% Version: 0.1.0


% Check number of input arguments.
%error(nargchk(1,5,nargin));

image=varargin{1};
d_image=double(image);

if nargin==1
    spoints=[-1 -1; -1 0; -1 1; 0 -1; -0 1; 1 -1; 1 0; 1 1];
    neighbors=8;
    lims=0;
    mode='i';
end


if (nargin > 2) 
    radius=varargin{2};
    neighbors=varargin{3};
    spoints=zeros(sum(neighbors),2);
    lims=0;
    % Angle step.
    id=1;
    for k=1:length(neighbors)
	    a = 2*pi/neighbors(k);    
	    for i = 1:neighbors(k)
		spoints(id,1) = -radius(k)*sin((i-1)*a);
		spoints(id,2) = radius(k)*cos((i-1)*a);
		id=id+1;
	    end
    end
    
  
    if(nargin >= 4)
        weight=varargin{4};
    end
    if nargin==5
	mode=varargin{5};
    else
	mode='mean';
    end
    
end

if (nargin == 2) && ischar(varargin{2})
    mode=varargin{2};
    spoints=[-1 -1; -1 0; -1 1; 0 -1; -0 1; 1 -1; 1 0; 1 1];
    neighbors=8;

end




% Determine the dimensions of the input image.
[ysize xsize] = size(image);

miny=min(spoints(:,1));
maxy=max(spoints(:,1));
minx=min(spoints(:,2));
maxx=max(spoints(:,2));

% Block size, each LBP code is computed within a block of size bsizey*bsizex
bsizey=ceil(max(maxy,0))-floor(min(miny,0))+1;
bsizex=ceil(max(maxx,0))-floor(min(minx,0))+1;

% Coordinates of origin (0,0) in the block
origy=1-floor(min(miny,0));
origx=1-floor(min(minx,0));

% Minimum allowed size for the input image depends
% on the radius of the used LBP operator.
if(xsize < bsizex || ysize < bsizey)
    error('Too small input image. Should be at least (2*radius+1) x (2*radius+1)');
end

% Calculate dx and dy;
dx = xsize - bsizex;
dy = ysize - bsizey;


% Fill the center pixel matrix C.
C = image(origy:origy+dy,origx:origx+dx);
d_C = double(C);


%Compute the local contrast

	for i = 1:sum(neighbors)+1
	    if i<=sum(neighbors)
		    y = spoints(i,1)+origy;
		    x = spoints(i,2)+origx;
		    % Calculate floors and ceils for the x and y.
		    fy = floor(y); cy = ceil(y);
		    fx = floor(x); cx = ceil(x);
		    
		    % Use double type images
		    ty = y - fy;
		    tx = x - fx;
		    
		    % Calculate the interpolation weights.
		    w1 = (1 - tx) * (1 - ty);
		    w2 =      tx  * (1 - ty);
		    w3 = (1 - tx) *      ty ;
		    w4 =      tx  *      ty ;
		    % Compute interpolated pixel values
		    N = w1*d_image(fy:fy+dy,fx:fx+dx) + w2*d_image(fy:fy+dy,cx:cx+dx) + ...
			w3*d_image(cy:cy+dy,fx:fx+dx) + w4*d_image(cy:cy+dy,cx:cx+dx);
		    % Compute the variance using on-line algorithm
		    % ( http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#On-line_algorithm ).
	    else
		  N=d_C;
	    end

   	    if i == 1
			M1=zeros(size(N));
			DELTA=zeros(size(N));
			M2=zeros(size(N));
			M3=zeros(size(N));
			M4=zeros(size(N));
	    end   
	    DELTA=N-M1;
	    DELTA_i= DELTA/i;
	    DELTA_i2=DELTA_i.*DELTA_i;	 
	    M1=M1+DELTA_i;
	    term=DELTA.*DELTA_i*(i-1);	
	    M4=M2+term.*DELTA_i2*(i*i-3*i+3)+6*DELTA_i2.*M2-4.*DELTA_i.*M3;	
	    M3=M3+term.*DELTA_i*(i-2)-3*DELTA_i.*M2;
	    M2=M2+term;
	end



