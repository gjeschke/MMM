function [fi,dfx,d2fx,d2fxy] = triparinterp(x,y,z,f,nodes)

% TRIPARINTERP 3-D piecewise parabolic interpolation/extrapolation of  
%    tabulated data and calculation of first and second derivatives. 
%    [FI,DFX,D2FX,D2FXY] = TRIPARINTERP(X,Y,Z,F,NODES) interpolates 
%    function values, first and second derivatives using the vectors 
%    X, Y and Z and the matrix F, at the nodes specified in the 
%    three-column matrix NODES. X, Y and Z can be row or column vectors, 
%    but the code expects that 
%
%              size( F ) = [length(X) length(Y) length(Z)].
% 
%    Each output is returned as a column vector. 
%    X and Y can be non-uniformly spaced, but repeated values will 
%    lead to NaN output. 
%    This function also extrapolates the input data (no NaNs are 
%    returned); use extrapolated values at your own risk. 
%
%    Requires bindex.m (available at Matlab FEX).  
% 
% Example: nx = 101; x = linspace(0,2*pi,nx); 
%          ny =  91; y = linspace(0,2*pi,ny);
%          nz =  81; z = linspace(0,2*pi,nz);  
%          [X,Y,Z] = ndgrid(x,y,z);
%          F  = ( sin( X ).*cos( Y ) ).*sin( Z ); 
%          m = 21; 
%          p = sort( 2*pi*rand(m,3) );
%          Fi  =  ( sin( p(:,1) ).*cos( p(:,2) ) ).*sin( p(:,3) );
%          Fx  =  ( cos( p(:,1) ).*cos( p(:,2) ) ).*sin( p(:,3) );
%          Fxx = -( sin( p(:,1) ).*cos( p(:,2) ) ).*sin( p(:,3) );
%          Fy  = -( sin( p(:,1) ).*sin( p(:,2) ) ).*sin( p(:,3) );
%          Fyy = -( sin( p(:,1) ).*cos( p(:,2) ) ).*sin( p(:,3) );
%          Fz  =  ( sin( p(:,1) ).*cos( p(:,2) ) ).*cos( p(:,3) );
%          Fzz = -( sin( p(:,1) ).*cos( p(:,2) ) ).*sin( p(:,3) );
%          Fxy = -( cos( p(:,1) ).*sin( p(:,2) ) ).*sin( p(:,3) );
%          Fxz =  ( cos( p(:,1) ).*cos( p(:,2) ) ).*cos( p(:,3) );
%          Fyz = -( sin( p(:,1) ).*sin( p(:,2) ) ).*cos( p(:,3) );
%          [fi,dfx,d2fx,d2fxy] = triparinterp(x,y,z,F,p);
%          [Fi(:),fi(:)]
%          [Fx(:),dfx(:,1)]
%          [Fy(:),dfx(:,2)]
%          [Fz(:),dfx(:,3)]
%          [Fxx(:),d2fx(:,1)]
%          [Fyy(:),d2fx(:,2)]
%          [Fzz(:),d2fx(:,3)]
%          [Fxy(:),d2fxy(:,1)]
%          [Fxz(:),d2fxy(:,2)]
%          [Fyz(:),d2fxy(:,3)]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First  version: 11/11/2007
% Second version: 16/11/2007
% Third  version: 23/11/2007
% 
% Contact: orodrig@ualg.pt
% 
% Any suggestions to improve the performance of this 
% code will be greatly appreciated. 
%
% Algorithm implemented: given (xi,yi,zi) find the three points along X, Y and Z, 
% given in the tabulated data, such that 
%                      x(i) <= xi < x(i+1) < x(i+2)  
%                      y(j) <= yi < y(j+1) < y(j+2)   
%                      z(k) <= zi < z(k+1) < z(k+2). 
% Let us call them (x1,x2,x3), (y1,y2,y3) and (z1,z2,z3).  
% Corresponding function values will be 
% (x1,y1,z1) -> f111, (x2,y1,z1) -> f211, (x3,y1,z1) -> f311, 
% (x1,y2,z1) -> f121, (x2,y2,z1) -> f221, (x3,y2,z1) -> f321, 
% (x1,y3,z1) -> f131, (x2,y3,z1) -> f231, (x3,y3,z1) -> f331,
%                 ... (more 18 points) ...   
% 
% Between those points the function can be interpolated as 
% 
% f(xi,yi,zi) = a111( xi - x2 )( xi - x3 )( yi - y2 )( yi - y3 )( zi - z2 )( zi - z3 ) 
%            +  a211( xi - x1 )( xi - x3 )( yi - y2 )( yi - y3 )( zi - z2 )( zi - z3 )
%            +  a311( xi - x1 )( xi - x2 )( yi - y2 )( yi - y3 )( zi - z2 )( zi - z3 )
%            + ... (more 24 terms) ...
% 
% From the conditions f(x1,y1,z1) = f111, f(x2,y1,z1) = f211, f(x3,y1,z1) = f311 ...
% it follows that 
%   
%                                 f111
% a111 = ------------------------------------------------------   
%        (x1 - x2)(x1 - x3)(y1 - y2)(y1 - y3)(z1 - z2)(z1 - z3)   
%
%                                 f211
% a211 = ------------------------------------------------------   
%        (x2 - x1)(x2 - x3)(y1 - y2)(y1 - y3)(z1 - z2)(z1 - z3)   
% 
%                                 ....
% 
% The first derivatives become  
% 
%  fx(xi,yi,zi) = a111( 2xi - x2 - x3 )( yi - y2 )( yi - y3 )( zi - z2 )( zi - z3 )
%               + a211( 2xi - x1 - x3 )( yi - y2 )( yi - y3 )( zi - z2 )( zi - z3 )
%               + ...
%
%  fy(xi,yi,zi) = a111( xi - x2 )( xi - x3 )( 2yi - y2 - y3 )( zi - z2 )( zi - z3 )
%               + a211( xi - x1 )( xi - x3 )( 2yi - y2 - y3 )( zi - z2 )( zi - z3 )
%               + ...
%
%  fz(xi,yi,zi) = a111( xi - x2 )( xi - x3 )( yi - y2 )( yi - y3 )( 2zi - z2 - z3 )
%               + a211( xi - x1 )( xi - x3 )( yi - y2 )( yi - y3 )( 2zi - z2 - z3 )
%               + ...
% 
% The second derivatives become 
% 
%  fxx(xi,yi,zi) = 2*a111( yi - y2 )( yi - y3 )( zi - z2 )( zi - z3 )
%                + 2*a211( yi - y2 )( yi - y3 )( zi - z2 )( zi - z3 )
%                + ...
%
%  fyy(xi,yi,zi) = 2*a111( xi - x2 )( xi - x3 )( zi - z2 )( zi - z3 )
%                + 2*a211( xi - x1 )( xi - x3 )( zi - z2 )( zi - z3 )
%                + ...
%
%  fzz(xi,yi,zi) = 2*a111( xi - x2 )( xi - x3 )( yi - y2 )( yi - y3 )
%                + 2*a211( xi - x1 )( xi - x3 )( yi - y2 )( yi - y3 )
%                + ...
%  fxy(xi,yi,zi) = a111( 2xi - x2 - x3 )( 2yi - y2 - y3 )( zi - z2 )( zi - z3 )
%                + a211( 2xi - x1 - x3 )( 2yi - y2 - y3 )( zi - z2 )( zi - z3 )
%                + ...
%
%  fxz(xi,yi,zi) = a111( 2xi - x2 - x3 )( yi - y2 )( yi - y3 )( 2zi - z2 - z3 )
%                + a211( 2xi - x1 - x3 )( yi - y2 )( yi - y3 )( 2zi - z2 - z3 )
%                + ...
%
%  fyz(xi,yi,zi) = a111( xi - x2 )( xi - x3 )( 2yi - y2 - y3 )( 2zi - z2 - z3 )
%                + a211( xi - x1 )( xi - x3 )( 2yi - y2 - y3 )( 2zi - z2 - z3 )
%                + ...
% 
% Extrapolation is included in order to avoid NaNs when 
% xi is outside the interval [x(1),x(nx)] 
% yi is outside the interval [y(1),y(ny)] 
% zi is outside the interval [z(1),z(nz)]. 
% However, it is up to the user to decide it the extrapolated values are of any use.  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Let's go: 
 
% Declare an empty output just in case calculations won't be 
% performed: 

   fi = [];
  dfx = [];
 d2fx = [];
d2fxy = [];  

% Get the coordinates of the interpolation points: 

xi = nodes(:,1)'; 
yi = nodes(:,2)'; 
zi = nodes(:,3)';
 
% Be sure that x, y and z are row vectors:   

x  =  x(:)'; 
y  =  y(:)';
z  =  z(:)';

nx = length(x);
ny = length(y);
nz = length(z); 

sx = size(x);
sy = size(y);
sz = size(z); 
sf = size(f); 

m  = length(xi);

% Error checking: 

if nargin ~= 5 
error('This function requires five input arguments...');
end

if nx ~= sf(1) 
error('The length of x must match the first dimension of f...')
end 

if ny ~= sf(2) 
error('The length of y must match the second dimension of f...')
end 

if nz ~= sf(3) 
error('The length of z must match the third dimension of f...')
end 

if min(sx) ~= 1 
error('Input x should be a row or column vector...')
end

if min(sy) ~= 1 
error('Input y should be a row or column vector...')
end 

if min(sz) ~= 1 
error('Input z should be a row or column vector...')
end 

% Proceed with calculations only when the # of function points 
% is greater than 2:
 
if ( nx > 2 )&( ny > 2 )&( nz > 2 ) 
    
% Preliminary allocation: 

      fi = zeros(m,1);
     dfx = zeros(m,3);  fxi = zeros(m,1);  fyi = zeros(m,1);  fzi = zeros(m,1);
    d2fx = zeros(m,3); fxxi = zeros(m,1); fyyi = zeros(m,1); fzzi = zeros(m,1);
   d2fxy = zeros(m,3); fxyi = zeros(m,1); fxzi = zeros(m,1); fyzi = zeros(m,1);
     
     fp = zeros(3,3,3,m);
      q = zeros(3,3,3,m);
      p = zeros(3,3,3,m);      
     Mx = zeros(3,3,3,m);
     My = zeros(3,3,3,m);
     Mz = zeros(3,3,3,m);
    Mxx = zeros(3,3,3,m);
    Myy = zeros(3,3,3,m);
    Mzz = zeros(3,3,3,m);     
    Mxy = zeros(3,3,3,m);
    Mxz = zeros(3,3,3,m);
    Myz = zeros(3,3,3,m);
    
% Find i,j and k: 
     
       i = bindex(xi,x,0);
       
       i( ( i == 0 ) ) = 1;
       
       i( ( i == nx - 1 )|( i == nx ) ) = nx - 2;
       
       ix(1,:) =   i  ; 
       ix(2,:) = i + 1;
       ix(3,:) = i + 2;
       
       j = bindex(yi,y,0);
       
       j( ( j == 0 ) ) = 1;
       
       j( ( j == ny - 1 )|( j == ny ) ) = ny - 2;
       
       jy(1,:) =   j  ; 
       jy(2,:) = j + 1;
       jy(3,:) = j + 2;
       
       k = bindex(zi,z,0);
       
       k( ( k == nz - 1 )|( k == nz ) ) = nz - 2;
       
       kz(1,:) =   k  ; 
       kz(2,:) = k + 1;
       kz(3,:) = k + 2;
       
% Pick up function values and calculate some auxiliary
% four-dimensional arrays:

       for ii = 1:3
           
	   i2 = [1:ii-1 ii+1:3];
	   
	   xii= x(ix(ii   ,:));
	   xb = x(ix(i2(1),:));
	   xa = x(ix(i2(2),:));
	   
	   px = ( xii - xb ).*( xii - xa );
	   pxi= ( xi  - xb ).*( xi  - xa );
	   sxi= 2*xi  - xb - xa;
	   
           for jj = 1:3
	   
	       j2 = [1:jj-1 jj+1:3];
	       
	       yjj= y(jy(jj   ,:));
	       yb = y(jy(j2(1),:));
	       ya = y(jy(j2(2),:));
	       
	       py = ( yjj - yb ).*( yjj - ya );
	       pyi= ( yi  - yb ).*( yi  - ya );
	       syi= 2*yi  - yb - ya;
	       
	       for kk = 1:3 
	       
		   k2 = [1:kk-1 kk+1:3];
	
		   zkk= z(kz(kk   ,:));
		   zb = z(kz(k2(1),:));
		   za = z(kz(k2(2),:));
	
		   pz = ( zkk - zb ).*( zkk - za );
		   pzi= ( zi  - zb ).*( zi  - za );
		   szi= 2*zi  - zb - za;
		  
		   ll = nx*ny*( kz(kk,:) - 1 ) + nx*( jy(jj,:) - 1 ) + ix(ii,:); 
		   		   
		   fp(ii,jj,kk,:) = f(ll);
		    q(ii,jj,kk,:) = ( px .*py  ).*pz ;
		    p(ii,jj,kk,:) = ( pxi.*pyi ).*pzi;
 		   Mx(ii,jj,kk,:) = ( sxi.*pyi ).*pzi;
		   My(ii,jj,kk,:) = ( syi.*pxi ).*pzi;
		   Mz(ii,jj,kk,:) = ( szi.*pxi ).*pyi;
		  Mxx(ii,jj,kk,:) =        2*pyi.*pzi;
		  Myy(ii,jj,kk,:) =        2*pxi.*pzi;
		  Mzz(ii,jj,kk,:) =        2*pxi.*pyi;
		  Mxy(ii,jj,kk,:) = ( sxi.*syi ).*pzi;
		  Mxz(ii,jj,kk,:) = ( sxi.*szi ).*pyi;
		  Myz(ii,jj,kk,:) = ( syi.*szi ).*pxi;
	      
	      end
	      
	   end
       
       end

% Let's get the coefficients: 
        
	a = fp./q;
	
% Now we can interpolate the function and its derivatives 
% at each node:
	    
       fi = squeeze( sum( sum( sum( a.*p  , 1 ), 2 ), 3 ) ); 
      fxi = squeeze( sum( sum( sum( a.*Mx , 1 ), 2 ), 3 ) );
      fyi = squeeze( sum( sum( sum( a.*My , 1 ), 2 ), 3 ) );
      fzi = squeeze( sum( sum( sum( a.*Mz , 1 ), 2 ), 3 ) );
     fxxi = squeeze( sum( sum( sum( a.*Mxx, 1 ), 2 ), 3 ) );
     fyyi = squeeze( sum( sum( sum( a.*Myy, 1 ), 2 ), 3 ) );
     fzzi = squeeze( sum( sum( sum( a.*Mzz, 1 ), 2 ), 3 ) );
     fxyi = squeeze( sum( sum( sum( a.*Mxy, 1 ), 2 ), 3 ) );
     fxzi = squeeze( sum( sum( sum( a.*Mxz, 1 ), 2 ), 3 ) );
     fyzi = squeeze( sum( sum( sum( a.*Myz, 1 ), 2 ), 3 ) );
         
      dfx = [ fxi  fyi  fzi]; 
     d2fx = [fxxi fyyi fzzi];
    d2fxy = [fxyi fxzi fyzi]; 

end 
