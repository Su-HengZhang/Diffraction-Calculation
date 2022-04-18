function g=fresnelas(f)
% Calculation fresnel diffraction integration by Angular Spectrum method
%   g=fresneldft(f)
    %   f stands for the matrix of the object light field in the calculation
    %     window on the back surface of the diffraction screen. The diffraction
    %     window is at the center of the calculation window. Suppose the number
    %     of rows and columns of f are both M , this function always
    %     think that the origin of the coordinates is at the center
    %     pixel of the calculation window
    %                       (floor(M/2)+1,floor(M/2)+1)
    %   g stands for the matrix of the Fresnel diffraction light field on
    %     the observation screen, the same size as f, the origin of the
    %     coordinate is also at the center pixel of the matrix

% The Fresnel diffraction integral formula ( ignor the constant phase
% factor exp(jkz) )
%
%       g(x,y)=f(x,y)(convolution)h(x,y)
%
% where f(x,y) stands for the object light field on the back surface of
% the diffraction window, g(x,y) stands for the Fresnel diffraction
% light field on the observation window, h(x,y) represents
%       h(x,y)=(1/(j\lambda z))exp(j\pi\frac{x^2+y^2}{\lambda z}
%
% In the Fourier spectral domain, the Fresnel integral formula can be
% expressed as
%       G(u,v)=F(u,v)H(u,v)
% where G(u,v) is the spectrum of g(x,y), F(u,v) is the spectrum of f(x,y),
% H(u,v) is the transfer function
%       H(u,v)=exp(-j\pi\lambda z(u^2+v^2))


%---Get the number of rows and colummns of f---%
[M,~]=size(f);


%---Calculate the spectrum of the object light field----%
% Compute the Fourier Transform using fft2
f=ifftshift(f);
F=fft2(f);

%---Generate the spectrum coordinate grid---%
% Set up range of variables
u=0:(M-1);

% Compute the indices for use in meshgrid
idu=find(u>=M/2);
u(idu)=u(idu)-M;
v=u;

% Compute the meshgrid arrays
[V,U]=meshgrid(v,u);

%---Generate frequency domain quadratic phase factor---%
H=exp(-1i*pi*(U.^2+V.^2)/M);

%----Compute the spectrum of the diffracted light field---%
% frequency domain filtering
G=F.*H;

%----Compute the diffracted light field-----%
g=ifft2(G);

%---Center the diffracted light field---%
g=fftshift(g);
