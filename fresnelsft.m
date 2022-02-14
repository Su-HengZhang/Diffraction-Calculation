function g=fresnelsft(f)
    % Calculation Fresnel diffraction integration by Single discrete Fourier
    % Transform
    %   g=fresnelsft(f)
    %   f stands for the matrix of the object light field in the calculation
    %     window on the back surface of the diffraction screen. The diffraction
    %     window is at the center of the calculation window. Suppose the number
    %     of rows and columns of f are M and N, respectively, this function
    %     always think that the origin of the coordinates is at the center
    %     pixel of the calculation window
    %                       (floor(M/2)+1,floor(N/2)+1)
    %   g stands for the matrix of the Fresnel diffraction light field on
    %     the observation screen, the same size as f, the origin of the
    %     coordinate is also at the center pixel of the matrix

    % The Fresnel diffraction integral formula ( ignor the constant phase
    % factor exp(jkz) )
    %
    % g(x,y)=1/(i\lambda z)exp[i*pi\frac{x^2+y^2}{\lambda z}]
    %         \iint f(x0,y0)exp[i*pi\frac{x0^2+y0^2}{\lambda z}]
    %               exp[-j2pi\frac{xx0+yy0}{\lambda z}]dx0dy0
    %
    % where f(x,y) stands for the object light field on the back surface of
    % the diffraction screen, g(x,y) stands for the Fresnel diffraction
    % light field on the observation screen

    %%
    %------------Generate quadratic phase factor function------%
    %Set up range of variables
    [M,N]=size(f);
    x=0:(M-1);
    y=0:(N-1);
    %Compute the indices for use in meshgrid
    idx=find(x>=M/2);x(idx)=x(idx)-M;
    idy=find(y>=N/2);y(idy)=y(idy)-N;
    %Compute the meshgrid arrays
    [Y,X]=meshgrid(y,x);
    %Quadractic phase factor funtion on diffraction and observation screen
    Q=exp(1i*pi*(X.^2/M+Y.^2/N));

    %%
    %----Calculate the spectrum of the product of the object matrix-------%
    %----and the quadratic phase factor matrix on the diffraction screen--%
    fc=ifftshift(f);
    FQc=fft2(fc.*Q);
    %-----Multiply the spectrum by the QPF on the observation screen---%
    gc=FQc.*Q;
    %-----Centered the origin of the coordinate-----%
    g=(1/(1i*sqrt(M*N)))*fftshift(gc);
