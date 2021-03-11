 clc,clear all
%savefile = 'eta_max1.mat';
% Solve ostrovski eqn on [-L,L] by  FFT with integrating factor 
nframes = 160;
kn = 723/745;    % less humps        
%kn = 871/745;
%mx = 235/67;
mx = 1;
chi = 2*(kn^2)/(12*kn^4+3);
a = (-1+sqrt(1+2*chi*mx))/(2*chi);
w0 = 1/kn - kn^3;
w2 = (2*kn^3)/(12*kn^4 +3);
w0kk = 2/(kn^3) - 6*kn;
cg = -3*kn^2 - 1/(kn^2);
K = a*sqrt(-w2/w0kk);
sigma = w2*(a^2)/2;
l = @(x) 2.*a.*sech(K.*x).*cos(kn.*x)+...
    2*a^2*chi.*((sech(K*x).^2).*cos(2*kn*x));
  xmax = 160; % 160
  xmin = -160;   
  N = 1024; 		% Number of x points: best is power of two
  xfac = (xmax-xmin)/(2*pi);
  xi = (2*pi/N)*(-N/2:N/2-1)';
  x = 0.5*(xmax+xmin) + xfac*xi;
  Ost = 1;
  
  kapwid = 20;
  kap0 = xmin + 3*kapwid;
  kapamp = 0.001; 
  
  
  kappa = kapamp*exp(-((x-kap0)/kapwid).^2) + ...
      kapamp*exp(-((x+kap0)/kapwid).^2); %Damping
  kap = [0:N/2-1 0 -N/2+1:-1]'; 
  k = kap/xfac;
  iom = 1i*k*cg + 1i*k.^3 - 1i*Ost./(k+eps);
  u = l(x);
  
  tmax = 15; %15
  
  dt = .4/N^2;
  dt = dt * xfac^2;
  E = exp(dt*iom/2); 
  E2 = E.^2;
 
  npltint = floor((tmax/nframes)/dt); 
  nmax = round(tmax/dt); 
  udata = u;
  tdata = 0; 
  v = fft(u); 
      
      
  for n = 1:nmax
    t = n*dt; 
    g = -.5i*dt*k;   
    u = real( ifft(     v    ) );
    a = g.*fft(u.^2) - fft(kappa .* u);
    u = real( ifft(E.*(v+a/2)) );
    b = g.*fft(u.^2) - fft(kappa .* u);  
    u = real( ifft(E.*v + b/2) );		% 4th-order
    c = g.*fft(u.^2) - fft(kappa .* u);     	% Runge-Kutta
    u = real( ifft(E2.*v+E.*c) );
    d = g.*fft(u.^2) - fft(kappa .* u);
    v = E2.*v + (E2.*a + 2*E.*(b+c) + d)/6;
    v(1) = 0;
    v(N/2+1) = 0;
   
    
    
    if mod(n,npltint) == 0 
      u = real(ifft(v)); 
      udata = [udata u]; 
      tdata = [tdata t];
    end
  end
 %%
 
 % Taking measurements
  udata = udata(:,end);
  tdata = 0;
  nframes = 320;
  tmax = 4;
  dt = .4/N^2;
  dt = dt * xfac^2;
  
  n_it = 1:250; % number of iterations
  cs = zeros(1,length(n_it)+1);
  cs(1) = cg;
  T = zeros(1,length(n_it));
  mc = zeros(1,length(n_it));
  npltint = floor((tmax/nframes)/dt); 
  nmax = round(tmax/dt);
  v = fft(udata); 
 
  h = waitbar(0, 'please wait . . . ' ) ;
  
  
  for j=1:length(n_it)
      
      tdata = 0;
      
      iom = 1i*k*cs(j) + 1i*k.^3 - 1i*Ost./(k+eps);
      E = exp(dt*iom/2); 
      E2 = E.^2;
      
      for n = 1:nmax
          t = n*dt;
          g = -.5i*dt*k;  
          u = real( ifft(     v    ) );
          a = g.*fft(u.^2) - fft(kappa .* u);
          u = real( ifft(E.*(v+a/2)) );
          b = g.*fft(u.^2) - fft(kappa .* u);   % 4th-order
          u = real( ifft(E.*v + b/2) );
          c = g.*fft(u.^2) - fft(kappa .* u);   	% Runge-Kutta
          u = real( ifft(E2.*v+E.*c) );
          d = g.*fft(u.^2) - fft(kappa .* u);
          v = E2.*v + (E2.*a + 2*E.*(b+c) + d)/6;
          v(1) = 0;
          v(N/2+1) = 0;
    
          if mod(n,npltint) == 0 
            u = real(ifft(v)); 
            udata = [udata u]; 
            tdata = [tdata t];
          end
      end
      
      
      sz = size(udata);
      max_eta = zeros(sz(2),1);
      max_x = zeros(sz(2),1);
      
      for i = 1:sz(2)
          [y2,p_x2] = max(udata(:,i));
          x2 = xmin + (x(2)-x(1))*(p_x2-1);
          x1 = xmin + (x(2)-x(1))*(p_x2-1-1);
          y1 = udata(p_x2-1,i);
          x3 = xmin + (x(2)-x(1))*(p_x2);
          y3 = udata(p_x2 +1,i);
          P = [x1^2, x1, 1;
              x2^2, x2, 1;
              x3^2, x3, 1];
          ys = [y1;y2;y3];
          cons = P\ys;
          a_0 =cons(1);
          b_0 = cons(2);
          c_0 = cons(3);
          df = @(x) a_0*x^2 + b_0*x + c_0;
          max_eta(i) = df(-b_0/(2*a_0));
          max_x(i) = -b_0/(2*a_0);
      end

      y_peaks = zeros(2,3); % row represents # of peak
      t_peaks = zeros(2,3);
      pos_onetwo = zeros(2,1);

      [on , pos_onetwo(1)] = min(abs(tdata-1));
      [tw , pos_onetwo(2)] = min(abs(tdata-2)); 
      
      [y_peaks(1,2) , pos_p2] = max(max_eta(pos_onetwo(1):pos_onetwo(2)));
      t_peaks(1,2) = tdata(pos_onetwo(1) + pos_p2 -1);
      t_peaks(1,1) = tdata(pos_onetwo(1) +pos_p2 - 2);
      t_peaks(1,3) = tdata(pos_onetwo(1) +pos_p2);
      y_peaks(1,1) = max_eta(pos_onetwo(1) + pos_p2 -2);
      y_peaks(1,3) = max_eta(pos_onetwo(1) + pos_p2);
      
      [y_peaks(2,2) , pos_p3] = max(max_eta(pos_onetwo(2):end));
      t_peaks(2,2) = tdata(pos_onetwo(2) + pos_p3 -1);
      t_peaks(2,1) = tdata(pos_onetwo(2) + pos_p3 -2);
      t_peaks(2,3) = tdata(pos_onetwo(2) + pos_p3);
      y_peaks(2,1) = max_eta(pos_onetwo(2) + pos_p3 -2);
      y_peaks(2,3) = max_eta(pos_onetwo(2) + pos_p3);
      

      P_p2 = [t_peaks(end-1,1)^2, t_peaks(end-1,1), 1;
         t_peaks(end-1,2)^2, t_peaks(end-1,2), 1;
         t_peaks(end-1,3)^2, t_peaks(end-1,3), 1];
      ys_p2 = [y_peaks(end-1,1); y_peaks(end-1,2); y_peaks(end-1,3)];
      cons_p2 = P_p2\ys_p2;
      a_p2 = cons_p2(1);
      b_p2 = cons_p2(2);
      c_p2 = cons_p2(3);
      t_max_p2 = -b_p2/(2*a_p2);
      y_max_p2 = a_p2*(-b_p2/(2*a_p2))^2 + b_p2*(-b_p2/(2*a_p2)) + c_p2;
     
      P_p3 = [t_peaks(end,1)^2, t_peaks(end,1), 1;
         t_peaks(end,2)^2, t_peaks(end,2), 1;
         t_peaks(end,3)^2, t_peaks(end,3), 1];
      ys_p3 = [y_peaks(end,1); y_peaks(end,2); y_peaks(end,3)];
      cons_p3 = P_p3\ys_p3;
      a_p3 = cons_p3(1);
      b_p3 = cons_p3(2);
      c_p3 = cons_p3(3);
      t_max_p3 = -b_p3/(2*a_p3);
      y_max_p3 = a_p3*(-b_p3/(2*a_p3))^2 + b_p3*(-b_p3/(2*a_p3)) + c_p3;
      
      mc(j) = y_max_p3;

      t_range = zeros(2,2);
      x_range = zeros(2,2);
     
      for i = 1:length(tdata)-2
          if (t_max_p2 >= tdata(i)) && (t_max_p2 <= tdata(i+1))
              t_range(1,1) = tdata(i);
              x_range(1,1) = max_x(i);
              t_range(1,2) = tdata(i+1);
              x_range(1,2) = max_x(i+1);
          end
          if (t_max_p3 >= tdata(i)) && (t_max_p3 <= tdata(i+1))
              t_range(2,1) = tdata(i);
              x_range(2,1) = max_x(i);
              t_range(2,2) = tdata(i+1);
              x_range(2,2) = max_x(i+1);
          end
      end
     
      P_x_0 = [t_range(1,1), 1;
         t_range(1,2), 1];
      x_0s = [x_range(1,1); x_range(1,2)];
      cons_x_0 = P_x_0\x_0s;
      d_0 = cons_x_0(1);
      e_0 = cons_x_0(2);
     
      x_m0 = @(t) d_0*t + e_0;
      x_0 = x_m0(t_max_p2);
     
      P_x_1 = [t_range(2,1), 1;
         t_range(2,2), 1];
      x_1s = [x_range(2,1); x_range(2,2)];
      cons_x_1 = P_x_1\x_1s;
      d_1 = cons_x_1(1);
      e_1 = cons_x_1(2);
     
      x_m1 = @(t) d_1*t + e_1;
      x_1 = x_m1(t_max_p3);
     
      T(j) = t_max_p3 - t_max_p2;
     
      cs(j+1) = (x_1 - x_0)/T(j) + cs(j);
      
      for i = 1:sz(2)
         if tdata(i) == t_peaks(2,2)
             udata = udata(:,i);
         end
      end
      
      v=fft(udata);          waitbar(j/(length(n_it)))
  end
  
    close(h)
    
figure(1)  
plot(n_it,T)
grid on
xlabel('Number of iterations')
ylabel('$T$','Interpreter','latex')
figure(2)
plot([n_it,n_it(end)+1],cs)
grid on
xlabel('Number of iterations')
ylabel('$s$','Interpreter','latex')
figure(3)
plot(n_it,mc)
grid on
ylabel('max $\eta$','Interpreter','latex')
xlabel('Number of iterations')

%save (savefile, 'udata')
%%
plot(tdata,max_x)
xlabel('$t$','Interpreter','latex')
ylabel('$x$','Interpreter','latex')