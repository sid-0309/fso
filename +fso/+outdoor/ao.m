classdef ao < handle
 
    properties

        % default values : change when necessary
        f_G = 60;
        f_3db = 60;
        r_0 = 6e-3;
        D = 12e-3;
        x = linspace(-6e-3,6e-3,101);
        y = linspace(-6e-3,6e-3,101);
        w_0 = 8e-6;
        lambda = 1.55e-6;
        f = 8.57e-5;
        w_alpha = (1.55e-6)*(8.57e-5)/(pi*8e-6);

    end
    
    methods

        

        function [coefficients, zernikeModes] = fitZernikeModes(obj,wavefront, N)
            if N > 40
                error('N cannot be greater than 40.');
            end
            
            [rows, cols] = size(wavefront);
            [Zx, Zy] = meshgrid(linspace(-1, 1, cols), linspace(-1, 1, rows));
            r = sqrt(Zx.^2 + Zy.^2);
            theta = atan2(Zy, Zx);
            
            mask = r <= 1;
            
            zernikeModes = cell(1, N);
            coefficients = zeros(1, N);
            
            for k = 1:N
                zernikeModes{k} = obj.zernikePolynomial(k, r, theta) .* mask;
            end
            
            wf_flat = wavefront(mask);
            zernike_flat = zeros(numel(wf_flat), N);
            
            for k = 1:N
                zernike_flat(:, k) = zernikeModes{k}(mask);
            end
            
            coefficients = zernike_flat \ wf_flat;
        end

        function zern = zernike(obj,r,t,n,m)
    
            if mod(n-m,2) == 1
                error('n-m must be even');
            end
            if n < 0
                error('n must both be positive')
            end
            if floor(n) ~= n || floor(m) ~= m
                error('n and m must both be integers')
            end
        
            if m < 0
                zern = -obj.zernike_radial(r,n,-m).*sin(m*t);
            else
                zern = obj.zernike_radial(r,n,m).*cos(m*t);
            end
        end

        function radial = zernike_radial(obj,r,n,m)
            if mod(n-m,2) == 1
                error('n-m must be even');
            end
            if n < 0 || m < 0
                error('n and m must both be positive in radial function')
            end
            if floor(n) ~= n || floor(m) ~= m
                error('n and m must both be integers')
            end
            if n == m
                radial = r.^n;
            elseif n - m == 2
                radial = n*obj.zernike_radial(r,n,n)-(n-1)*obj.zernike_radial(r,n-2,n-2);
            else
                H3 = (-4*((m+4)-2)*((m+4)-3)) / ((n+(m+4)-2)*(n-(m+4)+4));
                H2 = (H3*(n+(m+4))*(n-(m+4)+2)) / (4*((m+4)-1))  +  ((m+4)-2);
                H1 = ((m+4)*((m+4)-1) / 2)  -  (m+4)*H2  +  (H3*(n+(m+4)+2)*(n-(m+4))) / (8);
                radial = H1*obj.zernike_radial(r,n,m+4) + (H2+H3 ./ r.^2).*obj.zernike_radial(r,n,m+2);
                
                % Fill in NaN values that may have resulted from DIV/0 in prior
                % line. Evaluate these points directly (non-recursively) as they
                % are scarce if present.
                
                if sum(sum(isnan(radial))) > 0
                    [row, col] = find(isnan(radial));
                    c=1;
                    while c<=length(row)
                        x = 0;
                        for k = 0:(n-m)/2
                            ((-1)^k*factorial(n-k))/(factorial(k)*factorial((n+m)/2-k)*factorial((n-m)/2-k))*0^(n-2*k);
                            x = x + ((-1)^k*factorial(n-k))/(factorial(k)*factorial((n+m)/2-k)*factorial((n-m)/2-k))*0^(n-2*k);
                        end
                        radial(row(c),col(c)) = x;
                        c=c+1;
                    end
                end
        
            end
        
        end

        function Z = zernikePolynomial(obj,n,r,theta)

            switch n
                case 1, Z = obj.zernike(r,theta,0,0);
                case 2, Z = obj.zernike(r,theta,1,-1);
                case 3, Z = obj.zernike(r,theta,1,1);
                case 4, Z = obj.zernike(r,theta,2,0);
                case 5, Z = obj.zernike(r,theta,2,-2);
                case 6, Z = obj.zernike(r,theta,2,2);
                case 7, Z = obj.zernike(r,theta,3,-1);
                case 8, Z = obj.zernike(r,theta,3,1);
                case 9, Z = obj.zernike(r,theta,3,-3);
                case 10, Z = obj.zernike(r,theta,3,3);
                case 11, Z = obj.zernike(r,theta,4,0);
                case 12, Z = obj.zernike(r,theta,4,-2);
                case 13, Z = obj.zernike(r,theta,4,2);
                case 14, Z = obj.zernike(r,theta,4,-4);
                case 15, Z = obj.zernike(r,theta,4,4);
                case 16, Z = obj.zernike(r,theta,5,-1);
                case 17, Z = obj.zernike(r,theta,5,1);
                case 18, Z = obj.zernike(r,theta,5,-3);
                case 19, Z = obj.zernike(r,theta,5,3);
                case 20, Z = obj.zernike(r,theta,5,-5);
                case 21, Z = obj.zernike(r,theta,5,5);
                case 22, Z = obj.zernike(r,theta,6,0);
                case 23, Z = obj.zernike(r,theta,6,-2);
                case 24, Z = obj.zernike(r,theta,6,2);
                case 25, Z = obj.zernike(r,theta,6,-4);
                case 26, Z = obj.zernike(r,theta,6,4);
                case 27, Z = obj.zernike(r,theta,6,-6);
                case 28, Z = obj.zernike(r,theta,6,6);
                case 29, Z = obj.zernike(r,theta,7,-1);
                case 30, Z = obj.zernike(r,theta,7,1);
                case 31, Z = obj.zernike(r,theta,7,-3);
                case 32, Z = obj.zernike(r,theta,7,3);
                case 33, Z = obj.zernike(r,theta,7,-5);
                case 34, Z = obj.zernike(r,theta,7,5);
                case 35, Z = obj.zernike(r,theta,7,-7);
                case 36, Z = obj.zernike(r,theta,7,7);
                case 37, Z = obj.zernike(r,theta,8,0);
                case 38, Z = obj.zernike(r,theta,8,-2);
                case 39, Z = obj.zernike(r,theta,8,2);
                case 40, Z = obj.zernike(r,theta,8,-4);
                otherwise
                    error("Number of Zernike Modes > 40")
            end
        end



        function [x,y] =  receiver(obj, rx, ry, div_x, div_y)
            obj.x = linspace(-rx,rx,div_x);
            obj.y = linspace(-ry,ry,div_y);
            x = obj.x;
            y = obj.y;

        end

        function phi = received_wf(obj, div_x,div_y,fn,varargin)

            arguments
                obj 
                div_x (1,1) double = max(size(obj.x))
                div_y (1,1) double = max(size(obj.y))
                fn double = 0
            end
            arguments(Repeating)
                varargin
            end
            if nargin==5
                use_rand = varargin{1};
            else
                use_rand = true;
            end

            if size(fn) == [1 1]
                fn = zeros(div_x,div_y)+fn;
            end
            
            if use_rand
                phi = rand(div_x,div_y)*3e-2+fn;
            else
                phi = zeros(div_x,div_y)+fn;
            end
        end

        function E_A = receiver_field(obj,phi, x,y,D)

            arguments
                obj
                phi double
                x double = obj.x
                y double = obj.y
                D double = obj.D
            end

            E_A = zeros(size(x,2),size(y,2));
            for a = 1:length(x)
                for b = 1:length(y)
                    if (x(a)^2+y(b)^2 < (D/2)^2)
                        E_A(a,b) = exp(-1j*phi(a,b));
                    end
                end
            end
        end

        function F_A = of_mode_field(obj, x,y,w_alpha)

            arguments
                obj
                x double = obj.x
                y double = obj.y
                w_alpha double = obj.w_alpha
            end
            
            F_A = zeros(size(x,2),size(y,2));
            for a = 1:length(x)
                for b = 1:length(y)
                    F_A(a,b) = sqrt(2/(pi*w_alpha^2))*exp(-(x(a)^2+y(b)^2)/w_alpha^2);
                end
            end
        end

        % function [Pf,Pa,eta,int_ea,int_fa] = coupling_efficiency(obj, e_a, f_a,x,y)
        function eta = coupling_efficiency(obj, e_a, f_a,x,y)

            arguments
                obj
                e_a double
                f_a double
                x double = obj.x
                y double = obj.y
            end

            Pf = 0;
            int_ea = 0;
            int_fa = 0;

            ds = (x(2)-x(1)) * (y(2)-y(1));

            for a = 1:length(x)
                for b = 1:length(y)
                    Pf = Pf + conj(e_a(a,b))*f_a(a,b)*ds;
                    int_ea = int_ea + (abs(e_a(a,b))^2 *ds);
                    int_fa = int_fa + ds* (abs(f_a(a,b))^2);
                end
            end
            
            Pf = abs(Pf)^2;
            Pa = int_ea*int_fa;
            eta = (Pf)/(Pa);

        end

            
        function eta = ce(obj,phi,x,y,w_alpha)

            arguments
                obj
                phi double
                x double = obj.x
                y double = obj.y
                w_alpha double = obj.w_alpha
            end

            eta = 0;
            ds = (x(2)-x(1)) * (y(2)-y(1));
            for a = 1:size(x,2)
                for b = 1:size(y,2)
                    eta = eta + exp(-(x(a)^2 + y(b)^2)/w_alpha^2)*ds;
                end
            end
            eta = abs(eta*exp(-0.5*rms(rms(phi))^2))^2 * 8/(pi*obj.D*w_alpha)^2;
        end


        function eta_bar = average_ce(obj,Cj,x,y,w_alpha)

            arguments
                obj
                Cj double
                x double = obj.x
                y double = obj.y
                w_alpha double = obj.w_alpha
            end

            eta = 0;
            ds = (x(2)-x(1)) * (y(2)-y(1));
            for a = 1:size(x,2)
                for b = 1:size(y,2)
                    eta = eta + exp(-(x(a)^2 + y(b)^2)/w_alpha^2)*ds;
                end
            end
            eta = abs(eta * exp(-0.5*(Cj*((obj.D/obj.r_0)^(5/3)) + (obj.f_G/obj.f_3db)^(5/3))))^2;

            eta_bar = eta*8/(pi*obj.D*w_alpha)^2;
        end

        
    end

end

