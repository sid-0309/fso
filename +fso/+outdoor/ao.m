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
            if N > 180
                error('N cannot be greater than 180.');
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
                case 38, Z = obj.zernike(r,theta,8,2);
                case 39, Z = obj.zernike(r,theta,8,-2);
                case 40, Z = obj.zernike(r,theta,8,4);
                case 41, Z = obj.zernike(r,theta,8,-4);
                case 42, Z = obj.zernike(r,theta,8,6);
                case 43, Z = obj.zernike(r,theta,8,-6);
                case 44, Z = obj.zernike(r,theta,8,8);
                case 45, Z = obj.zernike(r,theta,8,-8);
                case 46, Z = obj.zernike(r,theta,9,1);
                case 47, Z = obj.zernike(r,theta,9,-1);
                case 48, Z = obj.zernike(r,theta,9,3);
                case 49, Z = obj.zernike(r,theta,9,-3);
                case 50, Z = obj.zernike(r,theta,9,5);
                case 51, Z = obj.zernike(r,theta,9,-5);
                case 52, Z = obj.zernike(r,theta,9,7);
                case 53, Z = obj.zernike(r,theta,9,-7);
                case 54, Z = obj.zernike(r,theta,9,9);
                case 55, Z = obj.zernike(r,theta,9,-9);
                case 56, Z = obj.zernike(r,theta,10,0);
                case 57, Z = obj.zernike(r,theta,10,2);
                case 58, Z = obj.zernike(r,theta,10,-2);
                case 59, Z = obj.zernike(r,theta,10,4);
                case 60, Z = obj.zernike(r,theta,10,-4);
                case 61, Z = obj.zernike(r,theta,10,6);
                case 62, Z = obj.zernike(r,theta,10,-6);
                case 63, Z = obj.zernike(r,theta,10,8);
                case 64, Z = obj.zernike(r,theta,10,-8);
                case 65, Z = obj.zernike(r,theta,10,10);
                case 66, Z = obj.zernike(r,theta,10,-10);
                case 67, Z = obj.zernike(r,theta,11,1);
                case 68, Z = obj.zernike(r,theta,11,-1);
                case 69, Z = obj.zernike(r,theta,11,3);
                case 70, Z = obj.zernike(r,theta,11,-3);
                case 71, Z = obj.zernike(r,theta,11,5);
                case 72, Z = obj.zernike(r,theta,11,-5);
                case 73, Z = obj.zernike(r,theta,11,7);
                case 74, Z = obj.zernike(r,theta,11,-7);
                case 75, Z = obj.zernike(r,theta,11,9);
                case 76, Z = obj.zernike(r,theta,11,-9);
                case 77, Z = obj.zernike(r,theta,11,11);
                case 78, Z = obj.zernike(r,theta,11,-11);
                case 79, Z = obj.zernike(r,theta,12,0);
                case 80, Z = obj.zernike(r,theta,12,2);
                case 81, Z = obj.zernike(r,theta,12,-2);
                case 82, Z = obj.zernike(r,theta,12,4);
                case 83, Z = obj.zernike(r,theta,12,-4);
                case 84, Z = obj.zernike(r,theta,12,6);
                case 85, Z = obj.zernike(r,theta,12,-6);
                case 86, Z = obj.zernike(r,theta,12,8);
                case 87, Z = obj.zernike(r,theta,12,-8);
                case 88, Z = obj.zernike(r,theta,12,10);
                case 89, Z = obj.zernike(r,theta,12,-10);
                case 90, Z = obj.zernike(r,theta,12,12);
                case 91, Z = obj.zernike(r,theta,12,-12);
                case 92, Z = obj.zernike(r,theta,13,1);
                case 93, Z = obj.zernike(r,theta,13,-1);
                case 94, Z = obj.zernike(r,theta,13,3);
                case 95, Z = obj.zernike(r,theta,13,-3);
                case 96, Z = obj.zernike(r,theta,13,5);
                case 97, Z = obj.zernike(r,theta,13,-5);
                case 98, Z = obj.zernike(r,theta,13,7);
                case 99, Z = obj.zernike(r,theta,13,-7);
                case 100, Z = obj.zernike(r,theta,13,9);
                case 101, Z = obj.zernike(r,theta,13,-9);
                case 102, Z = obj.zernike(r,theta,13,11);
                case 103, Z = obj.zernike(r,theta,13,-11);
                case 104, Z = obj.zernike(r,theta,13,13);
                case 105, Z = obj.zernike(r,theta,13,-13);
                case 106, Z = obj.zernike(r,theta,14,0);
                case 107, Z = obj.zernike(r,theta,14,2);
                case 108, Z = obj.zernike(r,theta,14,-2);
                case 109, Z = obj.zernike(r,theta,14,4);
                case 110, Z = obj.zernike(r,theta,14,-4);
                case 111, Z = obj.zernike(r,theta,14,6);
                case 112, Z = obj.zernike(r,theta,14,-6);
                case 113, Z = obj.zernike(r,theta,14,8);
                case 114, Z = obj.zernike(r,theta,14,-8);
                case 115, Z = obj.zernike(r,theta,14,10);
                case 116, Z = obj.zernike(r,theta,14,-10);
                case 117, Z = obj.zernike(r,theta,14,12);
                case 118, Z = obj.zernike(r,theta,14,-12);
                case 119, Z = obj.zernike(r,theta,14,14);
                case 120, Z = obj.zernike(r,theta,14,-14);
                case 121, Z = obj.zernike(r,theta,15,1);
                case 122, Z = obj.zernike(r,theta,15,-1);
                case 123, Z = obj.zernike(r,theta,15,3);
                case 124, Z = obj.zernike(r,theta,15,-3);
                case 125, Z = obj.zernike(r,theta,15,5);
                case 126, Z = obj.zernike(r,theta,15,-5);
                case 127, Z = obj.zernike(r,theta,15,7);
                case 128, Z = obj.zernike(r,theta,15,-7);
                case 129, Z = obj.zernike(r,theta,15,9);
                case 130, Z = obj.zernike(r,theta,15,-9);
                case 131, Z = obj.zernike(r,theta,15,11);
                case 132, Z = obj.zernike(r,theta,15,-11);
                case 133, Z = obj.zernike(r,theta,15,13);
                case 134, Z = obj.zernike(r,theta,15,-13);
                case 135, Z = obj.zernike(r,theta,15,15);
                case 136, Z = obj.zernike(r,theta,15,-15);
                case 137, Z = obj.zernike(r,theta,16,0);
                case 138, Z = obj.zernike(r,theta,16,2);
                case 139, Z = obj.zernike(r,theta,16,-2);
                case 140, Z = obj.zernike(r,theta,16,4);
                case 141, Z = obj.zernike(r,theta,16,-4);
                case 142, Z = obj.zernike(r,theta,16,6);
                case 143, Z = obj.zernike(r,theta,16,-6);
                case 144, Z = obj.zernike(r,theta,16,8);
                case 145, Z = obj.zernike(r,theta,16,-8);
                case 146, Z = obj.zernike(r,theta,16,10);
                case 147, Z = obj.zernike(r,theta,16,-10);
                case 148, Z = obj.zernike(r,theta,16,12);
                case 149, Z = obj.zernike(r,theta,16,-12);
                case 150, Z = obj.zernike(r,theta,16,14);
                case 151, Z = obj.zernike(r,theta,16,-14);
                case 152, Z = obj.zernike(r,theta,16,16);
                case 153, Z = obj.zernike(r,theta,16,-16);
                case 154, Z = obj.zernike(r,theta,17,1);
                case 155, Z = obj.zernike(r,theta,17,-1);
                case 156, Z = obj.zernike(r,theta,17,3);
                case 157, Z = obj.zernike(r,theta,17,-3);
                case 158, Z = obj.zernike(r,theta,17,5);
                case 159, Z = obj.zernike(r,theta,17,-5);
                case 160, Z = obj.zernike(r,theta,17,7);
                case 161, Z = obj.zernike(r,theta,17,-7);
                case 162, Z = obj.zernike(r,theta,17,9);
                case 163, Z = obj.zernike(r,theta,17,-9);
                case 164, Z = obj.zernike(r,theta,17,11);
                case 165, Z = obj.zernike(r,theta,17,-11);
                case 166, Z = obj.zernike(r,theta,17,13);
                case 167, Z = obj.zernike(r,theta,17,-13);
                case 168, Z = obj.zernike(r,theta,17,15);
                case 169, Z = obj.zernike(r,theta,17,-15);
                case 170, Z = obj.zernike(r,theta,17,17);
                case 171, Z = obj.zernike(r,theta,17,-17);
                case 172, Z = obj.zernike(r,theta,18,0);
                case 173, Z = obj.zernike(r,theta,18,2);
                case 174, Z = obj.zernike(r,theta,18,-2);
                case 175, Z = obj.zernike(r,theta,18,4);
                case 176, Z = obj.zernike(r,theta,18,-4);
                case 177, Z = obj.zernike(r,theta,18,6);
                case 178, Z = obj.zernike(r,theta,18,-6);
                case 179, Z = obj.zernike(r,theta,18,8);
                case 180, Z = obj.zernike(r,theta,18,-8);
                case 181, Z = obj.zernike(r,theta,18,10);
                case 182, Z = obj.zernike(r,theta,18,-10);
                case 183, Z = obj.zernike(r,theta,18,12);
                case 184, Z = obj.zernike(r,theta,18,-12);
                case 185, Z = obj.zernike(r,theta,18,14);
                case 186, Z = obj.zernike(r,theta,18,-14);
                case 187, Z = obj.zernike(r,theta,18,16);
                case 188, Z = obj.zernike(r,theta,18,-16);
                case 189, Z = obj.zernike(r,theta,18,18);
                case 190, Z = obj.zernike(r,theta,18,-18);

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

