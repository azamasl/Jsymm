
#######################################################
#           Newton with elimination
#######################################################
function Elimination_Newton_analytic_center_backtrackingLS(x,y,w,obj,fixed_stepsz,sp, itNum,prt,Ftol)
    val, ngx, ngy, ngw =0, 0,0,0 #norm of the gradinet for the primal, slack and dual var
    println("Newton method with backtracking line-search")
    n = length(x);m = length(w); H = I(n+2*m)
    F = [obj.∇xL(x,y,w);
         obj.∇yL(x,y,w);
        -obj.∇wL(x,y,w)]
    F_old = F
    #normF=LinearAlgebra.norm(F,Inf)
    normF=LinearAlgebra.norm(F)
    normFAll=[normF];z = [x;y;w];
    Zs=z;iter=0
    ALPHA = 0.1; BETA  = 0.5; #backtracking ls params
     val =obj.L(x,y,w)
     ngx = LinearAlgebra.norm(obj.∇xL(x,y,w))
     ngy = LinearAlgebra.norm(obj.∇yL(x,y,w))
     ngw = LinearAlgebra.norm(obj.∇wL(x,y,w))
     if prt==1
         println("L = $val")
         println("|∇xL| = $ngx")
         println("|∇yL| = $ngy")
         println("|∇wL| = $ngw")
         println("||F||_Inf = $normF")
         println("#################################End of iter $iter")
     end
    A, b= sp.A,sp.b
    while ( normF> Ftol ) && iter < itNum
        g   = obj.∇yL(x,y,w);
        Lyy = Diagonal(1 ./ (y .^2));
        rd  = g + w;
        rp  = A*x + y - b;
        res = [A'*w;rd;rp];
        dx  = (A'*Lyy*A)\(A'*(g - Lyy*rp));
        dy  = - A*dx - rp;
        dw  = - Lyy*dy - rd;
        ###backtracking ls
        F_old = F
        gama = 1;
        while(minimum(y+gama*dy) <= 0)
            gama = BETA*gama;
        end

        F = [obj.∇xL(x+gama*dx, y+gama*dy, w+gama*dw);
             obj.∇yL(x+gama*dx, y+gama*dy, w+gama*dw);
            -obj.∇wL(x+gama*dx, y+gama*dy, w+gama*dw)]

        while (norm(F) > (1-ALPHA*gama)*normF)
            gama = BETA*gama;
            F = [obj.∇xL(x+gama*dx, y+gama*dy, w+gama*dw);
                 obj.∇yL(x+gama*dx, y+gama*dy, w+gama*dw);
                -obj.∇wL(x+gama*dx, y+gama*dy, w+gama*dw)]
        end
        ###updating z
        p = [dx;dy;dw];
        s = gama*p;z = z+s;
        x = z[1:n];y = z[n+1:n+m];w = z[n+m+1:n+2*m];

        normF = norm(F);
        append!(normFAll, normF);Zs = hcat(Zs,z)
        val =obj.L(x,y,w)
        an_f = obj.f(b-A*x)
        ngx = LinearAlgebra.norm(obj.∇xL(x,y,w))
        ngy = LinearAlgebra.norm(obj.∇yL(x,y,w))
        ngw = LinearAlgebra.norm(obj.∇wL(x,y,w))
        iter=iter+1
        if prt==1
            println("t = $gama")
            println("L = $val")
            println("|∇xL| = $ngx")
            println("|∇yL| = $ngy")
            println("|∇wL| = $ngw")
            println("||F||_Inf = $normF")
            println("f(x) = $an_f")
            println("#################################End of iter $iter")
        end
    end
    append!(normFAll, normF)

    return x,y,w,iter,normFAll, val, normF,Zs
end


#######################################################
#           Newton with inversion
#######################################################
function Newton_analytic_center_backtrackingLS(x,y,w,obj,fixed_stepsz,sp, itNum,prt,Ftol)
    val, ngx, ngy, ngw =0, 0,0,0 #norm of the gradinet for the primal, slack and dual var
    println("Newton method with backtracking line-search")
    n = length(x);m = length(w); H = I(n+2*m)
    F = [obj.∇xL(x,y,w);
         obj.∇yL(x,y,w);
        -obj.∇wL(x,y,w)]
    F_old = F
    normF=LinearAlgebra.norm(F)
    normFAll=[normF];z = [x;y;w];
    Zs=z;iter=0
    ALPHA = 0.1; BETA  = 0.5; #backtracking ls params
    rho1 = 1e-15;
    Reg   = [rho1*I(n) zeros(n,2*m);
                      zeros(m,n+2m);
         zeros(m,n+m)   -rho1*I(m)];

         val =obj.L(x,y,w)
         ngx = LinearAlgebra.norm(obj.∇xL(x,y,w))
         ngy = LinearAlgebra.norm(obj.∇yL(x,y,w))
         ngw = LinearAlgebra.norm(obj.∇wL(x,y,w))
         if prt==1
             println("L = $val")
             println("|∇xL| = $ngx")
             println("|∇yL| = $ngy")
             println("|∇wL| = $ngw")
             println("||F||_Inf = $normF")
             println("#################################End of iter $iter")
         end

    while ( normF> Ftol ) && iter < itNum
        Hnt = (obj.∇F(x,y,w) +  Reg)\I(n+2*m)
        p = -Hnt*F
        ###backtracking ls
        F_old = F
        gama = 1;
        while(minimum(y+gama*p[n+1:n+m]) <= 0)
            gama = BETA*gama;
        end
        F = [obj.∇xL(x+gama*p[1:n], y+gama*p[n+1:n+m], w+gama*p[n+m+1:n+2*m]);
             obj.∇yL(x+gama*p[1:n], y+gama*p[n+1:n+m], w+gama*p[n+m+1:n+2*m]);
            -obj.∇wL(x+gama*p[1:n], y+gama*p[n+1:n+m], w+gama*p[n+m+1:n+2*m])]

        while (norm(F) > (1-ALPHA*gama)*normF)
            gama = BETA*gama;
            F = [obj.∇xL(x+gama*p[1:n], y+gama*p[n+1:n+m], w+gama*p[n+m+1:n+2*m]);
                 obj.∇yL(x+gama*p[1:n], y+gama*p[n+1:n+m], w+gama*p[n+m+1:n+2*m]);
                -obj.∇wL(x+gama*p[1:n], y+gama*p[n+1:n+m], w+gama*p[n+m+1:n+2*m])]
        end
        ###updating z
        s = gama*p;z = z+s;
        x = z[1:n];y = z[n+1:n+m];w = z[n+m+1:n+2*m];
        #normF = maximum(abs.(F))#
        normF = norm(F);
        append!(normFAll, normF);Zs = hcat(Zs,z)

        val =obj.L(x,y,w)
        ngx = LinearAlgebra.norm(obj.∇xL(x,y,w))
        ngy = LinearAlgebra.norm(obj.∇yL(x,y,w))
        ngw = LinearAlgebra.norm(obj.∇wL(x,y,w))
        iter=iter+1
        if prt==1
            println("t = $gama")
            println("L = $val")
            println("|∇xL| = $ngx")
            println("|∇yL| = $ngy")
            println("|∇wL| = $ngw")
            println("||F||_Inf = $normF")
            println("#################################End of iter $iter")
        end
    end
    append!(normFAll, normF)

    return x,y,w,iter,normFAll, val, normF,Zs
end

#######################################################
#         Good Broyden method with fixed-stepsize
#######################################################
function broyden_analytic_center_fixedstep(x,y,w,obj,fixed_stepsz,sp, itNum,prt,Ftol)
    println("Broyden method  for analytic center  with fixed step size")
    val, ngx, ngy, ngw =0, 0,0,0; #norm of the gradinet for the primal, slack and dual var

    n = length(x);m = length(w); H = I(n+2*m)
    F = [obj.∇xL(x,y,w);
         obj.∇yL(x,y,w);
        -obj.∇wL(x,y,w)]
    F_old = F;normF=LinearAlgebra.norm(F)
    normFAll=[normF];z = [x;y;w];Zs=z;iter=0
    ALPHA = 0.1; BETA  = 0.5; #backtracking ls params
    H = I(n+2*m);
    while ( normF> Ftol ) && iter < itNum
        p = -H*F;
        F_old = F;
        gama = fixed_stepsz;
        while(minimum(y+gama*p[n+1:n+m]) <= 0)
            gama = BETA*gama;
        end
        F = [obj.∇xL(x+gama*p[1:n], y+gama*p[n+1:n+m], w+gama*p[n+m+1:n+2*m]);
             obj.∇yL(x+gama*p[1:n], y+gama*p[n+1:n+m], w+gama*p[n+m+1:n+2*m]);
            -obj.∇wL(x+gama*p[1:n], y+gama*p[n+1:n+m], w+gama*p[n+m+1:n+2*m])]
        if prt==1
            println("t = $gama")
        end
        s = gama*p;z = z+s;
        x = z[1:n];y = z[n+1:n+m];w = z[n+m+1:n+2*m];
        normF = norm(F);append!(normFAll, normF);Zs = hcat(Zs,z)

        ### updating H
        Y = F - F_old # Y is the gradient differences
        stH = s'*H
        stHy = (stH*Y)[1]
        sty = s'*Y
        H = H + ((s-H*Y)*stH)/stHy

        val =obj.L(x,y,w)
        ngx = LinearAlgebra.norm(obj.∇xL(x,y,w))
        ngy = LinearAlgebra.norm(obj.∇yL(x,y,w))
        ngw = LinearAlgebra.norm(obj.∇wL(x,y,w))
        iter=iter+1
        if prt==1
            println("L = $val")
            println("|∇xL| = $ngx")
            println("|∇yL| = $ngy")
            println("|∇wL| = $ngw")
            println("||F|| = $normF")
            println("#################################End of iter $iter")
        end
    end
    append!(normFAll, normF)

    return x,y,w,iter,normFAll, val, normF,Zs
end

#######################################################
#              Jsymm with fixed-stepsize
#######################################################
function secant_inv_analytic_center_fixedstep(x,y,w,obj,fixed_stepsz,sp, itNum,prt,Ftol)
    println("Secant method  for analytic center  with fixed step size")
    val, ngx, ngy, ngw =0, 0,0,0; #norm of the gradinet for the primal, slack and dual var
    n = length(x);m = length(w);
    F = [obj.∇xL(x,y,w);
         obj.∇yL(x,y,w);
        -obj.∇wL(x,y,w)]
    F_old = F;normF=LinearAlgebra.norm(F)
    normFAll=[normF];z = [x;y;w];Zs=z;iter=0
    ALPHA = 0.1; BETA  = 0.5; #backtracking ls params

     H = I(n+2*m);
     val =obj.L(x,y,w)
     ngx = LinearAlgebra.norm(obj.∇xL(x,y,w))
     ngy = LinearAlgebra.norm(obj.∇yL(x,y,w))
     ngw = LinearAlgebra.norm(obj.∇wL(x,y,w))
     if prt==1
         println("L = $val")
         println("|∇xL| = $ngx")
         println("|∇yL| = $ngy")
         println("|∇wL| = $ngw")
         println("||F|| = $normF")
         println("#################################End of iter $iter")
     end

    while ( normF> Ftol ) && iter < itNum
        p = -H*F
        F_old = F
        gama = fixed_stepsz;
        while(minimum(y+gama*p[n+1:n+m]) <= 0)
            gama = BETA*gama;
        end
        F = [obj.∇xL(x+gama*p[1:n], y+gama*p[n+1:n+m], w+gama*p[n+m+1:n+2*m]);
             obj.∇yL(x+gama*p[1:n], y+gama*p[n+1:n+m], w+gama*p[n+m+1:n+2*m]);
            -obj.∇wL(x+gama*p[1:n], y+gama*p[n+1:n+m], w+gama*p[n+m+1:n+2*m])]

        ###updating z
        s = gama*p;z = z+s;
        x = z[1:n];y = z[n+1:n+m];w = z[n+m+1:n+2*m];
        normF = norm(F);append!(normFAll, normF);Zs = hcat(Zs,z)

        ### updating H
        Y = F - F_old # Y is y in my notes

        # Update:
        r = Y +gama*F_old #since Ds=-stepsize*F_old and  r = Y - D*s
        ns2 = (s'*s)[1]
        Js = [s[1:n+m];-s[n+m+1:n+2*m]];#J*s
        Jr = [r[1:n+m];-r[n+m+1:n+2*m]];#J*r
        α = (s'*Jr)[1]#α is supposed to be an scaler
        a = r .- α*Js/ns2
        Ha = H*a
        denom1 =ns2 + (s'*Ha)[1]
        Ainv = H - (Ha*s'*H)/denom1
        AinvJs = Ainv*Js
        denom3 = (Jr'*AinvJs)[1]
        denom2 =ns2 + denom3
        H = Ainv - (AinvJs*Jr'*Ainv)/denom2
        val =obj.L(x,y,w)
        ngx = LinearAlgebra.norm(obj.∇xL(x,y,w))
        ngy = LinearAlgebra.norm(obj.∇yL(x,y,w))
        ngw = LinearAlgebra.norm(obj.∇wL(x,y,w))
        iter=iter+1
        if prt==1
            println("t = $gama")
            println("L = $val")
            println("|∇xL| = $ngx")
            println("|∇yL| = $ngy")
            println("|∇wL| = $ngw")
            println("||F|| = $normF")
            println("#################################End of iter $iter")
        end
    end
    append!(normFAll, normF)

    return x,y,w,iter,normFAll, val, normF,Zs
end


#######################################################
# J-symm with line-search
#######################################################
# For the unconstrained analytic center y is the slack variable and w is the dual variable
function secant_inv_analytic_center_backtrackingLS(x,y,w,obj,fixed_stepsz,sp, itNum,prt,Ftol)
    println("Secant method for analytic center with backtracking line-search")
    val, ngx, ngy, ngw =0, 0,0,0; #norm of the gradinet for the primal, slack and dual var

    n = length(x);m = length(w); H = I(n+2*m)
    F = [obj.∇xL(x,y,w);
         obj.∇yL(x,y,w);
        -obj.∇wL(x,y,w)]
    F_old = F;normF=LinearAlgebra.norm(F)
    normFAll=[normF];z = [x;y;w];Zs=z;iter=0
    ALPHA = 0.1; BETA  = 0.5; #backtracking ls params
     H = I(n+2*m);
     val =obj.L(x,y,w)
     ngx = LinearAlgebra.norm(obj.∇xL(x,y,w))
     ngy = LinearAlgebra.norm(obj.∇yL(x,y,w))
     ngw = LinearAlgebra.norm(obj.∇wL(x,y,w))
     if prt==1
         println("L = $val")
         println("|∇xL| = $ngx")
         println("|∇yL| = $ngy")
         println("|∇wL| = $ngw")
         println("||F|| = $normF")
         println("#################################End of iter $iter")
     end

     A, b= sp.A,sp.b
    while ( normF> Ftol ) && iter < itNum
        p = -H*F
        ###backtracking ls
        gama = 1;
        while(minimum(y+gama*p[n+1:n+m]) <= 0)
            gama = BETA*gama;
        end
        F_temp = [obj.∇xL(x+gama*p[1:n], y+gama*p[n+1:n+m], w+gama*p[n+m+1:n+2*m]);
             obj.∇yL(x+gama*p[1:n], y+gama*p[n+1:n+m], w+gama*p[n+m+1:n+2*m]);
            -obj.∇wL(x+gama*p[1:n], y+gama*p[n+1:n+m], w+gama*p[n+m+1:n+2*m])]
        normF_temp=norm(F_temp);
        while (normF_temp > (1-ALPHA*gama)*normF && gama >= stepsize_tol)
            gama = BETA*gama;
            F_temp = [obj.∇xL(x+gama*p[1:n], y+gama*p[n+1:n+m], w+gama*p[n+m+1:n+2*m]);
                 obj.∇yL(x+gama*p[1:n], y+gama*p[n+1:n+m], w+gama*p[n+m+1:n+2*m]);
                -obj.∇wL(x+gama*p[1:n], y+gama*p[n+1:n+m], w+gama*p[n+m+1:n+2*m])];
            normF_temp=norm(F_temp);
        end
        s = gama*p;
        Y = F_temp -F;
        if  gama < stepsize_tol #reject the update
            append!(normFAll, normF);Zs = hcat(Zs,z)
        else                    #accept the update and move to the new point
            z = z+s; x = z[1:n];y = z[n+1:n+m];w = z[n+m+1:n+2*m];
            F_old = F; F = F_temp;normF = normF_temp;
            append!(normFAll, normF);Zs = hcat(Zs,z)
        end

        ### updating H via SWM
        r = Y +gama*F_old
        ns2 = LinearAlgebra.norm(s)^2
        Js = [s[1:n+m];-s[n+m+1:n+2*m]];#J*s
        Jr = [r[1:n+m];-r[n+m+1:n+2*m]];#J*r
        α = (s'*Jr)[1]#α is supposed to be an scaler
        a = r .- α*Js/ns2
        Ha = H*a
        denom1 =ns2 + (s'*Ha)[1]
        Ainv = H - (Ha*s'*H)/denom1
        AinvJs = Ainv*Js
        denom2 =ns2 + (Jr'*AinvJs)[1]
        H = Ainv - (AinvJs*Jr'*Ainv)/denom2

        val =obj.L(x,y,w)
        ngx = LinearAlgebra.norm(obj.∇xL(x,y,w))
        ngy = LinearAlgebra.norm(obj.∇yL(x,y,w))
        ngw = LinearAlgebra.norm(obj.∇wL(x,y,w))
        iter=iter+1
        if prt==1
            println("t = $gama")
            println("L = $val")
            println("|∇xL| = $ngx")
            println("|∇yL| = $ngy")
            println("|∇wL| = $ngw")
            println("||F|| = $normF")
            println("#################################End of iter $iter")
        end
    end
    append!(normFAll, normF)
    return x,y,w,iter,normFAll, val, normF,Zs
end
